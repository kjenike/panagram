"""Shared fixtures for the panagram test suite. See tests/README.md."""

import os

# Pin thread pools before numpy/pandas/numba/umap are imported anywhere.
for _var in (
	"OMP_NUM_THREADS",
	"OPENBLAS_NUM_THREADS",
	"MKL_NUM_THREADS",
	"NUMBA_NUM_THREADS",
):
	os.environ.setdefault(_var, "4")

import dataclasses
import importlib.metadata
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

import synthetic

REPO_ROOT = Path(__file__).resolve().parent.parent
K = 15  # small k for the synthetic genomes (must be >= synthetic.K_MIN)
ANCHOR_CORES = 2
EXTRA_BINARIES = ("kmc", "kmc_tools", "kmc_dump", "mash")


def pytest_configure(config):
	import panagram

	here = Path(panagram.__file__).resolve()
	if REPO_ROOT not in here.parents:
		pytest.exit(
			f"panagram imported from {here}, expected a path under {REPO_ROOT}. "
			f"Run with PYTHONPATH={REPO_ROOT}.",
			returncode=3,
		)


@pytest.fixture(autouse=True)
def _restore_cwd():
	"""panagram.index.Index chdirs into the index directory and never returns."""
	cwd = os.getcwd()
	try:
		yield
	finally:
		os.chdir(cwd)


def _find_extra_dir():
	"""Directory holding kmc, kmc_tools, mash and py_kmc_api, or None.

	The clone's panagram/extra has no binaries (they are built by setup.py). Look
	at $PANAGRAM_EXTRA_DIR first, then at an installed panagram distribution.
	"""
	candidates = []
	if os.environ.get("PANAGRAM_EXTRA_DIR"):
		candidates.append(Path(os.environ["PANAGRAM_EXTRA_DIR"]))
	try:
		dist = importlib.metadata.distribution("panagram")
		candidates.append(Path(str(dist.locate_file("panagram/extra"))))
	except importlib.metadata.PackageNotFoundError:
		pass
	for c in candidates:
		if all((c / b).is_file() for b in EXTRA_BINARIES) and list(c.glob("py_kmc_api*")):
			return c.resolve()
	return None


@pytest.fixture(scope="session")
def extra_dir():
	"""Path of a directory with the KMC and mash binaries; skips if none is installed."""
	d = _find_extra_dir()
	if d is None:
		pytest.skip("no installed panagram extra/ (kmc, mash, py_kmc_api); set PANAGRAM_EXTRA_DIR")
	for tool in ("samtools", "bgzip"):
		if shutil.which(tool) is None:
			pytest.skip(f"{tool} is not on PATH")
	return d


@pytest.fixture(scope="session")
def panagram_overlay(tmp_path_factory, extra_dir):
	"""A copy of the clone's panagram package whose extra/ links to the built binaries.

	Snakemake rules import panagram.index in fresh processes and read
	panagram.index.EXTRA_DIR, so monkeypatching in the test process is not enough.
	Returns an env dict (PYTHONPATH pointing at the overlay) for subprocesses. The
	code under test is the clone's, copied at session start.
	"""
	root = tmp_path_factory.mktemp("overlay")
	src = REPO_ROOT / "panagram"
	dst = root / "panagram"
	shutil.copytree(src, dst, ignore=shutil.ignore_patterns("__pycache__", "assets"))
	for f in extra_dir.iterdir():
		target = dst / "extra" / f.name
		if f.name == "__init__.py" or target.exists():
			continue
		target.symlink_to(f)
	env = dict(os.environ)
	env["PYTHONPATH"] = str(root)
	for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
		env[var] = os.environ[var]

	# py_kmc_api is linked against the system libstdc++; a conda python may load an
	# older one first. If the import fails, retry preloading a system libstdc++.
	probe = [sys.executable, "-c", "from panagram.extra import py_kmc_api"]
	if subprocess.run(probe, env=env, cwd=root, capture_output=True).returncode != 0:
		for lib in _system_libstdcxx():
			trial = dict(env, LD_PRELOAD=lib)
			if subprocess.run(probe, env=trial, cwd=root, capture_output=True).returncode == 0:
				env = trial
				break
		else:
			pytest.skip("py_kmc_api cannot be imported in a subprocess (libstdc++ mismatch?)")
	return env


def _system_libstdcxx():
	found = []
	for d in ("/usr/lib/x86_64-linux-gnu", "/usr/lib64", "/usr/lib", "/lib/x86_64-linux-gnu"):
		p = Path(d) / "libstdc++.so.6"
		if p.exists():
			found.append(str(p))
	return found


@dataclasses.dataclass
class Genomes:
	"""Files and ground truth of the synthetic genome set."""

	root: Path  # directory holding samples.tsv, FASTAS/, GFFS/
	samples_tsv: Path
	k: int
	names: list  # genome names in samples.tsv order (= bit order)
	data: synthetic.SyntheticSet
	fasta: dict  # genome -> Path
	gff: dict  # genome -> Path (no entry for the genome without annotation)

	def owners(self, genome, chrom):
		"""Per k-mer start: set of genome names sharing the k-mer by construction, or None if the layout does not determine it (boundary-straddling k-mers)."""
		return self.data.kmer_owners(genome, chrom, self.k)

	def expected_bits(self, genome, chrom):
		"""Matrix (k-mer starts x genomes) from the construction: 1/0 where determined, -1 where not."""
		import numpy as np

		own = self.owners(genome, chrom)
		out = np.full((len(own), len(self.names)), -1, dtype="int8")
		for p, o in enumerate(own):
			if o is not None:
				out[p] = [int(n in o) for n in self.names]
		return out

	def oracle_bits(self, genome, chrom):
		"""Brute-force membership matrix: bit j is 1 iff the canonical k-mer occurs anywhere in genome j.

		Independent of KMC and of the layout; fully determined for every position.
		"""
		import numpy as np

		if not hasattr(self, "_kmer_sets"):
			self._kmer_sets = {
				n: set().union(*(canonical_kmers(s, self.k) for s in self.data.seqs[n].values()))
				for n in self.names
			}
		return np.array(
			[
				[int(c in self._kmer_sets[n]) for n in self.names]
				for c in canonical_kmers(self.data.seqs[genome][chrom], self.k)
			],
			dtype="uint8",
		)

	@staticmethod
	def assert_matches_expected(got, expected):
		"""Assert `got` equals `expected` wherever expected is determined (>= 0)."""
		import numpy as np

		got, expected = np.asarray(got), np.asarray(expected)
		mask = expected >= 0
		np.testing.assert_array_equal(got[mask], expected[mask])


def canonical_kmers(seq, k):
	"""Canonical (lexicographically smaller of forward and reverse complement) k-mers of seq."""
	comp = str.maketrans("ACGT", "TGCA")
	out = []
	for i in range(len(seq) - k + 1):
		s = seq[i : i + k]
		out.append(min(s, s.translate(comp)[::-1]))
	return out


def _write_genomes(root, k):
	data = synthetic.build()
	(root / "FASTAS").mkdir(parents=True)
	(root / "GFFS").mkdir()
	fasta, gff, rows = {}, {}, []
	for name in synthetic.GENOMES:
		fasta[name] = root / "FASTAS" / f"{name}.fasta"
		synthetic.write_fasta(fasta[name], data.seqs[name])
		gff_rel = ""
		if name in data.genes:
			gff[name] = root / "GFFS" / f"{name}.gff"
			synthetic.write_gff3(gff[name], data.genes[name])
			gff_rel = f"GFFS/{name}.gff"
		rows.append((name, f"FASTAS/{name}.fasta", gff_rel))
	tsv = root / "samples.tsv"
	with open(tsv, "w") as f:
		f.write("name\tfasta\tgff\n")
		for r in rows:
			f.write("\t".join(r) + "\n")
	return Genomes(root, tsv, k, list(synthetic.GENOMES), data, fasta, gff)


@pytest.fixture(scope="session")
def synthetic_genomes(tmp_path_factory):
	"""Five synthetic genomes (FASTA + GFF3 + samples.tsv); epsilon has no GFF. Session-scoped, read-only."""
	return _write_genomes(tmp_path_factory.mktemp("genomes"), K)


def _run(cmd, cwd, env, timeout=600):
	r = subprocess.run(
		cmd, cwd=cwd, env=env, capture_output=True, text=True, timeout=timeout
	)
	if r.returncode != 0:
		raise RuntimeError(
			f"{' '.join(map(str, cmd))} failed ({r.returncode})\n"
			f"stdout:\n{r.stdout[-3000:]}\nstderr:\n{r.stderr[-3000:]}"
		)
	return r


# panagram/__main__.py has no `if __name__ == "__main__"` guard, so `python -m panagram`
# does nothing; call main() the way the console script does.
PANAGRAM_CLI = [sys.executable, "-c", "import sys; from panagram.__main__ import main; sys.argv[0] = 'panagram'; main()"]


def panagram_cli(*args):
	"""Command line (list) that runs `panagram <args>` from the code in the overlay."""
	return [*PANAGRAM_CLI, *map(str, args)]


def _prepare(genomes, dest, env):
	"""Copy the genomes into `dest` and run `panagram index samples.tsv -k K --prepare`."""
	shutil.copytree(genomes.root, dest)
	_run(
		panagram_cli("index", "samples.tsv", "-k", genomes.k, "--prepare"),
		dest,
		env,
	)
	return dest


@pytest.fixture(scope="session")
def prepared_index(synthetic_genomes, panagram_overlay, tmp_path_factory):
	"""Index directory after `panagram index samples.tsv -k 15 --prepare` (config.yaml, samples.tsv, Snakefile; not anchored)."""
	dest = tmp_path_factory.mktemp("prepared") / "index"
	return _prepare(synthetic_genomes, dest, panagram_overlay)


@pytest.fixture
def fresh_prepared_index(synthetic_genomes, panagram_overlay, tmp_path):
	"""A private, mutable prepared index (function-scoped) for tests that change files."""
	return _prepare(synthetic_genomes, tmp_path / "index", panagram_overlay)


@pytest.fixture(scope="session")
def anchored_index(synthetic_genomes, panagram_overlay, tmp_path_factory):
	"""Complete index directory (Path) built with `snakemake --cores 2 all`. Built once; treat as read-only.

	Open it with `panagram.index.Index(str(path))` (read mode). Index chdirs into
	the directory; the autouse fixture restores the working directory afterwards.
	"""
	dest = tmp_path_factory.mktemp("anchored") / "index"
	_prepare(synthetic_genomes, dest, panagram_overlay)
	_run(
		[sys.executable, "-m", "snakemake", "--cores", str(ANCHOR_CORES), "all"],
		dest,
		panagram_overlay,
	)
	return dest


def open_index(path):
	"""Open an anchored index directory in read mode and return the Index.

	Index.__post_init__ calls os.chdir(path); this restores the working directory.
	Open file handles stay valid, relative paths inside the Index do not.
	"""
	from panagram.index import Index

	cwd = os.getcwd()
	try:
		return Index(str(path))
	finally:
		os.chdir(cwd)


@pytest.fixture(scope="session")
def _session_reader(anchored_index):
	idx = open_index(anchored_index)
	yield idx
	idx.close()


@pytest.fixture
def anchored_reader(_session_reader, anchored_index, monkeypatch):
	"""A read-mode panagram.index.Index over `anchored_index` (shared object; do not close it).

	The working directory is the index directory for the test: pysam reopens the
	tabix files by relative path when genes or annotations are queried.
	"""
	monkeypatch.chdir(anchored_index)
	return _session_reader


@pytest.fixture(scope="session")
def run_panagram(panagram_overlay):
	"""Callable run_panagram(*args, cwd=None) -> CompletedProcess for the panagram CLI (raises on failure)."""

	def run(*args, cwd=None):
		return _run(panagram_cli(*args), cwd or os.getcwd(), panagram_overlay)

	return run

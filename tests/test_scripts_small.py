"""Small scripts: pairwise_comp, make_bins_bits, write_umaps, query_index, run_umaps.sh, preprocess.sh.

Expected values come from the brute-force k-mer oracle (`Genomes.oracle_bits`) summed or
binned by a few lines in the test, or are worked by hand.
"""

import io
import os
import re
import shutil
import stat
import subprocess
import types

os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
import pandas as pd
import pytest

from helpers_scripts import SCRIPTS, exec_script, gene_rows, oracle_by_chrom, run_script


@pytest.fixture(scope="module")
def oracles(synthetic_genomes):
	return {g: oracle_by_chrom(synthetic_genomes, g) for g in ("alpha", "gamma")}


def chrom_kmers(genomes, genome):
	return {c: len(s) - genomes.k + 1 for c, s in genomes.data.seqs[genome].items()}


class _RecordingIndex:
	"""Stand-in for panagram.index.Index: two genomes, one chromosome of `size` k-mers, bitmap rows all (1, 1); records query_bitmap calls."""

	size = 250
	calls = []

	def __init__(self, path):
		self.genomes = {"a": None, "b": None}
		self.genome_names = ["a", "b"]

	def __getitem__(self, genome):
		return types.SimpleNamespace(chrs=pd.DataFrame({"size": [self.size]}, index=["c1"]))

	def query_bitmap(self, genome, chrom, start, end, step):
		type(self).calls.append((genome, chrom, start, end, step))
		return pd.DataFrame(np.ones((len(range(start, end, step)), 2), dtype=int), columns=["a", "b"])


@pytest.fixture
def recording_index(monkeypatch):
	import panagram.index as pi

	_RecordingIndex.calls = []
	monkeypatch.setattr(pi, "Index", _RecordingIndex)
	return _RecordingIndex


# ------------------------------------------------------------------- pairwise_comp.py


@pytest.mark.parametrize("anchor", ["alpha", "gamma"])
def test_pairwise_comp_counts_shared_kmers_at_step_100(
	panagram_overlay, anchored_index, synthetic_genomes, oracles, tmp_path, anchor
):
	"""Per genome: number of the anchor's step-100 k-mers (all chromosomes) present in that genome, and percent of the anchor's own."""
	r = run_script(panagram_overlay, "pairwise_comp.py", anchor, anchored_index, cwd=tmp_path)
	rows = [line.split(",") for line in r.stdout.splitlines()]
	assert [x[0] for x in rows] == synthetic_genomes.names
	sums = sum(m[::100].sum(axis=0) for m in oracles[anchor].values())
	assert [int(x[1]) for x in rows] == sums.tolist()
	own = sums[synthetic_genomes.names.index(anchor)]
	np.testing.assert_allclose([float(x[2]) for x in rows], sums / own * 100, rtol=1e-12)
	assert float(rows[synthetic_genomes.names.index(anchor)][2]) == 100.0


def test_pairwise_comp_orders_genomes_by_shared_content(panagram_overlay, anchored_index, tmp_path):
	"""By construction alpha shares more with beta (CORE, PAIR, DIV) than with epsilon (CORE only)."""
	r = run_script(panagram_overlay, "pairwise_comp.py", "alpha", anchored_index, cwd=tmp_path)
	pct = {l.split(",")[0]: float(l.split(",")[2]) for l in r.stdout.splitlines()}
	assert pct["alpha"] > pct["beta"] > pct["gamma"] > pct["epsilon"] > 0


def test_pairwise_comp_queries_the_whole_chromosome_at_step_100(monkeypatch, capsys, recording_index):
	"""The bitmap is read from k-mer 0 up to the chromosome's size (here 201, so the row at 200 is included)."""
	recording_index.size = 201
	exec_script("pairwise_comp", ["a", "somewhere"], monkeypatch)
	assert recording_index.calls == [("a", "c1", 0, 201, 100)]
	assert capsys.readouterr().out.splitlines() == ["a,3,100.0", "b,3,100.0"]


# -------------------------------------------------------------------- make_bins_bits.py


def test_make_bars_by_hand(monkeypatch, anchored_index, tmp_path):
	"""Counts 1 2 1 0 2 1 at positions 0,100,...,500 with 300-wide bins: bin 0 holds [1 2 1], bin 1 holds [0 2 1].

	Unique (count 1) is 2 and 1; universal (count == num_samples == 2) is 1 and 1.
	"""
	ns = _run_make_bins_bits(monkeypatch, anchored_index, tmp_path, [("chr1", 2650)], 500)
	ns["num_bins"] = 300
	x, y_univ, y_uniq, z_unique, z_univ = ns["make_bars"]([1, 2, 1, 0, 2, 1], 7)
	assert x == [300]
	assert y_univ == ["chr7_univ"] and y_uniq == ["chr7_unique"]
	assert z_unique == [2, 1]
	assert z_univ == [1, 1]
	# no counts of 1 or 2 anywhere
	assert ns["make_bars"]([0, 0, 3, 3], 1)[3:] == ([0, 0], [0, 0])


def _run_make_bins_bits(monkeypatch, index_dir, tmp_path, chroms, num_bins):
	"""Run make_bins_bits.py on `index_dir`; anchor "Spri1" (hard coded in the script) is mapped to alpha.

	The script bins in 200000 bp steps; `num_bins` replaces that constant so the 2.6 kb synthetic
	chromosomes span several bins.
	"""
	import panagram.index as pi

	real = pi.Index

	class Proxy:
		def __init__(self, path):
			self.index = real(path)

		def query_bitmap(self, genome, *args):
			assert genome == "Spri1"
			return self.index.query_bitmap("alpha", *args)

	monkeypatch.setattr(pi, "Index", Proxy)
	fai = tmp_path / "alpha.fasta.fai"
	fai.write_text("".join(f"{c}\t{n}\t0\t60\t61\n" for c, n in chroms))
	cwd = os.getcwd()
	try:
		ns, _ = exec_script(
			"make_bins_bits",
			[fai, index_dir],
			monkeypatch,
			edits=[("num_bins = 200000", f"num_bins = {num_bins}")],
		)
	finally:
		os.chdir(cwd)
	return ns


def test_make_bins_bits_prints_bin_edges_and_unique_and_universal_counts(
	monkeypatch, capsys, anchored_index, synthetic_genomes, oracles, tmp_path
):
	"""Three lines per chromosome: bin starts, count of 2-genome k-mers, count of 1-genome (unique) k-mers, per 500 bp bin."""
	fai_lengths = [("chr1", 2650), ("chr2", 1700)]
	_run_make_bins_bits(monkeypatch, anchored_index, tmp_path, fai_lengths, 500)
	lines = capsys.readouterr().out.splitlines()
	assert len(lines) == 6
	alpha = oracles["alpha"]
	for i, (chrom, length) in enumerate(fai_lengths):
		rows = (length - 1) // 100 + 1
		counts = alpha[chrom][::100].sum(axis=1)[:rows]
		# k-mers past the end of chr1 are not part of the chromosome: use the bitmap rows queried
		nx = (rows * 100 - 1) // 500
		xs = [500 * (j + 1) for j in range(nx)]
		uniq = [0] * (nx + 1)
		univ = [0] * (nx + 1)
		for pos, c in enumerate(counts):
			b = (pos * 100) // 500
			if c == 1:
				uniq[b] += 1
			elif c == 2:
				univ[b] += 1
		xline, univ_line, uniq_line = lines[3 * i : 3 * i + 3]
		# the last, partial bin of z is not printed by the script (reported as a bug, not asserted)
		assert xline.startswith("".join(f"{v}," for v in xs))
		assert univ_line.startswith("".join(f"{v}," for v in univ[:nx]))
		assert uniq_line.startswith("".join(f"{v}," for v in uniq[:nx]))
	assert any(int(v) > 0 for v in lines[2].strip(",").split(","))


def test_make_bins_bits_queries_each_chromosome_up_to_its_fai_length(monkeypatch, tmp_path, recording_index):
	"""Chromosome i of the .fai is read as chr<i> from 0 to its length at step 100."""
	fai = tmp_path / "x.fai"
	fai.write_text("c1\t250\t0\t60\t61\nc2\t999\t0\t60\t61\n")
	exec_script("make_bins_bits", [fai, "somewhere"], monkeypatch)
	assert recording_index.calls == [("Spri1", "chr1", 0, 250, 100), ("Spri1", "chr2", 0, 999, 100)]


# ------------------------------------------------------------------ write_umaps.py


@pytest.fixture
def fake_umap(monkeypatch):
	"""Replace umap.UMAP by a recorder whose embedding of row i is (10 i, 0)."""
	umap = pytest.importorskip("umap")
	record = {}

	class Fake:
		def __init__(self, **kwargs):
			record["kwargs"] = kwargs

		def fit_transform(self, x):
			record["input"] = np.array(x)
			return np.column_stack([10.0 * np.arange(len(x)), np.zeros(len(x))])

	monkeypatch.setattr(umap, "UMAP", Fake)
	return record


def _expected_paircounts(oracle_alpha, chrom_sizes, bin_size, step):
	"""Rows = bins of each chromosome in order; columns = genomes: shared k-mers in the bin over the largest count in it."""
	rows, index = [], []
	for chrom in chrom_sizes:
		m = oracle_alpha[chrom][::step]
		bins = (np.arange(len(m)) * step) // bin_size
		for b in range(bins.max() + 1):
			counts = m[bins == b].sum(axis=0).astype(float)
			with np.errstate(invalid="ignore"):
				rows.append(np.nan_to_num(counts / counts.max()))
			index.append((chrom, b * bin_size))
	return np.array(rows), index


def _run_write_umaps(monkeypatch, anchored_index, capsys, *argv):
	monkeypatch.chdir(anchored_index)
	exec_script("write_umaps", argv, monkeypatch)
	return pd.read_csv(io.StringIO(capsys.readouterr().out), sep="\t")


def test_write_umaps_feeds_umap_the_binned_pair_counts(
	monkeypatch, capsys, anchored_index, synthetic_genomes, oracles, fake_umap
):
	"""The matrix given to UMAP is, per bin, each genome's shared k-mer count divided by the largest in the bin."""
	out = _run_write_umaps(monkeypatch, anchored_index, capsys, "alpha", "ALL", 500, 1, 4, 0.25, 100)
	expected, index = _expected_paircounts(
		oracles["alpha"], chrom_kmers(synthetic_genomes, "alpha"), 500, 100
	)
	np.testing.assert_allclose(fake_umap["input"], expected, rtol=1e-12)
	assert fake_umap["kwargs"]["n_neighbors"] == 4
	assert fake_umap["kwargs"]["min_dist"] == 0.25
	assert fake_umap["kwargs"]["n_components"] == 2
	assert list(out.columns) == ["chrom", "start", "end", "umap1", "umap2", "cluster"]
	assert list(zip(out["chrom"], out["start"])) == index
	assert (out["end"] - out["start"] == 500).all()
	# the anchor's own column is the maximum of every bin
	assert (expected[:, 0] == 1.0).all()
	# the fake embedding is (10 i, 0); with eps 1 every point is its own cluster
	np.testing.assert_allclose(out["umap1"], 10.0 * np.arange(len(out)))
	assert list(out["cluster"]) == list(range(len(out)))


def test_write_umaps_clusters_with_dbscan_eps(monkeypatch, capsys, anchored_index, fake_umap):
	"""Points 10 apart fall in one DBSCAN cluster when eps is 15 (min_samples 1)."""
	out = _run_write_umaps(monkeypatch, anchored_index, capsys, "alpha", "ALL", 500, 15, 4, 0.0, 100)
	assert set(out["cluster"]) == {0}
	assert len(out) == 10


def test_write_umaps_step_argument_thins_the_bitmap_rows(monkeypatch, capsys, anchored_index, synthetic_genomes, oracles, fake_umap):
	"""Step 200 reads every second row of the stored step-100 bitmap, so a 400 bp bin holds two rows."""
	out = _run_write_umaps(monkeypatch, anchored_index, capsys, "alpha", "ALL", 400, 1, 4, 0.0, 200)
	expected, index = _expected_paircounts(
		oracles["alpha"], chrom_kmers(synthetic_genomes, "alpha"), 400, 200
	)
	np.testing.assert_allclose(fake_umap["input"], expected, rtol=1e-12)
	assert list(zip(out["chrom"], out["start"])) == index


def test_write_umaps_treats_empty_bins_as_zero_pair_counts(monkeypatch, capsys, anchored_index, fake_umap):
	"""A bin in which no genome has any k-mer (0/0) reaches UMAP as a row of zeros, not NaN."""
	import panagram.index as pi

	real = pi.Index

	class Emptied:
		def __init__(self, path):
			self._index = real(path)

		def __getitem__(self, genome):
			return self._index[genome]

		def __getattr__(self, name):
			return getattr(self._index, name)

		def query_bitmap(self, *args, **kwargs):
			bitmap = self._index.query_bitmap(*args, **kwargs).copy()
			bitmap.iloc[:3] = 0  # the first 300 bp bin of the first chromosome
			return bitmap

	monkeypatch.setattr(pi, "Index", Emptied)
	_run_write_umaps(monkeypatch, anchored_index, capsys, "alpha", "ALL", 300, 1, 4, 0.0, 100)
	x = fake_umap["input"]
	assert not np.isnan(x).any()
	np.testing.assert_array_equal(x[0], np.zeros(x.shape[1]))
	assert (x[1] > 0).any()


def test_write_umaps_with_real_umap(panagram_overlay, anchored_index, synthetic_genomes, tmp_path):
	"""The real UMAP and DBSCAN give one row per bin and, with a huge eps, a single cluster."""
	pytest.importorskip("umap")
	r = run_script(
		panagram_overlay, "write_umaps.py", "alpha", "ALL", 200, 1e9, 4, 0.0, 100, cwd=anchored_index
	)
	out = pd.read_csv(io.StringIO(r.stdout), sep="\t")
	expected, index = _expected_paircounts(
		oracle_by_chrom(synthetic_genomes, "alpha"), chrom_kmers(synthetic_genomes, "alpha"), 200, 100
	)
	assert list(zip(out["chrom"], out["start"])) == index
	assert (out["end"] - out["start"] == 200).all()
	assert set(out["cluster"]) == {0}
	assert np.isfinite(out[["umap1", "umap2"]].to_numpy()).all()
	assert len(out) == 23  # 14 bins on chr1 + 9 on chr2


# ------------------------------------------------------------------- query_index.py


@pytest.fixture
def query_index_ns(monkeypatch, anchored_index):
	"""Namespace of query_index.py (run for the gene table) with the functions it defines.

	The module body runs at import, so the namespace is taken even if the body stops early.
	"""
	ns, _ = exec_script("query_index", [anchored_index, "alpha", "gene"], monkeypatch, tolerate=(Exception,))
	assert "find_gene_avgs" in ns and "find_cpg_avgs" in ns
	return ns


def _parse_avgs(text):
	rows = []
	for line in text.splitlines():
		f = line.split("\t")
		rows.append((f[0], int(f[1]), int(f[2]), float(f[3]), float(f[4]), *f[5:]))
	return rows


def test_find_gene_avgs_prints_mean_and_sd_of_genomes_sharing_each_kmer(query_index_ns, capsys, oracles, synthetic_genomes):
	"""Per element: mean and (n - 1) standard deviation over its positions of the number of genomes holding the k-mer."""
	genes = gene_rows(synthetic_genomes, "alpha")
	df = pd.DataFrame(
		[(c, s, e, n) for c, s, e, n in genes] + [("chr1", 700, 700, "empty")],
		columns=["chr", "start", "end", "attr"],
	)
	query_index_ns["find_gene_avgs"]([df.iloc[:3], df.iloc[3:]])
	rows = _parse_avgs(capsys.readouterr().out)
	assert [r[5] for r in rows] == [n for _, _, _, n in genes]  # zero-length element skipped
	for chrom, start, end, avg, sd, name in rows:
		counts = oracles["alpha"][chrom][start:end].sum(axis=1)
		assert avg == pytest.approx(counts.mean(), rel=1e-12)
		assert sd == pytest.approx(counts.std(ddof=1), rel=1e-9)
	# the CORE segment is carried by all five genomes: a mean near 5 (boundary k-mers lower it)
	core = {r[5]: r for r in rows}["alpha_chr1_CORE_1"]
	assert 4.5 < core[3] <= 5.0


def test_find_cpg_avgs_reads_regions_from_a_table(query_index_ns, capsys, tmp_path, oracles):
	"""Lines with tab-separated fields: label, chromosome, start, stop; blank lines are ignored."""
	f = tmp_path / "cpg.tsv"
	f.write_text("1\tchr1\t300\t900\n\n2\tchr2\t0\t500\n3\tchr1\t1200\t1500\n")
	query_index_ns["find_cpg_avgs"](f)
	rows = _parse_avgs(capsys.readouterr().out)
	assert [(r[0], r[1], r[2]) for r in rows] == [("chr1", 300, 900), ("chr2", 0, 500), ("chr1", 1200, 1500)]
	for chrom, start, end, avg, sd in [r[:5] for r in rows]:
		counts = oracles["alpha"][chrom][start:end].sum(axis=1)
		assert avg == pytest.approx(counts.mean(), rel=1e-12)
		assert sd == pytest.approx(counts.std(ddof=1), rel=1e-9)
	# region 1200-1500 lies in alpha's PAIR segment, shared with beta only (plus boundary effects)
	assert 1.9 < rows[2][3] < 2.1


# --------------------------------------------------------------------- run_umaps.sh

needs_bash = pytest.mark.skipif(shutil.which("bash") is None, reason="bash not available")


def _stub(path, body):
	path.write_text("#!/bin/sh\n" + body)
	path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


@needs_bash
def test_run_umaps_sh_runs_one_job_per_chromosome_then_all_then_pairwise(tmp_path):
	"""python is stubbed to echo its arguments: each chromosome of the .fai, then ALL, then each sample get their jobs and files."""
	stubs = tmp_path / "stubs"
	stubs.mkdir()
	_stub(stubs / "python", 'echo "$@"\n')
	work = tmp_path / "work"
	for name in ("alpha", "beta"):
		(work / "anchor" / name).mkdir(parents=True)
	fasta = work / "genome.fa"
	fasta.write_text(">c1\nA\n")
	(work / "genome.fa.fai").write_text("c1\t10\t0\t60\t61\nc2\t20\t0\t60\t61\n")
	(work / "samples.tsv").write_text("alpha\tFASTAS/a.fa\nbeta\tFASTAS/b.fa\textra\n")
	env = dict(os.environ, PATH=f"{stubs}{os.pathsep}{os.environ['PATH']}")
	r = subprocess.run(
		["bash", str(SCRIPTS / "run_umaps.sh"), "alpha", str(fasta), str(work)],
		cwd=work, env=env, capture_output=True, text=True, timeout=60,
	)
	assert r.returncode == 0
	assert r.stdout.split() == ["c1", "c2"]  # the script echoes each chromosome name
	tail = "100000 2 3 0 100"
	for chrom in ("c1", "c2", "ALL"):
		f = work / "anchor" / "alpha" / f"umap_clusters_{chrom}_100000_2_3_0_100.txt"
		toks = f.read_text().split()
		assert toks[1:] == ["alpha", chrom, *tail.split()], chrom
	for name in ("alpha", "beta"):
		f = work / "anchor" / name / f"perc_shared.{name}.txt"
		toks = f.read_text().split()
		assert toks[1:] == [name, str(work)]
	assert toks[0] == "pairwise_comp.py"


# ----------------------------------------------------------------------- preprocess.sh


def _time_v_works():
	try:
		r = subprocess.run(["/usr/bin/time", "-v", "true"], capture_output=True, timeout=10)
	except (OSError, subprocess.SubprocessError):
		return False
	return r.returncode == 0


needs_time_v = pytest.mark.skipif(
	shutil.which("bash") is None or not _time_v_works(), reason="needs bash and GNU /usr/bin/time -v"
)


def _run_preprocess(tmp_path, n_samples, k=21):
	"""Run preprocess.sh with its hard-coded KMC_DIR pointed at stubs that log their command lines."""
	kmc_dir = tmp_path / "kmc_root"
	(kmc_dir / "bin").mkdir(parents=True)
	(kmc_dir / "py_kmc_api").mkdir()
	(kmc_dir / "py_kmc_api" / "set_path.sh").write_text(":\n")
	log = tmp_path / "calls.log"
	for tool in ("kmc", "kmc_tools"):
		_stub(kmc_dir / "bin" / tool, f'echo "{tool} $*" >> "$STUB_LOG"\n')
	bin_dir = tmp_path / "bin"
	bin_dir.mkdir()
	for tool in ("panagram", "bgzip"):
		_stub(bin_dir / tool, f'echo "{tool} $*" >> "$STUB_LOG"\n')
	text = (SCRIPTS / "preprocess.sh").read_text()
	text, n = re.subn(r"(?m)^KMC_DIR=.*$", f'KMC_DIR="{kmc_dir}"', text)
	assert n == 1
	script = tmp_path / "preprocess.sh"
	script.write_text(text)
	samples = tmp_path / "samples.txt"
	samples.write_text("".join(f"g{i}\tFASTAS/g{i}.fa\n" for i in range(1, n_samples + 1)))
	out = tmp_path / "out"
	out.mkdir()
	env = dict(os.environ, PATH=f"{bin_dir}{os.pathsep}{os.environ['PATH']}", STUB_LOG=str(log))
	subprocess.run(
		["bash", str(script), str(samples), str(k), str(out)], env=env, cwd=tmp_path,
		capture_output=True, text=True, timeout=120, check=True,
	)
	return out, log.read_text().splitlines()


@needs_time_v
def test_preprocess_counts_each_genome_and_merges_into_one_bitvector(tmp_path):
	"""Three genomes: kmc per FASTA, set_counts 1, 2, 4, one opdef listing all three, one `complex` merge."""
	out, calls = _run_preprocess(tmp_path, 3, k=21)
	o = f"{out}/"
	kmc_calls = [c for c in calls if c.startswith("kmc ")]
	assert kmc_calls == [
		f"kmc -k21 -t12 -m64 -ci1 -cs10000000 -fm FASTAS/g{i}.fa {o}g{i}.cnt {o}tmp" for i in (1, 2, 3)
	]
	tool_calls = [c for c in calls if c.startswith("kmc_tools")]
	assert tool_calls[:3] == [
		f"kmc_tools -t12 transform {o}g{i}.cnt set_counts {2 ** (i - 1)} {o}g{i}.id" for i in (1, 2, 3)
	]
	assert tool_calls[3:] == [f"kmc_tools complex {o}opdef.0.txt"]
	assert (out / "opdef.0.txt").read_text() == (
		"INPUT:\n"
		f"g1 = {o}g1.id\ng2 = {o}g2.id\ng3 = {o}g3.id\n"
		"OUTPUT:\n"
		f"{o}bitvecs.0 = g1 + g2 + g3\n"
		"-ocsum\n"
	)
	assert (out / "tmp").is_dir()
	assert [c for c in calls if c.startswith(("panagram", "bgzip"))] == [
		f"panagram index {out.parent}/samples.txt {o}bitvecs.0  -k 21 -o {o}anchor".replace("  ", " "),
		f"bgzip -r {o}anchor.pank -I {o}anchor.panx",
		f"bgzip -r {o}anchor.100nt.pank -I {o}anchor.100nt.panx",
	]


@needs_time_v
def test_preprocess_starts_a_second_bitvector_after_32_genomes(tmp_path):
	"""35 genomes: counts tags 1..2^31 then restart at 1; opdef.0 merges g1..g32, opdef.1 merges g33..g35."""
	out, calls = _run_preprocess(tmp_path, 35)
	tags = [int(c.split()[-2]) for c in calls if "set_counts" in c]
	assert tags == [2**i for i in range(32)] + [1, 2, 4]
	o = f"{out}/"
	op0 = (out / "opdef.0.txt").read_text().splitlines()
	op1 = (out / "opdef.1.txt").read_text().splitlines()
	assert op0[1:33] == [f"g{i} = {o}g{i}.id" for i in range(1, 33)]
	assert op0[34] == f"{o}bitvecs.0 = " + " + ".join(f"g{i}" for i in range(1, 33))
	assert op1 == [
		"INPUT:", f"g33 = {o}g33.id", f"g34 = {o}g34.id", f"g35 = {o}g35.id",
		"OUTPUT:", f"{o}bitvecs.1 = g33 + g34 + g35", "-ocsum",
	]
	merges = [c for c in calls if c.startswith("kmc_tools complex")]
	assert merges == [f"kmc_tools complex {o}opdef.0.txt", f"kmc_tools complex {o}opdef.1.txt"]

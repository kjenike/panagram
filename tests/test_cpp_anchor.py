"""cpp/anchor.cpp: the C++ anchoring program `run_anchor`.

`run_anchor <ngenomes> <index dir> <name> <fasta> [<name> <fasta> ...]` reads the KMC bit vectors
in `<index dir>/kmc/bitvec*` (made by the Python workflow) and writes, per genome,
`anchor/<name>/{bitmap.1.gz, bitmap.1.gzi, bitmap.100.gz, bitmap.100.gzi, bitsum.bins.tsv, chrs.tsv}`.
cpp/Snakefile calls it that way after the KMC steps.

The program is built from cpp/anchor.cpp with the system C++ compiler and htslib when both are
present, and the committed cpp/run_anchor binary is also tested when it can run here; each test
is skipped for a binary that is not available. Expected bitmaps come from the brute-force
canonical k-mer oracle (`Genomes.oracle_bits`).
"""

import gzip
import os
import platform
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

from helpers_scripts import CPP, THREAD_VARS

NGENOMES = 5


def _compile(out_dir):
	"""Build run_anchor from cpp/anchor.cpp and KMC's API sources; returns the path or a skip reason string."""
	cxx = shutil.which("g++") or shutil.which("c++")
	if cxx is None:
		return "no C++ compiler"
	flags, libs = [], ["-lhts"]
	if shutil.which("pkg-config"):
		c = subprocess.run(["pkg-config", "--cflags", "htslib"], capture_output=True, text=True)
		l = subprocess.run(["pkg-config", "--libs", "htslib"], capture_output=True, text=True)
		if c.returncode == 0 and l.returncode == 0:
			flags, libs = c.stdout.split(), l.stdout.split()
	libs = libs + ["-lz", "-lbz2", "-llzma", "-lm", "-ldl", "-pthread"]
	out = out_dir / "run_anchor"
	cmd = [
		cxx, "-std=c++17", "-O1", "-fopenmp", *flags,
		"-o", str(out), str(CPP / "anchor.cpp"),
		*(str(CPP.parent / "KMC" / "kmc_api" / f) for f in ("mmer.cpp", "kmc_file.cpp", "kmer_api.cpp")),
		*libs,
	]
	try:
		r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
	except (OSError, subprocess.SubprocessError) as e:
		return f"compiler failed to run: {e}"
	if r.returncode != 0:
		return "anchor.cpp does not build here (htslib or OpenMP missing?)"
	return out


def _committed_binary_usable():
	exe = CPP / "run_anchor"
	if not exe.is_file():
		return "cpp/run_anchor is not present"
	if sys.platform != "linux" or platform.machine() not in ("x86_64", "AMD64"):
		return "the committed cpp/run_anchor is a Linux x86-64 executable"
	ldd = shutil.which("ldd")
	if ldd:
		r = subprocess.run([ldd, str(exe)], capture_output=True, text=True)
		if "not found" in r.stdout or r.returncode != 0:
			return "cpp/run_anchor needs shared libraries that are not installed"
	return exe


@pytest.fixture(scope="session", params=["source", "committed"])
def anchor_binary(request, tmp_path_factory):
	"""Path of a run_anchor binary: built from anchor.cpp, or the committed one."""
	if request.param == "source":
		got = _compile(tmp_path_factory.mktemp("anchor_build"))
	else:
		got = _committed_binary_usable()
	if isinstance(got, str):
		pytest.skip(got)
	return got


def run_anchor(binary, env, kmc_source, workdir, fastas, threads=2, check=True, ngenomes=NGENOMES):
	"""Run the binary on `fastas` ({name: path}) with the KMC database copied from `kmc_source` (an index dir)."""
	root = workdir / "root"
	if not (root / "kmc").exists():
		shutil.copytree(kmc_source / "kmc", root / "kmc")
	args = []
	for name, fasta in fastas.items():
		(root / "anchor" / name).mkdir(parents=True, exist_ok=True)
		args += [name, str(fasta)]
	env = dict(env, **{v: str(threads) for v in THREAD_VARS})
	r = subprocess.run(
		[str(binary), str(ngenomes), str(root), *args], env=env, capture_output=True, text=True, timeout=300
	)
	if check and r.returncode != 0:
		raise RuntimeError(f"run_anchor failed ({r.returncode})\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}")
	return root, r


@pytest.fixture(scope="session")
def cpp_run(anchor_binary, panagram_overlay, anchored_index, synthetic_genomes, tmp_path_factory):
	"""Root directory after run_anchor ran on the five synthetic genomes with two threads."""
	root, _ = run_anchor(
		anchor_binary, panagram_overlay, anchored_index, tmp_path_factory.mktemp("cpp_run"), synthetic_genomes.fasta
	)
	return root


def read_rows(path, nbytes=1, ngenomes=NGENOMES):
	"""Decompressed bitmap as an array of genome-bit rows (bit j of a row is genome j)."""
	raw = np.frombuffer(gzip.open(path, "rb").read(), dtype="uint8").reshape(-1, nbytes)
	return np.unpackbits(raw, axis=1, bitorder="little")[:, :ngenomes]


def chrom_order(genomes, genome):
	return list(genomes.data.seqs[genome])


# ------------------------------------------------------------------------- outputs


def test_chrs_tsv_lists_chromosomes_in_file_order_with_kmer_counts(cpp_run, synthetic_genomes):
	"""Columns name, id, size, gene_count; size is the number of k-mers (length - k + 1); gene_count is 0."""
	for genome in synthetic_genomes.names:
		df = pd.read_csv(cpp_run / "anchor" / genome / "chrs.tsv", sep="\t")
		assert list(df.columns) == ["name", "id", "size", "gene_count"]
		seqs = synthetic_genomes.data.seqs[genome]
		assert list(df["name"]) == list(seqs)
		assert list(df["id"]) == list(range(len(seqs)))
		assert list(df["size"]) == [len(s) - synthetic_genomes.k + 1 for s in seqs.values()]
		assert list(df["gene_count"]) == [0] * len(seqs)


def test_full_resolution_bitmap_equals_the_kmer_oracle(cpp_run, synthetic_genomes):
	"""Row i of a chromosome has bit j set iff the canonical k-mer at i occurs in genome j (also where the layout fixes it)."""
	for genome in synthetic_genomes.names:
		rows = read_rows(cpp_run / "anchor" / genome / "bitmap.1.gz")
		start = 0
		for chrom in chrom_order(synthetic_genomes, genome):
			oracle = synthetic_genomes.oracle_bits(genome, chrom)
			got = rows[start : start + len(oracle)]
			np.testing.assert_array_equal(got, oracle, err_msg=f"{genome} {chrom}")
			synthetic_genomes.assert_matches_expected(got, synthetic_genomes.expected_bits(genome, chrom))
			start += len(oracle)
		assert start == len(rows)


def test_step_100_bitmap_holds_every_hundredth_row(cpp_run, synthetic_genomes):
	"""bitmap.100 has rows at k-mer positions 0, 100, 200, ... of each chromosome (restarting at 0 per chromosome)."""
	for genome in synthetic_genomes.names:
		rows = read_rows(cpp_run / "anchor" / genome / "bitmap.100.gz")
		expected = np.vstack(
			[synthetic_genomes.oracle_bits(genome, c)[::100] for c in chrom_order(synthetic_genomes, genome)]
		)
		np.testing.assert_array_equal(rows, expected, err_msg=genome)


def test_bitsum_bins_count_kmers_by_number_of_genomes(cpp_run, synthetic_genomes):
	"""Rows chromosome index, bin start, then counts of k-mers present in 0..5 genomes.

	A chromosome with fewer than 100 * 200000 k-mers is cut into bins of nkmers // 100 k-mers
	(the last bin may be shorter), so each chromosome has 100 or 101 bins here.
	"""
	for genome in synthetic_genomes.names:
		df = pd.read_csv(cpp_run / "anchor" / genome / "bitsum.bins.tsv", sep="\t")
		assert list(df.columns) == ["chr", "start"] + [str(i) for i in range(NGENOMES + 1)]
		expected = []
		for ci, chrom in enumerate(chrom_order(synthetic_genomes, genome)):
			sums = synthetic_genomes.oracle_bits(genome, chrom).sum(axis=1)
			binlen = len(sums) // 100
			for start in range(0, len(sums), binlen):
				hist = np.bincount(sums[start : start + binlen], minlength=NGENOMES + 1)
				expected.append([ci, start, *hist])
		np.testing.assert_array_equal(df.to_numpy(), np.array(expected), err_msg=genome)
		total = sum(len(s) - synthetic_genomes.k + 1 for s in synthetic_genomes.data.seqs[genome].values())
		assert df[[str(i) for i in range(NGENOMES + 1)]].to_numpy().sum() == total


def test_gzi_index_gives_random_access_through_panagram(cpp_run, anchored_index, synthetic_genomes, tmp_path):
	"""panagram's reader, using the .gzi files written by the C++ program, returns the oracle rows of a window."""
	from conftest import open_index

	dest = tmp_path / "idx"
	shutil.copytree(anchored_index, dest)
	for genome in ("gamma", "delta"):
		for f in ("bitmap.1.gz", "bitmap.1.gzi", "bitmap.100.gz", "bitmap.100.gzi", "bitsum.bins.tsv", "chrs.tsv"):
			shutil.copy(cpp_run / "anchor" / genome / f, dest / "anchor" / genome / f)
	index = open_index(dest)
	try:
		for genome, chrom, lo, hi in [("gamma", "chr2", 100, 400), ("delta", "chr1", 900, 1300)]:
			oracle = synthetic_genomes.oracle_bits(genome, chrom)
			got = index.query_bitmap(genome, chrom, lo, hi, 1).to_numpy()
			np.testing.assert_array_equal(got, oracle[lo:hi], err_msg=f"{genome} {chrom} step 1")
			got100 = index.query_bitmap(genome, chrom, 0, len(oracle), 100).to_numpy()
			np.testing.assert_array_equal(got100, oracle[::100], err_msg=f"{genome} {chrom} step 100")
	finally:
		index.close()


# ----------------------------------------------------------------- inputs and options


def test_threads_do_not_change_the_output(anchor_binary, panagram_overlay, anchored_index, synthetic_genomes, tmp_path):
	"""One and three OpenMP threads write the same decompressed bitmaps, bins and chromosome tables."""
	roots = [
		run_anchor(anchor_binary, panagram_overlay, anchored_index, tmp_path / f"t{n}", synthetic_genomes.fasta, threads=n)[0]
		for n in (1, 3)
	]
	for genome in synthetic_genomes.names:
		for f in ("bitmap.1.gz", "bitmap.100.gz"):
			a = gzip.open(roots[0] / "anchor" / genome / f, "rb").read()
			b = gzip.open(roots[1] / "anchor" / genome / f, "rb").read()
			assert a == b, (genome, f)
		for f in ("bitsum.bins.tsv", "chrs.tsv"):
			assert (roots[0] / "anchor" / genome / f).read_text() == (roots[1] / "anchor" / genome / f).read_text()


def _alpha_chr1(genomes):
	return genomes.data.seqs["alpha"]["chr1"]


def test_fasta_header_description_and_unwrapped_sequence(
	anchor_binary, panagram_overlay, anchored_index, synthetic_genomes, tmp_path
):
	"""The chromosome name stops at the first space; sequence on one line, with no final newline, is read whole."""
	seq = _alpha_chr1(synthetic_genomes)[:1000]
	fa = tmp_path / "custom.fa"
	fa.write_text(">chr1 some description here\n" + seq[:500] + "\n" + seq[500:])
	root, _ = run_anchor(anchor_binary, panagram_overlay, anchored_index, tmp_path, {"custom": fa})
	df = pd.read_csv(root / "anchor" / "custom" / "chrs.tsv", sep="\t")
	assert df["name"].tolist() == ["chr1"]
	assert df["size"].tolist() == [1000 - synthetic_genomes.k + 1]
	rows = read_rows(root / "anchor" / "custom" / "bitmap.1.gz")
	np.testing.assert_array_equal(rows, synthetic_genomes.oracle_bits("alpha", "chr1")[: len(rows)])


def test_every_header_stops_at_its_first_space_and_chromosomes_stay_separate(
	anchor_binary, panagram_overlay, anchored_index, synthetic_genomes, tmp_path
):
	"""Two records with descriptions: names are cut at the first space, ids count 0, 1, and the rows keep the file order."""
	seq = _alpha_chr1(synthetic_genomes)
	fa = tmp_path / "two.fa"
	fa.write_text(f">first desc one\n{seq[:700]}\n>second desc two\n{seq[700:1500]}\n>third x\n{seq[1500:2200]}\n")
	root, _ = run_anchor(anchor_binary, panagram_overlay, anchored_index, tmp_path, {"two": fa})
	df = pd.read_csv(root / "anchor" / "two" / "chrs.tsv", sep="\t")
	k = synthetic_genomes.k
	assert df["name"].tolist() == ["first", "second", "third"]
	assert df["id"].tolist() == [0, 1, 2]
	assert df["size"].tolist() == [700 - k + 1, 800 - k + 1, 700 - k + 1]
	rows = read_rows(root / "anchor" / "two" / "bitmap.1.gz")
	assert len(rows) == sum(df["size"])


def test_a_base_n_clears_every_kmer_that_covers_it(anchor_binary, panagram_overlay, anchored_index, synthetic_genomes, tmp_path):
	"""N at position 600 removes the 15-mers starting at 586..600; all other rows keep the oracle value."""
	k = synthetic_genomes.k
	seq = _alpha_chr1(synthetic_genomes)[:1000]
	masked = seq[:600] + "N" + seq[601:]
	fa = tmp_path / "n.fa"
	fa.write_text(">chr1\n" + masked + "\n")
	root, _ = run_anchor(anchor_binary, panagram_overlay, anchored_index, tmp_path, {"withn": fa})
	rows = read_rows(root / "anchor" / "withn" / "bitmap.1.gz")
	expected = synthetic_genomes.oracle_bits("alpha", "chr1")[: 1000 - k + 1].copy()
	expected[600 - k + 1 : 601] = 0
	np.testing.assert_array_equal(rows, expected)


def test_too_many_arguments_is_an_error_that_writes_nothing(anchor_binary, panagram_overlay, anchored_index, synthetic_genomes, tmp_path):
	"""More than 2 * ngenomes + 1 trailing words exits with status 1 and an error message."""
	root = tmp_path / "root"
	shutil.copytree(anchored_index / "kmc", root / "kmc")
	(root / "anchor").mkdir()
	env = dict(panagram_overlay, OMP_NUM_THREADS="1")
	args = [str(anchor_binary), "1", str(root), "a", str(synthetic_genomes.fasta["alpha"]), "extra", "words", "here"]
	r = subprocess.run(args, env=env, capture_output=True, text=True, timeout=60)
	assert r.returncode == 1
	assert "Error" in r.stderr
	assert list((root / "anchor").iterdir()) == []


# ------------------------------------------------- more than 8 genomes: several bytes per row


def _many_genomes(synthetic_genomes, n):
	"""n genomes as {name: {chrom: seq}}: the five synthetic ones, then slices of them (every other one reverse-complemented)."""
	from synthetic import revcomp

	genomes = {name: dict(synthetic_genomes.data.seqs[name]) for name in synthetic_genomes.names}
	bases = [synthetic_genomes.data.seqs[g]["chr1"] for g in ("alpha", "gamma", "delta", "beta")]
	for i in range(len(genomes), n):
		base = bases[i % len(bases)]
		start = (i * 37) % 300
		seq = base[start : start + 700 + 11 * i]
		genomes[f"g{i}"] = {"chr1": revcomp(seq) if i % 2 else seq}
	return genomes


@pytest.fixture(scope="module", params=[12, 35])
def many_genomes_db(request, extra_dir, synthetic_genomes, tmp_path_factory):
	"""(genomes, fastas, kmc dir) for `request.param` genomes: kmc/bitvec<i> built the way the Snakefile builds them.

	32 genomes go into one database; genome j of a database is tagged 2**j (set_counts) and the
	tagged databases are summed with `kmc_tools complex ... -ocsum`.
	"""
	from concurrent.futures import ThreadPoolExecutor

	from synthetic import write_fasta

	n, k = request.param, synthetic_genomes.k
	genomes = _many_genomes(synthetic_genomes, n)
	names = list(genomes)
	workdir = tmp_path_factory.mktemp(f"many{n}")
	kmc = workdir / "kmc"
	kmc.mkdir()
	(workdir / "tmp").mkdir()
	fastas = {name: workdir / f"{name}.fa" for name in names}

	def count(j):
		name = names[j]
		write_fasta(fastas[name], genomes[name])
		scratch = workdir / "tmp" / name
		scratch.mkdir()
		for cmd in (
			[extra_dir / "kmc", f"-k{k}", "-t1", "-m2", "-ci1", "-cs1000", "-fm", fastas[name], kmc / f"{name}.count", scratch],
			[extra_dir / "kmc_tools", "-t1", "transform", kmc / f"{name}.count", "set_counts", 1 << (j % 32), kmc / f"{name}.onehot"],
		):
			subprocess.run(list(map(str, cmd)), check=True, capture_output=True, timeout=120)

	with ThreadPoolExecutor(2) as pool:
		list(pool.map(count, range(n)))
	for i in range(0, n, 32):
		block = names[i : i + 32]
		opdef = kmc / f"opdef{i // 32}.txt"
		opdef.write_text(
			"INPUT:\n"
			+ "".join(f"{b} = {kmc}/{b}.onehot\n" for b in block)
			+ f"OUTPUT:\n{kmc}/bitvec{i // 32} = {' + '.join(block)}\n-ocsum\n"
		)
		subprocess.run([str(extra_dir / "kmc_tools"), "-t1", "complex", str(opdef)], check=True, capture_output=True, timeout=120)
	return genomes, fastas, kmc


def test_more_than_eight_genomes_use_several_bytes_per_row(anchor_binary, many_genomes_db, synthetic_genomes, tmp_path):
	"""12 genomes need 2 bytes per row; 35 need 5 bytes read from two k-mer databases of 32 and 3 genomes.

	Row i of an anchored genome has bit j (little-endian across the bytes) set iff its canonical
	k-mer occurs in genome j; the step-100 bitmap and the per-bin count histogram follow from it.
	"""
	from conftest import canonical_kmers

	k = synthetic_genomes.k
	genomes, fastas, kmc = many_genomes_db
	names = list(genomes)
	n = len(names)
	kmer_sets = [set().union(*(canonical_kmers(s, k) for s in genomes[name].values())) for name in names]
	root = tmp_path / "root"
	shutil.copytree(kmc, root / "kmc", ignore=shutil.ignore_patterns("*.count*"))
	anchored = [names[0], names[n // 2], names[-1]]
	args = []
	for name in anchored:
		(root / "anchor" / name).mkdir(parents=True)
		args += [name, str(fastas[name])]
	env = dict(os.environ, **{v: "2" for v in THREAD_VARS})
	r = subprocess.run([str(anchor_binary), str(n), str(root), *args], env=env, capture_output=True, text=True, timeout=300)
	assert r.returncode == 0, r.stderr[-2000:]
	nbytes = (n + 7) // 8
	for name in anchored:
		chroms = genomes[name]
		oracle = np.array(
			[[int(c in ks) for ks in kmer_sets] for seq in chroms.values() for c in canonical_kmers(seq, k)],
			dtype="uint8",
		)
		out = root / "anchor" / name
		rows = read_rows(out / "bitmap.1.gz", nbytes, n)
		np.testing.assert_array_equal(rows, oracle, err_msg=name)
		assert oracle[:, names.index(name)].all()
		starts = np.cumsum([0] + [len(s) - k + 1 for s in chroms.values()])[:-1]
		oracle100 = np.vstack([oracle[a : b][::100] for a, b in zip(starts, [*starts[1:], len(oracle)])])
		np.testing.assert_array_equal(read_rows(out / "bitmap.100.gz", nbytes, n), oracle100, err_msg=f"{name} step 100")
		bins = pd.read_csv(out / "bitsum.bins.tsv", sep="\t")
		assert list(bins.columns) == ["chr", "start"] + [str(i) for i in range(n + 1)]
		expected = []
		for ci, (a, b) in enumerate(zip(starts, [*starts[1:], len(oracle)])):
			sums = oracle[a:b].sum(axis=1)
			binlen = len(sums) // 100
			expected += [[ci, s, *np.bincount(sums[s : s + binlen], minlength=n + 1)] for s in range(0, len(sums), binlen)]
		np.testing.assert_array_equal(bins.to_numpy(), np.array(expected), err_msg=name)

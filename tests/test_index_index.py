"""panagram.index.Index: configuration, samples table, properties, aggregate tables, bitmap-table helpers.

Write-mode tests build index directories in-process (no KMC). Read-mode tests use the anchored index and
compare aggregate tables with brute-force counts from the synthetic genomes.
"""

import dataclasses
import logging
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

import helpers_index as H
import synthetic
from conftest import K, open_index
from panagram import index as pidx


def _small_samples(tmp_path, rows, header="name\tfasta\tgff", fasta_len=60):
	"""Write samples.tsv in tmp_path (rows: tuples) and a tiny FASTA for every fasta cell; return its path."""
	lines = [header]
	for r in rows:
		lines.append("\t".join("" if c is None else str(c) for c in r))
		fasta = r[1] if len(r) > 1 else None
		if fasta:
			synthetic.write_fasta(tmp_path / fasta, {"c1": "ACGT" * (fasta_len // 4)})
	tsv = tmp_path / "samples.tsv"
	tsv.write_text("\n".join(lines) + "\n")
	return tsv


def _write_index(tsv, **kw):
	"""Index(tsv, ...) in write mode, restoring the working directory afterwards."""
	cwd = os.getcwd()
	try:
		return pidx.Index(str(tsv), **kw)
	finally:
		os.chdir(cwd)


# ----------------------------------------------------------------- module constants


def test_name_regex_accepts_letters_digits_underscore_dash():
	"""Genome names are one or more of A-Za-z0-9_-."""
	import re

	for ok in ["a", "Hap1_chr-2", "0", "A-b_C9"]:
		assert re.fullmatch(pidx.NAME_REGEX, ok), ok
	for bad in ["", "a b", "a.b", "a/b", "a:b", "é"]:
		assert not re.fullmatch(pidx.NAME_REGEX, bad), repr(bad)


def test_module_paths_point_into_the_package():
	"""SRC_DIR is the package directory, SNAKEFILE the workflow Snakefile (which exists), EXTRA_DIR the extra/ folder."""
	import panagram

	assert pidx.SRC_DIR == os.path.dirname(os.path.realpath(panagram.__file__))
	assert pidx.SNAKEFILE == os.path.join(pidx.SRC_DIR, "workflow", "Snakefile")
	assert os.path.isfile(pidx.SNAKEFILE)
	assert pidx.EXTRA_DIR == os.path.join(pidx.SRC_DIR, "extra")
	assert (pidx.BGZ_SUFFIX, pidx.IDX_SUFFIX, pidx.ANCHOR_DIR) == ("gz", "gzi", "anchor")
	assert pidx.MODES == {"r", "w"}


def test_init_logger_writes_formatted_records_to_the_file(tmp_path):
	"""init_logger sends INFO records to the given file as '[date time LEVEL] message'."""
	root = logging.getLogger()
	saved, level = root.handlers[:], root.level
	root.handlers[:] = []
	try:
		log = tmp_path / "run.log"
		pidx.init_logger(str(log))
		logging.getLogger("panagram.test").info("hello %d", 7)
		logging.getLogger("panagram.test").debug("hidden")
		for h in root.handlers:
			h.flush()
		text = log.read_text().splitlines()
	finally:
		for h in root.handlers[:]:
			h.close()
		root.handlers[:] = saved
		root.setLevel(level)
	import re

	assert len(text) == 1
	assert re.fullmatch(r"\[\d{4}-\d\d-\d\d \d\d:\d\d:\d\d INFO\] hello 7", text[0])


# ------------------------------------------------------------- parameter dataclasses


def test_kmc_and_umap_defaults():
	"""KMC and UMAP parameter objects carry their documented defaults."""
	assert dataclasses.asdict(pidx.KMC()) == {"memory": 8, "threads": 1, "use_existing": False}
	assert dataclasses.asdict(pidx.UMAP()) == {
		"neighbors": 4,
		"dist": 0,
		"eps": 1,
		"samples": 1,
		"bin_size": 100000,
	}


def test_index_defaults(tmp_path):
	"""A new write-mode Index has the documented parameter defaults."""
	tsv = _small_samples(tmp_path, [("a", "a.fa", None)])
	idx = _write_index(tsv)
	p = idx.params
	assert (p["k"], p["cores"], p["lowres_step"], p["max_bin_kbp"], p["min_bin_count"]) == (21, 1, 100, 200, 100)
	assert p["max_view_chrs"] == 50
	assert p["gff_gene_types"] == ["gene"] and p["gff_anno_types"] is None and p["gff_name"] == "Name"
	assert p["prepare"] is False and p["mode"] is None
	assert p["genome_umap"] == p["chrom_umap"] == dataclasses.asdict(pidx.UMAP())
	assert idx.steps == (1, 100)


# --------------------------------------------------------------------- init_config


def test_init_config_writes_samples_with_ids_and_anchors(tmp_path):
	"""samples.tsv gains an id (row order) and an anchor column; every genome with a FASTA is an anchor by default."""
	tsv = _small_samples(tmp_path, [("x", "x.fa", "x.gff"), ("y", "y.fa", None), ("z", "z.fa", None)])
	_write_index(tsv)
	df = pd.read_table(tmp_path / "samples.tsv")
	assert list(df.columns) == ["name", "fasta", "gff", "id", "anchor"]
	assert df["name"].tolist() == ["x", "y", "z"]
	assert df["id"].tolist() == [0, 1, 2]
	assert df["anchor"].tolist() == [True, True, True]
	assert df["gff"].isna().tolist() == [False, True, True]


def test_init_config_anchor_genomes_option_selects_anchors(tmp_path):
	"""anchor_genomes restricts the anchors; others get chrs=None and are not anchored but keep their bit id."""
	tsv = _small_samples(tmp_path, [("x", "x.fa", None), ("y", "y.fa", None), ("z", "z.fa", None)])
	idx = _write_index(tsv, anchor_genomes=["y"])
	df = pd.read_table(tmp_path / "samples.tsv")
	assert df["anchor"].tolist() == [False, True, False]
	assert [idx[n].anchored for n in "xyz"] == [False, True, False]
	assert idx["z"].id == 2
	assert idx["y"].chrs is not None and idx["x"].chrs is None


def test_init_config_without_gff_column_adds_empty_gff(tmp_path):
	"""A samples file with only name and fasta columns is accepted; every gff is missing."""
	tsv = _small_samples(tmp_path, [("x", "x.fa"), ("y", "y.fa")], header="name\tfasta")
	idx = _write_index(tsv)
	assert pd.read_table(tmp_path / "samples.tsv")["gff"].isna().all()
	assert not idx["x"].annotated


def test_init_config_drops_rows_without_fasta_and_gff(tmp_path):
	"""A row with neither FASTA nor GFF is dropped and the remaining ids stay consecutive."""
	tsv = _small_samples(tmp_path, [("x", "x.fa", None), ("empty", None, None), ("y", "y.fa", None)])
	_write_index(tsv)
	df = pd.read_table(tmp_path / "samples.tsv")
	assert df["name"].tolist() == ["x", "y"]
	assert df["id"].tolist() == [0, 1]


@pytest.mark.parametrize("header", ["name\tgff", "fasta\tgff", "genome\tfile"])
def test_init_config_requires_name_and_fasta_columns(tmp_path, header):
	"""Input without a name or a fasta column header raises ValueError naming both columns."""
	tsv = tmp_path / "samples.tsv"
	tsv.write_text(header + "\nx\ty\n")
	with pytest.raises(ValueError, match="'name' and 'fasta'"):
		_write_index(tsv)


@pytest.mark.parametrize("bad", ["a b", "a.b", "a/b", "x:y"])
def test_init_config_rejects_invalid_genome_names(tmp_path, bad):
	"""Names outside NAME_REGEX raise ValueError that lists the offending name."""
	tsv = tmp_path / "samples.tsv"
	tsv.write_text(f"name\tfasta\tgff\nok\tok.fa\t\n{bad}\tb.fa\t\n")
	with pytest.raises(ValueError, match="Invalid genome names") as e:
		_write_index(tsv)
	assert bad in str(e.value) and "'ok'" not in str(e.value)


# ----------------------------------------------------------- config.yaml round trip


def test_write_config_omits_prefix_and_load_config_restores_values(tmp_path, monkeypatch):
	"""config.yaml has every parameter except prefix; load_config copies nested values back into the dataclasses."""
	tsv = _small_samples(tmp_path, [("x", "x.fa", None)])
	idx = _write_index(tsv, k=13, max_bin_kbp=7, gff_gene_types=["gene", "mRNA"])
	monkeypatch.chdir(tmp_path)  # the Index keeps relative paths, valid inside the index directory
	idx.chrom_umap.neighbors = 9
	idx.kmc.memory = 3
	idx.write_config()
	cfg = yaml.safe_load((tmp_path / "config.yaml").read_text())
	assert "prefix" not in cfg
	assert cfg["k"] == 13 and cfg["max_bin_kbp"] == 7 and cfg["gff_gene_types"] == ["gene", "mRNA"]
	assert cfg["chrom_umap"]["neighbors"] == 9 and cfg["kmc"]["memory"] == 3
	idx.write_config(exclude=["prefix", "input"])
	assert "input" not in yaml.safe_load((tmp_path / "config.yaml").read_text())

	reopened = H.open_write(tmp_path)
	assert reopened.k == 13 and reopened.max_bin_kbp == 7
	assert reopened.chrom_umap.neighbors == 9 and reopened.kmc.memory == 3
	assert reopened.gff_gene_types == ["gene", "mRNA"]
	assert reopened.genome_umap.neighbors == 4  # untouched value keeps its default
	assert list(reopened.genome_names) == ["x"]


def test_load_dict_updates_nested_dataclasses_in_place(tmp_path):
	"""_load_dict sets plain values, recurses into dicts for dataclass fields and accepts dataclass instances."""
	tsv = _small_samples(tmp_path, [("x", "x.fa", None)])
	idx = _write_index(tsv)
	kmc_obj = idx.kmc
	idx._load_dict(idx, {"k": 9, "kmc": {"memory": 2}, "genome_umap": {"dist": 0.3, "bin_size": 5}})
	assert idx.k == 9
	assert idx.kmc is kmc_obj and idx.kmc.memory == 2 and idx.kmc.threads == 1
	assert idx.genome_umap.dist == 0.3 and idx.genome_umap.bin_size == 5 and idx.genome_umap.neighbors == 4
	new = pidx.UMAP(neighbors=11)
	idx._load_dict(idx, {"chrom_umap": new})
	assert idx.chrom_umap is new
	with pytest.raises(ValueError, match="kmc"):
		idx._load_dict(idx, {"kmc": 5})


def test_params_is_a_plain_dict_of_the_dataclass(tmp_path):
	"""Index.params returns nested plain dicts, not dataclass objects."""
	idx = _write_index(_small_samples(tmp_path, [("x", "x.fa", None)]), k=17)
	assert idx.params["k"] == 17
	assert isinstance(idx.params["kmc"], dict) and isinstance(idx.params["genome_umap"], dict)


# ----------------------------------------------------------------- __post_init__


def test_post_init_rejects_an_invalid_mode(tmp_path):
	"""mode other than None, 'r' or 'w' raises ValueError."""
	tsv = _small_samples(tmp_path, [("x", "x.fa", None)])
	with pytest.raises(ValueError):
		_write_index(tsv, mode="x")


def test_post_init_input_checks(tmp_path):
	"""Reading needs a directory; writing needs a samples file or an initialized directory."""
	with pytest.raises(ValueError, match="directory"):
		pidx.Index(str(tmp_path / "missing"), mode="r")
	tsv = _small_samples(tmp_path, [("x", "x.fa", None)])
	with pytest.raises(ValueError, match="directory"):
		pidx.Index(str(tsv), mode="r")
	with pytest.raises(ValueError, match="TSV or initialized directory"):
		pidx.Index(str(tmp_path / "missing"), mode="w")
	empty = tmp_path / "empty"
	empty.mkdir()
	with pytest.raises(ValueError, match="not initialized"):
		pidx.Index(str(empty), mode="w")


def test_post_init_detects_mode_from_input_type(tmp_path):
	"""A file input means write mode; the prefix is its directory and the cwd becomes the prefix."""
	tsv = _small_samples(tmp_path, [("x", "x.fa", None)])
	cwd = os.getcwd()
	idx = pidx.Index(str(tsv))
	here = os.getcwd()
	os.chdir(cwd)
	assert idx.write_mode is True
	assert Path(here).resolve() == tmp_path.resolve()
	assert idx.prefix == ""
	assert idx.config_fname == "config.yaml" and idx.samples_fname == "samples.tsv"
	assert idx.genome_dist_fname == "genome_dist.tsv"


def test_post_init_prefix_option_sets_the_output_directory(tmp_path):
	"""With -o/--prefix the config and samples are written to that directory (absolute FASTA paths keep working)."""
	out = tmp_path / "out"
	out.mkdir()
	synthetic.write_fasta(tmp_path / "x.fa", {"c1": "ACGT" * 15})
	tsv = tmp_path / "in.tsv"
	tsv.write_text(f"name\tfasta\nx\t{tmp_path / 'x.fa'}\n")
	idx = _write_index(tsv, prefix=str(out))
	assert (out / "config.yaml").is_file() and (out / "samples.tsv").is_file()
	assert not (tmp_path / "config.yaml").exists()
	assert idx["x"].chrs["size"].tolist() == [60 - 21 + 1]


def test_samples_without_required_columns_give_a_helpful_error(tmp_path):
	"""A directory whose samples.tsv lacks the id column raises ValueError telling the user to prepare first."""
	(tmp_path / "config.yaml").write_text("k: 15\n")
	(tmp_path / "samples.tsv").write_text("name\tfasta\tgff\nx\tx.fa\t\n")
	with pytest.raises(ValueError, match="missing required column") as e:
		pidx.Index(str(tmp_path), mode="w")
	assert "id" in str(e.value) and "--prepare" in str(e.value)
	assert "name" not in str(e.value).split("column(s):")[1].split(".")[0]


def test_reading_a_prepared_but_unanchored_index_fails(prepared_copy_dir):
	"""Read mode needs the anchor files; a prepared-only directory raises OSError instead of returning a half-built Index."""
	with pytest.raises(OSError):
		open_index(prepared_copy_dir)


@pytest.fixture
def prepared_copy_dir(synthetic_genomes, tmp_path):
	return H.prepared_copy(synthetic_genomes, tmp_path / "idx")


# ---------------------------------------------------------------------------- run


def test_run_prepare_writes_files_and_prints_the_next_step(synthetic_genomes, tmp_path, capsys, monkeypatch):
	"""run() with prepare writes config.yaml, samples.tsv and the Snakefile and does not call snakemake."""
	import snakemake.cli

	called = []
	monkeypatch.setattr(snakemake.cli, "main", lambda args: called.append(args))
	dest = H.prepared_copy(synthetic_genomes, tmp_path / "idx")
	out = capsys.readouterr().out
	assert called == []
	assert "Wrote config.yaml and samples.tsv" in out
	assert "Wrote Snakefile" in out
	assert "snakemake --cores 1 all" in out
	assert (dest / "Snakefile").read_text() == Path(pidx.SNAKEFILE).read_text()
	assert yaml.safe_load((dest / "config.yaml").read_text())["k"] == K


def test_run_keeps_an_existing_snakefile(synthetic_genomes, tmp_path, capsys):
	"""A second run does not overwrite a Snakefile the user edited."""
	dest = H.prepared_copy(synthetic_genomes, tmp_path / "idx")
	(dest / "Snakefile").write_text("# edited\n")
	cwd = os.getcwd()
	try:
		pidx.Index(str(dest / "samples.tsv"), k=K, prepare=True).run()
	finally:
		os.chdir(cwd)
	assert (dest / "Snakefile").read_text() == "# edited\n"
	assert "Using existing" in capsys.readouterr().out


def test_run_without_prepare_calls_snakemake_with_cores(synthetic_genomes, tmp_path, monkeypatch, capsys):
	"""Without prepare, run() invokes snakemake's CLI with --cores <cores> all."""
	import shutil
	import snakemake.cli

	called = []
	monkeypatch.setattr(snakemake.cli, "main", lambda args: called.append(list(args)))
	shutil.copytree(synthetic_genomes.root, tmp_path / "idx")
	cwd = os.getcwd()
	try:
		pidx.Index(str(tmp_path / "idx" / "samples.tsv"), k=K, cores=3).run()
	finally:
		os.chdir(cwd)
	assert called == [["--cores", "3", "all"]]
	assert "Running 'snakemake --cores 3 all'" in capsys.readouterr().out


# ----------------------------------------------------------- write-mode properties


@pytest.fixture
def wide_index(tmp_path, synthetic_genomes):
	"""Write-mode index with 33 genomes (the same FASTA 33 times) so two KMC bit vectors are needed."""
	fa = synthetic_genomes.fasta["alpha"]
	rows = [(f"g{i}", str(fa), None) for i in range(33)]
	tsv = _small_samples(tmp_path, [], header="name\tfasta\tgff")
	tsv.write_text("name\tfasta\tgff\n" + "".join(f"{n}\t{f}\t\n" for n, f, _ in rows))
	return _write_index(tsv, k=15)


def test_kmc_bitvec_count_is_genomes_over_32_rounded_up(wide_index, tmp_path):
	"""33 genomes need 2 bit vectors, with matching opdef and bitvec file names under kmc/."""
	assert wide_index.kmc_bitvec_count == 2
	assert wide_index.opdef_filenames == ["kmc/opdef0.txt", "kmc/opdef1.txt"]
	assert wide_index.bitvec_prefixes == ["kmc/bitvec0", "kmc/bitvec1"]
	assert wide_index.kmc_prefix("a", "onehot") == "kmc/a.onehot"
	assert wide_index.get_subdir("kmc") == "kmc"
	assert wide_index.snakefile == "Snakefile"


def test_kmc_bitvec_count_five_genomes_is_one(synthetic_genomes, tmp_path):
	"""The five-genome set needs one bit vector."""
	idx = _write_index(H.prepared_copy(synthetic_genomes, tmp_path / "i") / "samples.tsv", k=K)
	assert idx.kmc_bitvec_count == 1


def test_init_opdefs_writes_one_complex_definition_per_bit_vector(wide_index, tmp_path, monkeypatch):
	"""Each opdef file lists its genomes' onehot databases and ORs them into one bit vector with -ocsum.

	Genomes 0-31 go to opdef0, genome 32 alone to opdef1.
	"""
	monkeypatch.chdir(tmp_path)
	os.makedirs("kmc")
	wide_index.init_opdefs()
	first = Path("kmc/opdef0.txt").read_text()
	lines = first.splitlines()
	assert lines[0] == "INPUT:"
	assert lines[1:33] == [f"g{i} = kmc/g{i}.onehot" for i in range(32)]
	assert lines[33] == "OUTPUT:"
	assert lines[34] == "kmc/bitvec0 = g0" + "".join(f" + g{i}" for i in range(1, 32))
	assert lines[35] == "-ocsum" and len(lines) == 36
	assert Path("kmc/opdef1.txt").read_text() == "INPUT:\ng32 = kmc/g32.onehot\nOUTPUT:\nkmc/bitvec1 = g32\n-ocsum\n"


# ------------------------------------------------------------ read-mode aggregates


def test_index_read_mode_shape_and_chromosome_table(anchored_reader, synthetic_genomes):
	"""ngenomes, the (genome, chr) table and the genome size table reflect the FASTA files."""
	idx, g = anchored_reader, synthetic_genomes
	assert idx.write_mode is False and idx.ngenomes == 5
	assert list(idx.bitsum_index) == [0, 1, 2, 3, 4, 5]
	assert list(idx.genome_names) == synthetic.GENOMES
	assert list(idx.chrs.index.names) == ["genome", "chr"]
	for genome in g.names:
		for chrom, seq in g.data.seqs[genome].items():
			assert idx.chrs.loc[(genome, chrom), "size"] == len(seq) - K + 1
	want = {
		n: (sum(len(s) - K + 1 for s in g.data.seqs[n].values()), len(g.data.seqs[n])) for n in g.names
	}
	for n in g.names:
		assert idx.genome_sizes.loc[n, "length"] == want[n][0]
		assert idx.genome_sizes.loc[n, "chr_count"] == want[n][1]
	assert idx[synthetic.NO_GFF] is idx.genomes[synthetic.NO_GFF]


def test_index_totals_are_sums_of_the_brute_force_histograms(anchored_reader, synthetic_genomes):
	"""bitsum_totals / bitfreq_totals / the averages equal what the occupancy of every position gives."""
	idx, g = anchored_reader, synthetic_genomes
	totals = {}
	for genome in g.names:
		occ = np.concatenate([g.oracle_bits(genome, c).sum(axis=1) for c in g.data.seqs[genome]])
		totals[genome] = np.bincount(occ, minlength=6)
		np.testing.assert_array_equal(idx.bitsum_totals.loc[genome].to_numpy(), totals[genome])
		freq = totals[genome] / totals[genome].sum()
		np.testing.assert_allclose(idx.bitfreq_totals.loc[genome].to_numpy(), freq)
		np.testing.assert_allclose(idx.bitsum_totals_avg[genome], (freq * np.arange(6)).sum())
	assert list(idx.bitsum_totals_avg.index) == sorted(totals, key=lambda n: (totals[n] * np.arange(6)).sum() / totals[n].sum())
	assert list(idx.bitsum_totals_avg.values) == sorted(idx.bitsum_totals_avg.values)


def test_index_chromosome_tables_cover_every_anchor_chromosome(anchored_reader, synthetic_genomes):
	"""bitsum_chrs, bitfreq_chrs, bitsum_bins and their averages exist for every (genome, chr)."""
	idx, g = anchored_reader, synthetic_genomes
	for genome in g.names:
		for chrom in g.data.seqs[genome]:
			occ = g.oracle_bits(genome, chrom).sum(axis=1)
			hist = np.bincount(occ, minlength=6)
			np.testing.assert_array_equal(idx.bitsum_chrs.loc[(genome, chrom)].to_numpy(), hist)
			np.testing.assert_allclose(idx.bitfreq_chrs.loc[(genome, chrom)].to_numpy(), hist / hist.sum())
			np.testing.assert_allclose(idx.bitsum_chrs_avg[(genome, chrom)], (hist / hist.sum() * np.arange(6)).sum())
			np.testing.assert_array_equal(idx.bitsum_bins.loc[(genome, chrom)].sum().to_numpy(), hist)
	assert len(idx.bitsum_chrs_avg) == len(idx.chrs)
	assert idx.bitsum_bins.index.is_monotonic_increasing


def test_index_bitsum_count_histogram(anchored_reader):
	"""bitsum_count turns occupancies 1..n into an array of n counts (index = occupancy - 1)."""
	got = anchored_reader.bitsum_count(np.array([1, 1, 3, 5, 5, 5]))
	assert got.dtype == np.uint32
	assert got.tolist() == [2, 0, 1, 0, 3]


def test_index_query_methods_delegate_to_the_genome(anchored_reader):
	"""query_bitmap, query_genes and query_anno return the named genome's answers."""
	idx = anchored_reader
	a = idx["gamma"]
	pd.testing.assert_frame_equal(idx.query_bitmap("gamma", "chr1", 10, 40, 3), a.query("chr1", 10, 40, 3))
	pd.testing.assert_frame_equal(idx.query_genes("gamma", "chr1"), a.query_genes("chr1"))
	pd.testing.assert_frame_equal(idx.query_anno("gamma", "chr1", 0, 5000), a.query_anno("chr1", 0, 5000))
	assert len(idx.query_genes("gamma", "chr1")) == 4
	with pytest.raises(KeyError):
		idx.query_bitmap("nope", "chr1")


def test_close_closes_every_genome(anchored_index):
	"""Index.close closes the bitmap handles of all genomes (opened privately here)."""
	idx = open_index(anchored_index)
	handles = [h for g in idx.genomes.values() for h in g.bitmaps.values()]
	assert len(handles) == 10
	idx.close()
	assert all(h._handle.closed for h in handles)


# --------------------------------------------------- bitmap-table helper functions

# a 6-position bitmap over the 5 genomes; row sums are 1, 2, 5, 1, 1, 2
BITMAP = pd.DataFrame(
	[
		[1, 0, 0, 0, 0],
		[1, 1, 0, 0, 0],
		[1, 1, 1, 1, 1],
		[1, 0, 0, 0, 0],
		[1, 0, 0, 0, 0],
		[1, 1, 0, 0, 0],
	],
	columns=list("abcde"),
)


def test_bitmap_to_pancount_is_the_row_sum(anchored_reader):
	"""The pan-count of a position is how many genomes carry its k-mer."""
	got = anchored_reader.bitmap_to_pancount(BITMAP)
	assert got.tolist() == [1, 2, 5, 1, 1, 2]
	assert got.index.tolist() == list(range(6))


def test_pancount_to_bins_counts_positions_per_bin(anchored_reader):
	"""pancount_to_bins gives rows pan-count 0..5 and columns bins; each cell counts positions."""
	got = anchored_reader.pancount_to_bins(pd.Series([1, 2, 5, 1, 1, 2]), 3)
	assert got.index.tolist() == [0, 1, 2, 3, 4, 5]
	assert got.columns.tolist() == [0, 1]
	assert got.to_numpy().tolist() == [[0, 0], [1, 2], [1, 1], [0, 0], [0, 0], [1, 0]]


def test_bitmap_to_paircount_bins_normalises_each_bin_by_its_maximum(anchored_reader):
	"""Per bin, the number of positions shared with each genome divided by the largest of those numbers."""
	got = anchored_reader.bitmap_to_paircount_bins(BITMAP, 3)
	assert got.index.tolist() == list("abcde")
	assert got.columns.tolist() == [0, 3]
	np.testing.assert_allclose(
		got.to_numpy(),
		[[1, 1], [2 / 3, 1 / 3], [1 / 3, 0], [1 / 3, 0], [1 / 3, 0]],
	)


def test_bitmap_to_bins_returns_both_tables(anchored_reader):
	"""bitmap_to_bins is the pair (pancount bins, paircount bins) of the two functions above."""
	pan, pair = anchored_reader.bitmap_to_bins(BITMAP, 3)
	assert pan.to_numpy().tolist() == [[0, 0], [1, 2], [1, 1], [0, 0], [0, 0], [1, 0]]
	np.testing.assert_allclose(pair.to_numpy(), anchored_reader.bitmap_to_paircount_bins(BITMAP, 3).to_numpy())
	assert pair.columns.tolist() == [0, 3]


def test_bitmap_to_paircount_bins_of_a_real_window(anchored_reader, synthetic_genomes):
	"""For a real lowres query the first-bin column counts the sampled positions that carry each genome's k-mer."""
	g = synthetic_genomes
	bm = anchored_reader.query_bitmap("alpha", "chr1", step=100)
	got = anchored_reader.bitmap_to_paircount_bins(bm, 1000)
	full = g.oracle_bits("alpha", "chr1")[::100]
	first = full[:10].sum(axis=0)
	np.testing.assert_allclose(got[0].to_numpy(), first / first.max())

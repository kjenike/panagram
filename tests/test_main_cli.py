"""panagram.__main__: argument parsing for every subcommand, dispatch, helper functions, error paths.

The CLI is driven in-process by setting sys.argv and calling main(), with `run` methods replaced by recorders where a
command would start a server or touch an anchored index.
"""

import dataclasses
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

import helpers_index as H
import synthetic
from conftest import K
from panagram import __main__ as cli


def run_main(monkeypatch, *args):
	monkeypatch.setattr(sys, "argv", ["panagram", *map(str, args)])
	return cli.main()


@pytest.fixture
def recorded(monkeypatch):
	"""Replace View/Bitdump/Annotate/Intros .run by recorders; returns the list of (class name, instance)."""
	calls = []
	for klass in (cli.View, cli.Bitdump, cli.Annotate, cli.Intros):
		monkeypatch.setattr(klass, "run", lambda self, _n=klass.__name__: calls.append((_n, self)))
	return calls


# --------------------------------------------------------------------- main()


def test_main_without_arguments_prints_help_and_runs_nothing(monkeypatch, capsys, recorded):
	"""`panagram` alone prints the usage with all five subcommands and the descriptions."""
	run_main(monkeypatch)
	out = capsys.readouterr().out
	assert "usage: panagram" in out
	for cmd in ("view", "index", "bitdump", "annotate", "intros"):
		assert cmd in out
	assert "Alignment-free pan-genome viewer" in out
	assert recorded == []


def test_main_help_hides_the_profiling_option(monkeypatch, capsys):
	"""--cprof is a hidden option: the help does not list it."""
	with pytest.raises(SystemExit) as e:
		run_main(monkeypatch, "--help")
	assert e.value.code == 0
	assert "cprof" not in capsys.readouterr().out


@pytest.mark.parametrize(
	"args",
	[["bogus"], ["view"], ["bitdump", "idx"], ["bitdump", "idx", "g"], ["annotate", "idx", "g"], ["view", "idx", "--port"]],
)
def test_main_rejects_bad_command_lines(monkeypatch, capsys, recorded, args):
	"""Unknown subcommands and missing or incomplete arguments exit with status 2 and run nothing."""
	with pytest.raises(SystemExit) as e:
		run_main(monkeypatch, *args)
	assert e.value.code == 2
	assert recorded == []


def test_main_rejects_non_integer_coordinates(monkeypatch, capsys, recorded):
	"""Coordinates are integers: `bitdump idx g chr1 a b` is a usage error."""
	with pytest.raises(SystemExit) as e:
		run_main(monkeypatch, "bitdump", "idx", "g", "chr1", "a", "b")
	assert e.value.code == 2
	assert recorded == []


def test_main_dataclass_run_returns_the_command_result():
	"""Main.run calls the chosen command's run and returns its result."""

	class Cmd:
		def run(self):
			return 42

	assert cli.Main(cmd=Cmd()).run() == 42
	assert cli.Main(cmd=Cmd()).cprof is None


# ------------------------------------------------------------------------ view


def test_view_positionals_and_defaults(monkeypatch, recorded):
	"""`view index_dir` leaves the optional positionals None and options at their defaults."""
	run_main(monkeypatch, "view", "idx")
	(name, v), = recorded
	assert name == "View"
	assert dataclasses.asdict(v) == {
		"index_dir": "idx",
		"genome": None,
		"chrom": None,
		"start": None,
		"end": None,
		"order": None,
		"ndebug": False,
		"port": "8050",
		"host": "127.0.0.1",
		"url_base": "/",
		"max_chr_bins": 350,
		"bookmarks": None,
	}


def test_view_all_options(monkeypatch, recorded):
	"""Every view positional and option is parsed to its field, integers as int."""
	run_main(
		monkeypatch, "view", "idx", "g1", "chr2", "5", "900", "--order", "a", "b", "c", "--ndebug",
		"--port", "9000", "--host", "0.0.0.0", "--url_base", "/pan/", "--max_chr_bins", "12", "--bookmarks", "m.bed",
	)
	(_, v), = recorded
	assert (v.index_dir, v.genome, v.chrom, v.start, v.end) == ("idx", "g1", "chr2", 5, 900)
	assert v.order == ["a", "b", "c"]
	assert v.ndebug is True
	assert (v.port, v.host, v.url_base, v.max_chr_bins, v.bookmarks) == ("9000", "0.0.0.0", "/pan/", 12, "m.bed")


def test_view_run_calls_the_viewer_with_itself(monkeypatch):
	"""View.run imports panagram.view.view and calls it with the View object."""
	pytest.importorskip("dash")
	import panagram.view

	seen = []
	monkeypatch.setattr(panagram.view, "view", lambda args: seen.append(args))
	v = cli.View(index_dir="idx", genome=None, chrom=None, start=None, end=None)
	v.run()
	assert seen == [v]


# --------------------------------------------------------------------- bitdump


def test_bitdump_parsing(monkeypatch, recorded):
	"""bitdump takes index_dir genome chromosome [start [end [step]]] and -v."""
	run_main(monkeypatch, "bitdump", "idx", "alpha", "chr1")
	run_main(monkeypatch, "bitdump", "idx", "alpha", "chr1", "10", "20", "-v")
	run_main(monkeypatch, "bitdump", "idx", "alpha", "chr1", "10", "2000", "100")
	(_, a), (_, b), (_, c) = recorded
	assert dataclasses.asdict(a) == dict(index_dir="idx", genome="alpha", chrom="chr1", start=None, end=None, step=1, verbose=False)
	assert (b.start, b.end, b.step, b.verbose) == (10, 20, 1, True)
	assert (c.start, c.end, c.step, c.verbose) == (10, 2000, 100, False)


def _expected_frame(genomes, genome, chrom, start, end, step=1):
	oracle = genomes.oracle_bits(genome, chrom)
	idx = pd.RangeIndex(start, end, step)
	return pd.DataFrame(oracle[start:end:step], index=idx, columns=pd.Index(genomes.names, name="name"))


def test_bitdump_prints_the_window_as_a_table(monkeypatch, anchored_index, synthetic_genomes, capsys):
	"""Bitdump.run prints the DataFrame of the requested window, equal to the brute-force bitmap printed the same way."""
	run_main(monkeypatch, "bitdump", anchored_index, "delta", "chr1", "100", "130")
	assert capsys.readouterr().out == str(_expected_frame(synthetic_genomes, "delta", "chr1", 100, 130)) + "\n"


def test_bitdump_whole_chromosome_and_step(monkeypatch, anchored_index, synthetic_genomes, capsys):
	"""Without start/end the whole chromosome is dumped; a step of 100 prints every hundredth row."""
	run_main(monkeypatch, "bitdump", anchored_index, "gamma", "chr2")
	n = len(synthetic_genomes.oracle_bits("gamma", "chr2"))
	assert capsys.readouterr().out == str(_expected_frame(synthetic_genomes, "gamma", "chr2", 0, n)) + "\n"
	run_main(monkeypatch, "bitdump", anchored_index, "gamma", "chr1", "0", "1000", "100")
	assert capsys.readouterr().out == str(_expected_frame(synthetic_genomes, "gamma", "chr1", 0, 1000, 100)) + "\n"


def test_bitdump_with_start_only_runs_to_the_chromosome_end(monkeypatch, anchored_index, synthetic_genomes, capsys):
	"""`bitdump idx g chr 600` (no end) dumps positions 600 to the chromosome's last k-mer."""
	run_main(monkeypatch, "bitdump", anchored_index, "gamma", "chr2", "600")
	n = len(synthetic_genomes.oracle_bits("gamma", "chr2"))
	assert capsys.readouterr().out == str(_expected_frame(synthetic_genomes, "gamma", "chr2", 600, n)) + "\n"


@pytest.mark.parametrize("genome,chrom", [("nope", "chr1"), ("alpha", "chrZ")])
def test_bitdump_unknown_genome_or_chromosome_raises_keyerror(monkeypatch, anchored_index, genome, chrom):
	"""Querying a genome or chromosome the index does not have raises KeyError."""
	with pytest.raises(KeyError):
		run_main(monkeypatch, "bitdump", anchored_index, genome, chrom, "0", "10")


def test_bitdump_on_a_missing_directory_fails(monkeypatch, tmp_path):
	"""A directory that is no index raises ValueError (not a directory) or OSError (not an anchored index)."""
	with pytest.raises((ValueError, OSError)):
		run_main(monkeypatch, "bitdump", tmp_path / "missing", "alpha", "chr1")


# -------------------------------------------------------------------- annotate


def test_annotate_parsing(monkeypatch, recorded):
	"""annotate takes index_dir genome_name gff_file and --nogene."""
	run_main(monkeypatch, "annotate", "idx", "alpha", "new.gff")
	run_main(monkeypatch, "annotate", "idx", "alpha", "new.gff", "--nogene")
	(_, a), (_, b) = recorded
	assert dataclasses.asdict(a) == dict(index_dir="idx", genome="alpha", gff_file="new.gff", nogene=False)
	assert b.nogene is True


def test_annotate_run_calls_run_annotate_on_the_named_genome(monkeypatch, anchored_index):
	"""Annotate.run opens the index and calls run_annotate(gff, nogene=...) on the chosen genome only."""
	from panagram.index import Genome

	calls = []
	monkeypatch.setattr(Genome, "run_annotate", lambda self, gff_file=None, logfile=None, nogene=False: calls.append((self.name, gff_file, nogene)))
	run_main(monkeypatch, "annotate", anchored_index, "beta", "x.gff", "--nogene")
	run_main(monkeypatch, "annotate", anchored_index, "gamma", "y.gff")
	assert calls == [("beta", "x.gff", True), ("gamma", "y.gff", False)]


def test_annotate_unknown_genome_raises_keyerror(monkeypatch, anchored_index):
	"""An unknown genome name raises KeyError."""
	with pytest.raises(KeyError):
		run_main(monkeypatch, "annotate", anchored_index, "nope", "x.gff")


# ----------------------------------------------------------------------- intros


def test_intros_sweep_flag_is_parsed(monkeypatch, recorded):
	"""intros accepts --sweep; without it the flag is False."""
	run_main(monkeypatch, "intros", "--sweep")
	run_main(monkeypatch, "intros")
	(_, a), (_, b) = recorded
	assert a.sweep is True and b.sweep is False


def test_intros_helper_dispatch_through_main_passes_the_remaining_arguments(monkeypatch):
	"""`panagram intros heatmap -i a.txt --flag` hands everything after the helper name to the script."""
	seen = []
	monkeypatch.setattr(cli.subprocess, "call", lambda cmd: seen.append(cmd) or 0)
	with pytest.raises(SystemExit) as e:
		run_main(monkeypatch, "intros", "heatmap", "-i", "a.txt", "--flag")
	assert e.value.code == 0
	assert seen[0][2:] == ["-i", "a.txt", "--flag"]
	assert seen[0][1].endswith("create_heatmap.py")


@pytest.mark.parametrize(
	"name,script",
	[("heatmap", "create_heatmap.py"), ("bed2txt", "bed_to_txt.py"), ("simulate", "simulate_introgressions.py")],
)
def test_intros_helper_scripts_run_in_a_subprocess_and_exit_with_its_status(monkeypatch, name, script):
	"""`intros heatmap|bed2txt|simulate ARGS` runs the helper script with the same interpreter and exits with its status."""
	seen = []
	monkeypatch.setattr(cli.subprocess, "call", lambda cmd: seen.append(cmd) or 7)
	with pytest.raises(SystemExit) as e:
		cli.Intros(extra_args=[name, "--x", "1"]).run()
	assert e.value.code == 7
	(cmd,) = seen
	assert cmd[0] == sys.executable
	assert Path(cmd[1]) == Path(cli.__file__).resolve().parent / "introgressions" / script
	assert Path(cmd[1]).is_file()
	assert cmd[2:] == ["--x", "1"]


def test_intros_reads_the_target_from_sys_argv_when_no_extra_args(monkeypatch):
	"""Without extra_args the helper name and its arguments are taken from sys.argv after 'intros'."""
	seen = []
	monkeypatch.setattr(cli.subprocess, "call", lambda cmd: seen.append(cmd) or 0)
	monkeypatch.setattr(sys, "argv", ["panagram", "intros", "heatmap", "a", "b"])
	with pytest.raises(SystemExit) as e:
		cli.Intros().run()
	assert e.value.code == 0
	assert seen[0][2:] == ["a", "b"]


def test_intros_without_target_prints_usage(monkeypatch, capsys):
	"""`panagram intros` with nothing else prints the usage lines and returns."""
	monkeypatch.setattr(sys, "argv", ["panagram", "intros"])
	assert cli.Intros().run() is None
	out = capsys.readouterr().out
	assert "usage: panagram intros <config.yaml> [--sweep]" in out
	assert "heatmap" in out and "simulate" in out


@pytest.mark.parametrize("sweep", [False, True])
def test_intros_runs_the_introgression_runner_on_a_config(monkeypatch, sweep):
	"""Any other first argument is a config path handed to introgression_runner.main with the sweep flag."""
	import types

	seen = []
	fake = types.ModuleType("panagram.introgressions.introgression_runner")
	fake.main = lambda path, sweep=False: seen.append((path, sweep))
	monkeypatch.setitem(sys.modules, "panagram.introgressions.introgression_runner", fake)
	cli.Intros(sweep=sweep, extra_args=["conf.yaml"]).run()
	assert seen == [(Path("conf.yaml"), sweep)]


# ------------------------------------------------------------------- index


def _samples(tmp_path, genomes=("alpha", "beta", "gamma")):
	"""samples.tsv with absolute FASTA paths (the index may live in another directory)."""
	for n in genomes:
		synthetic.write_fasta(tmp_path / f"{n}.fa", {"c1": "ACGT" * 20})
	tsv = tmp_path / "samples.tsv"
	tsv.write_text("name\tfasta\n" + "".join(f"{n}\t{tmp_path / (n + '.fa')}\n" for n in genomes))
	return tsv


def _index_config(monkeypatch, tmp_path, *opts, prefix=None):
	tsv = _samples(tmp_path)
	cwd = os.getcwd()
	try:
		run_main(monkeypatch, "index", tsv, "--prepare", *opts)
	finally:
		os.chdir(cwd)
	return yaml.safe_load(((Path(prefix) if prefix else tmp_path) / "config.yaml").read_text())


def test_index_defaults_in_config(monkeypatch, tmp_path, capsys):
	"""`index samples.tsv --prepare` writes config.yaml with the documented defaults and prints the next step."""
	cfg = _index_config(monkeypatch, tmp_path)
	out = capsys.readouterr().out
	assert "Prepared. Run 'snakemake --cores 1 all' to build index" in out
	assert (tmp_path / "Snakefile").is_file() and (tmp_path / "samples.tsv").is_file()
	assert cfg["k"] == 21 and cfg["cores"] == 1 and cfg["lowres_step"] == 100
	assert cfg["max_bin_kbp"] == 200 and cfg["min_bin_count"] == 100 and cfg["max_view_chrs"] == 50
	assert cfg["gff_gene_types"] == ["gene"] and cfg["gff_anno_types"] is None and cfg["gff_name"] == "Name"
	assert cfg["anchor_genomes"] == ["alpha", "beta", "gamma"] and cfg["prepare"] is True


@pytest.mark.parametrize(
	"opts,path,expected",
	[
		(["-k", "11"], ("k",), 11),
		(["-c", "3"], ("cores",), 3),
		(["--lowres_step", "50"], ("lowres_step",), 50),
		(["--max_bin_kbp", "7"], ("max_bin_kbp",), 7),
		(["--min_bin_count", "9"], ("min_bin_count",), 9),
		(["--max_view_chrs", "4"], ("max_view_chrs",), 4),
		(["--gff_gene_types", "gene", "mRNA"], ("gff_gene_types",), ["gene", "mRNA"]),
		(["--gff_anno_types", "exon", "CDS"], ("gff_anno_types",), ["exon", "CDS"]),
		(["--gff_name", "ID"], ("gff_name",), "ID"),
		(["--anchor_genomes", "alpha", "gamma"], ("anchor_genomes",), ["alpha", "gamma"]),
		(["--kmc.use_existing"], ("kmc", "use_existing"), True),
		(["--genome_umap.neighbors", "7"], ("genome_umap", "neighbors"), 7),
		(["--genome_umap.dist", "0.5"], ("genome_umap", "dist"), 0.5),
		(["--genome_umap.eps", "2"], ("genome_umap", "eps"), 2),
		(["--genome_umap.samples", "3"], ("genome_umap", "samples"), 3),
		(["--genome_umap.bin_size", "5000"], ("genome_umap", "bin_size"), 5000),
		(["--chrom_umap.neighbors", "6"], ("chrom_umap", "neighbors"), 6),
		(["--chrom_umap.dist", "0.25"], ("chrom_umap", "dist"), 0.25),
		(["--chrom_umap.eps", "3"], ("chrom_umap", "eps"), 3),
		(["--chrom_umap.samples", "2"], ("chrom_umap", "samples"), 2),
		(["--chrom_umap.bin_size", "50"], ("chrom_umap", "bin_size"), 50),
	],
)
def test_index_option_reaches_config_yaml(monkeypatch, tmp_path, opts, path, expected):
	"""Each `panagram index` option is written under its name in config.yaml (nested for kmc.*, genome_umap.*, chrom_umap.*)."""
	cfg = _index_config(monkeypatch, tmp_path, *opts)
	for p in path:
		cfg = cfg[p]
	assert cfg == expected


def test_index_prefix_option_sends_config_to_another_directory(monkeypatch, tmp_path):
	"""-o/--prefix writes config.yaml and samples.tsv there, not next to the input."""
	out = tmp_path / "out"
	out.mkdir()
	cfg = _index_config(monkeypatch, tmp_path, "-o", out, prefix=out)
	assert cfg["k"] == 21
	assert (out / "samples.tsv").is_file() and (out / "Snakefile").is_file()
	assert not (tmp_path / "config.yaml").exists()
	cfg = _index_config(monkeypatch, tmp_path, "--prefix", out, "-k", "9", prefix=out)
	assert cfg["k"] == 9


def test_index_options_combine_and_anchor_genomes_mark_anchors(monkeypatch, tmp_path):
	"""Several options at once; --anchor_genomes limits the anchor column of samples.tsv."""
	cfg = _index_config(monkeypatch, tmp_path, "-k", "13", "--anchor_genomes", "beta", "--chrom_umap.bin_size", "500")
	assert cfg["k"] == 13 and cfg["chrom_umap"]["bin_size"] == 500 and cfg["genome_umap"]["bin_size"] == 100000
	df = pd.read_table(tmp_path / "samples.tsv")
	assert df["anchor"].tolist() == [False, True, False]


def test_index_command_with_invalid_genome_name_fails_before_writing(monkeypatch, tmp_path):
	"""A genome name that does not match the name pattern aborts `panagram index` with ValueError and writes no config."""
	tsv = tmp_path / "samples.tsv"
	tsv.write_text("name\tfasta\nbad name\tx.fa\n")
	with pytest.raises(ValueError, match="Invalid genome names"):
		run_main(monkeypatch, "index", tsv, "--prepare")
	assert not (tmp_path / "config.yaml").exists()


def test_index_without_input_is_a_usage_error(monkeypatch, capsys):
	"""`panagram index` needs the config/samples file argument."""
	with pytest.raises(SystemExit) as e:
		run_main(monkeypatch, "index")
	assert e.value.code == 2


# --------------------------------------------------------------- helper functions


def test_parse_coords_splits_genome_chrom_and_range():
	"""'genome:chrom:start-end' gives (genome, chrom, int start, int end)."""
	assert cli.parse_coords("alpha:chr1:100-2000") == ("alpha", "chr1", 100, 2000)
	assert cli.parse_coords("g:c:0-1") == ("g", "c", 0, 1)


@pytest.mark.parametrize("bad", ["alpha:chr1", "alpha:chr1:100", "alpha:chr1:a-b", "a:b:c:1-2"])
def test_parse_coords_rejects_malformed_text(bad):
	"""Anything but genome:chrom:start-end raises ValueError."""
	with pytest.raises(ValueError):
		cli.parse_coords(bad)


def test_comma_split():
	"""comma_split splits on commas and keeps empty fields."""
	assert cli.comma_split("a,b,c") == ["a", "b", "c"]
	assert cli.comma_split("a") == ["a"]
	assert cli.comma_split("a,,b") == ["a", "", "b"]

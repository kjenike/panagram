"""Tests of the postprocess_introgressions command line (--bed --idx --out --act --ref --map --threads
--paf --bin --min --gap) on the anchored synthetic index.

Chromosome lengths in the index are k-mer counts: alpha chr1 2636, alpha chr2 1686, beta chr1 2886
(k = 15). The tool is called in-process with sys.argv set, as the console entry would.
"""

import sys
from pathlib import Path

import pytest

from helpers_introgression_post import fake_tools, index_copy, write_bed  # noqa: F401
from panagram.introgressions import postprocess_introgressions as pp


@pytest.fixture
def run_pp(anchored_index, monkeypatch, tmp_path):
	"""run_pp(bed_dir_or_file, *args, idx=anchored_index) -> output directory; raises what the tool raises."""

	def run(bed, *args, idx=None, out=None):
		out = out or tmp_path / "out"
		argv = ["postprocess", "--bed", str(bed), "--idx", str(idx or anchored_index), "--out", str(out), *args]
		monkeypatch.setattr(sys, "argv", argv)
		pp.postprocess_introgressions()
		return out

	return run


def _beds(tmp_path, **files):
	d = tmp_path / "beds"
	d.mkdir(exist_ok=True)
	for name, rows in files.items():
		write_bed(d / f"{name}.bed", rows)
	return d


def _lines(path):
	return Path(path).read_text().splitlines()


TWO_REGIONS = [("chr1", 200, 400, "x"), ("chr1", 600, 800, "x")]  # bins 2,3 and 6,7 at --bin 100


def test_fgap_merges_regions_separated_by_at_most_gap_bins(run_pp, tmp_path):
	beds = _beds(tmp_path, alpha_chr1_WT=TWO_REGIONS)
	out = run_pp(beds, "--act", "fgap", "--bin", "100", "--gap", "2")
	assert _lines(out / "alpha_chr1_WT.bed") == ["chr1\t200\t799\tWT_intro"]


def test_fgap_keeps_regions_separated_by_more_than_gap_bins(run_pp, tmp_path):
	beds = _beds(tmp_path, alpha_chr1_WT=TWO_REGIONS)
	out = run_pp(beds, "--act", "fgap", "--bin", "100", "--gap", "1")
	assert _lines(out / "alpha_chr1_WT.bed") == ["chr1\t200\t399\tWT_intro", "chr1\t600\t799\tWT_intro"]


def test_default_gap_is_one_bin(run_pp, tmp_path):
	# bins 2,3 and 5: one empty bin between them
	beds = _beds(tmp_path, alpha_chr1_WT=[("chr1", 200, 400, "x"), ("chr1", 500, 600, "x")])
	out = run_pp(beds, "--act", "fgap", "--bin", "100")
	assert _lines(out / "alpha_chr1_WT.bed") == ["chr1\t200\t599\tWT_intro"]


def test_rmbn_removes_runs_shorter_than_min_bins(run_pp, tmp_path):
	beds = _beds(tmp_path, alpha_chr1_WT=[("chr1", 200, 400, "x"), ("chr1", 600, 1000, "x")])
	out = run_pp(beds, "--act", "rmbn", "--bin", "100", "--min", "3")
	assert _lines(out / "alpha_chr1_WT.bed") == ["chr1\t600\t999\tWT_intro"]


def test_rmbn_writes_an_empty_file_when_everything_is_removed(run_pp, tmp_path):
	beds = _beds(tmp_path, alpha_chr1_WT=TWO_REGIONS)
	out = run_pp(beds, "--act", "rmbn", "--bin", "100", "--min", "3")
	assert (out / "alpha_chr1_WT.bed").read_text() == ""


def test_default_min_is_four_bins(run_pp, tmp_path):
	beds = _beds(
		tmp_path,
		alpha_chr1_WT=[("chr1", 200, 500, "x")],  # 3 bins: removed
		beta_chr1_WT=[("chr1", 200, 600, "x")],  # 4 bins: kept
	)
	out = run_pp(beds, "--act", "rmbn", "--bin", "100")
	assert (out / "alpha_chr1_WT.bed").read_text() == ""
	assert _lines(out / "beta_chr1_WT.bed") == ["chr1\t200\t599\tWT_intro"]


def test_actions_are_applied_in_the_order_given(run_pp, tmp_path):
	beds = _beds(tmp_path, alpha_chr1_WT=TWO_REGIONS)
	# fill the gap (one 6-bin region), then the region is long enough to survive --min 5
	out = run_pp(beds, "--act", "fgap", "rmbn", "--bin", "100", "--gap", "2", "--min", "5")
	assert _lines(out / "alpha_chr1_WT.bed") == ["chr1\t200\t799\tWT_intro"]
	# removing small regions first deletes both 2-bin regions, so the later gap fill has nothing to join
	out2 = run_pp(beds, "--act", "rmbn", "fgap", "--bin", "100", "--gap", "2", "--min", "5", out=tmp_path / "out2")
	assert (out2 / "alpha_chr1_WT.bed").read_text() == ""


def test_a_folder_of_bed_files_is_processed_and_only_bed_files_are_picked_up(run_pp, tmp_path):
	beds = _beds(tmp_path, alpha_chr1_WT=TWO_REGIONS, beta_chr1_WT=[("chr1", 1000, 1200, "x")])
	(beds / "notes.txt").write_text("not a bed file")
	out = run_pp(beds, "--act", "fgap", "--bin", "100", "--gap", "2")
	assert sorted(p.name for p in out.iterdir()) == ["alpha_chr1_WT.bed", "beta_chr1_WT.bed"]
	assert _lines(out / "beta_chr1_WT.bed") == ["chr1\t1000\t1199\tWT_intro"]


def test_a_single_bed_file_can_be_given_instead_of_a_folder(run_pp, tmp_path):
	beds = _beds(tmp_path, alpha_chr1_WT=TWO_REGIONS, beta_chr1_WT=TWO_REGIONS)
	out = run_pp(beds / "alpha_chr1_WT.bed", "--act", "fgap", "--bin", "100", "--gap", "2")
	assert [p.name for p in out.iterdir()] == ["alpha_chr1_WT.bed"]


def test_an_empty_bed_file_gives_an_empty_output_file(run_pp, tmp_path):
	beds = _beds(tmp_path)
	(beds / "alpha_chr1_WT.bed").write_text("")
	out = run_pp(beds, "--act", "fgap", "rmbn", "--bin", "100")
	assert (out / "alpha_chr1_WT.bed").read_text() == ""


def test_output_folder_is_created_with_parents(run_pp, tmp_path):
	beds = _beds(tmp_path, alpha_chr1_WT=TWO_REGIONS)
	out = run_pp(beds, "--act", "fgap", "--bin", "100", out=tmp_path / "a" / "b" / "c")
	assert (out / "alpha_chr1_WT.bed").is_file()


def test_ref_files_take_their_chromosome_length_from_the_reference(run_pp, tmp_path):
	# beta has no chr2; alpha chr2 has 1686 positions = 17 bins of 100, so bins 15 and 16 fit
	beds = _beds(tmp_path, beta_chr2_REF=[("chr2", 1500, 1700, "x")])
	out = run_pp(beds, "--act", "fgap", "--bin", "100", "--ref", "alpha")
	assert _lines(out / "beta_chr2_REF.bed") == ["chr2\t1500\t1699\tREF_intro"]


def test_ref_files_without_ref_are_rejected(run_pp, tmp_path):
	beds = _beds(tmp_path, alpha_chr1_REF=TWO_REGIONS)
	with pytest.raises(ValueError, match="--ref"):
		run_pp(beds, "--act", "fgap", "--bin", "100")


@pytest.mark.parametrize(
	"extra,match",
	[
		(["--act", "bogus"], "Unrecognized action"),
		(["--act", "fgap", "nope"], "Unrecognized action"),
	],
)
def test_unknown_actions_are_rejected(run_pp, tmp_path, extra, match):
	beds = _beds(tmp_path, alpha_chr1_WT=TWO_REGIONS)
	with pytest.raises(ValueError, match=match):
		run_pp(beds, *extra)


def test_missing_index_and_missing_bed_are_rejected(run_pp, tmp_path):
	beds = _beds(tmp_path, alpha_chr1_WT=TWO_REGIONS)
	with pytest.raises(ValueError, match="Index directory not found"):
		run_pp(beds, "--act", "fgap", idx=tmp_path / "no_such_index")
	with pytest.raises(ValueError, match="Bed file/folder not found"):
		run_pp(tmp_path / "no_such_beds", "--act", "fgap")


def test_act_and_bed_and_idx_and_out_are_required(anchored_index, monkeypatch, capsys):
	for argv in (
		["pp", "--idx", str(anchored_index), "--out", "o", "--act", "fgap"],
		["pp", "--bed", "b", "--out", "o", "--act", "fgap"],
		["pp", "--bed", "b", "--idx", str(anchored_index), "--act", "fgap"],
		["pp", "--bed", "b", "--idx", str(anchored_index), "--out", "o"],
	):
		monkeypatch.setattr(sys, "argv", argv)
		with pytest.raises(SystemExit) as e:
			pp.postprocess_introgressions()
		assert e.value.code == 2


def test_lift_requires_ref(run_pp, tmp_path):
	beds = _beds(tmp_path, beta_chr1_WT=TWO_REGIONS)
	with pytest.raises(ValueError, match="--ref must be specified for liftover"):
		run_pp(beds, "--act", "lift", "--paf", str(tmp_path))


def test_lift_with_a_paf_folder_lifts_coordinates_through_paftools(run_pp, tmp_path, fake_tools):
	beds = _beds(tmp_path, beta_chr1_WT=[("chr1", 1, 2, "x")])
	paf = tmp_path / "paf"
	paf.mkdir()
	fake_tools.set_lift("chr1\t200\t400\tq\nchr1\t600\t800\tq\n")
	out = run_pp(beds, "--act", "lift", "--ref", "alpha", "--paf", str(paf), "--bin", "100")
	assert _lines(out / "beta_chr1_WT.bed") == ["chr1\t200\t400\tWT_intro", "chr1\t600\t800\tWT_intro"]
	assert fake_tools.args_file.read_text().split()[:2] == ["liftover", str(paf / "beta_alpha.paf")]
	assert not (Path(paf).parent / "alignments").exists()


def test_lift_then_fgap_processes_the_lifted_coordinates(run_pp, tmp_path, fake_tools):
	beds = _beds(tmp_path, beta_chr1_WT=[("chr1", 1, 2, "x")])
	paf = tmp_path / "paf"
	paf.mkdir()
	fake_tools.set_lift("chr1\t200\t400\tq\nchr1\t600\t800\tq\n")
	out = run_pp(beds, "--act", "lift", "fgap", "--ref", "alpha", "--paf", str(paf), "--bin", "100", "--gap", "2")
	assert _lines(out / "beta_chr1_WT.bed") == ["chr1\t200\t799\tWT_intro"]


def test_lift_without_paf_aligns_with_the_default_minimap2_flags(run_pp, tmp_path, fake_tools, index_copy):
	beds = _beds(tmp_path, beta_chr1_WT=[("chr1", 1, 2, "x")])
	fake_tools.set_lift("chr1\t200\t400\tq\n")
	run_pp(beds, "--act", "lift", "--ref", "alpha", idx=index_copy)
	paf = index_copy / "alignments" / "beta_alpha.paf"
	assert paf.read_text().split() == [
		"-x", "asm20", "-c", "-t", "1",
		str(index_copy / "FASTAS/alpha.fasta"),
		str(index_copy / "FASTAS/beta.fasta"),
	]  # fmt: skip


def test_map_option_replaces_the_minimap2_flags_and_threads_runs_accessions_in_parallel(
	run_pp, tmp_path, fake_tools, index_copy
):
	beds = _beds(
		tmp_path,
		beta_chr1_WT=[("chr1", 1, 2, "x")],
		gamma_chr1_WT=[("chr1", 1, 2, "x")],
	)
	fake_tools.set_lift("chr1\t200\t400\tq\n")
	run_pp(beds, "--act", "lift", "--ref", "alpha", "--map", "-x map-ont", "--threads", "2", idx=index_copy)
	for acc in ("beta", "gamma"):
		got = (index_copy / "alignments" / f"{acc}_alpha.paf").read_text().split()
		assert got == [
			"-x", "map-ont",
			str(index_copy / "FASTAS/alpha.fasta"),
			str(index_copy / f"FASTAS/{acc}.fasta"),
		]  # fmt: skip

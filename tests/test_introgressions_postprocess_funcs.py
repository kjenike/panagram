"""Unit tests for the functions of panagram/introgressions/postprocess_introgressions.py.

Expected values are worked by hand on tiny inputs. minimap2 and paftools.js are
replaced by shell stand-ins (helpers_introgression_post.fake_tools).
"""

import io
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from helpers_introgression_post import fake_tools, write_bed  # noqa: F401
from panagram.introgressions import postprocess_introgressions as pp


# --- FASTA reading -----------------------------------------------------------------


def test_read_fasta_generator_joins_wrapped_lines_and_keeps_header_text():
	lines = io.StringIO(">a some description\nACGT\nTTGG\n>b\nCC\n")
	assert list(pp.read_fasta_generator(lines)) == [(">a some description", "ACGTTTGG"), (">b", "CC")]


def test_read_fasta_generator_empty_record_and_empty_input():
	assert list(pp.read_fasta_generator(io.StringIO(">e\n>f\nA\n"))) == [(">e", ""), (">f", "A")]
	assert list(pp.read_fasta_generator(io.StringIO(""))) == []


def test_read_fasta_returns_name_and_sequence_columns(tmp_path):
	fa = tmp_path / "x.fa"
	fa.write_text(">chr1\nACGT\nAC\n>chr2\nggnn\n")
	df = pp.read_fasta(fa)
	assert list(df.columns) == ["Name", "Sequence"]
	assert df["Name"].tolist() == ["chr1", "chr2"]
	assert df["Sequence"].tolist() == ["ACGTAC", "ggnn"]


def test_extract_bed_region_is_zero_based_half_open_and_picks_the_chromosome():
	fasta_df = pd.DataFrame({"Name": ["c1", "c2"], "Sequence": ["ACGTACGT", "TTTTGGGG"]})
	row = {"Chromosome": "c1", "Start": 2, "End": 5}
	assert pp.extract_bed_region(row, fasta_df) == "GTA"
	assert pp.extract_bed_region({"Chromosome": "c2", "Start": 3, "End": 5}, fasta_df) == "TG"
	assert pp.extract_bed_region({"Chromosome": "c2", "Start": 0, "End": 8}, fasta_df) == "TTTTGGGG"
	# floats (from a shifted column) are accepted
	assert pp.extract_bed_region({"Chromosome": "c1", "Start": 2.0, "End": 5.0}, fasta_df) == "GTA"


# --- BED reading -----------------------------------------------------------------------


def test_bed_file_is_empty(tmp_path):
	empty = tmp_path / "empty.bed"
	empty.write_text("")
	full = write_bed(tmp_path / "full.bed", [("chr1", 0, 10, "x")])
	assert pp.bed_file_is_empty(empty) is True
	assert pp.bed_file_is_empty(str(full)) is False


def test_read_bed_file_keeps_first_four_columns_and_names_them(tmp_path):
	bed = write_bed(tmp_path / "a.bed", [("chr1", 5, 20, "n1", 0, "+"), ("chr2", 30, 40, "n2", 1, "-")])
	df = pp.read_bed_file(bed)
	assert list(df.columns) == ["Chromosome", "Start", "End", "Notes", "Sequence"]
	assert df["Chromosome"].tolist() == ["chr1", "chr2"]
	assert df["Start"].tolist() == [5, 30]
	assert df["End"].tolist() == [20, 40]
	assert df["Notes"].tolist() == ["n1", "n2"]
	assert df["Sequence"].isna().all()


def test_read_bed_file_returns_none_for_an_empty_file(tmp_path):
	empty = tmp_path / "empty.bed"
	empty.write_text("")
	assert pp.read_bed_file(empty) is None


# --- centromere merging -------------------------------------------------------------------


def _bed(rows):
	return pd.DataFrame(rows, columns=["Chromosome", "Start", "End", "Notes"]).assign(Sequence=None)


def _fasta(seq):
	return pd.DataFrame({"Name": ["c"], "Sequence": [seq]})


def test_merge_centromere_regions_joins_two_regions_separated_by_a_run_of_n():
	# bin 100: regions [0,100) and [300,400); the 200 bp between them hold 100 N
	seq = "A" * 100 + "C" * 50 + "N" * 100 + "C" * 50 + "G" * 100 + "T" * 100
	out = pp.merge_centromere_regions(_bed([("c", 0, 100, "a"), ("c", 300, 400, "b")]), _fasta(seq), 100)
	assert out["Start"].tolist() == [0]
	assert out["End"].tolist() == [400]
	assert out["Chromosome"].tolist() == ["c"]
	assert out["Sequence"].tolist() == [seq[0:400]]
	assert out["Start"].dtype.kind == "i" and out["End"].dtype.kind == "i"


def test_merge_centromere_regions_requires_50_n_in_a_row():
	seq = "A" * 100 + ("N" * 49 + "A") * 4 + "G" * 100
	out = pp.merge_centromere_regions(_bed([("c", 0, 100, "a"), ("c", 300, 400, "b")]), _fasta(seq), 100)
	assert out["Start"].tolist() == [0, 300]
	assert out["End"].tolist() == [100, 400]


def test_merge_centromere_regions_ignores_regions_not_exactly_two_bins_apart():
	seq = "N" * 800
	bed = _bed([("c", 0, 100, "a"), ("c", 200, 300, "b"), ("c", 600, 700, "c")])
	out = pp.merge_centromere_regions(bed, _fasta(seq), 100)
	assert out["Start"].tolist() == [0, 200, 600]
	assert out["End"].tolist() == [100, 300, 700]


def test_merge_centromere_regions_merges_only_the_qualifying_pair():
	# gap [100,300) is N, gap [400,600) is not, so only the first pair merges
	seq = "A" * 100 + "N" * 200 + "A" * 100 + "C" * 200 + "G" * 100
	bed = _bed([("c", 0, 100, "a"), ("c", 300, 400, "b"), ("c", 600, 700, "c")])
	out = pp.merge_centromere_regions(bed, _fasta(seq), 100)
	assert sorted(zip(out["Start"], out["End"])) == [(0, 400), (600, 700)]


# --- calling minimap2 / paftools -------------------------------------------------------------


def test_align_assemblies_query_to_reference_puts_reference_first(tmp_path, fake_tools):
	out = tmp_path / "q_r.paf"
	pp.align_assemblies("REF.fa", "QRY.fa", ["-x", "asm20", "-c"], out)
	assert out.read_text().split() == ["-x", "asm20", "-c", "REF.fa", "QRY.fa"]
	assert (tmp_path / "q_r.log").read_text().strip() == "fake minimap2 log"


def test_align_assemblies_reference_to_query_swaps_the_files(tmp_path, fake_tools):
	out = tmp_path / "r_q.paf"
	pp.align_assemblies("REF.fa", "QRY.fa", ["-t", "1"], out, direction="reference_to_query")
	assert out.read_text().split() == ["-t", "1", "QRY.fa", "REF.fa"]


def test_align_assemblies_rejects_unknown_direction(tmp_path, fake_tools):
	with pytest.raises(ValueError, match="Direction"):
		pp.align_assemblies("R", "Q", [], tmp_path / "o.paf", direction="sideways")
	assert not (tmp_path / "o.paf").exists()


def test_align_assemblies_exits_with_status_1_when_minimap2_fails(tmp_path, fake_tools, monkeypatch):
	monkeypatch.setenv("FAKE_MINIMAP_FAIL", "1")
	with pytest.raises(SystemExit) as e:
		pp.align_assemblies("R", "Q", [], tmp_path / "o.paf")
	assert e.value.code == 1
	assert not (tmp_path / "o.paf").exists()


def test_liftover_coordinates_passes_files_and_parses_four_columns(tmp_path, fake_tools):
	fake_tools.set_lift("chr1\t10\t20\tnote1\textra\nchr2\t5\t9\tnote2\textra\n")
	df = pp.liftover_coordinates("in.bed", "map.paf")
	assert fake_tools.args_file.read_text().split() == ["liftover", "map.paf", "in.bed", "output.liftover"]
	assert list(df.columns) == ["Chromosome", "Start", "End", "Notes"]
	assert df["Chromosome"].tolist() == ["chr1", "chr2"]
	assert df["Start"].tolist() == [10, 5]
	assert df["End"].tolist() == [20, 9]
	assert df["Notes"].tolist() == ["note1", "note2"]


def test_liftover_coordinates_exits_with_status_1_when_paftools_fails(fake_tools, monkeypatch):
	monkeypatch.setenv("FAKE_PAFTOOLS_FAIL", "1")
	with pytest.raises(SystemExit) as e:
		pp.liftover_coordinates("in.bed", "map.paf")
	assert e.value.code == 1


def _stub_index(names):
	return SimpleNamespace(genomes={n: SimpleNamespace(fasta=f"FASTAS/{n}.fasta") for n in names})


def test_run_alignments_with_paf_folder_runs_nothing_and_returns_it(tmp_path, fake_tools):
	paf_dir = tmp_path / "given"
	got = pp.run_alignments(["a"], str(paf_dir), tmp_path, _stub_index(["a"]), "ref", tmp_path / "ref.fa", "-x y", 1)
	assert got == paf_dir
	assert not (tmp_path / "alignments").exists()


def test_run_alignments_runs_one_minimap2_per_accession(tmp_path, fake_tools):
	index_dir = tmp_path / "idx"
	index_dir.mkdir()
	ref_file = index_dir / "FASTAS" / "ref.fasta"
	got = pp.run_alignments(
		["a", "b"], None, index_dir, _stub_index(["a", "b"]), "ref", ref_file, "-x asm20 -c -t 1", 2
	)
	assert got == index_dir / "alignments"
	assert sorted(p.name for p in got.glob("*.paf")) == ["a_ref.paf", "b_ref.paf"]
	for acc in ("a", "b"):
		args = (got / f"{acc}_ref.paf").read_text().split()
		assert args == ["-x", "asm20", "-c", "-t", "1", str(ref_file), str(index_dir / f"FASTAS/{acc}.fasta")]


# --- liftover of BED files ---------------------------------------------------------------------


def _intro_lines(path):
	return Path(path).read_text().splitlines()


def test_run_liftovers_sorts_ref_files_lifts_others_and_writes_empty_placeholders(tmp_path, fake_tools):
	beds = tmp_path / "beds"
	beds.mkdir()
	out = tmp_path / "out"
	out.mkdir()
	paf_dir = tmp_path / "paf"
	paf_dir.mkdir()
	b_ref = write_bed(beds / "alpha_chr1_REF.bed", [("chr1", 500, 900, "n1"), ("chr1", 100, 300, "n2")])
	b_wt = write_bed(beds / "beta_chr1_WT.bed", [("chr1", 1, 2, "n")])
	b_empty = beds / "gamma_chr2_WT.bed"
	b_empty.write_text("")
	fake_tools.set_lift("chr1\t10\t20\tx\nchr1\t5\t8\ty\n")

	got = pp.run_liftovers([b_ref, b_wt, b_empty], ["alpha", "beta", "gamma"], paf_dir, "alpha", out)

	assert got == [out / "alpha_chr1_REF.bed", out / "beta_chr1_WT.bed", out / "gamma_chr2_WT.bed"]
	assert _intro_lines(got[0]) == ["chr1\t100\t300\tREF_intro", "chr1\t500\t900\tREF_intro"]
	assert _intro_lines(got[1]) == ["chr1\t5\t8\tWT_intro", "chr1\t10\t20\tWT_intro"]
	assert got[2].read_text() == ""
	# only the non-REF, non-empty file went through paftools, against the paf named accession_reference
	assert fake_tools.args_file.read_text().split() == [
		"liftover",
		str(paf_dir / "beta_alpha.paf"),
		str(b_wt),
		"output.liftover",
	]


def test_run_liftovers_treats_refa_files_as_lifted_ref_files(tmp_path, fake_tools):
	beds = tmp_path / "beds"
	beds.mkdir()
	out = tmp_path / "out"
	out.mkdir()
	paf_dir = tmp_path / "paf"
	paf_dir.mkdir()
	b = write_bed(beds / "beta_chr1_REFA.bed", [("chr1", 1, 2, "n")])
	fake_tools.set_lift("chr1\t40\t60\tx\n")
	got = pp.run_liftovers([b], ["beta"], paf_dir, "alpha", out)
	assert got == [out / "beta_chr1_REF.bed"]
	assert _intro_lines(got[0]) == ["chr1\t40\t60\tREF_intro"]


def test_run_liftovers_requires_the_paf_directory(tmp_path):
	b = write_bed(tmp_path / "beta_chr1_WT.bed", [("chr1", 1, 2, "n")])
	with pytest.raises(ValueError, match="paf"):
		pp.run_liftovers([b], ["beta"], tmp_path / "nope", "alpha", tmp_path)


# --- binning ---------------------------------------------------------------------------------------


def test_get_intro_df_template_has_ceil_chr_length_over_bin_columns_of_zero():
	df = pp.get_intro_df_template(100, 250)
	assert list(df.columns) == [0, 100, 200]
	assert df.shape == (1, 3)
	assert (df.values == 0).all()
	assert list(pp.get_intro_df_template(100, 300).columns) == [0, 100, 200]
	assert list(pp.get_intro_df_template(100, 301).columns) == [0, 100, 200, 300]


def _bins(rows, bin_size, chr_length):
	bed = pd.DataFrame(rows, columns=["Chromosome", "Start", "End", "Notes"])
	return pp.bed_to_bins(bed, bin_size, chr_length)


def test_bed_to_bins_marks_bins_from_the_rounded_start_up_to_the_rounded_end():
	out = _bins([("c", 100, 400, "x")], 100, 1000)
	assert list(out.columns) == ["introgression"]
	assert out.index.tolist() == list(range(0, 1000, 100))
	assert out["introgression"].tolist() == [0, 1, 1, 1, 0, 0, 0, 0, 0, 0]


def test_bed_to_bins_rounds_coordinates_to_the_nearest_bin():
	# 130 rounds to 100, 380 rounds to 400: bins 100, 200, 300
	out = _bins([("c", 130, 380, "x")], 100, 1000)
	assert out["introgression"].tolist() == [0, 1, 1, 1, 0, 0, 0, 0, 0, 0]
	# 160 rounds to 200, 440 rounds to 400: bins 200, 300
	out = _bins([("c", 160, 440, "x")], 100, 1000)
	assert out["introgression"].tolist() == [0, 0, 1, 1, 0, 0, 0, 0, 0, 0]


def test_bed_to_bins_marks_a_region_inside_one_bin_only_if_at_least_a_quarter_bin_long():
	# both ends round to bin 100; 30 bp >= 25 marks bin 100, 20 bp < 25 marks nothing
	assert _bins([("c", 110, 140, "x")], 100, 300)["introgression"].tolist() == [0, 1, 0]
	assert _bins([("c", 110, 130, "x")], 100, 300)["introgression"].tolist() == [0, 0, 0]


def test_bed_to_bins_accumulates_several_regions_and_clips_past_the_last_bin():
	out = _bins([("c", 0, 100, "x"), ("c", 200, 300, "y"), ("c", 900, 2000, "z")], 100, 450)
	assert out.index.tolist() == [0, 100, 200, 300, 400]
	assert out["introgression"].tolist() == [1, 0, 1, 0, 0]


# --- gap filling and small-region removal -----------------------------------------------------------


def test_fill_gaps_fills_only_gaps_between_introgressions_up_to_gap_size():
	row = pd.Series([0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1])
	assert list(pp.fill_gaps(row, 2)) == [0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1]
	assert list(pp.fill_gaps(row, 3)) == [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
	assert list(pp.fill_gaps(row, 0)) == list(row)


def test_fill_gaps_does_not_modify_its_input_and_leaves_leading_zeros():
	row = pd.Series([0, 0, 1, 0, 1])
	out = pp.fill_gaps(row, 1)
	assert list(out) == [0, 0, 1, 1, 1]
	assert list(row) == [0, 0, 1, 0, 1]


def test_remove_small_regions_drops_runs_shorter_than_min_size():
	row = pd.Series([1, 1, 0, 1, 0, 1, 1, 1, 0, 1])
	assert list(pp.remove_small_regions(row, 2)) == [1, 1, 0, 0, 0, 1, 1, 1, 0, 0]
	assert list(pp.remove_small_regions(row, 3)) == [0, 0, 0, 0, 0, 1, 1, 1, 0, 0]
	assert list(pp.remove_small_regions(row, 4)) == [0] * 10
	assert list(pp.remove_small_regions(row, 1)) == list(row)


def test_remove_small_regions_keeps_a_run_that_reaches_the_end_if_long_enough():
	row = pd.Series([0, 1, 1, 1])
	assert list(pp.remove_small_regions(row, 3)) == [0, 1, 1, 1]
	assert list(pp.remove_small_regions(row, 4)) == [0, 0, 0, 0]
	assert list(row) == [0, 1, 1, 1]


# --- file name parsing ---------------------------------------------------------------------------------


def test_get_bed_pieces_splits_accession_chromosome_and_type():
	assert pp.get_bed_pieces("alpha_chr1_REF.bed", ["alpha", "beta"]) == ("chr1", "alpha", "REF")
	assert pp.get_bed_pieces(Path("/x/y/beta_chr2_WT.bed"), ["alpha", "beta"]) == ("chr2", "beta", "WT")


def test_get_bed_pieces_allows_underscores_in_chromosome_and_accession():
	assert pp.get_bed_pieces("alpha_chr_1_WT.bed", ["alpha"]) == ("chr_1", "alpha", "WT")
	assert pp.get_bed_pieces("my_acc_chr1_WT.bed", ["my_acc"]) == ("chr1", "my_acc", "WT")


def test_get_bed_pieces_picks_the_longest_matching_accession():
	cands = ["my", "my_acc"]
	assert pp.get_bed_pieces("my_acc_chr1_WT.bed", cands) == ("chr1", "my_acc", "WT")
	assert pp.get_bed_pieces("my_chr1_WT.bed", cands) == ("chr1", "my", "WT")


@pytest.mark.parametrize("name", ["alpha.bed", "zeta_chr1_WT.bed", "alpha_WT.bed", "alphabet_chr1_WT.bed"])
def test_get_bed_pieces_rejects_names_it_cannot_parse(name):
	with pytest.raises(ValueError, match="Unable to parse"):
		pp.get_bed_pieces(name, ["alpha"])

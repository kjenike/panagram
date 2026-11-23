"""panagram.introgressions.score_introgressions: reading, merging, thresholding, confusion metrics, CLI.

Truth is worked by hand on tiny tables: each case lists the predicted and ground-truth bins,
and the expected counts follow from comparing them cell by cell.

CLI chromosome used: alpha chr1 of the synthetic index has 2636 k-mer positions, i.e. six bins
of 500 (0, 500, ..., 2500); alpha chr2 has 1686 positions (four bins); gamma chr2 has 686 (two).
"""

import sys

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import pytest

import helpers_introgressions as H
from panagram.introgressions import score_introgressions as sc

# ------------------------------------------------------------------ file helpers


def write_gt(path, rows, bins):
	"""Jaccard text file: header Sample + bin starts, one row per accession."""
	with open(path, "w") as f:
		f.write("\t".join(["Sample", *map(str, bins)]) + "\n")
		for name, vals in rows.items():
			f.write("\t".join([name, *("" if v is None else repr(v) for v in vals)]) + "\n")
	return path


def write_bed(path, rows):
	with open(path, "w") as f:
		for r in rows:
			f.write("\t".join(map(str, r)) + "\n")
	return path


# ------------------------------------------------------------------ read_text_file / merge_text_files


def test_read_text_file_uses_first_column_as_index_ints_as_columns_and_fills_nan_with_zero(tmp_path):
	p = write_gt(tmp_path / "t.txt", {"a": [0.1, None, 0.3], "b": [1.0, 0.5, 0.0]}, [0, 500, 1000])
	df = sc.read_text_file(p)
	assert list(df.index) == ["a", "b"]
	assert list(df.columns) == [0, 500, 1000]
	assert all(isinstance(c, (int, np.integer)) for c in df.columns)
	assert df.loc["a"].tolist() == [0.1, 0.0, 0.3]
	assert df.loc["b"].tolist() == [1.0, 0.5, 0.0]


def test_merge_text_files_takes_elementwise_maximum_and_unions_accessions(tmp_path):
	f1 = write_gt(tmp_path / "1.txt", {"a": [0.1, 0.9, 0.3], "b": [0.5, 0.5, 0.5]}, [0, 100, 200])
	f2 = write_gt(tmp_path / "2.txt", {"a": [0.4, 0.2, 0.3], "c": [0.7, 0.0, 0.1]}, [0, 100, 200])
	out = sc.merge_text_files([f1, f2])
	assert sorted(out.index) == ["a", "b", "c"]
	assert out.loc["a"].tolist() == [0.4, 0.9, 0.3]
	assert out.loc["b"].tolist() == [0.5, 0.5, 0.5]
	assert out.loc["c"].tolist() == [0.7, 0.0, 0.1]


def test_merge_text_files_single_file_is_that_file(tmp_path):
	f1 = write_gt(tmp_path / "1.txt", {"a": [0.1, 0.9]}, [0, 100])
	assert sc.merge_text_files([f1]).loc["a"].tolist() == [0.1, 0.9]


# ------------------------------------------------------------------ merge_bed_files


def test_merge_bed_files_rounds_to_nearest_bin_and_sorts_accessions(tmp_path):
	"""Bins of 100 over 500 bases. 100-299 -> start bin round(1.0)=1, end bin round(2.99)=3: bins 100 and 200."""
	b1 = write_bed(tmp_path / "zed_chr1_DON.bed", [("chr1", 100, 299, "x")])
	b2 = write_bed(tmp_path / "abe_chr1_DON.bed", [("chr1", 0, 99, "x"), ("chr1", 300, 499, "x")])
	out = sc.merge_bed_files([b1, b2], 100, 500)
	assert list(out.index) == ["abe", "zed"]
	assert out.index.name == "Sample"
	assert list(out.columns) == [0, 100, 200, 300, 400]
	assert out.loc["zed"].tolist() == [0, 1, 1, 0, 0]
	assert out.loc["abe"].tolist() == [1, 0, 0, 1, 1]


def test_merge_bed_files_empty_bed_gives_an_all_zero_row_and_keeps_underscored_names(tmp_path):
	b1 = write_bed(tmp_path / "my_acc_chr2_REF.bed", [])
	b2 = write_bed(tmp_path / "other_chr2_REF.bed", [("chr2", 200, 399, "x")])
	out = sc.merge_bed_files([b1, b2], 100, 400)
	assert list(out.index) == ["my_acc", "other"]
	assert out.loc["my_acc"].tolist() == [0, 0, 0, 0]
	assert out.loc["other"].tolist() == [0, 0, 1, 1]


def test_merge_bed_files_short_entry_within_a_bin_counts_from_a_quarter_bin(tmp_path):
	"""110-130 stays in bin 100 and spans 20 < 25 bases: not marked; 110-140 spans 30 >= 25: marked."""
	short = write_bed(tmp_path / "a_c_G.bed", [("c", 110, 130, "x")])
	long_ = write_bed(tmp_path / "b_c_G.bed", [("c", 110, 140, "x")])
	out = sc.merge_bed_files([short, long_], 100, 300)
	assert out.loc["a"].tolist() == [0, 0, 0]
	assert out.loc["b"].tolist() == [0, 1, 0]


def test_merge_bed_files_chromosome_length_sets_number_of_bins(tmp_path):
	b = write_bed(tmp_path / "a_c_G.bed", [])
	assert sc.merge_bed_files([b], 100, 250).shape == (1, 3)
	assert sc.merge_bed_files([b], 100, 300).shape == (1, 3)
	assert sc.merge_bed_files([b], 100, 301).shape == (1, 4)


# ------------------------------------------------------------------ threshold helpers


def test_threshold_introgressions_helper_cut_is_inclusive_above():
	df = pd.DataFrame({0: [0.2, 0.5], 1: [0.9, 0.0]}, index=["a", "b"])
	out = sc.threshold_introgressions_helper(df, 0.5)
	assert out.values.tolist() == [[0, 1], [1, 0]]
	assert out.dtypes.tolist() == [np.dtype(int)] * 2
	assert df.loc["a", 0] == 0.2  # input untouched


def test_threshold_introgressions_gt_uses_threshold_pred_needs_exact_ones():
	gt = pd.DataFrame([[0.2, 0.6, 0.5]], index=["a"], columns=[0, 1, 2])
	pred = pd.DataFrame([[0, 1, 1]], index=["a"], columns=[0, 1, 2])
	p, g = sc.threshold_introgressions(pred, gt, 0.5)
	assert g.values.tolist() == [[0, 1, 1]]
	assert p.values.tolist() == [[0, 1, 1]]
	# predictions below 1 are never counted even if above the ground-truth threshold
	p2, _ = sc.threshold_introgressions(pd.DataFrame([[0.9, 1.0]], index=["a"]), gt.iloc[:, :2], 0.1)
	assert p2.values.tolist() == [[0, 1]]


# ------------------------------------------------------------------ score_introgressions


def _df(rows, bins=None):
	cols = bins or list(range(len(next(iter(rows.values())))))
	return pd.DataFrame(rows, index=cols).T


def test_score_introgressions_confusion_counts_and_rates():
	"""Two accessions x four bins. a: pred 1100, gt 1010 -> TP 1, FP 1, FN 1, TN 1. b: pred 0001, gt 0001 -> TP 1, TN 3."""
	pred = _df({"a": [1, 1, 0, 0], "b": [0, 0, 0, 1]})
	gt = _df({"a": [1, 0, 1, 0], "b": [0, 0, 0, 1]})
	m = sc.score_introgressions(pred, gt)
	assert m.shape == (1, 8)
	r = m.iloc[0]
	assert (r["True Positive"], r["False Positive"], r["False Negative"], r["True Negative"]) == (2, 1, 1, 4)
	assert r["Accuracy"] == pytest.approx(6 / 8)
	assert r["Precision"] == pytest.approx(2 / 3)
	assert r["Recall"] == pytest.approx(2 / 3)
	assert r["FPR"] == pytest.approx(1 / 5)
	assert list(m.columns) == [
		"True Positive", "True Negative", "False Positive", "False Negative",
		"Accuracy", "Precision", "Recall", "FPR",
	]


def test_score_introgressions_uses_only_shared_accessions():
	pred = _df({"a": [1, 0], "only_pred": [1, 1]})
	gt = _df({"a": [1, 1], "only_gt": [0, 0]})
	r = sc.score_introgressions(pred, gt).iloc[0]
	assert (r["True Positive"], r["False Negative"], r["False Positive"], r["True Negative"]) == (1, 1, 0, 0)
	assert r["Accuracy"] == pytest.approx(0.5)


def test_score_introgressions_undefined_rates_are_nan():
	"""No predictions and no ground truth: precision, recall and FPR have zero denominators except FPR (all negatives)."""
	zero = _df({"a": [0, 0, 0]})
	r = sc.score_introgressions(zero, zero).iloc[0]
	assert r["Accuracy"] == 1.0
	assert np.isnan(r["Precision"]) and np.isnan(r["Recall"])
	assert r["FPR"] == 0.0
	ones = _df({"a": [1, 1]})
	r = sc.score_introgressions(ones, ones).iloc[0]
	assert r["Precision"] == 1.0 and r["Recall"] == 1.0
	assert np.isnan(r["FPR"])


def test_score_introgressions_all_wrong():
	r = sc.score_introgressions(_df({"a": [1, 0]}), _df({"a": [0, 1]})).iloc[0]
	assert (r["True Positive"], r["True Negative"], r["False Positive"], r["False Negative"]) == (0, 0, 1, 1)
	assert r["Accuracy"] == 0.0 and r["Precision"] == 0.0 and r["Recall"] == 0.0 and r["FPR"] == 1.0


# ------------------------------------------------------------------ rescale_introgressions_helper


def test_rescale_introgressions_helper_hand_worked():
	"""Bins of 100 -> 200 over 800 bases: bins 200-400 (bases 200-499) round to start bin 200, end bin 400: one bin of 200 at 200."""
	row = pd.Series([0, 0, 1, 1, 1, 0, 0, 0], index=[i * 100 for i in range(8)])
	out = sc.rescale_introgressions_helper(row, 100, 200, 800)
	assert list(out.index) == [0, 200, 400, 600]
	assert out.tolist() == [0, 1, 0, 0]


def test_rescale_introgressions_helper_to_finer_bins_and_empty_row():
	row = pd.Series([0, 1, 1, 0], index=[0, 200, 400, 600])
	out = sc.rescale_introgressions_helper(row, 200, 100, 800)
	assert out.tolist() == [0, 0, 1, 1, 1, 1, 0, 0]
	assert list(out.index) == [i * 100 for i in range(8)]
	empty = sc.rescale_introgressions_helper(pd.Series([0, 0, 0, 0], index=[0, 200, 400, 600]), 200, 100, 800)
	assert empty.tolist() == [0] * 8


# ------------------------------------------------------------------ merge_centromere_regions_helper


def test_merge_centromere_regions_helper_round_trips_rows_without_centromere_gaps():
	"""A run, two runs three bins apart, and no run: the bins -> BED -> bins conversion gives the row back."""
	fasta = pd.DataFrame({"Name": ["nan"], "Sequence": ["A" * 800]})
	for vals in ([0, 1, 1, 0, 0, 0, 0, 0], [1, 0, 0, 0, 1, 1, 0, 0], [0, 0, 1, 0, 0, 0, 1, 0]):
		row = pd.Series(vals, index=[i * 100 for i in range(8)])
		out = sc.merge_centromere_regions_helper(row, 100, 800, fasta)
		assert out.tolist() == vals


# ------------------------------------------------------------------ create_scored_heatmap


@pytest.fixture
def figures(monkeypatch):
	seen = []

	def fake(self, file, *args, **kwargs):
		seen.append((str(file), self))

	monkeypatch.setattr(go.Figure, "write_image", fake)
	return seen


def test_create_scored_heatmap_codes_and_csv(figures, tmp_path):
	"""TP 5, TN 4, FP 3, FN 2; accessions not in both inputs are dropped; rows sorted by name."""
	pred = _df({"b": [1, 1, 0, 0], "a": [0, 0, 1, 1], "extra": [1, 1, 1, 1]}, [0, 100, 200, 300])
	gt = _df({"b": [1, 0, 1, 0], "a": [0, 0, 1, 0]}, [0, 100, 200, 300])
	out = tmp_path / "chr1_DON.png"
	sc.create_scored_heatmap(pred, gt, out)
	csv = pd.read_csv(out.with_suffix(".csv"), sep="\t", index_col=0)
	assert list(csv.index) == ["a", "b"]
	assert csv.loc["a"].tolist() == [4, 4, 5, 3]
	assert csv.loc["b"].tolist() == [5, 3, 2, 4]
	assert list(csv.columns) == ["0", "100", "200", "300"]
	path, fig = figures[0]
	assert path == str(out)
	assert fig.layout.title.text == "chr1_DON"
	assert fig.layout.xaxis.dtick == 2000000
	assert np.asarray(fig.data[0].z).tolist() == [[4, 4, 5, 3], [5, 3, 2, 4]]
	assert pred.loc["b"].tolist() == [1, 1, 0, 0]  # caller's frame unchanged


def test_create_scored_heatmap_groups_order_rows_and_dtick(figures, tmp_path):
	pred = _df({"a": [1, 0], "b": [0, 0], "c": [0, 1]})
	gt = _df({"a": [1, 0], "b": [0, 1], "c": [0, 1]})
	groups = pd.DataFrame({"group": ["x", "y", "z", "w"]}, index=["c", "zzz", "a", "b"])
	out = tmp_path / "g.png"
	sc.create_scored_heatmap(pred, gt, out, groups=groups, xaxis_dtick=7)
	csv = pd.read_csv(out.with_suffix(".csv"), sep="\t", index_col=0)
	assert list(csv.index) == ["c", "a", "b"]  # unknown names in groups are skipped
	assert figures[0][1].layout.xaxis.dtick == 7


def test_create_scored_heatmap_renders_when_a_renderer_exists(tmp_path):
	pytest.importorskip("kaleido")
	if not H.can_render_images(tmp_path):
		pytest.skip("plotly cannot write images here (Chrome missing)")
	pred = _df({"a": [1, 0]})
	out = tmp_path / "r.svg"
	sc.create_scored_heatmap(pred, pred, out)
	assert out.stat().st_size > 0


# ------------------------------------------------------------------ main()


BINS6 = [0, 500, 1000, 1500, 2000, 2500]


@pytest.fixture
def score_env(tmp_path):
	"""Beds for alpha (bins 500, 1000) and beta (none) on chr1, and a ground truth with Jaccard values.

	Ground truth at threshold 0.5: alpha 011000, beta 001011.
	"""
	pre, gdt = tmp_path / "pre", tmp_path / "gdt"
	pre.mkdir()
	gdt.mkdir()
	write_bed(pre / "alpha_chr1_DON.bed", [("chr1", 500, 1499, "x")])
	write_bed(pre / "beta_chr1_DON.bed", [])
	write_gt(
		gdt / "chr1_DON.txt",
		{"alpha": [0.1, 0.9, 0.8, 0.2, 0.0, 0.0], "beta": [0.0, 0.0, 0.7, 0.1, 0.6, 0.9]},
		BINS6,
	)
	return pre, gdt


@pytest.fixture
def run_score(monkeypatch, anchored_index, tmp_path):
	def run(*args, idx=None, out=None):
		out = out or (tmp_path / "out")
		argv = [
			"score_introgressions",
			"--idx", str(idx or anchored_index),
			"--out", str(out),
			*[str(a) for a in args],
		]
		monkeypatch.setattr(sys, "argv", argv)
		sc.main()
		return out

	return run


def _metrics(out, kind):
	return pd.read_csv(out / f"metrics_{kind}.tsv", sep="\t", index_col=0)


def test_main_scores_a_folder_against_ground_truth(run_score, score_env):
	"""alpha pred 011000 = gt: TP 2. beta pred none, gt 001011: FN 3. TN = 12 - 2 - 3 = 7."""
	pre, gdt = score_env
	out = run_score("--pre", pre, "--gdt", gdt, "--ref", "alpha", "--bin", 500, "--thr", 0.5)
	m = _metrics(out, "DON")
	assert list(m.index) == ["chr1"]
	r = m.iloc[0]
	assert (r["True Positive"], r["True Negative"], r["False Positive"], r["False Negative"]) == (2, 7, 0, 3)
	assert r["Accuracy"] == pytest.approx(9 / 12)
	assert r["Precision"] == 1.0
	assert r["Recall"] == pytest.approx(0.4)
	assert r["FPR"] == 0.0
	pred_txt = pd.read_csv(out / "pred" / "chr1_DON.txt", sep="\t", index_col=0)
	assert pred_txt.loc["alpha"].tolist() == [0, 1, 1, 0, 0, 0]
	assert pred_txt.loc["beta"].tolist() == [0, 0, 0, 0, 0, 0]
	assert list((out / "gt_postprocessed").iterdir()) == []  # no --act: nothing postprocessed
	assert not (out / "heatmaps").exists()


def test_main_accepts_a_single_bed_file(run_score, score_env):
	pre, gdt = score_env
	out = run_score("--pre", pre / "alpha_chr1_DON.bed", "--gdt", gdt, "--ref", "alpha", "--bin", 500)
	r = _metrics(out, "DON").iloc[0]
	assert (r["True Positive"], r["False Negative"], r["True Negative"]) == (2, 0, 4)  # alpha only
	assert r["Recall"] == 1.0


def test_main_thr_changes_the_ground_truth_cut(run_score, score_env):
	"""At --thr 0.85 the ground truth is alpha bin 500 (0.9) and beta bin 2500 (0.9): alpha pred has bin 1000 extra -> FP 1; beta bin missed -> FN 1."""
	pre, gdt = score_env
	out = run_score("--pre", pre, "--gdt", gdt, "--ref", "alpha", "--bin", 500, "--thr", 0.85)
	r = _metrics(out, "DON").iloc[0]
	assert (r["True Positive"], r["False Positive"], r["False Negative"], r["True Negative"]) == (1, 1, 1, 9)


def test_main_act_fgap_fills_ground_truth_gaps_only(run_score, score_env):
	"""beta gt 001011 -> with --gap 1 the one-bin gap fills: 001111, so FN rises to 4 and the postprocessed gt is saved."""
	pre, gdt = score_env
	out = run_score(
		"--pre", pre, "--gdt", gdt, "--ref", "alpha", "--bin", 500, "--act", "fgap", "--gap", 1
	)
	r = _metrics(out, "DON").iloc[0]
	assert (r["True Positive"], r["False Negative"], r["False Positive"]) == (2, 4, 0)
	saved = pd.read_csv(out / "gt_postprocessed" / "chr1_DON.txt", sep="\t", index_col=0)
	assert saved.loc["beta"].tolist() == [0, 0, 1, 1, 1, 1]
	assert saved.loc["alpha"].tolist() == [0, 1, 1, 0, 0, 0]
	# prediction file is unaffected by actions
	assert pd.read_csv(out / "pred" / "chr1_DON.txt", sep="\t", index_col=0).loc["beta"].sum() == 0


def test_main_act_rmbn_removes_small_ground_truth_regions(run_score, score_env):
	"""--min 2: beta's single-bin region at 1000 (001011) vanishes, its two-bin region and alpha's stay."""
	pre, gdt = score_env
	out = run_score(
		"--pre", pre, "--gdt", gdt, "--ref", "alpha", "--bin", 500, "--act", "rmbn", "--min", 2
	)
	r = _metrics(out, "DON").iloc[0]
	assert (r["True Positive"], r["False Negative"], r["False Positive"], r["True Negative"]) == (2, 2, 0, 8)
	saved = pd.read_csv(out / "gt_postprocessed" / "chr1_DON.txt", sep="\t", index_col=0)
	assert saved.loc["beta"].tolist() == [0, 0, 0, 0, 1, 1]


def test_main_actions_run_in_the_order_given(run_score, score_env, tmp_path):
	"""fgap then rmbn (--gap 1, --min 5): beta becomes 001111 (4 bins) and is removed; alpha (2 bins) is removed too."""
	pre, gdt = score_env
	out = run_score(
		"--pre", pre, "--gdt", gdt, "--ref", "alpha", "--bin", 500,
		"--act", "fgap", "rmbn", "--gap", 1, "--min", 5,
	)
	saved = pd.read_csv(out / "gt_postprocessed" / "chr1_DON.txt", sep="\t", index_col=0)
	assert saved.values.sum() == 0
	r = _metrics(out, "DON").iloc[0]
	assert (r["True Positive"], r["False Positive"], r["False Negative"]) == (0, 2, 0)
	assert r["Precision"] == 0.0 and np.isnan(r["Recall"])
	# rmbn (--min 2) then fgap: beta's single bin at 1000 is dropped first, so no gap is filled
	out2 = run_score(
		"--pre", pre, "--gdt", gdt, "--ref", "alpha", "--bin", 500,
		"--act", "rmbn", "fgap", "--gap", 1, "--min", 2, out=tmp_path / "out2",
	)
	saved2 = pd.read_csv(out2 / "gt_postprocessed" / "chr1_DON.txt", sep="\t", index_col=0)
	assert saved2.loc["beta"].tolist() == [0, 0, 0, 0, 1, 1]
	assert saved2.loc["alpha"].tolist() == [0, 1, 1, 0, 0, 0]


def test_main_rescales_predictions_to_the_ground_truth_bin_size(run_score, tmp_path):
	"""Beds at --bin 250, ground truth in bins of 500: the same alpha region 500-1499 scores as before."""
	pre, gdt = tmp_path / "pre", tmp_path / "gdt"
	pre.mkdir()
	gdt.mkdir()
	write_bed(pre / "alpha_chr1_DON.bed", [("chr1", 500, 1499, "x")])
	write_gt(gdt / "chr1_DON.txt", {"alpha": [0.1, 0.9, 0.8, 0.2, 0.0, 0.0]}, BINS6)
	out = run_score("--pre", pre, "--gdt", gdt, "--ref", "alpha", "--bin", 250)
	r = _metrics(out, "DON").iloc[0]
	assert (r["True Positive"], r["False Positive"], r["False Negative"], r["True Negative"]) == (2, 0, 0, 4)
	saved = pd.read_csv(out / "pred" / "chr1_DON.txt", sep="\t", index_col=0)
	assert saved.shape[1] == 11  # the raw prediction table is saved at the --bin size (ceil(2636 / 250))


def test_main_act_fcen_after_rescaling_uses_the_ground_truth_bin_size(run_score, tmp_path):
	"""Beds at --bin 250, ground truth in bins of 500, --act fcen: the reference has no run of 50 N, so the rows are
	unchanged and the score is the same as without the action (TP 2, TN 4). The bins of the postprocessed ground truth
	are the 500-base ones."""
	pre, gdt = tmp_path / "pre", tmp_path / "gdt"
	pre.mkdir()
	gdt.mkdir()
	write_bed(pre / "alpha_chr1_DON.bed", [("chr1", 500, 1499, "x")])
	write_gt(gdt / "chr1_DON.txt", {"alpha": [0.1, 0.9, 0.8, 0.2, 0.0, 0.0]}, BINS6)
	out = run_score("--pre", pre, "--gdt", gdt, "--ref", "alpha", "--bin", 250, "--act", "fcen")
	r = _metrics(out, "DON").iloc[0]
	assert (r["True Positive"], r["False Positive"], r["False Negative"], r["True Negative"]) == (2, 0, 0, 4)
	saved = pd.read_csv(out / "gt_postprocessed" / "chr1_DON.txt", sep="\t", index_col=0)
	assert list(saved.columns) == [str(b) for b in BINS6]
	assert saved.loc["alpha"].tolist() == [0, 1, 1, 0, 0, 0]


def test_main_vis_writes_scored_heatmap_tables(run_score, score_env, figures, tmp_path):
	"""Codes per bin: alpha TN TP TP TN TN TN, beta TN TN FN TN FN FN; --grp lists beta first."""
	pre, gdt = score_env
	groups = H.write_groups(tmp_path / "g.tsv", {"beta": "x", "alpha": "y"})
	out = run_score(
		"--pre", pre, "--gdt", gdt, "--ref", "alpha", "--bin", 500, "--vis", "--grp", groups
	)
	assert [p for p, _ in figures] == [str(out / "heatmaps" / "chr1_DON.png")]
	csv = pd.read_csv(out / "heatmaps" / "chr1_DON.csv", sep="\t", index_col=0)
	assert list(csv.index) == ["beta", "alpha"]
	assert csv.loc["alpha"].tolist() == [4, 5, 5, 4, 4, 4]
	assert csv.loc["beta"].tolist() == [4, 4, 2, 4, 2, 2]


def test_main_vis_without_groups_sorts_accessions(run_score, score_env, figures):
	pre, gdt = score_env
	out = run_score("--pre", pre, "--gdt", gdt, "--ref", "alpha", "--bin", 500, "--vis")
	csv = pd.read_csv(out / "heatmaps" / "chr1_DON.csv", sep="\t", index_col=0)
	assert list(csv.index) == ["alpha", "beta"]


def test_main_ref_and_merged_types_use_cmp_groups_and_take_the_max_of_their_ground_truths(
	run_score, tmp_path
):
	"""REF beds are scored against max(chr1_DON, chr1_OTH): alpha gt 0.1/0.9 and 0.6/0.2 -> max 0.6/0.9."""
	pre, gdt = tmp_path / "pre", tmp_path / "gdt"
	pre.mkdir()
	gdt.mkdir()
	write_bed(pre / "alpha_chr1_REF.bed", [("chr1", 0, 999, "x")])
	write_bed(pre / "alpha_chr1_merged.bed", [("chr1", 500, 999, "x")])
	bins = BINS6
	write_gt(gdt / "chr1_DON.txt", {"alpha": [0.1, 0.9, 0.0, 0.0, 0.0, 0.0]}, bins)
	write_gt(gdt / "chr1_OTH.txt", {"alpha": [0.6, 0.2, 0.0, 0.0, 0.0, 0.0]}, bins)
	out = run_score(
		"--pre", pre, "--gdt", gdt, "--ref", "alpha", "--bin", 500, "--cmp", "DON", "OTH"
	)
	for kind, tp, fn in (("REF", 2, 0), ("merged", 1, 1)):
		r = _metrics(out, kind).iloc[0]
		assert (r["True Positive"], r["False Negative"], r["False Positive"]) == (tp, fn, 0), kind


def test_main_scores_several_chromosomes_with_the_reference_chromosome_lengths(run_score, tmp_path):
	"""Rows of the metrics file are chromosomes. alpha chr2 has four bins, gamma chr2 two (--ref sets the lengths)."""
	pre, gdt = tmp_path / "pre", tmp_path / "gdt"
	pre.mkdir()
	gdt.mkdir()
	write_bed(pre / "alpha_chr1_DON.bed", [("chr1", 0, 999, "x")])
	write_bed(pre / "alpha_chr2_DON.bed", [("chr2", 500, 1499, "x")])
	write_gt(gdt / "chr1_DON.txt", {"alpha": [0.9, 0.9, 0.0, 0.0, 0.0, 0.0]}, BINS6)
	write_gt(gdt / "chr2_DON.txt", {"alpha": [0.0, 0.9, 0.0, 0.0]}, [0, 500, 1000, 1500])
	out = run_score("--pre", pre, "--gdt", gdt, "--ref", "alpha", "--bin", 500)
	m = _metrics(out, "DON")
	assert list(m.index) == ["chr1", "chr2"]
	assert m["True Positive"].tolist() == [2, 1]
	assert m["False Positive"].tolist() == [0, 1]
	assert m["True Negative"].tolist() == [4, 2]
	assert pd.read_csv(out / "pred" / "chr2_DON.txt", sep="\t", index_col=0).shape == (1, 4)
	# gamma chr2 is only two bins long
	pre2, gdt2 = tmp_path / "pre2", tmp_path / "gdt2"
	pre2.mkdir()
	gdt2.mkdir()
	write_bed(pre2 / "alpha_chr2_DON.bed", [("chr2", 0, 499, "x")])
	write_gt(gdt2 / "chr2_DON.txt", {"alpha": [0.9, 0.0]}, [0, 500])
	out2 = run_score("--pre", pre2, "--gdt", gdt2, "--ref", "gamma", "--bin", 500, out=tmp_path / "o2")
	assert pd.read_csv(out2 / "pred" / "chr2_DON.txt", sep="\t", index_col=0).shape == (1, 2)
	assert _metrics(out2, "DON").iloc[0]["True Positive"] == 1


# ------------------------------------------------------------------ main(): argument errors


def test_main_requires_pre_gdt_idx_ref_out(monkeypatch):
	full = ["--pre", "p", "--gdt", "g", "--idx", "i", "--ref", "r", "--out", "o"]
	for drop in range(0, len(full), 2):
		argv = ["x", *full[:drop], *full[drop + 2 :]]
		monkeypatch.setattr(sys, "argv", argv)
		with pytest.raises(SystemExit) as e:
			sc.main()
		assert e.value.code == 2


@pytest.mark.parametrize(
	"extra,fragment",
	[
		(["--act", "foo"], "Unrecognized action"),
		(["--act", "fgap", "bogus"], "Unrecognized action"),
	],
)
def test_main_rejects_unknown_actions(run_score, score_env, extra, fragment):
	pre, gdt = score_env
	with pytest.raises(ValueError, match=fragment):
		run_score("--pre", pre, "--gdt", gdt, "--ref", "alpha", *extra)


def test_main_rejects_missing_inputs(run_score, score_env, tmp_path):
	pre, gdt = score_env
	with pytest.raises(ValueError, match="Ground truth"):
		run_score("--pre", pre, "--gdt", tmp_path / "nope", "--ref", "alpha")
	with pytest.raises(ValueError, match="Index directory"):
		run_score("--pre", pre, "--gdt", gdt, "--ref", "alpha", idx=tmp_path / "noindex")
	with pytest.raises(ValueError, match="Bed file"):
		run_score("--pre", tmp_path / "nobeds", "--gdt", gdt, "--ref", "alpha")


def test_main_requires_cmp_for_ref_and_merged_and_all_group_files(run_score, tmp_path):
	pre, gdt = tmp_path / "pre", tmp_path / "gdt"
	pre.mkdir()
	gdt.mkdir()
	write_bed(pre / "alpha_chr1_REF.bed", [])
	write_gt(gdt / "chr1_DON.txt", {"alpha": [0.0] * 6}, BINS6)
	with pytest.raises(ValueError, match="--cmp"):
		run_score("--pre", pre, "--gdt", gdt, "--ref", "alpha", "--bin", 500)
	with pytest.raises(ValueError, match="all files"):
		run_score("--pre", pre, "--gdt", gdt, "--ref", "alpha", "--bin", 500, "--cmp", "DON", "OTH")


def test_main_missing_ground_truth_file_for_a_type(run_score, tmp_path):
	pre, gdt = tmp_path / "pre", tmp_path / "gdt"
	pre.mkdir()
	gdt.mkdir()
	write_bed(pre / "alpha_chr1_DON.bed", [])
	with pytest.raises(ValueError, match="chr1_DON.txt"):
		run_score("--pre", pre, "--gdt", gdt, "--ref", "alpha", "--bin", 500)

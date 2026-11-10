"""Tests of panagram/introgressions/visualize_introgressions.py.

Figures are captured through a write_image stub (no Chrome needed) and checked as plotly figure
objects against values worked by hand. One test renders for real and skips without kaleido/Chrome.
"""

import math
import re
import sys

import numpy as np
import pandas as pd
import pytest

from helpers_introgression_post import figure_log  # noqa: F401
from panagram.introgressions import visualize_introgressions as vis


def _metrics_file(root, thr, intro, rows):
	"""root/<root.name>_<thr>/scored/metrics_<intro>.tsv with (chromosome, TP, FP, FN, TN[, precision, recall]) rows."""
	d = root / f"{root.name}_{thr}" / "scored"
	d.mkdir(parents=True, exist_ok=True)
	cols = ["True Positive", "False Positive", "False Negative", "True Negative", "Precision", "Recall"]
	df = pd.DataFrame(
		[r[1:] + (float("nan"),) * (len(cols) - len(r) + 1) for r in rows], columns=cols, index=[r[0] for r in rows]
	)
	df.to_csv(d / f"metrics_{intro}.tsv", sep="\t")


def _mcc(tp, fp, fn, tn):
	return (tp * tn - fp * fn) / math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))


@pytest.fixture
def res(tmp_path):
	"""Results folder 'res' with two thresholds, intro type WT and two chromosomes each.

	threshold 0.2: TP 5+3=8, FP 1+1=2, FN 1+1=2, TN 40+48=88 -> precision .8, recall .8
	threshold 0.4: TP 3+1=4, FP 0+0=0, FN 3+3=6, TN 40+50=90 -> precision 1,  recall .4
	"""
	root = tmp_path / "res"
	root.mkdir()
	_metrics_file(root, 0.2, "WT", [("chr1", 5, 1, 1, 40), ("chr2", 3, 1, 1, 48)])
	_metrics_file(root, 0.4, "WT", [("chr1", 3, 0, 3, 40), ("chr2", 1, 0, 3, 50)])
	return root


# --- calculate_pr_auc ------------------------------------------------------------------------


def test_calculate_pr_auc_sums_counts_over_chromosomes_and_computes_metrics(res):
	df = vis.calculate_pr_auc(res, "WT", [0.2, 0.4])
	rows = df.set_index("threshold")
	assert rows.loc[0.2, ["TP", "FP", "FN", "TN"]].tolist() == [8, 2, 2, 88]
	assert rows.loc[0.4, ["TP", "FP", "FN", "TN"]].tolist() == [4, 0, 6, 90]
	assert rows.loc[0.2, "precision"] == pytest.approx(8 / 10)
	assert rows.loc[0.2, "recall"] == pytest.approx(8 / 10)
	assert rows.loc[0.4, "precision"] == pytest.approx(1.0)
	assert rows.loc[0.4, "recall"] == pytest.approx(4 / 10)
	assert rows.loc[0.2, "MCC"] == pytest.approx(_mcc(8, 2, 2, 88))
	assert rows.loc[0.4, "MCC"] == pytest.approx(_mcc(4, 0, 6, 90))


def test_calculate_pr_auc_adds_the_zero_recall_point_below_when_recall_rises_with_threshold(res):
	# recall .4 at 0.2 and .8 at 0.4 would rise with threshold: build that case by swapping the files
	root = res.parent / "up"
	root.mkdir()
	_metrics_file(root, 0.2, "WT", [("chr1", 3, 0, 3, 40), ("chr2", 1, 0, 3, 50)])
	_metrics_file(root, 0.4, "WT", [("chr1", 5, 1, 1, 40), ("chr2", 3, 1, 1, 48)])
	df = vis.calculate_pr_auc(root, "WT", [0.2, 0.4])
	assert df["threshold"].tolist() == pytest.approx([-0.8, 0.2, 0.4])
	first = df.iloc[0]
	assert (first["precision"], first["recall"], first["MCC"], first["TN"], first["TP"]) == (1, 0, 0, 1, 0)
	assert df["recall"].tolist() == pytest.approx([0, 0.4, 0.8])


def test_calculate_pr_auc_adds_the_zero_recall_point_above_when_recall_falls_with_threshold(res):
	df = vis.calculate_pr_auc(res, "WT", [0.2, 0.4])
	assert df["threshold"].tolist() == pytest.approx([0.2, 0.4, 1.4])
	last = df.iloc[-1]
	assert (last["precision"], last["recall"], last["MCC"]) == (1, 0, 0)


def test_calculate_pr_auc_prints_the_trapezoid_area_under_precision_against_recall(tmp_path, capsys):
	root = tmp_path / "up"
	root.mkdir()
	_metrics_file(root, 0.2, "WT", [("chr1", 3, 0, 3, 40), ("chr2", 1, 0, 3, 50)])  # precision 1, recall .4
	_metrics_file(root, 0.4, "WT", [("chr1", 5, 1, 1, 40), ("chr2", 3, 1, 1, 48)])  # precision .8, recall .8
	vis.calculate_pr_auc(root, "WT", [0.2, 0.4])
	# points (recall, precision): (0, 1), (.4, 1), (.8, .8): .4 * 1 + .4 * (1 + .8) / 2 = .76
	out = capsys.readouterr().out
	area = float(re.search(r"PR AUC WT:\s*(\S+)", out).group(1))
	assert area == pytest.approx(0.76)
	max_mcc = float(re.search(r"Max MCC WT:\s*(\S+)", out).group(1))
	assert max_mcc == pytest.approx(_mcc(8, 2, 2, 88))


def test_calculate_pr_auc_zero_denominators_give_zero_not_nan(tmp_path):
	root = tmp_path / "z"
	root.mkdir()
	_metrics_file(root, 0.1, "WT", [("chr1", 0, 0, 0, 10)])  # nothing called, nothing to find
	_metrics_file(root, 0.2, "WT", [("chr1", 0, 5, 0, 5)])  # only false positives
	df = vis.calculate_pr_auc(root, "WT", [0.1, 0.2]).set_index("threshold")
	assert df.loc[0.1, ["precision", "recall", "MCC"]].tolist() == [0, 0, 0]
	assert df.loc[0.2, "precision"] == 0
	assert df.loc[0.2, "recall"] == 0
	assert not df.isna().any().any()


def test_calculate_pr_auc_reports_a_missing_metrics_file(res):
	with pytest.raises(ValueError, match="Metrics file not found"):
		vis.calculate_pr_auc(res, "WT", [0.2, 0.9])
	with pytest.raises(ValueError, match="Metrics file not found"):
		vis.calculate_pr_auc(res, "REF", [0.2])


# --- curves ----------------------------------------------------------------------------------------


def test_create_mcc_curve_plots_mcc_against_threshold_and_names_the_file(res, figure_log):
	vis.create_mcc_curve(res, "WT", [0.2, 0.4])
	fig = figure_log.by_name("WT_mcc.svg")
	assert figure_log.paths()[0].parent == res
	assert len(fig.data) == 1
	assert list(fig.data[0].x) == pytest.approx([0.2, 0.4, 1.4])
	assert list(fig.data[0].y) == pytest.approx([_mcc(8, 2, 2, 88), _mcc(4, 0, 6, 90), 0])
	assert fig.layout.xaxis.title.text == "Threshold"
	assert fig.layout.yaxis.title.text == "Matthews Correlation Coefficient (MCC)"
	assert (fig.layout.width, fig.layout.height) == (500, 500)


def test_create_pr_curve_plots_precision_against_recall(res, figure_log):
	vis.create_pr_curve(res, "WT", [0.2, 0.4])
	fig = figure_log.by_name("WT_prc.svg")
	assert list(fig.data[0].x) == pytest.approx([0.8, 0.4, 0.0])
	assert list(fig.data[0].y) == pytest.approx([0.8, 1.0, 1.0])
	assert fig.layout.xaxis.title.text == "Recall"
	assert fig.layout.yaxis.title.text == "Precision"


def test_create_pr_curve_drops_points_with_zero_precision_and_zero_recall(tmp_path, figure_log):
	root = tmp_path / "z"
	root.mkdir()
	_metrics_file(root, 0.2, "WT", [("chr1", 0, 5, 5, 5)])  # TP 0 -> precision 0, recall 0
	_metrics_file(root, 0.4, "WT", [("chr1", 5, 5, 5, 5)])  # precision .5, recall .5
	vis.create_pr_curve(root, "WT", [0.2, 0.4])
	fig = figure_log.by_name("WT_prc.svg")
	# recall rises with threshold, so the (recall 0, precision 1) end point comes first; the
	# threshold 0.2 point (precision 0, recall 0) is dropped, leaving two points
	assert list(fig.data[0].x) == pytest.approx([0.0, 0.5])
	assert list(fig.data[0].y) == pytest.approx([1.0, 0.5])


def test_create_pr_curve_chromosomes_has_one_line_per_chromosome(res, tmp_path, figure_log):
	# the per-chromosome file carries Precision and Recall columns
	root = tmp_path / "pc"
	root.mkdir()
	for thr, vals in ((0.2, [(0.9, 0.5), (0.7, 0.6)]), (0.4, [(1.0, 0.2), (float("nan"), 0.0)])):
		d = root / f"pc_{thr}" / "scored"
		d.mkdir(parents=True)
		pd.DataFrame(vals, columns=["Precision", "Recall"], index=["chr1", "chr2"]).to_csv(
			d / "metrics_WT.tsv", sep="\t"
		)
	vis.create_pr_curve_chromosomes(root, "WT", [0.2, 0.4])
	fig = figure_log.by_name("WT_prcc.svg")
	by_name = {t.name: t for t in fig.data}
	assert set(by_name) == {"chr1", "chr2"}
	assert list(by_name["chr1"].x) == pytest.approx([0.5, 0.2])
	assert list(by_name["chr1"].y) == pytest.approx([0.9, 1.0])
	# the NaN precision of chr2 becomes 0
	assert list(by_name["chr2"].x) == pytest.approx([0.6, 0.0])
	assert list(by_name["chr2"].y) == pytest.approx([0.7, 0.0])
	assert fig.layout.title.text == "Precision-Recall Curves per Chromosome"


def _heatmap_csv(root, thr, chrom, intro, rows):
	"""root/<root.name>_<thr>/scored/heatmaps/<chrom>_<intro>.csv, codes 5 TP, 4 TN, 3 FP, 2 FN per bin."""
	d = root / f"{root.name}_{thr}" / "scored" / "heatmaps"
	d.mkdir(parents=True, exist_ok=True)
	pd.DataFrame(rows).T.to_csv(d / f"{chrom}_{intro}.csv", sep="\t")
	return d


def test_create_pr_curve_accessions_has_one_line_per_sample_from_summed_codes(tmp_path, figure_log):
	root = tmp_path / "acc"
	root.mkdir()
	# rows are samples, columns are bins; chr1 + chr2 are summed per sample
	# threshold 0.2: s1 chr1 TP2 FP1 FN1 TN2, chr2 TP1 TN3 FN1 -> TP3 FP1 FN2 TN5
	#                s2 chr1 TN4 FN2, chr2 TN3 -> TP0 FP0 FN2 TN7
	_heatmap_csv(root, 0.2, "chr1", "WT", {"s1": [5, 5, 3, 2, 4, 4], "s2": [4, 4, 4, 4, 2, 2]})
	_heatmap_csv(root, 0.2, "chr2", "WT", {"s1": [5, 4, 4, 4, 2, 4], "s2": [4, 4, 4, 4, 4, 4]})
	# threshold 0.4: s1 all TP/TN, s2 TP2 FP1
	_heatmap_csv(root, 0.4, "chr1", "WT", {"s1": [5, 5, 4, 4], "s2": [5, 5, 3, 4]})
	_heatmap_csv(root, 0.4, "chr2", "WT", {"s1": [4, 4, 4, 4], "s2": [4, 4, 4, 4]})
	vis.create_pr_curve_accessions(root, "WT", [0.2, 0.4])
	fig = figure_log.by_name("WT_prca.svg")
	by_name = {t.name: t for t in fig.data}
	assert set(by_name) == {"s1", "s2"}
	# s1: precision 3/4, recall 3/5 at 0.2; precision 1, recall 1 at 0.4
	assert list(by_name["s1"].x) == pytest.approx([3 / 5, 1.0])
	assert list(by_name["s1"].y) == pytest.approx([3 / 4, 1.0])
	# s2: no TP at 0.2 gives precision 0 and recall 0, which is plotted as precision 1;
	# at 0.4 TP 2, FP 1, FN 0 gives precision 2/3, recall 1
	assert list(by_name["s2"].x) == pytest.approx([0.0, 1.0])
	assert list(by_name["s2"].y) == pytest.approx([1.0, 2 / 3])


def test_create_pr_curve_accessions_requires_the_heatmaps_folder(res):
	with pytest.raises(ValueError, match="Heatmaps directory not found"):
		vis.create_pr_curve_accessions(res, "WT", [0.2])


# --- heatmap collage -----------------------------------------------------------------------------


def _png(path, color):
	from PIL import Image

	Image.new("RGB", (4, 4), color).save(path)


def _page_count(pdf_bytes):
	return len(re.findall(rb"/Type\s*/Page(?![A-Za-z])", pdf_bytes))


def test_collage_with_other_than_18_thresholds_warns_and_writes_nothing(res, capsys):
	vis.create_scored_heatmap_collage(res, "WT", [0.2, 0.4])
	assert "Expected 18 thresholds" in capsys.readouterr().out
	assert not (res / "WT_scored_heatmaps.pdf").exists()


def test_collage_writes_two_landscape_pages_per_chromosome(tmp_path):
	pytest.importorskip("reportlab")
	pytest.importorskip("PIL")
	root = tmp_path / "col"
	root.mkdir()
	thresholds = [round(0.1 + 0.05 * i, 2) for i in range(18)]
	for thr in thresholds:
		d = root / f"col_{thr}" / "scored" / "heatmaps"
		d.mkdir(parents=True)
		for chrom in ("chr1", "chr2", "chr3"):
			_png(d / f"{chrom}_WT.png", (200, 100, 50))
		_png(d / "chr1_REF.png", (0, 0, 0))  # another intro type, must be ignored
	vis.create_scored_heatmap_collage(root, "WT", thresholds)
	pdf = (root / "WT_scored_heatmaps.pdf").read_bytes()
	assert pdf.startswith(b"%PDF")
	assert _page_count(pdf) == 6
	box = re.search(rb"/MediaBox\s*\[\s*0\s+0\s+([\d.]+)\s+([\d.]+)\s*\]", pdf)
	width, height = float(box.group(1)), float(box.group(2))
	assert (width, height) == pytest.approx((841.89, 595.27), abs=0.01)


def test_collage_requires_the_heatmaps_folder_of_every_threshold(tmp_path):
	pytest.importorskip("reportlab")
	root = tmp_path / "col"
	root.mkdir()
	thresholds = [round(0.1 + 0.05 * i, 2) for i in range(18)]
	with pytest.raises(ValueError, match="Heatmaps directory not found"):
		vis.create_scored_heatmap_collage(root, "WT", thresholds)


# --- distance heatmaps -----------------------------------------------------------------------------


def _distances(path, index, bins, values):
	pd.DataFrame(values, index=index, columns=bins).rename_axis("Sample").to_csv(path, sep="\t")


def test_create_heatmap_plots_the_sorted_matrix_with_nan_as_zero(tmp_path, figure_log):
	f = tmp_path / "chr1_WT.txt"
	_distances(f, ["b", "a"], [0, 100, 200], [[0.5, np.nan, 1.0], [0.0, 0.25, 0.75]])
	vis.create_heatmap(f)
	fig = figure_log.by_name("chr1_WT.svg")
	assert figure_log.paths()[0].parent == tmp_path
	assert list(fig.data[0].y) == ["a", "b"]
	assert list(fig.data[0].x) == [0, 100, 200]
	np.testing.assert_allclose(fig.data[0].z, [[0.0, 0.25, 0.75], [0.5, 0.0, 1.0]])
	assert fig.layout.title.text == "chr1_WT"
	assert (fig.layout.coloraxis.cmin, fig.layout.coloraxis.cmax) == (0, 1)


def test_create_heatmap_orders_rows_by_groups_and_drops_unknown_names(tmp_path, figure_log):
	f = tmp_path / "chr1_WT.txt"
	_distances(f, ["a", "b", "c"], [0, 100], [[0.1, 0.2], [0.3, 0.4], [0.5, 0.6]])
	groups = pd.DataFrame({"group": ["X", "X", "Y", "Y"]}, index=["c", "ghost", "a", "b"])
	vis.create_heatmap(f, groups=groups)
	fig = figure_log.by_name("chr1_WT.svg")
	assert list(fig.data[0].y) == ["c", "a", "b"]
	np.testing.assert_allclose(fig.data[0].z, [[0.5, 0.6], [0.1, 0.2], [0.3, 0.4]])


def test_create_heatmap_runner_visualizes_chr_files_but_not_max_species_files(tmp_path, figure_log):
	for name in ("chr1_WT.txt", "chr2_WT.txt", "chr1_max_species.txt", "notchr.txt"):
		_distances(tmp_path / name, ["a"], [0, 100], [[0.1, 0.2]])
	vis.create_heatmap_runner(tmp_path)
	assert sorted(p.name for p in figure_log.paths()) == ["chr1_WT.svg", "chr2_WT.svg"]


def test_create_heatmap_renders_an_svg_file(tmp_path):
	pytest.importorskip("kaleido")
	f = tmp_path / "chr1_WT.txt"
	_distances(f, ["a", "b"], [0, 100], [[0.1, 0.2], [0.3, 0.4]])
	try:
		vis.create_heatmap(f)
	except Exception as e:  # kaleido needs Chrome, which may be absent
		pytest.skip(f"cannot render images here: {type(e).__name__}")
	assert "<svg" in (tmp_path / "chr1_WT.svg").read_text()[:500]


# --- command line -------------------------------------------------------------------------------------


@pytest.fixture
def run_vis(monkeypatch):
	def run(*args):
		monkeypatch.setattr(sys, "argv", ["visualize", *map(str, args)])
		vis.main()

	return run


def test_main_creates_every_requested_curve_for_every_intro_type(res, run_vis, figure_log, capsys):
	# a second intro type in the scored folder
	for thr, rows in ((0.2, [("chr1", 5, 1, 1, 40)]), (0.4, [("chr1", 3, 0, 3, 40)])):
		_metrics_file(res, thr, "REF", rows)
	run_vis("-v", "prc", "prcc", "mcc", "--dir", res, "--thresholds", 0.2, 0.4)
	names = sorted(p.name for p in figure_log.paths())
	assert names == ["REF_mcc.svg", "REF_prc.svg", "REF_prcc.svg", "WT_mcc.svg", "WT_prc.svg", "WT_prcc.svg"]
	assert "Done." in capsys.readouterr().out


def test_main_prca_and_shtmp(res, run_vis, figure_log, capsys):
	_heatmap_csv(res, 0.2, "chr1", "WT", {"s1": [5, 4]})
	_heatmap_csv(res, 0.4, "chr1", "WT", {"s1": [5, 4]})
	run_vis("-v", "prca", "shtmp", "--dir", res, "--thresholds", 0.2, 0.4)
	assert [p.name for p in figure_log.paths()] == ["WT_prca.svg"]
	out = capsys.readouterr().out
	assert "Expected 18 thresholds" in out
	assert not (res / "WT_scored_heatmaps.pdf").exists()


def test_main_htmp_uses_the_groups_file_for_row_order(tmp_path, run_vis, figure_log):
	d = tmp_path / "dist"
	d.mkdir()
	_distances(d / "chr1_WT.txt", ["a", "b"], [0, 100], [[0.1, 0.2], [0.3, 0.4]])
	grp = tmp_path / "groups.tsv"
	grp.write_text("name\tgroup\nb\tX\na\tY\n")
	run_vis("-v", "htmp", "--dir", d, "--grp", grp)
	assert list(figure_log.by_name("chr1_WT.svg").data[0].y) == ["b", "a"]


def test_main_htmp_without_groups_sorts_rows_by_name(tmp_path, run_vis, figure_log):
	d = tmp_path / "dist"
	d.mkdir()
	_distances(d / "chr1_WT.txt", ["b", "a"], [0, 100], [[0.1, 0.2], [0.3, 0.4]])
	run_vis("-v", "htmp", "--dir", d)
	assert list(figure_log.by_name("chr1_WT.svg").data[0].y) == ["a", "b"]


def test_main_default_thresholds_are_0_1_to_0_95_in_steps_of_0_05(tmp_path, run_vis, monkeypatch):
	root = tmp_path / "def"
	root.mkdir()
	seen = {}
	monkeypatch.setattr(vis, "create_pr_curve", lambda d, t, thr: seen.update(thr=thr, t=t))
	_metrics_file(root, 0.1, "WT", [("chr1", 1, 1, 1, 1)])
	run_vis("-v", "prc", "--dir", root)
	assert seen["t"] == "WT"
	assert seen["thr"] == pytest.approx([round(0.1 + 0.05 * i, 2) for i in range(18)])


def test_main_rejects_unknown_visual_missing_directory_and_missing_scored_folder(tmp_path, res, run_vis):
	with pytest.raises(ValueError, match="Unknown visualization"):
		run_vis("-v", "bogus", "--dir", res)
	with pytest.raises(ValueError, match="Input directory not found"):
		run_vis("-v", "prc", "--dir", tmp_path / "nowhere")
	with pytest.raises(ValueError, match="Scored directory not found"):
		run_vis("-v", "prc", "--dir", res, "--thresholds", 0.9)


def test_main_requires_v_and_dir(res, monkeypatch):
	for argv in (["visualize", "--dir", str(res)], ["visualize", "-v", "prc"]):
		monkeypatch.setattr(sys, "argv", argv)
		with pytest.raises(SystemExit) as e:
			vis.main()
		assert e.value.code == 2

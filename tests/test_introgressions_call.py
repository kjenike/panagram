"""panagram.introgressions.call_introgressions: binning, preprocessing, thresholding, BED output, CLI.

Unit tests use hand-worked tiny tables. The CLI tests run `main()` on the synthetic anchored
index and compare the BED files with calls derived from the brute-force k-mer oracle of the
synthetic genomes (helpers_introgressions.reference_bins), never from panagram's bitmap code.

Chromosome layout used throughout: bins of 500 over k-mer positions, five sampled rows per bin
(step 100). Similarity of alpha chr1 to epsilon per bin (from the oracle) is
[0.4, 0.8, 0.0, 0.0, 0.0, 0.0]; to beta [0.4, 0.8, 0.6, 0.6, 0.6, 0.0].
"""

import sys
import types

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import pytest
from scipy.ndimage import median_filter, uniform_filter1d

import helpers_introgressions as H
from panagram.introgressions import call_introgressions as call

GROUPS = {"alpha": "ING", "beta": "DON", "gamma": "OTH", "delta": "OTH", "epsilon": "REF"}


# ------------------------------------------------------------------ fill_snps


def test_fill_snps_fills_only_gaps_of_exactly_gap_size():
	"""Hand-worked: ones at 0,2,3,6; the one-zero gap is filled for gap_size 1, the two-zero gap for gap_size 2."""
	row = types.SimpleNamespace(values=np.array([1, 0, 1, 1, 0, 0, 1]))
	assert call.fill_snps(row, 1).tolist() == [1, 1, 1, 1, 0, 0, 1]
	row = types.SimpleNamespace(values=np.array([1, 0, 1, 1, 0, 0, 1]))
	assert call.fill_snps(row, 2).tolist() == [1, 0, 1, 1, 1, 1, 1]


def test_fill_snps_on_series_row():
	"""A pandas row works too (pandas 3 hands out a read-only .values, so that case is skipped)."""
	row = pd.Series([1, 0, 1, 1, 0, 0, 1])
	if not row.values.flags.writeable:
		pytest.skip("pandas returns read-only .values (copy-on-write)")
	assert call.fill_snps(row, 1).tolist() == [1, 1, 1, 1, 0, 0, 1]


def test_fill_snps_fills_a_leading_or_trailing_gap_of_exactly_gap_size():
	"""Hand-worked: a run of gap_size zeros at either end of the row is filled; a shorter run at the end is not."""
	row = lambda v: types.SimpleNamespace(values=np.array(v))
	assert call.fill_snps(row([0, 1, 1, 1]), 1).tolist() == [1, 1, 1, 1]
	assert call.fill_snps(row([1, 1, 1, 0]), 1).tolist() == [1, 1, 1, 1]
	assert call.fill_snps(row([0, 0, 1, 1, 0, 0]), 2).tolist() == [1, 1, 1, 1, 1, 1]
	# one leading zero is shorter than gap_size 2
	assert call.fill_snps(row([0, 1, 1, 1]), 2).tolist() == [0, 1, 1, 1]
	assert call.fill_snps(row([1, 1, 1, 0]), 2).tolist() == [1, 1, 1, 0]


# ------------------------------------------------------------------ bitmap_to_bins


def _toy_bitmap():
	# anchor A is 1 everywhere; B is 1 in 2 of the first 4 rows, 0 in the next 4, 1 in the last 4
	idx = list(range(0, 1200, 100))
	return pd.DataFrame(
		{"A": [1] * 12, "B": [1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1]},
		index=idx,
	)


def test_bitmap_to_bins_hand_worked():
	"""Bins of 400 over 12 rows: B has 2/4, 0/4 and 4/4 ones, A is 4/4 everywhere."""
	out = call.bitmap_to_bins(_toy_bitmap(), 400)
	assert list(out.index) == ["A", "B"]
	assert list(out.columns) == [0, 400, 800]
	assert out.loc["A"].tolist() == [1.0, 1.0, 1.0]
	assert out.loc["B"].tolist() == [0.5, 0.0, 1.0]


def test_bitmap_to_bins_omit_fixed_kmers():
	"""Rows that are 1 in every accession are dropped: B becomes 0/2 in bin 0 and 0/4 in bin 400; the all-fixed bin 800 becomes 1."""
	out = call.bitmap_to_bins(_toy_bitmap(), 400, omit_fixed_kmers=True)
	assert out.loc["B"].tolist() == [0.0, 0.0, 1.0]
	assert out.loc["A"].tolist() == [1.0, 1.0, 1.0]


def test_bitmap_to_bins_uneven_last_bin_is_divided_by_its_own_maximum():
	"""Ten rows, bins of 400: the last bin has 2 rows; B has 1 of them."""
	bm = pd.DataFrame({"A": [1] * 10, "B": [1, 1, 1, 1, 0, 0, 0, 0, 1, 0]}, index=range(0, 1000, 100))
	out = call.bitmap_to_bins(bm, 400)
	assert out.loc["B"].tolist() == [1.0, 0.0, 0.5]


def test_bitmap_to_bins_omit_unique_kmers_gives_ref_and_outgroup_the_anchor_only_kmers():
	"""Rows with no 1 in {R, O} get 1 in both first. Rows (1,0,0), (1,1,0), (1,0,1) -> R and O are 2 of 3 each."""
	bm = pd.DataFrame({"A": [1, 1, 1], "R": [0, 1, 0], "O": [0, 0, 1]}, index=[0, 100, 200])
	outgroup = ["O"]
	out = call.bitmap_to_bins(bm, 300, False, True, "R", outgroup)
	assert out.loc["R", 0] == pytest.approx(2 / 3)
	assert out.loc["O", 0] == pytest.approx(2 / 3)
	assert outgroup == ["O"]
	plain = call.bitmap_to_bins(bm, 300)
	assert plain.loc["R", 0] == pytest.approx(1 / 3)


def test_bitmap_to_bins_omit_unique_then_fixed():
	"""The anchor-only row becomes (1,1,1) and is then removed as fixed: R and O are 1 of 2 each."""
	bm = pd.DataFrame({"A": [1, 1, 1], "R": [0, 1, 0], "O": [0, 0, 1]}, index=[0, 100, 200])
	out = call.bitmap_to_bins(bm, 300, True, True, "R", ["O"])
	assert out.loc["R", 0] == pytest.approx(0.5)
	assert out.loc["O", 0] == pytest.approx(0.5)


def test_bitmap_to_bins_matches_the_reference_on_synthetic_bitmaps(synthetic_genomes):
	"""Every genome and chromosome of the synthetic set, with and without omitting fixed k-mers."""
	for genome, chrom in [("alpha", "chr1"), ("alpha", "chr2"), ("gamma", "chr1"), ("beta", "chr1")]:
		bm = H.oracle_bitmap(synthetic_genomes, genome, chrom)
		for rmf in (False, True):
			got = call.bitmap_to_bins(bm, H.BIN, rmf)
			want = H.reference_bins(bm, H.BIN, rmf=rmf)
			np.testing.assert_allclose(got.values, want.values, err_msg=f"{genome} {chrom} {rmf}")
			assert list(got.columns) == list(want.columns)


# ------------------------------------------------------------------ row_trimmed_mean


def test_row_trimmed_mean_drops_the_outlier():
	"""Twenty 0s and one 10: mean 10/21, sample std 2.18; with 3 std the 10 lies outside (7.0) so the trimmed mean is 0."""
	row = pd.Series([0.0] * 20 + [10.0])
	assert call.row_trimmed_mean(row, 3) == 0.0
	assert call.row_trimmed_mean(row, -1) == pytest.approx(10 / 21)
	assert call.row_trimmed_mean(row, 100) == pytest.approx(10 / 21)


def test_row_trimmed_mean_uses_the_sample_standard_deviation():
	"""[0, 0, 0, 0, 5]: mean 1, squared deviations sum to 20, so the sample std is sqrt(20/4) = 2.236 and 1.9 std = 4.25.
	The 5 lies 4 from the mean and stays (mean 1.0); a population std (2.0, 1.9 std = 3.8) would drop it (mean 0.0)."""
	row = pd.Series([0.0, 0.0, 0.0, 0.0, 5.0])
	assert call.row_trimmed_mean(row, 1.9) == pytest.approx(1.0)
	assert call.row_trimmed_mean(row, 1.7) == pytest.approx(0.0)  # 1.7 * 2.236 = 3.80 < 4: the 5 is dropped


# ------------------------------------------------------------------ get_genome_similarities


def _expected_similarity(genomes, anchor, trim, rmf):
	cols = []
	for chrom in genomes.data.seqs[anchor]:
		cols.append(H.reference_bins(H.oracle_bitmap(genomes, anchor, chrom), H.BIN, rmf=rmf))
	allbins = pd.concat(cols, axis=1)
	out = {}
	for name, row in allbins.iterrows():
		v = row.values.astype(float)
		if trim != -1:
			m, s = v.mean(), v.std(ddof=1)
			v = v[(v >= m - trim * s) & (v <= m + trim * s)]
		out[name] = v.mean()
	return pd.Series(out)


@pytest.mark.parametrize("trim,rmf", [(-1, False), (1.0, False), (3.0, True)])
def test_get_genome_similarities_is_mean_over_all_chromosomes(
	anchored_reader, synthetic_genomes, trim, rmf
):
	"""Per accession, (trimmed) mean of the binned similarities over alpha's bins on chr1 and chr2."""
	genome = anchored_reader.genomes["alpha"]
	got = call.get_genome_similarities(genome, H.STEP, H.BIN, rmf, None, None, None, trim)
	want = _expected_similarity(synthetic_genomes, "alpha", trim, rmf)
	assert list(got.index) == list(want.index)
	np.testing.assert_allclose(got.values, want.values)


# ------------------------------------------------------------------ smooth_row


def test_smooth_row_mean_and_median_hand_worked():
	"""Size-3 filters with reflected edges: mean of [1,2,3,4,5] is [4/3, 2, 3, 4, 14/3]; median of [0,5,0,0,0] is all 0."""
	row = pd.Series([1.0, 2, 3, 4, 5], index=[0, 10, 20, 30, 40])
	out = call.smooth_row(row, "mean", 3)
	np.testing.assert_allclose(out.values, [4 / 3, 2, 3, 4, 14 / 3])
	assert list(out.index) == [0, 10, 20, 30, 40]
	spike = pd.Series([0.0, 5, 0, 0, 0])
	assert call.smooth_row(spike, "median", 3).tolist() == [0, 0, 0, 0, 0]
	assert call.smooth_row(pd.Series([1.0, 9, 2, 8, 3]), "median", 3).tolist() == [1, 2, 8, 3, 3]


# ------------------------------------------------------------------ edge_tapered_row_normalization


def test_edge_tapered_row_normalization_hand_worked():
	"""Three columns: boost = 0.1*[exp(-4), 1, exp(-4)]; value*(1+boost), clip to 1, minus 0.2 unless it is 1, clip to [0, 1]."""
	df = pd.DataFrame(
		[[0.5, 0.5, 0.5], [1.0, 1.0, 1.0], [0.95, 0.95, 0.95], [0.1, 0.1, 0.1]],
		index=["a", "b", "c", "d"],
	)
	out = call.edge_tapered_row_normalization(df)
	e = 0.1 * np.exp(-4)
	np.testing.assert_allclose(out.loc["a"], [0.5 * (1 + e) - 0.2, 0.5 * 1.1 - 0.2, 0.5 * (1 + e) - 0.2])
	assert out.loc["b"].tolist() == [1.0, 1.0, 1.0]
	# the centre reaches 1 after boosting and stays 1; the edges are lowered by 0.2
	np.testing.assert_allclose(out.loc["c"], [0.95 * (1 + e) - 0.2, 1.0, 0.95 * (1 + e) - 0.2])
	assert out.loc["d"].tolist() == [0.0, 0.0, 0.0]  # 0.1*1.0x - 0.2 < 0 clipped
	assert df.loc["a"].tolist() == [0.5, 0.5, 0.5]  # input untouched


def test_edge_tapered_row_normalization_intensity():
	"""A larger intensity raises the centre value: 0.5*(1+0.3)-0.2 = 0.45."""
	df = pd.DataFrame([[0.5, 0.5, 0.5]])
	out = call.edge_tapered_row_normalization(df, intensity=0.3)
	assert out.iloc[0, 1] == pytest.approx(0.45)


# ------------------------------------------------------------------ preprocess_binned_bitmap


def _prep_input():
	return pd.DataFrame(
		[[0.5, 0.99, 1.0], [0.3, 0.6, 1.0], [1.0, 1.0, 1.0]],
		index=["A", "B", "ANC"],
		columns=[0, 100, 200],
	)


def _prep(df, sims=None, gnm=None, sft=None, ssz=5, edg=False):
	return call.preprocess_binned_bitmap(df, sims, gnm, sft, ssz, edg)


def test_preprocess_rounds_and_leaves_input_alone():
	df = pd.DataFrame([[0.3333, 0.6666]], index=["A"], columns=[0, 100])
	out = _prep(df)
	assert out.loc["A"].tolist() == [0.33, 0.67]
	assert df.loc["A"].tolist() == [0.3333, 0.6666]


def test_preprocess_normalizes_to_given_mean_only_below_098():
	"""Shift = target - similarity per accession: A +0.05 (0.5 -> 0.55; 0.99 and 1.0 stay), B +0.15."""
	sims = pd.Series({"A": 0.9, "B": 0.8, "ANC": 1.0})
	out = _prep(_prep_input(), sims, gnm=0.95)
	np.testing.assert_allclose(out.loc["A"], [0.55, 0.99, 1.0])
	np.testing.assert_allclose(out.loc["B"], [0.45, 0.75, 1.0])
	assert out.loc["ANC"].tolist() == [1.0, 1.0, 1.0]


def test_preprocess_normalization_includes_bins_at_exactly_098():
	"""A bin at 0.98 moves with the shift (0.98 - 0.05 = 0.93); 0.99 does not."""
	df = pd.DataFrame([[0.98, 0.99]], index=["A"], columns=[0, 100])
	out = _prep(df, pd.Series({"A": 0.95}), gnm=0.9)
	np.testing.assert_allclose(out.loc["A"], [0.93, 0.99])


def test_preprocess_gnm_minus_one_uses_largest_mean_below_one():
	"""The target is max of the similarities that are not 1: 0.9 here, so A moves by 0 and B by +0.1."""
	sims = pd.Series({"A": 0.9, "B": 0.8, "ANC": 1.0})
	out = _prep(_prep_input(), sims, gnm=-1)
	np.testing.assert_allclose(out.loc["A"], [0.5, 0.99, 1.0])
	np.testing.assert_allclose(out.loc["B"], [0.4, 0.7, 1.0])


def test_preprocess_normalization_is_clipped_to_unit_interval():
	df = pd.DataFrame([[0.9, 0.1]], index=["A"], columns=[0, 100])
	up = _prep(df, pd.Series({"A": 0.2}), gnm=0.9)  # +0.7
	np.testing.assert_allclose(up.loc["A"], [1.0, 0.8])
	down = _prep(df, pd.Series({"A": 0.9}), gnm=0.5)  # -0.4
	np.testing.assert_allclose(down.loc["A"], [0.5, 0.0])


def test_preprocess_smoothing_and_edge_order():
	"""Edge normalization runs before smoothing: [0.5]*3 -> [0.3009, 0.35, 0.3009] -> size-3 mean."""
	df = pd.DataFrame([[0.5, 0.5, 0.5]], index=["A"], columns=[0, 100, 200])
	out = _prep(df, sft="mean", ssz=3, edg=True)
	e = 0.5 * (1 + 0.1 * np.exp(-4)) - 0.2
	c = 0.5 * 1.1 - 0.2
	want = [(2 * e + c) / 3, (2 * e + c) / 3, (2 * e + c) / 3]
	np.testing.assert_allclose(out.loc["A"], want)
	assert list(out.columns) == [0, 100, 200]
	only_smooth = _prep(pd.DataFrame([[0.0, 0.9, 0.0]], index=["A"], columns=[0, 1, 2]), sft="median", ssz=3)
	assert only_smooth.loc["A"].tolist() == [0.0, 0.0, 0.0]


# ------------------------------------------------------------------ threshold_introgressions


def _groups_table():
	return pd.DataFrame(
		{
			0: [0.5, 0.9, 0.9, 0.95, 0.0],
			1: [0.5, 0.2, 0.5, 0.9, 0.0],
			2: [0.5, 0.4, 0.97, 0.99, 0.0],
			"group": ["REC", "REC", "REF", "DON", "DON"],
		},
		index=["A", "A2", "R", "D", "D2"],
	)


def test_threshold_two_way_flags_bins_below_threshold_in_best_ref_member():
	"""Hand-worked: comp_sim is the REF row [0.9, 0.5, 0.97]; below 0.6 only the middle bin."""
	out = call.threshold_introgressions(_groups_table(), "A", "REF", 0.6)
	assert out["introgression"].tolist() == [0, 1, 0]
	np.testing.assert_allclose(out["comp_sim"], [0.9, 0.5, 0.97])
	# the anchor's own group without the anchor: only A2
	np.testing.assert_allclose(out["anchor_sim"], [0.9, 0.2, 0.4])
	assert list(out.index) == [0, 1, 2]


def test_threshold_two_way_equal_to_threshold_is_not_introgressed():
	tbl = pd.DataFrame({0: [0.5, 0.5], "group": ["REC", "REF"]}, index=["A", "R"])
	assert call.threshold_introgressions(tbl, "A", "REF", 0.5)["introgression"].tolist() == [0]
	assert call.threshold_introgressions(tbl, "A", "REF", 0.51)["introgression"].tolist() == [1]


def test_threshold_two_way_takes_max_over_ref_group_members():
	tbl = pd.DataFrame(
		{0: [0.5, 0.2, 0.9], 1: [0.5, 0.2, 0.1], "group": ["REC", "REF", "REF"]},
		index=["A", "R1", "R2"],
	)
	out = call.threshold_introgressions(tbl, "A", "REF", 0.5)
	assert out["comp_sim"].tolist() == [0.9, 0.2]
	assert out["introgression"].tolist() == [0, 1]


def test_threshold_three_way_hand_worked():
	"""comp = max(D, D2) = [0.95, 0.9, 0.99]; ref = [0.9, 0.5, 0.97]. Called when ref < 0.95 and comp >= ref + thr."""
	tbl = _groups_table()
	loose = call.threshold_introgressions(tbl, "A", "DON", 0.01)
	np.testing.assert_allclose(loose["ref_sim"], [0.9, 0.5, 0.97])
	np.testing.assert_allclose(loose["comp_sim"], [0.95, 0.9, 0.99])
	# bin 2: comp - ref = 0.02 but ref >= 0.95, so it is not called
	assert loose["introgression"].tolist() == [1, 1, 0]
	strict = call.threshold_introgressions(tbl, "A", "DON", 0.2)
	assert strict["introgression"].tolist() == [0, 1, 0]


def test_threshold_three_way_difference_equal_to_threshold_is_called():
	"""Exact binary fractions: bin 0 has comp 0.75 = ref 0.5 + 0.25 (called), bin 1 has comp 0.74 (not called)."""
	tbl = pd.DataFrame(
		{0: [0.0, 0.5, 0.75], 1: [0.0, 0.5, 0.74], "group": ["REC", "REF", "DON"]},
		index=["A", "R", "D"],
	)
	out = call.threshold_introgressions(tbl, "A", "DON", 0.25)
	assert out["introgression"].tolist() == [1, 0]


def test_threshold_introgressions_simple_uses_only_the_anchor_row():
	tbl = pd.DataFrame(
		{0: [0.5, 0.1], 1: [0.4, 0.9], 2: [0.6, 0.9], "group": ["REF", "REC"]}, index=["R", "A"]
	)
	out = call.threshold_introgressions_simple(tbl, "A", 0.5)
	assert out["introgression"].tolist() == [1, 0, 0]
	np.testing.assert_allclose(out["anchor_sim"], [0.1, 0.9, 0.9])
	assert out["comp_sim"].isna().all()
	# exactly at the threshold: not introgressed
	assert call.threshold_introgressions_simple(tbl, "R", 0.5)["introgression"].tolist() == [0, 1, 0]


# ------------------------------------------------------------------ bins_to_bed


def _bins(values, bin_size=500):
	return pd.DataFrame({"introgression": values}, index=[i * bin_size for i in range(len(values))])


def test_bins_to_bed_merges_adjacent_bins_with_inclusive_ends():
	out = call.bins_to_bed(_bins([0, 1, 1, 0, 1, 0]), 500, "chr1", "G")
	assert list(out.columns) == ["chr", "start", "end", "name"]
	assert [tuple(r) for r in out.values.tolist()] == [
		("chr1", 500, 1499, "G_intro"),
		("chr1", 2000, 2499, "G_intro"),
	]


def test_bins_to_bed_run_at_the_start_and_whole_chromosome():
	out = call.bins_to_bed(_bins([1, 1, 0], 100), 100, "c", "REFA")
	assert out.values.tolist() == [["c", 0, 199, "REFA_intro"]]
	out = call.bins_to_bed(_bins([1, 1, 1, 1], 100), 100, "c", "X")
	assert out.values.tolist() == [["c", 0, 399, "X_intro"]]


def test_bins_to_bed_counts_a_merged_value_above_one_as_introgressed():
	"""Merged calls of several groups add up to values of 2; any value > 0 is an introgression."""
	out = call.bins_to_bed(_bins([2, 1, 0, 2], 10), 10, "c", "merged")
	assert out[["start", "end"]].values.tolist() == [[0, 19], [30, 39]]


def test_bins_to_bed_no_introgressions_gives_an_empty_frame():
	out = call.bins_to_bed(_bins([0, 0, 0], 100), 100, "c", "G")
	assert out.empty
	assert list(out.columns) == ["chr", "start", "end", "name"]


def test_bins_to_bed_separated_by_one_bin_stay_separate():
	out = call.bins_to_bed(_bins([1, 0, 1], 100), 100, "c", "G")
	assert out[["start", "end"]].values.tolist() == [[0, 99], [200, 299]]


# ------------------------------------------------------------------ visualize


@pytest.fixture
def figures(monkeypatch):
	"""Capture (path, figure) instead of rendering; kaleido and Chrome are not needed."""
	seen = []

	def fake(self, file, *args, **kwargs):
		seen.append((str(file), self))

	monkeypatch.setattr(go.Figure, "write_image", fake)
	return seen


def test_visualize_orders_rows_by_groups_and_keeps_introgressions_row_last(figures, tmp_path):
	bm = pd.DataFrame(
		[[0.1, 0.2], [0.3, 0.4], [0.5, 0.6], [1.0, 0.0]],
		index=["b", "a", "c", "Introgressions"],
		columns=[0, 500],
	)
	groups = pd.DataFrame({"group": ["x", "x", "x"]}, index=["c", "a", "b"])
	call.visualize(bm, tmp_path / "v.svg", inverse=True, title="T", groups=groups)
	path, fig = figures[0]
	assert path == str(tmp_path / "v.svg")
	assert list(fig.data[0].y) == ["c", "a", "b", "Introgressions"]
	np.testing.assert_allclose(
		np.asarray(fig.data[0].z), [[0.5, 0.6], [0.3, 0.4], [0.1, 0.2], [1.0, 0.0]]
	)
	assert fig.layout.title.text == "T"
	assert fig.layout.xaxis.dtick == 2000000
	assert fig.layout.coloraxis.cmin == 0 and fig.layout.coloraxis.cmax == 1


def test_visualize_without_groups_keeps_order_and_scales_figure(figures, tmp_path):
	"""Without groups the rows stay as given; 40 rows -> height 12*40+200 = 680; last bin at 200 Mb -> width 1400."""
	idx = [f"r{i}" for i in range(40)]
	bm = pd.DataFrame(np.zeros((40, 2)), index=idx, columns=[0, 200_000_000])
	call.visualize(bm, tmp_path / "w.svg", inverse=False, title="plain", xaxis_dtick=5)
	_, fig = figures[0]
	assert list(fig.data[0].y) == idx
	assert fig.layout.height == 680
	assert fig.layout.width == 1400
	assert fig.layout.xaxis.dtick == 5
	small = pd.DataFrame(np.zeros((3, 2)), index=["a", "b", "c"], columns=[0, 500])
	call.visualize(small, tmp_path / "s.svg")
	assert (figures[1][1].layout.height, figures[1][1].layout.width) == (500, 500)


# ------------------------------------------------------------------ main(): fixtures and runners


@pytest.fixture(scope="module")
def groups_tsv(tmp_path_factory):
	return H.write_groups(tmp_path_factory.mktemp("grp") / "group.tsv", GROUPS)


@pytest.fixture
def run_main(monkeypatch, anchored_index, groups_tsv, tmp_path):
	"""Run call_introgressions.main() in this process. Returns the output directory.

	run_main(*args, inline=True, tsv=None, idx=None): --idx, --tsv and --out are added.
	"""

	def run(*args, inline=True, tsv=None, idx=None, out=None):
		out = out or (tmp_path / "out")
		args = [str(a) for a in args]
		if "--bin" not in args:
			args += ["--bin", str(H.BIN)]
		argv = [
			"call_introgressions",
			"--idx", str(idx or anchored_index),
			"--tsv", str(tsv or groups_tsv),
			"--out", str(out),
			*[str(a) for a in args],
		]
		monkeypatch.setattr(sys, "argv", argv)
		if inline:
			monkeypatch.setattr(call, "ProcessPoolExecutor", H.InlineExecutor)
		call.main()
		return out

	return run


def _bed(out, thr, anchor, chrom, group):
	return out / f"{out.name}_{thr}" / "raw" / f"{anchor}_{chrom}_{group}.bed"


def _expected_two_way(genomes, anchor, chrom, ref, thr, rmf=False, bin_size=H.BIN, transform=None, sims_of=None):
	sims = H.reference_bins(H.oracle_bitmap(genomes, sims_of or anchor, chrom), bin_size, rmf=rmf)
	if transform is not None:
		sims = transform(sims)
	name = "REFA" if sims_of is None else "REF"
	calls = [int(c) for c in sims.columns if sims.loc[anchor if sims_of else ref, c] < thr]
	return H.runs_to_bed(calls, bin_size, chrom, f"{name}_intro")


# ------------------------------------------------------------------ main(): two-way


def test_main_two_way_accession_space_matches_oracle_calls(run_main, synthetic_genomes):
	"""--grp picks anchors (alpha, gamma, delta of ING and OTH), --chr is absent so every chromosome is run; two thresholds."""
	out = run_main("--grp", "ING", "OTH", "--cmp", "REF", "--thr", "0.3", "0.7")
	seen = 0
	for anchor in ("alpha", "gamma", "delta"):
		for chrom in synthetic_genomes.data.seqs[anchor]:
			for thr in (0.3, 0.7):
				got = H.read_bed_rows(_bed(out, thr, anchor, chrom, "REFA"))
				assert got == _expected_two_way(synthetic_genomes, anchor, chrom, "epsilon", thr)
				seen += 1
	assert seen == 10
	# the one worked example from the module docstring: alpha chr1 vs epsilon [.4,.8,0,0,0,0]
	assert H.read_bed_rows(_bed(out, 0.3, "alpha", "chr1", "REFA")) == [
		("chr1", 1000, 2999, "REFA_intro")
	]
	assert H.read_bed_rows(_bed(out, 0.7, "alpha", "chr1", "REFA")) == [
		("chr1", 0, 499, "REFA_intro"),
		("chr1", 1000, 2999, "REFA_intro"),
	]


def test_main_anc_and_chr_restrict_what_is_written(run_main):
	out = run_main("--anc", "beta", "--chr", "chr1", "--cmp", "REF", "--thr", "0.3")
	files = sorted(p.name for p in (out / f"{out.name}_0.3" / "raw").iterdir())
	assert files == ["beta_chr1_REFA.bed"]


def test_main_process_pool_and_threads_option(run_main, synthetic_genomes):
	"""Same calls through a real two-worker process pool (--threads 2)."""
	out = run_main("--anc", "alpha", "--cmp", "REF", "--thr", "0.3", "--threads", "2", inline=False)
	for chrom in ("chr1", "chr2"):
		assert H.read_bed_rows(_bed(out, 0.3, "alpha", chrom, "REFA")) == _expected_two_way(
			synthetic_genomes, "alpha", chrom, "epsilon", 0.3
		)


def test_main_rmf_changes_values_and_matches_reference(run_main, synthetic_genomes):
	"""Without fixed k-mers alpha-vs-epsilon on chr1 is 0 in every bin, so the whole chromosome is called."""
	plain = _expected_two_way(synthetic_genomes, "alpha", "chr1", "epsilon", 0.3)
	withrmf = _expected_two_way(synthetic_genomes, "alpha", "chr1", "epsilon", 0.3, rmf=True)
	assert plain != withrmf
	out = run_main("--anc", "alpha", "--chr", "chr1", "--cmp", "REF", "--thr", "0.3", "--rmf")
	assert H.read_bed_rows(_bed(out, 0.3, "alpha", "chr1", "REFA")) == withrmf
	assert withrmf == [("chr1", 0, 2999, "REFA_intro")]


def test_main_bin_size_option(run_main, synthetic_genomes):
	"""--bin 1000: three bins; alpha vs epsilon is (0.4+0.8... computed from the oracle)."""
	out = run_main("--anc", "alpha", "--chr", "chr1", "--cmp", "REF", "--thr", "0.5", "--bin", "1000")
	want = _expected_two_way(synthetic_genomes, "alpha", "chr1", "epsilon", 0.5, bin_size=1000)
	assert H.read_bed_rows(_bed(out, 0.5, "alpha", "chr1", "REFA")) == want
	assert want != _expected_two_way(synthetic_genomes, "alpha", "chr1", "epsilon", 0.5)


def test_main_stp_option_samples_fewer_rows(run_main, synthetic_genomes):
	"""--stp 50 doubles the sampled rows; similarity is the fraction of sampled rows with a 1 in each bin."""
	out = run_main("--anc", "alpha", "--chr", "chr1", "--cmp", "REF", "--thr", "0.3", "--stp", "50")
	bm = H.oracle_bitmap(synthetic_genomes, "alpha", "chr1", step=50)
	sims = H.reference_bins(bm, H.BIN)
	want = H.runs_to_bed(H.reference_calls_2way(sims, "alpha", "epsilon", 0.3), H.BIN, "chr1", "REFA_intro")
	assert H.read_bed_rows(_bed(out, 0.3, "alpha", "chr1", "REFA")) == want


def test_main_gnm_and_trm_shift_each_genome_to_the_target_mean(run_main, synthetic_genomes):
	"""--gnm 0.7: each accession's rows move by (0.7 - its trimmed mean over both chromosomes), only where <= 0.98."""
	trim = 1.0
	sims = _expected_similarity(synthetic_genomes, "alpha", trim, False)

	def shift(bins):
		bins = bins.round(2)
		for name in bins.index:
			row = bins.loc[name]
			delta = 0.7 - sims[name]
			bins.loc[name] = np.where(row <= 0.98, row + delta, row).clip(0, 1)
		return bins

	want = _expected_two_way(synthetic_genomes, "alpha", "chr1", "epsilon", 0.3, transform=shift)
	assert want != _expected_two_way(synthetic_genomes, "alpha", "chr1", "epsilon", 0.3)
	out = run_main(
		"--anc", "alpha", "--chr", "chr1", "--cmp", "REF", "--thr", "0.3",
		"--gnm", "0.7", "--trm", trim,
	)
	assert H.read_bed_rows(_bed(out, 0.3, "alpha", "chr1", "REFA")) == want


def test_main_gnm_minus_one_uses_largest_non_one_mean(run_main, synthetic_genomes):
	sims = _expected_similarity(synthetic_genomes, "alpha", 3.0, False)
	target = sims[sims != 1].max()

	def shift(bins):
		bins = bins.round(2)
		for name in bins.index:
			row = bins.loc[name]
			bins.loc[name] = np.where(row <= 0.98, row + (target - sims[name]), row).clip(0, 1)
		return bins

	for thr in (0.3, 0.7):
		want = _expected_two_way(synthetic_genomes, "alpha", "chr1", "epsilon", thr, transform=shift)
		out = run_main(
			"--anc", "alpha", "--chr", "chr1", "--cmp", "REF", "--thr", thr, "--gnm", "-1",
		)
		assert H.read_bed_rows(_bed(out, thr, "alpha", "chr1", "REFA")) == want


@pytest.mark.parametrize("kind", ["mean", "median"])
def test_main_smoothing_filter_and_size(run_main, synthetic_genomes, kind):
	"""--sft with --ssz 3 smooths every row before thresholding; compared with scipy filters on the oracle bins."""
	filt = uniform_filter1d if kind == "mean" else median_filter

	def smooth(bins):
		bins = bins.round(2)
		return pd.DataFrame(
			[filt(r.values.astype(float), size=3) for _, r in bins.iterrows()],
			index=bins.index,
			columns=bins.columns,
		)

	thr = 0.45
	want = _expected_two_way(synthetic_genomes, "alpha", "chr1", "epsilon", thr, transform=smooth)
	assert want != _expected_two_way(synthetic_genomes, "alpha", "chr1", "epsilon", thr)
	out = run_main(
		"--anc", "alpha", "--chr", "chr1", "--cmp", "REF", "--thr", thr, "--sft", kind, "--ssz", 3
	)
	assert H.read_bed_rows(_bed(out, thr, "alpha", "chr1", "REFA")) == want


def test_main_edge_tapered_normalization(run_main, synthetic_genomes):
	"""--edg: value*(1+0.1*window) clipped to 1, minus 0.2 unless 1. Bin 0 of alpha vs epsilon goes from 0.4 to about 0.2."""

	def edge(bins):
		bins = bins.round(2)
		x = np.linspace(-1, 1, bins.shape[1])
		boost = 0.1 * np.exp(-4 * x**2) / np.exp(-4 * x**2).max()
		v = np.clip(bins.values * (1 + boost), 0, 1)
		v = np.where(v == 1, v, v - 0.2).clip(0, 1)
		return pd.DataFrame(v, index=bins.index, columns=bins.columns)

	want = _expected_two_way(synthetic_genomes, "alpha", "chr1", "epsilon", 0.3, transform=edge)
	assert want == [("chr1", 0, 499, "REFA_intro"), ("chr1", 1000, 2999, "REFA_intro")]
	assert want != _expected_two_way(synthetic_genomes, "alpha", "chr1", "epsilon", 0.3)
	out = run_main("--anc", "alpha", "--chr", "chr1", "--cmp", "REF", "--thr", "0.3", "--edg")
	assert H.read_bed_rows(_bed(out, 0.3, "alpha", "chr1", "REFA")) == want


def test_main_skips_anchor_inside_the_comparison_group(run_main):
	out = run_main("--anc", "epsilon", "--cmp", "REF", "--thr", "0.5")
	assert not any(out.rglob("*.bed")) if out.exists() else True


# ------------------------------------------------------------------ main(): reference view (--urf)


def test_main_urf_calls_in_reference_coordinates(run_main, synthetic_genomes):
	"""--urf --ref epsilon: bins of epsilon's chr1, alpha's row; written as alpha_chr1_REF.bed."""
	out = run_main(
		"--anc", "alpha", "--chr", "chr1", "--cmp", "REF", "--thr", "0.3", "0.7",
		"--urf", "--ref", "epsilon",
	)
	sims = H.reference_bins(H.oracle_bitmap(synthetic_genomes, "epsilon", "chr1"), H.BIN)
	for thr in (0.3, 0.7):
		calls = [int(c) for c in sims.columns if sims.loc["alpha", c] < thr]
		want = H.runs_to_bed(calls, H.BIN, "chr1", "REF_intro")
		assert H.read_bed_rows(_bed(out, thr, "alpha", "chr1", "REF")) == want
	# epsilon chr1 is 2400 positions long: five bins (the test is not vacuous)
	assert sims.shape[1] == 5
	assert (sims.loc["alpha"] < 0.7).any()


def test_main_urf_with_gnm(run_main, synthetic_genomes):
	"""In the reference view the target mean comes from the reference genome's own similarities."""
	sims_all = _expected_similarity(synthetic_genomes, "epsilon", 3.0, False)
	target = sims_all[sims_all != 1].max()
	bins = H.reference_bins(H.oracle_bitmap(synthetic_genomes, "epsilon", "chr1"), H.BIN).round(2)
	row = bins.loc["alpha"]
	shifted = np.where(row <= 0.98, row + (target - sims_all["alpha"]), row).clip(0, 1)
	thr = 0.4
	want = H.runs_to_bed(
		[int(c) for c, v in zip(bins.columns, shifted) if v < thr], H.BIN, "chr1", "REF_intro"
	)
	out = run_main(
		"--anc", "alpha", "--chr", "chr1", "--cmp", "REF", "--thr", thr,
		"--urf", "--ref", "epsilon", "--gnm", "-1",
	)
	assert H.read_bed_rows(_bed(out, thr, "alpha", "chr1", "REF")) == want


# ------------------------------------------------------------------ main(): three-way and merged


def test_main_three_way_single_group(run_main, synthetic_genomes):
	"""alpha vs REF epsilon and donor beta: called where epsilon < 0.95 and beta >= epsilon + 0.3 (bins 1000-2499)."""
	out = run_main("--anc", "alpha", "--chr", "chr1", "--cmp", "DON", "--thr", "0.3")
	bm = lambda g: H.reference_bins(H.oracle_bitmap(synthetic_genomes, "alpha", "chr1"), H.BIN).loc[g]
	calls = [int(c) for c in bm("beta").index if bm("epsilon")[c] < 0.95 and bm("beta")[c] >= bm("epsilon")[c] + 0.3]
	assert calls == [1000, 1500, 2000]
	assert H.read_bed_rows(_bed(out, 0.3, "alpha", "chr1", "DON")) == [("chr1", 1000, 2499, "DON_intro")]
	assert not _bed(out, 0.3, "alpha", "chr1", "merged").exists()


def test_main_three_way_merged_of_two_groups(run_main, synthetic_genomes):
	"""DON = beta and OTH = max(gamma, delta). Each group alone, plus a merged file holding their union."""
	bins = H.reference_bins(H.oracle_bitmap(synthetic_genomes, "alpha", "chr1"), H.BIN)
	ref = bins.loc["epsilon"]

	def called(rows):
		comp = bins.loc[rows].max(axis=0)
		return [int(c) for c in bins.columns if ref[c] < 0.95 and comp[c] >= ref[c] + 0.3]

	don, oth = called(["beta"]), called(["gamma", "delta"])
	assert don == [1000, 1500, 2000] and oth == [1500, 2000]
	out = run_main("--anc", "alpha", "--chr", "chr1", "--cmp", "DON", "OTH", "--thr", "0.3")
	got = lambda g: H.read_bed_rows(_bed(out, 0.3, "alpha", "chr1", g))
	assert got("DON") == H.runs_to_bed(don, H.BIN, "chr1", "DON_intro")
	assert got("OTH") == H.runs_to_bed(oth, H.BIN, "chr1", "OTH_intro")
	assert got("merged") == H.runs_to_bed(sorted(set(don) | set(oth)), H.BIN, "chr1", "merged_intro")


def test_main_duplicate_cmp_groups_are_collapsed(run_main):
	out = run_main("--anc", "alpha", "--chr", "chr1", "--cmp", "DON", "DON", "--thr", "0.3")
	files = sorted(p.name for p in (out / f"{out.name}_0.3" / "raw").iterdir())
	assert files == ["alpha_chr1_DON.bed"]


# ------------------------------------------------------------------ main(): rmu / ogrp


def test_main_rmu_outgroup_view_is_written_as_refa(run_main, synthetic_genomes):
	"""--rmu alpha --ogrp DON: k-mers absent from epsilon (REF) and beta (outgroup) count as shared with both. --urf is overridden for alpha."""
	out = run_main(
		"--anc", "alpha", "--chr", "chr1", "--cmp", "REF", "--thr", "0.3", "0.7",
		"--urf", "--ref", "epsilon", "--rmu", "alpha", "--ogrp", "DON",
	)
	bm = H.oracle_bitmap(synthetic_genomes, "alpha", "chr1")
	sims = H.reference_bins(bm, H.BIN, rmu_keep=["beta", "epsilon"])
	plain = H.reference_bins(bm, H.BIN)
	assert not sims.equals(plain)
	for thr in (0.3, 0.7):
		want = H.runs_to_bed(H.reference_calls_2way(sims, "alpha", "epsilon", thr), H.BIN, "chr1", "REFA_intro")
		assert H.read_bed_rows(_bed(out, thr, "alpha", "chr1", "REFA")) == want
		assert not _bed(out, thr, "alpha", "chr1", "REF").exists()


def test_main_rmu_true_applies_to_every_anchor(run_main, synthetic_genomes):
	out = run_main(
		"--anc", "alpha", "gamma", "--chr", "chr1", "--cmp", "REF", "--thr", "0.5",
		"--ref", "epsilon", "--rmu", "true", "--ogrp", "DON",
	)
	for anchor in ("alpha", "gamma"):
		bm = H.oracle_bitmap(synthetic_genomes, anchor, "chr1")
		sims = H.reference_bins(bm, H.BIN, rmu_keep=["beta", "epsilon"])
		want = H.runs_to_bed(H.reference_calls_2way(sims, anchor, "epsilon", 0.5), H.BIN, "chr1", "REFA_intro")
		assert H.read_bed_rows(_bed(out, 0.5, anchor, "chr1", "REFA")) == want


def test_main_rmu_combined_with_rmf(run_main, synthetic_genomes):
	out = run_main(
		"--anc", "alpha", "--chr", "chr1", "--cmp", "REF", "--thr", "0.5",
		"--ref", "epsilon", "--rmu", "alpha", "--ogrp", "DON", "--rmf",
	)
	bm = H.oracle_bitmap(synthetic_genomes, "alpha", "chr1")
	sims = H.reference_bins(bm, H.BIN, rmf=True, rmu_keep=["beta", "epsilon"])
	want = H.runs_to_bed(H.reference_calls_2way(sims, "alpha", "epsilon", 0.5), H.BIN, "chr1", "REFA_intro")
	assert H.read_bed_rows(_bed(out, 0.5, "alpha", "chr1", "REFA")) == want


# ------------------------------------------------------------------ main(): --vis


def test_main_vis_writes_heatmaps_with_inverted_introgression_row(run_main, figures):
	"""The last row of the figure is 0 where a bin was called, 1 elsewhere; one SVG per comparison group, plus merged."""
	out = run_main(
		"--anc", "alpha", "--chr", "chr1", "--cmp", "DON", "OTH", "--thr", "0.3", "--vis"
	)
	names = sorted(p.rsplit("/", 1)[1] for p, _ in figures)
	assert names == [
		"alpha_chr1_DON_heatmap.svg",
		"alpha_chr1_OTH_heatmap.svg",
		"alpha_chr1_merged_heatmap.svg",
	]
	for path, fig in figures:
		assert f"{out.name}_0.3/heatmaps/" in path.replace("\\", "/")
		assert list(fig.data[0].y)[-1] == "Introgressions"
	by_name = {p.rsplit("/", 1)[1]: f for p, f in figures}
	z = np.asarray(by_name["alpha_chr1_DON_heatmap.svg"].data[0].z, dtype=float)
	assert z[-1].tolist() == [1, 1, 0, 0, 0, 1]
	zm = np.asarray(by_name["alpha_chr1_merged_heatmap.svg"].data[0].z, dtype=float)
	# merged row is 1 - (number of groups calling the bin) / (largest such number): OTH misses bin 1000
	assert zm[-1].tolist() == [1, 1, 0.5, 0, 0, 1]
	assert by_name["alpha_chr1_DON_heatmap.svg"].layout.title.text == "alpha chr1 Introgressions Called with DON"
	# row order follows group.tsv
	assert list(by_name["alpha_chr1_DON_heatmap.svg"].data[0].y)[:5] == list(GROUPS)


def test_main_vis_renders_real_svgs_when_a_renderer_exists(run_main, tmp_path):
	pytest.importorskip("kaleido")
	if not H.can_render_images(tmp_path):
		pytest.skip("plotly cannot write images here (Chrome missing)")
	out = run_main("--anc", "alpha", "--chr", "chr1", "--cmp", "REF", "--thr", "0.3", "--vis")
	svg = out / f"{out.name}_0.3" / "heatmaps" / "alpha_chr1_REFA_heatmap.svg"
	assert svg.stat().st_size > 0 and b"<svg" in svg.read_bytes()[:2000]


# ------------------------------------------------------------------ main(): argument errors


def test_main_requires_cmp_idx_tsv_out(monkeypatch):
	for argv in (
		["x", "--idx", "i", "--tsv", "t", "--out", "o"],  # no --cmp
		["x", "--cmp", "REF", "--tsv", "t", "--out", "o"],  # no --idx
		["x", "--cmp", "REF", "--idx", "i", "--out", "o"],  # no --tsv
		["x", "--cmp", "REF", "--idx", "i", "--tsv", "t"],  # no --out
	):
		monkeypatch.setattr(sys, "argv", argv)
		with pytest.raises(SystemExit) as e:
			call.main()
		assert e.value.code == 2


@pytest.mark.parametrize(
	"args,fragment",
	[
		(["--anc", "alpha", "--cmp", "REF"], "threshold"),
		(["--anc", "alpha", "--cmp", "REF", "--thr", "0.5", "--sft", "max"], "smoothing"),
		(["--cmp", "REF", "--thr", "0.5"], "anchor"),
		(["--anc", "alpha", "--grp", "ING", "--cmp", "REF", "--thr", "0.5"], "both"),
		(["--anc", "alpha", "--cmp", "REF", "DON", "--thr", "0.5"], "REF must be the only"),
		(["--anc", "alpha", "--cmp", "DON", "--thr", "0.5", "--urf", "--ref", "epsilon"], "REF must be the only"),
		(["--anc", "alpha", "--cmp", "REF", "--thr", "0.5", "--rmu", "alpha"], "--ref"),
		(["--anc", "alpha", "--cmp", "REF", "--thr", "0.5", "--rmu", "alpha", "--ref", "epsilon"], "--ogrp"),
		(
			["--anc", "alpha", "--cmp", "REF", "--thr", "0.5", "--rmu", "alpha", "--ref", "epsilon", "--ogrp", "ING"],
			"outgroup",
		),
	],
)
def test_main_rejects_invalid_option_combinations(run_main, args, fragment):
	with pytest.raises(ValueError) as e:
		run_main(*args)
	assert fragment in str(e.value)


def test_main_rejects_missing_files_and_underscore_group_names(run_main, tmp_path):
	with pytest.raises(ValueError, match="not found"):
		run_main("--anc", "alpha", "--cmp", "REF", "--thr", "0.5", tsv=tmp_path / "missing.tsv")
	with pytest.raises(ValueError, match="not found"):
		run_main("--anc", "alpha", "--cmp", "REF", "--thr", "0.5", idx=tmp_path / "no_index")
	bad = H.write_groups(tmp_path / "bad.tsv", {**GROUPS, "beta": "DO_N"})
	with pytest.raises(ValueError, match="underscore"):
		run_main("--anc", "alpha", "--cmp", "REF", "--thr", "0.5", tsv=bad)

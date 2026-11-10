"""panagram.view: the Chromosome tab (plot_chr_whole, plot_gene_content, plot_interactive via update_all_figs, small helpers).

Bar heights, heatmap rows and percentages are compared with the brute-force k-mer oracle of the
synthetic genomes; gene tracks with the synthetic GFF layout.
"""

import numpy as np
import pandas as pd
import pytest

from helpers_view import UMAP_ROWS, popcount_bins, umaps, v, view_ns  # noqa: F401

NGENOMES = 5


def _region_bits(genomes, genome, chrom, start, end):
	return genomes.oracle_bits(genome, chrom)[start:end]


def _bar_counts(bits, start, bin_size):
	"""Counts per (popcount, bin) where bin b holds coordinates [b*bin_size, (b+1)*bin_size)."""
	counts = bits.sum(axis=1)
	bins = (np.arange(len(bits)) + start) // bin_size
	first = start // bin_size
	out = np.zeros((bits.shape[1] + 1, bins.max() - first + 1), dtype=int)
	for c, b in zip(counts, bins):
		out[c, b - first] += 1
	return out, (np.arange(out.shape[1]) + first) * bin_size


def _run(v, chrom, start, end, anchor="alpha", step=1):
	return v.update_all_figs(1, chrom, anchor, 0, start, end, step)


def _bars(fig):
	"""Stacked count bars by name; the unnamed region histogram at the top left is left out."""
	return {t.name: t for t in fig.data if t.type == "bar" and t.name is not None}


# ---------------------------------------------------------------- small helpers


def test_get_buffer_threshold(v):
	"""min distance = n_skips * max_chr_bins * 8; shorter regions ask for step 1 (returns 1)."""
	need = 1 * v.params.max_chr_bins * 8
	assert v.get_buffer(0, need, 1) == 0
	assert v.get_buffer(100, 100 + need - 1, 1) == 1
	assert v.get_buffer(0, 100 * v.params.max_chr_bins * 8, 100) == 0
	assert v.get_buffer(0, 100 * v.params.max_chr_bins * 8 - 1, 100) == 1


def test_set_coords_defaults_to_whole_first_chromosome(v):
	assert v.set_coords("gamma") == ("gamma", "chr1", 0, 2486)


def test_set_coords_known_chromosome_and_end(v):
	assert v.set_coords("gamma", "chr2") == ("gamma", "chr2", 0, 686)
	assert v.set_coords("gamma", "chr2", 5, 50) == ("gamma", "chr2", 5, 50)


def test_set_coords_unknown_chromosome_falls_back_to_first(v):
	assert v.set_coords("alpha", "chrZ") == ("alpha", "chr1", 0, 2636)


def test_update_text_helpers(v):
	assert v.update_output_div("chr1", 10, 20, "alpha") == "10-20"
	assert v.update_out_chr("chr2", "alpha") == "alpha.chr2:"
	assert v.update_chromosome_list("gamma") == [
		{"label": "chr1", "value": "chr1"},
		{"label": "chr2", "value": "chr2"},
	]


def _genes(names):
	return pd.DataFrame({"name": names})


def test_update_gene_locals_by_count(v):
	assert v.update_gene_locals(None, "chr1", 0, 10, "alpha") == ""
	assert v.update_gene_locals(_genes([]), "chr1", 0, 10, "alpha") == ""
	assert v.update_gene_locals(_genes(["g1"]), "chr1", 0, 10, "alpha") == "Genes: g1"
	few = v.update_gene_locals(_genes(["g1", "g2", "g3"]), "chr1", 0, 10, "alpha")
	assert few.startswith("Genes: g1, g2, g3")
	assert v.update_gene_locals(_genes([f"g{i}" for i in range(25)]), "chr1", 0, 10, "alpha").startswith(
		"Genes: g0, g1"
	)
	assert v.update_gene_locals(_genes([f"g{i}" for i in range(26)]), "chr1", 0, 10, "alpha") == (
		"Genes in this region: 26"
	)


def test_make_gene_whole_chr_counts_every_gene_once(v):
	x = pd.Index([0, 26, 52, 78])
	z = v.make_gene_whole_chr("alpha", x, np.array([1, 30, 31, 70, 1000]))
	assert z.sum() == 5
	assert z.dtype.kind == "i"
	assert v.make_gene_whole_chr("alpha", x, np.array([], dtype=int)).sum() == 0


def test_make_gene_whole_chr_edge_placement(v):
	"""A gene start on a bin edge belongs to the bin that starts there: index k holds starts in [x[k-1], x[k])."""
	x = pd.Index([0, 26, 52, 78])
	z = v.make_gene_whole_chr("alpha", x, np.array([0, 25, 26, 77, 78, 1000]))
	assert list(z) == [0, 2, 1, 1, 2, 0]


def test_get_local_info_just_one_ignores_zero_genome_bin(v):
	"""Entry 0 counts k-mers shared by no genome and is not part of the percentages (hand-worked: 6 of 8)."""
	bar = v.get_local_info_just_one(pd.Series([100, 6, 2, 0, 0, 0]), "alpha", "chr1")
	assert list(bar.y) == pytest.approx([75, 25, 0, 0, 0])


def test_get_local_info_just_one_percentages(v):
	"""bar_sum_regional = [0, 8, 2, 0, 0, 0] -> 80, 20, 0, 0, 0 percent for 1..5 genomes."""
	bar = v.get_local_info_just_one(pd.Series([0, 8, 2, 0, 0, 0]), "alpha", "chr1")
	assert list(bar.x) == [1, 2, 3, 4, 5]
	assert list(bar.y) == pytest.approx([80, 20, 0, 0, 0])
	assert list(bar.marker.color) == list(v.colors)
	again = v.get_local_info_just_one(np.array([0, 1, 1, 1, 1, 0]), "alpha", "chr1")
	assert list(again.y) == pytest.approx([25, 25, 25, 25, 0])


def test_get_local_info_region_panel_and_titles(v):
	"""The first panel is the regional distribution in percent over 1..5 genomes (hand-worked)."""
	fig = v.get_local_info(pd.Series([0, 6, 2, 2, 0, 0]), "alpha", "chr1")
	first = fig.data[0]
	assert list(first.x) == [1, 2, 3, 4, 5]
	assert list(first.y) == pytest.approx([60, 20, 20, 0, 0])
	assert [a.text for a in fig.layout.annotations] == ["Whole chromosome", "This region", "Genes"]
	assert len(fig.data) == 3


# ---------------------------------------------------------------- plot_chr_whole


def test_plot_chr_whole_density_rows(v, umaps, synthetic_genomes):
	"""Rows: clusters, universal fraction, private fraction, gene density; boxes mark the selection."""
	g = synthetic_genomes
	genes = v.index.query_genes("alpha", "chr1", 0, 2636)
	fig = v.plot_chr_whole(100, 900, "alpha", "chr1", genes)
	cluster, box1, universal, box2, private, box3, gene_row, box4 = fig.data
	assert list(cluster.x) == [0, 1000, 2000, 0, 1000, 2000]
	assert list(cluster.y) == [1, 1, 1, 2, 2, 2]
	assert list(cluster.z) == [1000, 2000, 3000, 1000, 2000, 3000]
	edges = [int(e) for e in universal.x]
	hist = popcount_bins(g.oracle_bits("alpha", "chr1"), edges)
	totals = hist.sum(axis=0)
	assert totals.sum() == 2636
	assert list(universal.z) == pytest.approx(list(hist[NGENOMES] / totals))
	assert list(private.z) == pytest.approx(list(hist[1] / totals))
	assert list(universal.y) == [1] * len(edges)
	for box in (box1, box2, box3, box4):
		assert list(box.x)[:3] == [100, 100, None]
		assert list(box.x)[3:5] == [900, 900]
	assert fig.layout.xaxis.range == (0, 2636)
	assert fig.layout.dragmode == "select" and fig.layout.height == 350
	assert sum(gene_row.z) == len(genes) == 4


def test_plot_chr_whole_without_genes_has_empty_gene_row(v, umaps):
	fig = v.plot_chr_whole(0, 2636, "alpha", "chr1", pd.DataFrame(columns=["start"]))
	gene_row = fig.data[6]
	assert list(gene_row.z) == [0] * len(gene_row.x)


# ---------------------------------------------------------------- plot_gene_content


def _gene_table(v):
	genes = v.index.query_genes("alpha", "chr1", 0, 2636)
	genes["local"] = genes["name"].str.contains("CORE|PAIR")
	return genes


def test_plot_gene_content_points_follow_the_table(v):
	"""Each gene is a point at its rank by universal fraction; fractions are k-mer count / (end - start)."""
	genes = _gene_table(v)
	fig = v.plot_gene_content(genes, "alpha", "chr1", 0, 2636)
	uniq, univ, hi_univ, hi_uniq = fig.data
	assert uniq.name == "% Unique" and univ.name == "% Universal"
	size = genes["end"] - genes["start"]
	expect_univ = dict(zip(genes["name"], genes[NGENOMES] / size))
	expect_uniq = dict(zip(genes["name"], genes[1] / size))
	assert list(uniq.x) == list(range(len(genes)))
	assert dict(zip(uniq.text, uniq.y)) == pytest.approx(expect_uniq)
	assert dict(zip(uniq.text, univ.y)) == pytest.approx(expect_univ)
	assert list(univ.y) == sorted(univ.y)
	by_x = dict(zip(uniq.x, uniq.text))
	local_names = {by_x[x] for x in hi_univ.x}
	assert local_names == set(genes.loc[genes["local"], "name"])
	assert {by_x[x] for x in hi_uniq.x} == local_names
	assert fig.layout.xaxis.title.text == "Unique k-mer content rank"
	assert fig.layout.yaxis.title.text == "% of k-mers"


def test_plot_gene_content_gene_classes(v):
	"""Private filler genes are unique, CORE genes universal (synthetic construction)."""
	genes = _gene_table(v)
	fig = v.plot_gene_content(genes, "alpha", "chr1", 0, 2636)
	uniq, univ = fig.data[:2]
	seen = {}
	for name, a, b in zip(uniq.text, uniq.y, univ.y):
		seen[name] = (a, b)
	assert seen["alpha_chr1_private"][0] > 0.9 and seen["alpha_chr1_private"][1] < 0.01
	assert seen["alpha_chr1_CORE_1"][1] > 0.9
	assert seen["alpha_chr1_PAIR_3"][1] == 0


# ---------------------------------------------------------------- update_all_figs / plot_interactive


@pytest.fixture
def whole_chr1(v, umaps):
	return _run(v, "chr1", 0, 2636)


def test_update_all_figs_returns_seven_outputs(whole_chr1):
	chr_fig, main, shared, genes_fig, umap, umap_hist, text = whole_chr1
	assert [f.__class__.__name__ for f in (chr_fig, main, shared, genes_fig, umap, umap_hist)] == ["Figure"] * 6
	assert text.startswith("Genes: alpha_chr1_private, alpha_chr1_CORE_1")


def test_main_figure_stacked_bars_are_oracle_counts(v, umaps, synthetic_genomes):
	"""Region 0-2636, bin size 2636 // 350 + 1 = 8: bar i holds the k-mers shared by i genomes per bin."""
	_, main, *_ = _run(v, "chr1", 0, 2636)
	bits = _region_bits(synthetic_genomes, "alpha", "chr1", 0, 2636)
	expected, xs = _bar_counts(bits, 0, 8)
	bars = _bars(main)
	assert sorted(bars) == ["0", "1", "2", "3", "4", "5"]
	for i in range(NGENOMES + 1):
		assert list(bars[str(i)].y) == list(expected[i]), i
		assert list(bars[str(i)].x) == list(xs)
	assert main.layout.barmode == "stack"
	assert main.layout.yaxis6.range == (0, 8)


def test_main_figure_subregion_bins_align_to_multiples(v, umaps, synthetic_genomes):
	"""Region 100-900: bin size 3, the first bin holds only the coordinates 100 and 101."""
	_, main, *_ = _run(v, "chr1", 100, 900)
	bits = _region_bits(synthetic_genomes, "alpha", "chr1", 100, 900)
	expected, xs = _bar_counts(bits, 100, 3)
	bars = _bars(main)
	for i in range(NGENOMES + 1):
		assert list(bars[str(i)].y) == list(expected[i]), i
	assert list(bars["1"].x) == list(xs)
	assert xs[0] == 99
	assert main.layout.xaxis2.range == (100, 900)


def test_main_figure_pair_heatmap_rows_per_genome(v, umaps, synthetic_genomes):
	"""Heatmap z[row g, bin] = fraction of the bin's k-mers that genome g also has."""
	g = synthetic_genomes
	_, main, *_ = _run(v, "chr1", 0, 2636)
	bits = _region_bits(g, "alpha", "chr1", 0, 2636)
	heat = [t for t in main.data if t.type == "heatmap" and t.coloraxis == "coloraxis2"]
	assert len(heat) == 1
	heat = heat[0]
	assert sorted(heat.y) == sorted(g.names)
	for row, name in zip(np.asarray(heat.z), heat.y):
		col = bits[:, g.names.index(name)]
		want = [col[b * 8 : (b + 1) * 8].mean() for b in range(int(np.ceil(2636 / 8)))]
		assert list(row) == pytest.approx(want), name


def test_main_figure_gene_tracks_use_synthetic_gene_ends(v, umaps, synthetic_genomes):
	_, main, *_ = _run(v, "chr1", 0, 2636)
	ends = sorted(e for c, s, e, n in synthetic_genomes.data.genes["alpha"] if c == "chr1")
	genes = [t for t in main.data if t.name == "gene"]
	assert len(genes) == 1
	assert list(genes[0].x)[1::3] == ends
	assert list(genes[0].x)[2::3] == [None] * len(ends)
	by_end = sorted((e, n) for c, s, e, n in synthetic_genomes.data.genes["alpha"] if c == "chr1")
	assert list(genes[0].text)[::3] == [n for e, n in by_end]
	exons = [t for t in main.data if t.name == "exon"]
	assert len(exons) == 1
	assert list(exons[0].x)[1::3] == ends


def test_main_figure_gene_track_limited_to_region(v, umaps):
	"""Region 100-900 overlaps the private gene (ends 300) and the CORE gene (ends 900) only."""
	_, main, *_ = _run(v, "chr1", 100, 900)
	genes = [t for t in main.data if t.name == "gene"][0]
	assert list(genes.x)[1::3] == [300, 900]


def test_main_figure_coordinate_ruler_and_clusters(v, umaps):
	_, main, *_ = _run(v, "chr1", 100, 900)
	ruler = [t for t in main.data if getattr(t, "mode", None) == "lines+markers+text"][0]
	assert list(ruler.x) == [int(round(x)) for x in np.linspace(100, 901, 10)]
	assert list(ruler.text) == [str(int(x)) for x in ruler.x]
	clusters = [t for t in main.data if t.type == "heatmap" and t.coloraxis is None][-1]
	assert list(clusters.x) == [0, 1000, 2000, 0, 1000, 2000]
	assert list(clusters.y) == [1, 1, 1, 2, 2, 2]


def test_main_figure_tree_shapes(v, umaps):
	_, main, *_ = _run(v, "chr1", 0, 2636)
	tree_nodes = [t for t in main.data if t.type == "scatter" and t.mode == "markers"]
	assert len(tree_nodes) == 1
	assert len(main.layout.shapes) == 2 * NGENOMES - 1 + NGENOMES - 1


def test_percent_shared_figure_difference_to_genome_wide(v, umaps, synthetic_genomes):
	"""fig2 = % of the region's k-mers shared (relative to the largest genome) - genome-wide % shared (oracle)."""
	g = synthetic_genomes
	_, _, shared, *_ = _run(v, "chr1", 0, 2636)
	region = _region_bits(g, "alpha", "chr1", 0, 2636).sum(axis=0)
	whole = np.vstack([g.oracle_bits("alpha", c) for c in g.data.seqs["alpha"]]).mean(axis=0)
	got = dict(zip(shared.data[0].x, shared.data[0].y))
	for j, name in enumerate(g.names):
		assert got[name] == pytest.approx(100.0 * region[j] / region.max() - 100.0 * whole[j])


def test_order_parameter_sets_heatmap_row_order(v, umaps, monkeypatch, synthetic_genomes):
	order = ["epsilon", "delta", "gamma", "beta", "alpha"]
	monkeypatch.setattr(v.params, "order", order)
	_, main, *_ = _run(v, "chr1", 0, 2636)
	heat = [t for t in main.data if t.type == "heatmap" and t.coloraxis == "coloraxis2"][0]
	assert list(heat.y) == order
	bits = _region_bits(synthetic_genomes, "alpha", "chr1", 0, 2636)
	assert list(heat.z[0]) == pytest.approx([bits[b * 8 : (b + 1) * 8, 4].mean() for b in range(330)])


def test_max_chr_bins_sets_bin_size(v, umaps, monkeypatch):
	monkeypatch.setattr(v.params, "max_chr_bins", 100)
	_, main, *_ = _run(v, "chr1", 0, 2636)
	bar = _bars(main)["1"]
	assert list(bar.x) == [27 * b for b in range(int(np.ceil(2636 / 27)))]
	assert main.layout.yaxis6.range == (0, 27)


def test_low_resolution_step_samples_every_hundredth_kmer(v, umaps, monkeypatch, synthetic_genomes):
	"""max_chr_bins=3 makes a 2636-wide region long enough to keep step 100: 27 sampled k-mers."""
	monkeypatch.setattr(v.params, "max_chr_bins", 3)
	_, main, shared, *_ = _run(v, "chr1", 0, 2636, step=100)
	bits = _region_bits(synthetic_genomes, "alpha", "chr1", 0, 2636)[::100]
	assert len(bits) == 27
	bars = _bars(main)
	totals = {i: int(sum(bars[str(i)].y)) for i in range(NGENOMES + 1)}
	want = np.bincount(bits.sum(axis=1), minlength=NGENOMES + 1)
	assert totals == {i: int(want[i]) for i in range(NGENOMES + 1)}
	assert list(bars["1"].x) == [0, 879, 1758]


def test_genome_without_annotation_has_no_gene_tracks(v, umaps, monkeypatch):
	"""epsilon has no GFF: no gene or exon trace and an empty gene text; the figures still build."""
	genome = v.index.genomes["epsilon"]
	table = pd.DataFrame(
		{"chrom": ["chr1", "chr1"], "start": [0, 1000], "end": [1000, 2000], "umap1": [0.0, 1.0], "umap2": [1.0, 0.0], "cluster": [0, 1]}
	)
	monkeypatch.setattr(genome, "genome_umap", table)
	monkeypatch.setattr(genome, "chrom_umaps", table.set_index("chrom"))
	chr_fig, main, shared, genes_fig, umap, hist, text = _run(v, "chr1", 0, 2486, anchor="epsilon")
	assert text == ""
	assert not [t for t in main.data if t.name in ("gene", "exon")]
	assert len(genes_fig.data) >= 2 and list(genes_fig.data[0].x) == []


def test_plot_interactive_direct_call_returns_figure_and_percent_shared(v, umaps, synthetic_genomes):
	"""Called with an oracle bitmap: the returned percentages are each genome's share of the largest genome's k-mers."""
	g = synthetic_genomes
	bits = g.oracle_bits("alpha", "chr1")
	bitmap = v.index.query_bitmap("alpha", "chr1", 0, 2636, 1)
	assert (bitmap.to_numpy() == bits).all()
	pan, pair = v.index.bitmap_to_bins(bitmap, 8)
	genes = v.index.query_genes("alpha", "chr1", 0, 2636)
	fig, shared = v.plot_interactive(pan.sum(axis=1), "alpha", "chr1", 0, 2636, 1, bitmap, pan, pair, genes)
	sums = bits.sum(axis=0)
	assert shared == {n: pytest.approx(100.0 * s / sums.max()) for n, s in zip(g.names, sums)}
	assert [t.name for t in fig.data if t.type == "bar" and t.name is not None] == ["0", "1", "2", "3", "4", "5"]
	assert fig.layout.barmode == "stack"

"""panagram.view: figures of the Pangenome and Anchor tabs (counts, sizes, averages, UMAPs, per-chromosome rows).

Expected values come from the synthetic ground truth (tests/synthetic.py layout and the
brute-force k-mer oracle) and from hand-written UMAP tables.
"""

import numpy as np
import pandas as pd
import pytest

import synthetic
from conftest import K
from helpers_view import UMAP_ROWS, popcount_bins, umaps, v, view_ns  # noqa: F401


def _xy(trace):
	return list(trace.x), list(trace.y)


def _vline_x(fig):
	return [s.x0 for s in fig.layout.shapes if s.type == "line" and s.x0 == s.x1]


def _positions(genomes, genome):
	"""Number of k-mer positions of a genome (the index reports sizes in k-mer starts)."""
	return sum(len(s) - K + 1 for s in genomes.data.seqs[genome].values())


# ---------------------------------------------------------------- pangenome tab


def test_genome_count_plot_is_scaffold_count_ascending(v, synthetic_genomes):
	fig = v.make_genome_count_plot("beta")
	x, y = _xy(fig.data[0])
	expected = {g: len(synthetic.LAYOUT[g]) for g in synthetic.GENOMES}
	assert dict(zip(x, y)) == expected
	assert y == sorted(y)
	assert fig.layout.yaxis.title.text == "# of scaffolds"


def test_genome_count_plot_marks_the_anchor(v):
	for anchor in ("alpha", "delta"):
		assert _vline_x(v.make_genome_count_plot(anchor)) == [anchor]


def test_genome_size_plot_is_total_kmer_positions_ascending(v, synthetic_genomes):
	fig = v.make_genome_size_plot("gamma")
	x, y = _xy(fig.data[0])
	assert dict(zip(x, y)) == {g: _positions(synthetic_genomes, g) for g in synthetic.GENOMES}
	assert y == sorted(y)
	assert _vline_x(fig) == ["gamma"]


def test_plot_avgs_is_mean_number_of_genomes_per_kmer(v, synthetic_genomes):
	"""Average k-mer content of a genome = mean popcount of its bitmap rows (oracle)."""
	g = synthetic_genomes
	fig = v.plot_avgs("alpha")
	x, y = _xy(fig.data[0])
	for name, value in zip(x, y):
		bits = np.vstack([g.oracle_bits(name, c) for c in g.data.seqs[name]])
		assert value == pytest.approx(bits.sum(axis=1).mean())
	assert sorted(x) == sorted(synthetic.GENOMES)
	assert y == sorted(y)
	assert _vline_x(fig) == ["alpha"]
	assert fig.layout.yaxis.title.text == "Average k-mer"


def test_read_genome_size_files_sorted_tables(v, synthetic_genomes):
	sizes, seqs = v.read_genome_size_files()
	assert list(sizes.columns) == ["Sample", "Size"]
	assert list(seqs.columns) == ["Sample", "Seq"]
	assert dict(zip(sizes["Sample"], sizes["Size"])) == {
		g: _positions(synthetic_genomes, g) for g in synthetic.GENOMES
	}
	assert dict(zip(seqs["Sample"], seqs["Seq"])) == {g: len(synthetic.LAYOUT[g]) for g in synthetic.GENOMES}
	assert list(sizes["Size"]) == sorted(sizes["Size"])
	assert list(seqs["Seq"]) == sorted(seqs["Seq"])


def test_avg_kmer_fig_per_chromosome_means(v, synthetic_genomes):
	g = synthetic_genomes
	fig = v.make_avg_kmer_fig("alpha")
	x, y = _xy(fig.data[0])
	assert sorted(x) == ["chr1", "chr2"]
	for chrom, value in zip(x, y):
		assert value == pytest.approx(g.oracle_bits("alpha", chrom).sum(axis=1).mean())
	buttons = fig.layout.updatemenus[0].buttons
	assert [b.label for b in buttons] == ["Linear", "Log"]
	assert [b.args[0]["yaxis.type"] for b in buttons] == ["linear", "log"]


def test_genes_per_chr_fig_counts_genes(v, synthetic_genomes):
	"""Gene count per chromosome equals the number of genes in the synthetic GFF."""
	fig = v.make_genes_per_chr_fig("alpha")
	x, y = _xy(fig.data[0])
	expected = {}
	for chrom, *_ in synthetic_genomes.data.genes["alpha"]:
		expected[chrom] = expected.get(chrom, 0) + 1
	assert dict(zip(x, y)) == expected
	assert [b.label for b in fig.layout.updatemenus[0].buttons] == ["Linear", "Log"]


def test_gene_per_genome_fig_unique_and_universal_fractions(v, synthetic_genomes):
	"""Private filler genes are unique and not universal; CORE genes universal; PAIR/DIV genes neither."""
	fig = v.make_gene_per_genome_fig("alpha")
	unique, universal = fig.data
	assert unique.name == "Proportion Unique" and universal.name == "Proportion Universal"
	expected_names = {f"{c}:{n}" for c, _, _, n in synthetic_genomes.data.genes["alpha"]}
	assert set(unique.text) == expected_names
	assert list(unique.x) == list(range(len(expected_names)))
	assert list(universal.y) == sorted(universal.y)
	by_name = {t: (u, w) for t, u, w in zip(unique.text, unique.y, universal.y)}
	for name, (u, w) in by_name.items():
		if name.endswith("_private"):
			assert u > 0.9 and w < 0.01, name
		elif "_CORE_" in name:
			assert w > 0.9 and u < 0.1, name
		else:  # PAIR (alpha, beta) and DIV (not epsilon): shared by some genomes only
			assert u < 0.1 and w == 0, name


def test_anchor_tab_legend_figure_has_hidden_axes(v):
	fig = v.make_anchor_tab_legends("alpha")
	assert fig.layout.xaxis.visible is False
	assert fig.layout.yaxis.visible is False
	assert len(fig.data) == 0


# ---------------------------------------------------------------- perc-shared bars


def test_plot_perc_shared_difference_and_colours(v, synthetic_genomes):
	"""Bar = % of k-mers shared in the region minus the genome-wide % (oracle); colour by sign."""
	g = synthetic_genomes
	genome = g.data.seqs["alpha"]
	bits = np.vstack([g.oracle_bits("alpha", c) for c in genome])
	frac = dict(zip(g.names, bits.mean(axis=0)))
	perc = {n: 100 * frac[n] for n in g.names}
	perc["beta"] += 7.5
	perc["gamma"] -= 3.0
	table = v.index.genomes["alpha"].total_paircounts
	perc["delta"] = table.loc["delta", "frac"] * 100  # exactly equal -> grey
	fig = v.plot_perc_shared(perc, "alpha")
	trace = fig.data[0]
	got = dict(zip(trace.x, trace.y))
	assert got["beta"] == pytest.approx(7.5)
	assert got["gamma"] == pytest.approx(-3.0)
	assert got["delta"] == 0
	assert got["epsilon"] == pytest.approx(0, abs=1e-9)
	assert got["alpha"] == pytest.approx(0, abs=1e-9)
	colors = dict(zip(trace.x, trace.marker.color))
	assert colors["beta"] == "#00b4d8"
	assert colors["gamma"] == "#c1121f"
	assert colors["delta"] == "grey"
	assert list(trace.x) == list(table.sort_values("frac").index)


def test_total_paircounts_match_the_oracle(v, synthetic_genomes):
	"""The fractions plot_perc_shared sorts by agree with the oracle share of alpha's k-mers per genome."""
	g = synthetic_genomes
	bits = np.vstack([g.oracle_bits("alpha", c) for c in g.data.seqs["alpha"]])
	table = v.index.genomes["alpha"].total_paircounts
	for j, name in enumerate(g.names):
		assert table.loc[name, "frac"] == pytest.approx(bits[:, j].mean())


# ---------------------------------------------------------------- UMAP figures (hand-written tables)


def test_anchor_umaps_per_chr(v, umaps):
	data, zmax = v.anchor_umaps_per_chr("alpha")
	assert zmax == 2
	assert sorted(data) == ["chr1", "chr2"]
	assert list(data["chr1"]["umap1"]) == [-1.0, 0.5, 2.0]
	assert list(data["chr1"]["cluster"]) == [0, 1, 1]
	assert list(data["chr2"]["umap2"]) == [2.0, 0.25]
	assert list(data["chr2"].columns) == ["umap1", "umap2", "cluster"]


def test_make_genome_umap_points_and_hover(v, umaps):
	fig = v.make_genome_umap("alpha")
	trace = fig.data[0]
	assert list(trace.x) == UMAP_ROWS["umap1"]
	assert list(trace.y) == UMAP_ROWS["umap2"]
	assert list(trace.marker.color) == UMAP_ROWS["cluster"]
	assert trace.hovertext[0] == "Cluster: 0<br>Chromosome: chr1:0-1000<br>X: -1.0<br>Y: 0.0"
	assert trace.hovertext[4] == "Cluster: 2<br>Chromosome: chr2:1000-2000<br>X: -2.5<br>Y: 0.25"
	assert len(trace.hovertext) == 5
	assert fig.layout.coloraxis.colorbar.title.text == "DBSCAN Cluster"
	# crosshair lines through the origin
	lines = {(s.x0, s.y0, s.x1, s.y1) for s in fig.layout.shapes}
	assert any(s[0] == s[2] == 0 for s in lines) and any(s[1] == s[3] == 0 for s in lines)


def test_get_genome_umap_hist_bp_per_cluster(v, umaps):
	"""Clusters 0, 1, 2 hold 2, 2 and 1 windows of 1000 bp."""
	fig = v.get_genome_umap_hist("alpha")
	bar = fig.data[0]
	assert list(bar.y) == [0, 1, 2]
	assert list(bar.x) == [2000, 2000, 1000]
	assert bar.orientation == "h"
	assert fig.layout.xaxis.type == "log"


def test_get_umap_chrom_hist(v, umaps):
	one = v.get_umap_chrom_hist("alpha", "chr1")  # clusters 0,1,1
	assert list(one.data[0].y) == [0, 1]
	assert list(one.data[0].x) == [1000, 2000]
	two = v.get_umap_chrom_hist("alpha", "chr2")  # clusters 0,2
	assert list(two.data[0].y) == [0, 1, 2]
	assert list(two.data[0].x) == [1000, 0, 1000]


def test_get_umap_from_file_one_chromosome(v, umaps):
	fig = v.get_umap_from_file("chr1", "alpha")
	trace = fig.data[0]
	assert list(trace.x) == [-1.0, 0.5, 2.0]
	assert list(trace.y) == [0.0, 1.0, -1.0]
	assert list(trace.marker.color) == [0, 1, 1]
	assert list(trace.hovertext) == [
		"Cluster: 0<br>Chromosome: chr1:0-1000",
		"Cluster: 1<br>Chromosome: chr1:1000-2000",
		"Cluster: 1<br>Chromosome: chr1:2000-3000",
	]
	assert fig.layout.xaxis.showticklabels is False


def test_read_clusters_two_rows_of_windows(v, umaps):
	data = v.read_clusters("alpha", "chr1")
	assert data == {
		"x": [0, 1000, 2000, 0, 1000, 2000],
		"y": [1, 1, 1, 2, 2, 2],
		"z": [1000, 2000, 3000, 1000, 2000, 3000],
	}


# ---------------------------------------------------------------- whole-genome rows


def test_plot_whole_genome_rows_per_chromosome(v, umaps, synthetic_genomes):
	"""Per chromosome: a cluster heatmap, a universal row and a private row, from the bin table and UMAP."""
	g = synthetic_genomes
	fig = v.plot_whole_genome("alpha")
	heatmaps = [t for t in fig.data if t.type == "heatmap"]
	assert len(heatmaps) == 6  # 2 chromosomes x (cluster, universal, private)
	for n, chrom in enumerate(["chr1", "chr2"]):
		cluster, universal, private = heatmaps[3 * n : 3 * n + 3]
		umap = umaps[umaps["chrom"] == chrom]
		assert list(cluster.x) == list(umap["umap1"])
		assert list(cluster.z) == list(umap["cluster"])
		edges = [int(e) for e in universal.x]
		assert edges[0] == 0 and edges == sorted(edges)
		hist = popcount_bins(g.oracle_bits("alpha", chrom), edges)
		assert list(universal.z) == list(hist[len(g.names)])
		assert list(private.z) == list(hist[1])
	assert fig.layout.height == 2 * 250
	assert [a.text for a in fig.layout.annotations] == ["chr1", "chr2"]

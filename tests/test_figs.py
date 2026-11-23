"""panagram.figs: genome colours and the three summary figures, checked against hand-worked and brute-force values."""

import itertools
import types

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("plotly")
pytest.importorskip("scipy")
pytest.importorskip("matplotlib")

import synthetic
from panagram import figs


def _ns(n):
	return types.SimpleNamespace(ngenomes=n)


# ------------------------------------------------------------------ genome_colors


def test_genome_colors_empty_and_single():
	"""No genomes give no colours; one genome gets the colormap's midpoint (viridis midpoint is #21908d)."""
	assert figs.genome_colors(_ns(0)) == []
	assert figs.genome_colors(_ns(1)) == ["#21908d"]


def test_genome_colors_spread_viridis_reversed():
	"""n genomes sample the reversed viridis evenly: first is its yellow end (#fde725), last its purple end (#440154)."""
	assert figs.genome_colors(_ns(2)) == ["#fde725", "#440154"]
	assert figs.genome_colors(_ns(3)) == ["#fde725", "#21908d", "#440154"]
	five = figs.genome_colors(_ns(5))
	assert len(five) == 5 and len(set(five)) == 5
	assert five[0] == "#fde725" and five[-1] == "#440154" and five[2] == "#21908d"
	assert all(len(c) == 7 and c[0] == "#" and c == c.lower() for c in five)


def test_genome_colors_custom_colormap():
	"""A named colormap argument replaces the default: gray runs black, mid-gray, white."""
	assert figs.genome_colors(_ns(3), cmap="gray") == ["#000000", "#808080", "#ffffff"]
	assert figs.genome_colors(_ns(1), cmap="gray") == ["#808080"]


# ------------------------------------------------------------- read_genome_comp


def test_read_genome_comp_has_one_bar_panel_per_chromosome(anchored_reader):
	"""One bar trace per chromosome of the anchor, bars at x = 0..n with y = percent of positions, genome colours, log y axes."""
	idx = anchored_reader
	fig = figs.read_genome_comp(idx, "alpha")
	totals = idx["alpha"].bitfreq_chrs
	assert len(fig.data) == len(totals) == 2
	for trace, (chrom, freqs) in zip(fig.data, totals.iterrows()):
		np.testing.assert_allclose(np.array(trace.y, dtype=float), freqs.to_numpy() * 100)
		assert list(trace.x) == [0, 1, 2, 3, 4, 5]
		assert list(trace.marker.color) == figs.genome_colors(idx)
		assert trace.showlegend is False
		assert trace.type == "bar"
	for ax in ("yaxis", "yaxis2"):
		assert fig.layout[ax].type == "log"
	assert fig.layout.paper_bgcolor == "rgba(0,0,0,0)"


def test_read_genome_comp_percentages_match_bruteforce_occupancy(anchored_reader, synthetic_genomes):
	"""Panel c of genome g plots, for sharing level n = 0..5, the percent of that chromosome's positions shared by n genomes (brute force)."""
	g = synthetic_genomes
	for genome in ("alpha", "gamma"):
		fig = figs.read_genome_comp(anchored_reader, genome)
		chroms = list(g.data.seqs[genome])
		assert len(fig.data) == len(chroms)
		for trace, chrom in zip(fig.data, chroms):
			occ = g.oracle_bits(genome, chrom).sum(axis=1)
			want = np.bincount(occ, minlength=6) / len(occ) * 100
			np.testing.assert_allclose(np.array(trace.y, dtype=float), want, atol=1e-9)


def test_read_genome_comp_percentages_sum_to_100(anchored_reader):
	"""Each panel's bars sum to 100 percent."""
	fig = figs.read_genome_comp(anchored_reader, "gamma")
	for trace in fig.data:
		assert abs(np.sum(np.array(trace.y, dtype=float)) - 100) < 1e-9


# ------------------------------------------------------------ read_pangenome_comp


def test_read_pangenome_comp_stacks_one_trace_per_sharing_level(anchored_reader):
	"""Trace i (named i+1) holds, per genome, the fraction of its positions shared by i+1 genomes; bars are horizontal and stacked."""
	idx = anchored_reader
	fig = figs.read_pangenome_comp(idx)
	assert len(fig.data) == 5
	colors = figs.genome_colors(idx)
	for i, trace in enumerate(fig.data):
		assert trace.name == str(i + 1)
		assert list(trace.y) == synthetic.GENOMES
		np.testing.assert_allclose(np.array(trace.x, dtype=float), idx.bitfreq_totals[i + 1].to_numpy())
		assert trace.orientation == "h"
		assert trace.marker.color == colors[i]
	assert fig.layout.barmode == "stack"
	total = sum(np.array(t.x, dtype=float) for t in fig.data)
	np.testing.assert_allclose(total, 1.0)  # position 0 never occurs, so the five levels cover everything


def test_read_pangenome_comp_fractions_match_bruteforce_occupancy(anchored_reader, synthetic_genomes):
	"""Trace n, bar of genome g: the fraction of g's positions whose k-mer is present in exactly n genomes (brute force)."""
	g = synthetic_genomes
	fig = figs.read_pangenome_comp(anchored_reader)
	for j, genome in enumerate(g.names):
		occ = np.concatenate([g.oracle_bits(genome, c).sum(axis=1) for c in g.data.seqs[genome]])
		want = np.bincount(occ, minlength=6) / len(occ)
		got = [float(fig.data[n - 1].x[j]) for n in range(1, 6)]
		np.testing.assert_allclose(got, want[1:], atol=1e-12)


def test_read_pangenome_comp_has_linear_log_buttons():
	"""The axis buttons switch the x axis between linear and log with relayout."""
	idx = types.SimpleNamespace(
		ngenomes=2,
		genome_names=pd.Index(["a", "b"], name="name"),
		bitfreq_totals=pd.DataFrame({0: [0.0, 0.0], 1: [0.5, 0.2], 2: [0.5, 0.8]}, index=["a", "b"]),
	)
	fig = figs.read_pangenome_comp(idx)
	(menu,) = fig.layout.updatemenus
	assert [b.label for b in menu.buttons] == ["Linear", "Log"]
	assert [b.method for b in menu.buttons] == ["relayout", "relayout"]
	assert [b.args[0]["xaxis.type"] for b in menu.buttons] == ["linear", "log"]
	assert list(fig.data[0].x) == [0.5, 0.2] and list(fig.data[1].x) == [0.5, 0.8]


# ----------------------------------------------------------- make_all_genome_dend

DIST = [("a", "b", 1.0), ("a", "c", 4.0), ("a", "d", 4.5), ("b", "c", 3.5), ("b", "d", 4.2), ("c", "d", 1.5)]


def _dist_index(tmp_path, rows=DIST, names="abcd"):
	f = tmp_path / "genome_dist.tsv"
	f.write_text("".join(f"{x}\t{y}\t{v}\t0\t1/2\n" for x, y, v in rows))
	return types.SimpleNamespace(
		ngenomes=len(names),
		genome_dist_fname=str(f),
		genomes={n: types.SimpleNamespace(id=i) for i, n in enumerate(names)},
	)


def _matrix(rows, names):
	m = np.zeros((len(names), len(names)))
	for x, y, v in rows:
		i, j = names.index(x), names.index(y)
		m[i, j] = m[j, i] = v
	return m


def test_dendrogram_figure_has_two_trees_and_a_heatmap(tmp_path):
	"""Top tree, side tree (n - 1 merges each) and one heatmap; leaves are the genome names, closest genomes adjacent."""
	idx = _dist_index(tmp_path)
	fig = figs.make_all_genome_dend(idx)
	kinds = [t.type for t in fig.data]
	ntree = (len(fig.data) - 1) // 2  # traces per tree: plotly groups merges by colour, so the count depends on the tree
	assert len(fig.data) == 2 * ntree + 1 and ntree >= 1
	assert kinds == ["scatter"] * (2 * ntree) + ["heatmap"]
	assert [t.yaxis for t in fig.data[:ntree]] == ["y2"] * ntree
	assert [t.xaxis for t in fig.data[ntree : 2 * ntree]] == ["x2"] * ntree
	labels = list(fig.layout.xaxis.ticktext)
	assert sorted(labels) == ["a", "b", "c", "d"]
	assert {frozenset(labels[:2]), frozenset(labels[2:])} == {frozenset("ab"), frozenset("cd")}
	assert fig.layout.showlegend is False
	assert fig.layout.xaxis.domain == (0.2, 1)
	assert fig.layout.yaxis2.domain == (0.79, 0.975)


def test_dendrogram_heatmap_holds_pairwise_distances_between_distance_rows(tmp_path):
	"""The heatmap is symmetric with a zero diagonal; its cells are Euclidean distances between rows of the genome distance matrix."""
	rows = DIST
	idx = _dist_index(tmp_path, rows)
	fig = figs.make_all_genome_dend(idx)
	z = np.array(fig.data[-1].z, dtype=float)
	assert z.shape == (4, 4)
	np.testing.assert_allclose(z, z.T)
	assert (np.diag(z) == 0).all()
	m = _matrix(rows, list("abcd"))
	want = sorted(np.linalg.norm(m[i] - m[j]) for i, j in itertools.combinations(range(4), 2))
	got = sorted(z[np.triu_indices(4, 1)])
	np.testing.assert_allclose(got, want)
	assert z[np.triu_indices(4, 1)].min() == pytest.approx(np.linalg.norm(m[0] - m[1]))  # a and b are the nearest rows
	# rows and columns follow the leaf order of the tree (or its mirror), not the input order
	full = np.linalg.norm(m[:, None] - m[None], axis=2)
	order = ["abcd".index(label) for label in fig.layout.xaxis.ticktext]
	assert any(
		np.allclose(z, full[np.ix_(o, o)]) for o in (order, order[::-1])
	), "heatmap is not the distance matrix in leaf order"


def test_dendrogram_heatmap_follows_leaf_order_when_it_differs_from_input_order(tmp_path):
	"""With a and c close and b and d close, the leaves are not in input order; the heatmap is still the distance matrix in leaf order."""
	rows = [("a", "b", 4.0), ("a", "c", 1.0), ("a", "d", 4.5), ("b", "c", 3.5), ("b", "d", 1.5), ("c", "d", 4.2)]
	fig = figs.make_all_genome_dend(_dist_index(tmp_path, rows))
	order = ["abcd".index(label) for label in fig.layout.xaxis.ticktext]
	assert order != [0, 1, 2, 3] and order != [3, 2, 1, 0]
	m = _matrix(rows, list("abcd"))
	full = np.linalg.norm(m[:, None] - m[None], axis=2)
	z = np.array(fig.data[-1].z, dtype=float)
	assert any(np.allclose(z, full[np.ix_(o, o)]) for o in (order, order[::-1]))


def test_dendrogram_from_the_real_index_has_five_leaves(anchored_reader):
	"""With the anchored index (mash distances written by the workflow) the figure has five leaves and a 5x5 heatmap."""
	fig = figs.make_all_genome_dend(anchored_reader)
	assert sorted(fig.layout.xaxis.ticktext) == sorted(synthetic.GENOMES)
	z = np.array(fig.data[-1].z, dtype=float)
	assert z.shape == (5, 5)
	np.testing.assert_allclose(z, z.T)
	assert (np.diag(z) == 0).all()
	assert len(fig.data) % 2 == 1 and fig.data[-1].type == "heatmap"

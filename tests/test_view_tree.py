"""panagram.view: the dendrogram helpers (newick, coordinates, branch shapes, create_tree).

Expected values are worked by hand on tiny trees or taken from scipy and seaborn.
"""

from io import StringIO

import numpy as np
import pandas as pd
import pytest
import seaborn as sns
from Bio import Phylo
from scipy.cluster import hierarchy

from helpers_view import v, view_ns  # noqa: F401


def _phylo(newick):
	return Phylo.read(StringIO(newick), "newick")


def _leaf_depths(tree):
	d = tree.depths()
	return {t.name: d[t] for t in tree.get_terminals()}


def test_get_newick_three_leaves(v):
	"""Linkage [[0,1,1],[2,3,3]]: a,b join at height 1, c joins at 3; branch lengths are height differences."""
	z = np.array([[0, 1, 1.0, 2], [2, 3, 3.0, 3]])
	root = hierarchy.to_tree(z, False)
	text = v.get_newick(root, root.dist, ["a", "b", "c"])
	assert text.endswith(";")
	tree = _phylo(text)
	lengths = {t.name: t.branch_length for t in tree.get_terminals()}
	assert lengths == {"a": pytest.approx(1.0), "b": pytest.approx(1.0), "c": pytest.approx(3.0)}
	inner = [c for c in tree.root.clades if not c.is_terminal()]
	assert len(inner) == 1
	assert inner[0].branch_length == pytest.approx(2.0)
	assert sorted(t.name for t in inner[0].get_terminals()) == ["a", "b"]


def test_get_newick_four_leaves_two_subtrees(v):
	"""Linkage heights 1, 2 and 6: leaf lengths 1/1 and 2/2, internal lengths 5 and 4, every leaf 6 below the root."""
	z = np.array([[0, 1, 1.0, 2], [2, 3, 2.0, 2], [4, 5, 6.0, 4]])
	root = hierarchy.to_tree(z, False)
	tree = _phylo(v.get_newick(root, root.dist, ["a", "b", "c", "d"]))
	depths = _leaf_depths(tree)
	assert depths == {k: pytest.approx(6.0) for k in "abcd"}
	groups = {}
	for clade in tree.root.clades:
		groups[tuple(sorted(t.name for t in clade.get_terminals()))] = clade.branch_length
	assert groups == {("a", "b"): pytest.approx(5.0), ("c", "d"): pytest.approx(4.0)}
	by_name = {t.name: t.branch_length for t in tree.get_terminals()}
	assert by_name == {"a": pytest.approx(1.0), "b": pytest.approx(1.0), "c": pytest.approx(2.0), "d": pytest.approx(2.0)}


def test_get_x_coordinates_uses_branch_lengths(v):
	tree = _phylo("((A:1,B:2):0.5,C:4);")
	x = v.get_x_coordinates(tree)
	assert {t.name: x[t] for t in tree.get_terminals()} == {"A": 1.5, "B": 2.5, "C": 4.0}
	assert x[tree.root] == 0


def test_get_x_coordinates_without_lengths_falls_back_to_unit_branches(v):
	tree = _phylo("((A,B),C);")
	x = v.get_x_coordinates(tree)
	assert {t.name: x[t] for t in tree.get_terminals()} == {"A": 2, "B": 2, "C": 1}


def test_get_y_coordinates_three_leaves(v):
	"""Leaves are spaced by dist from the top: with 3 leaves A,B,C -> 0.4, 1.7, 3.0; parents sit midway between first and last child."""
	tree = _phylo("(A,B,C);")
	y = v.get_y_coordinates(tree)
	assert [y[t] for t in tree.get_terminals()] == pytest.approx([0.4, 1.7, 3.0])
	assert y[tree.root] == pytest.approx(1.7)


def test_get_y_coordinates_nested_and_custom_spacing(v):
	tree = _phylo("((A,B),C);")
	y = v.get_y_coordinates(tree, dist=1.0)
	names = {t.name: y[t] for t in tree.get_terminals()}
	assert names == {"A": pytest.approx(1.0), "B": pytest.approx(2.0), "C": pytest.approx(3.0)}
	inner = tree.root.clades[0]
	assert y[inner] == pytest.approx(1.5)
	assert y[tree.root] == pytest.approx(2.25)


def test_get_clade_lines_horizontal(v):
	line = v.get_clade_lines(
		orientation="horizontal", y_curr=2.5, start_coord=1, x_curr=4, line_color="red", line_width=3
	)
	assert (line["x0"], line["y0"], line["x1"], line["y1"]) == (1, 2.5, 4, 2.5)
	assert line["type"] == "line"
	assert line["line"] == {"color": "red", "width": 3}


def test_get_clade_lines_vertical(v):
	line = v.get_clade_lines(orientation="vertical", x_curr=7, y_bot=1, y_top=9)
	assert (line["x0"], line["y0"], line["x1"], line["y1"]) == (7, 1, 7, 9)
	assert line["line"]["color"] == "rgb(25,25,25)"
	assert line["line"]["width"] == 0.5


def test_get_clade_lines_rejects_unknown_orientation(v):
	with pytest.raises(ValueError):
		v.get_clade_lines(orientation="diagonal")


def test_biggest_num_in_clade(v):
	tree = _phylo("((A,B),C);")
	kmer_num = {"A": 3, "B": 9, "C": 5}
	inner = tree.root.clades[0]
	assert v.biggest_num_in_clade(tree.root, kmer_num) == 9
	assert v.biggest_num_in_clade(inner, kmer_num) == 9
	assert v.biggest_num_in_clade(tree.root.clades[1], kmer_num) == 5
	assert v.biggest_num_in_clade(inner.clades[0], kmer_num) == 3


def test_draw_clade_shapes(v):
	"""Worked tree ((A:1,B:1):1,C:2): 3 leaves + inner + root give 5 horizontal and 2 vertical segments."""
	tree = _phylo("((A:1,B:1):1,C:2);")
	x = v.get_x_coordinates(tree)
	y = v.get_y_coordinates(tree)  # A 0.4, B 1.7, C 3.0, inner 1.05, root 2.025
	palette = [f"p{i}" for i in range(130)]
	color_code = {"A": "#a", "B": "#b", "C": "#c"}
	kmer_num = {"A": 10, "B": 40, "C": 20}
	shapes = []
	v.draw_clade(color_code, 40, kmer_num, palette, tree.root, 0, shapes, x_coords=x, y_coords=y)

	def seg(s):
		return (s["x0"], s["y0"], s["x1"], s["y1"], s["line"]["color"])

	expected = [
		(0, 2.025, 0, 2.025, "p50"),  # root horizontal, colour of the largest leaf + 10
		(0, 3.0, 0, 1.05, "p50"),  # root vertical from last child C to first child inner
		(0, 1.05, 1, 1.05, "p50"),  # inner horizontal
		(1, 1.7, 1, 0.4, "p50"),  # inner vertical from B to A
		(1, 0.4, 2, 0.4, "#a"),
		(1, 1.7, 2, 1.7, "#b"),
		(0, 3.0, 2, 3.0, "#c"),
	]
	assert len(shapes) == len(expected)
	for got, want in zip(shapes, expected):
		assert seg(got)[:4] == pytest.approx(want[:4])
		assert seg(got)[4] == want[4]
	assert all(s["type"] == "line" for s in shapes)


# ---------------------------------------------------------------- create_tree


def _bitmap_with_twins(v):
	"""40 k-mers x 5 genomes in which alpha/beta and gamma/delta have identical columns."""
	i = np.arange(40)
	p1 = (i % 2).astype(int)
	p2 = (i % 3 == 0).astype(int)
	p3 = (i < 10).astype(int)
	return pd.DataFrame(
		{"alpha": p1, "beta": p1, "gamma": p2, "delta": p2, "epsilon": p3},
		columns=list(v.index.genome_names),
	)


def test_create_tree_percent_of_largest_genome(v):
	"""Column sums are 20, 20, 14, 14, 10 -> 100, 100, 70, 70, 50 percent of the largest."""
	bitmap = _bitmap_with_twins(v)
	shapes, annotations, node, labels, kmer_num = v.create_tree(bitmap)
	assert kmer_num == {
		"alpha": pytest.approx(100.0),
		"beta": pytest.approx(100.0),
		"gamma": pytest.approx(70.0),
		"delta": pytest.approx(70.0),
		"epsilon": pytest.approx(50.0),
	}


def test_create_tree_groups_identical_genomes(v):
	"""Genomes with identical columns are at distance 0, so they are adjacent leaves."""
	_, _, _, labels, _ = v.create_tree(_bitmap_with_twins(v))
	assert sorted(labels) == sorted(v.index.genome_names)
	pos = {name: i for i, name in enumerate(labels)}
	assert abs(pos["alpha"] - pos["beta"]) == 1
	assert abs(pos["gamma"] - pos["delta"]) == 1


def test_create_tree_leaf_depth_is_ward_root_height(v):
	"""All leaves sit at the root height of the Ward linkage of the columns (scipy), within newick rounding."""
	bitmap = _bitmap_with_twins(v)
	_, _, node, _, _ = v.create_tree(bitmap)
	root_height = hierarchy.linkage(bitmap.T.to_numpy(), method="ward")[-1, 2]
	leaf_x = [x for x, t in zip(node["x"], node["text"]) if t is not None]
	assert len(leaf_x) == 5
	assert leaf_x == pytest.approx([root_height] * 5, abs=0.05)
	assert max(node["x"]) == pytest.approx(root_height, abs=0.05)


def test_create_tree_shapes_count_and_colours(v):
	"""5 leaves: 9 horizontal segments (one per clade) and 4 vertical ones (one per internal clade)."""
	shapes, annotations, node, labels, kmer_num = v.create_tree(_bitmap_with_twins(v))
	horizontal = [s for s in shapes if s["y0"] == s["y1"]]
	vertical = [s for s in shapes if s["x0"] == s["x1"] and s["y0"] != s["y1"]]
	assert len(horizontal) == 9
	assert len(vertical) == 4
	palette = sns.color_palette("RdPu", 130).as_hex()
	leaf_colors = {t: c for t, c in zip(node["text"], node["marker"]["color"]) if t is not None}
	assert leaf_colors == {g: palette[int(p) + 10] for g, p in kmer_num.items()}
	assert node["marker"]["size"] == [0] * len(node["x"])


def test_create_tree_annotations_label_each_leaf(v):
	_, annotations, node, labels, kmer_num = v.create_tree(_bitmap_with_twins(v))
	assert len(annotations) == 5
	texts = {a["text"] for a in annotations}
	assert "\tepsilon - 50.0% " in texts
	assert "\tgamma - 70.0% " in texts
	assert all(a["xref"] == "x" and a["showarrow"] is False for a in annotations)
	by_label = {a["text"].split(" - ")[0].strip(): (a["x"], a["y"]) for a in annotations}
	leaf_xy = {t: (x, y) for t, x, y in zip(node["text"], node["x"], node["y"]) if t is not None}
	assert by_label == leaf_xy


def test_create_tree_on_synthetic_bitmap_matches_oracle(v, synthetic_genomes):
	"""Percentages come from the brute-force k-mer oracle of alpha's chromosomes."""
	g = synthetic_genomes
	bits = np.vstack([g.oracle_bits("alpha", c) for c in g.data.seqs["alpha"]])
	bitmap = pd.DataFrame(bits, columns=list(v.index.genome_names))
	_, _, _, labels, kmer_num = v.create_tree(bitmap)
	sums = bits.sum(axis=0)
	for name, s in zip(g.names, sums):
		assert kmer_num[name] == pytest.approx(100.0 * s / sums.max())
	assert sorted(labels) == sorted(g.names)

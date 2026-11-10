"""scripts/query.py and scripts/plot_umaps.py.

Both scripts run at import time, import `mycolorpy` (not a declared dependency; its release 1.5.3
calls a removed matplotlib function) and write figures with kaleido. The tests run the script
bodies in process with a small `mycolorpy` replacement, `Figure.write_image` recording the figure
instead of rendering it, and (for plot_umaps) `umap.UMAP` replaced by a recorder. Expected
values come from the brute-force k-mer oracle or are worked by hand.
"""

import io
import os

os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("dash")
pytest.importorskip("plotly")
pytest.importorskip("seaborn")
pytest.importorskip("pysam")

import plotly.graph_objects as go
pytest.importorskip("Bio")
from Bio import Phylo

from helpers_scripts import exec_script, install_mycolorpy_shim, oracle_by_chrom


@pytest.fixture(scope="module")
def oracle(synthetic_genomes):
	return oracle_by_chrom(synthetic_genomes, "alpha")


@pytest.fixture
def figures(monkeypatch):
	"""List that receives (figure, path) instead of rendering for every Figure.write_image call."""
	install_mycolorpy_shim(monkeypatch)
	captured = []
	monkeypatch.setattr(go.Figure, "write_image", lambda self, path, *a, **k: captured.append((self, str(path))))
	return captured


def run_query(monkeypatch, index_dir, chrom, start, end):
	return exec_script("query", ["alpha", chrom, start, end, index_dir], monkeypatch)[0]


# ------------------------------------------------------------------------ query.py


@pytest.mark.parametrize("start, end", [(0, 2636), (300, 1700)])
def test_query_figure_bars_and_heatmap_hold_the_oracle_bitmap(
	monkeypatch, figures, anchored_index, oracle, synthetic_genomes, start, end
):
	"""One bar trace per pan-count 0..5 and a genome x bin heatmap, from the step-100 bitmap rows of the region.

	The region is shorter than 20 kb so each bitmap row is its own bin; the anchor holds every one
	of its k-mers so each bin's maximum is 1 and the heatmap shows the 0/1 bits.
	"""
	run_query(monkeypatch, anchored_index, "chr1", start, end)
	assert [os.path.basename(p) for _, p in figures] == [
		f"alpha.chr1.{start}.{end}.svg",
		f"alpha.chr1.{start}.{end}.png",
	]
	fig = figures[0][0]
	positions = list(range(start, end, 100))
	bits = oracle["chr1"][positions]
	sums = bits.sum(axis=1)
	bars = {t.name: np.asarray(t.y) for t in fig.data if t.type == "bar"}
	assert sorted(bars) == ["0", "1", "2", "3", "4", "5"]
	for i in range(6):
		np.testing.assert_array_equal(bars[str(i)], (sums == i).astype(int), err_msg=f"pan-count {i}")
	(hm,) = [t for t in fig.data if t.type == "heatmap"]
	z = np.asarray(hm.z)
	genome_rows = list(hm.y)
	assert sorted(genome_rows) == sorted(synthetic_genomes.names)
	for r, name in enumerate(genome_rows):
		np.testing.assert_array_equal(z[r], bits[:, synthetic_genomes.names.index(name)])
	bin_size = (end - start) // 200 + 1
	assert list(hm.x) == [(p // bin_size) * bin_size for p in positions]
	assert fig.layout.xaxis2.title.text == f"Sequence position (chr1:{start}-{end})"


def test_query_figure_has_two_rows_with_titled_panels(monkeypatch, figures, anchored_index):
	run_query(monkeypatch, anchored_index, "chr2", 0, 1686)
	fig = figures[0][0]
	assert sorted(a.text for a in fig.layout.annotations) == [
		"Pan-k-mer containment",
		"Pan-k-mer multiplicity",
	]
	assert fig.layout.barmode == "stack"
	assert (fig.layout.width, fig.layout.height) == (2000, 2000)


@pytest.fixture
def qns(monkeypatch, figures, anchored_index):
	return run_query(monkeypatch, anchored_index, "chr1", 0, 2636)


def test_create_tree_reports_each_genomes_share_of_the_largest_count(qns, oracle, synthetic_genomes):
	"""kmer_num is each genome's k-mer count over the largest count, in percent; the tree has all five genomes and 13 drawn lines."""
	names = synthetic_genomes.names
	bitmap = pd.DataFrame(oracle["chr1"][::100], columns=names)
	shapes, annotations, node, legend, kmer_num = qns["create_tree"](bitmap)
	sums = bitmap.sum(axis=0)
	expected = sums / sums.max() * 100
	assert sorted(legend) == sorted(names)
	for n in names:
		assert kmer_num[n] == pytest.approx(expected[n])
	# 9 branches (5 tips + 4 internal nodes incl. the root) and 4 vertical joins
	assert len(shapes) == 13
	assert sum(1 for s in shapes if s["y0"] != s["y1"]) == 4
	text = {a["text"].split(" - ")[0].strip("\t"): a["text"] for a in annotations}
	assert sorted(text) == sorted(names)
	for n in names:
		assert text[n] == f"\t{n} - {str(kmer_num[n])[:4]}% "
	assert node["type"] == "scatter" and len(node["x"]) == 9


def test_create_tree_groups_genomes_with_similar_bitmaps(qns):
	"""Genomes a and b share every k-mer, c and d share every other one: the tree pairs a with b and c with d."""
	names = qns["index"].genome_names.tolist()
	base = np.array([1, 1, 1, 0, 0, 0] * 5)
	cols = {names[0]: base, names[1]: base, names[2]: 1 - base, names[3]: 1 - base, names[4]: base}
	_, _, _, legend, _ = qns["create_tree"](pd.DataFrame(cols))
	# neighbours in the drawing order are the pairs
	order = list(legend)
	pos = {n: order.index(n) for n in names}
	assert abs(pos[names[2]] - pos[names[3]]) == 1
	assert abs(pos[names[0]] - pos[names[1]]) in (1, 2)


def test_get_newick_in_query_matches_the_dendrogram(qns):
	from scipy.cluster import hierarchy
	from scipy.cluster.hierarchy import linkage

	z = linkage(np.array([[0.0], [1.0], [10.0]]), method="average")
	root = hierarchy.to_tree(z, False)
	text = qns["get_newick"](root, root.dist, ["A", "B", "C"])
	tree = Phylo.read(io.StringIO(text), "newick")
	assert tree.distance("A", "B") == pytest.approx(2.0)
	assert tree.distance("A", "C") == pytest.approx(19.0)


def test_tree_coordinates_by_hand(qns):
	tree = Phylo.read(io.StringIO("((A:1,B:1):1,C:2);"), "newick")
	x = qns["get_x_coordinates"](tree)
	byname = {(c.name or "internal"): v for c, v in x.items()}
	assert byname["A"] == byname["B"] == byname["C"] == 2
	ab = tree.common_ancestor("A", "B")
	assert x[ab] == 1 and x[tree.root] == 0
	y = qns["get_y_coordinates"](tree)
	assert y[tree.find_any("C")] == pytest.approx(3.0)
	assert y[tree.find_any("B")] == pytest.approx(1.7)
	assert y[tree.find_any("A")] == pytest.approx(0.4)
	assert y[ab] == pytest.approx((0.4 + 1.7) / 2)
	assert y[tree.root] == pytest.approx((y[ab] + 3.0) / 2)


def test_x_coordinates_fall_back_to_unit_branch_lengths(qns):
	"""A tree without branch lengths is drawn with every branch of length 1."""
	tree = Phylo.read(io.StringIO("((A,B),C);"), "newick")
	x = qns["get_x_coordinates"](tree)
	assert {c.name: v for c, v in x.items() if c.name} == {"A": 2, "B": 2, "C": 1}


def test_get_clade_lines(qns):
	h = qns["get_clade_lines"](orientation="horizontal", y_curr=3, start_coord=1, x_curr=5)
	assert (h["x0"], h["y0"], h["x1"], h["y1"]) == (1, 3, 5, 3)
	assert h["type"] == "line" and h["line"]["width"] == 0.5
	v = qns["get_clade_lines"](orientation="vertical", x_curr=5, y_bot=2, y_top=4, line_color="red", line_width=3)
	assert (v["x0"], v["y0"], v["x1"], v["y1"]) == (5, 2, 5, 4)
	assert v["line"] == {"color": "red", "width": 3}
	with pytest.raises(ValueError):
		qns["get_clade_lines"](orientation="diagonal")


def test_biggest_num_in_clade_takes_the_largest_tip_value(qns):
	tree = Phylo.read(io.StringIO("((A:1,B:1):1,C:2);"), "newick")
	num = {"A": 10, "B": 30, "C": 20}
	assert qns["biggest_num_in_clade"](tree.root, num) == 30
	assert qns["biggest_num_in_clade"](tree.common_ancestor("A", "B"), num) == 30
	assert qns["biggest_num_in_clade"](tree.find_any("C"), num) == 20


def test_set_color(qns):
	assert qns["SetColor"]("P") == "#ff006e"
	assert qns["SetColor"]("A") == "#3a86ff"
	assert qns["SetColor"]("M") == "#ffbe0b"
	assert qns["SetColor"]("Q") is None


def test_update_biodigs_orders_rows_by_reversed_tree_order(qns):
	df = pd.DataFrame({"Accession_rename": ["x", "y", "z", "w"], "val": [1, 2, 3, 4]})
	out = qns["update_biodigs"](df, "v", ["x", "y", "z"])
	assert list(out["Accession_rename"])[:3] == ["z", "y", "x"]
	assert list(out["Accession_rename"])[3] == "w"  # not in the tree: missing order sorts last
	assert out.set_index("Accession_rename").loc["x", "var_order"] == 1
	assert out.set_index("Accession_rename").loc["z", "var_order"] == 3
	assert "v_order_reverse" in out.columns


def test_read_sample_order_and_clade_order_read_files_from_the_working_directory(qns, tmp_path, monkeypatch):
	monkeypatch.chdir(tmp_path)
	(tmp_path / "sample_order.txt").write_text("s1\ns2\ns3\n")
	assert qns["read_sample_order"]() == ["s1", "s2", "s3"]
	(tmp_path / "aves_v2.txt").write_text("bird1\nbird2\n")
	(tmp_path / "primates.txt").write_text("ape1\n")
	x, y, z = qns["read_clade_order"](["ape1", "bird2", "cow", "bird1"])
	assert x == [1, 1, 1, 1]
	assert y == [0, 1, 2, 3]
	assert z == ["P", "A", "M", "A"]


def test_read_gff_clips_genes_and_exons_to_the_region(qns, tmp_path):
	"""GTF-style gene lines inside the region are clipped to it; the script expects `highlight` and `gff_file` to be set."""
	gff = tmp_path / "x.gtf"
	gff.write_text(
		'chr1\ts\tgene\t100\t200\t.\t+\t.\tgene_id "G1"; \n'
		'chr1\ts\tgene\t280\t400\t.\t+\t.\tgene_id "G2"; \n'
		'chr1\ts\tgene\t500\t600\t.\t+\t.\tgene_id "G3"; \n'
		'chr2\ts\tgene\t150\t250\t.\t+\t.\tgene_id "G4"; \n'
		'chr1\ts\texon\t90\t160\t.\t+\t.\tgene_id "G1"; \n'
	)
	qns["highlight"] = 0
	qns["gff_file"] = str(gff)
	genes, exons, hl = qns["read_gff"]("chr1", 150, 300)
	assert genes == [150, 200, None, 280, 300, None]
	assert exons == [150, 160, None]
	assert hl == []


def test_read_clusters_returns_columns_for_the_queried_chromosome(qns, tmp_path):
	f = tmp_path / "clusters.tsv"
	f.write_text("chr1\t100\tx\t3\nchr2\t200\tx\t4\nchr1\t300\tx\t5\n")
	assert qns["read_clusters"](f) == ([100, 300], [1, 1], [3, 5])


# -------------------------------------------------------------------- plot_umaps.py


@pytest.fixture
def umap_run(monkeypatch, figures, anchored_index, tmp_path):
	"""Run plot_umaps.py for chr1 in tmp_path with a recording UMAP; returns a function of (bin size, eps)."""
	umap = pytest.importorskip("umap")
	record = {}

	class Fake:
		def __init__(self, **kwargs):
			record["kwargs"] = kwargs

		def fit_transform(self, x):
			record["input"] = np.array(x)
			return np.column_stack([10.0 * np.arange(len(x)), np.zeros(len(x))])

	monkeypatch.setattr(umap, "UMAP", Fake)
	import panagram.index as pi

	# the script opens an index at a path on the author's machine
	real_index = pi.Index

	def open_anchored(path):
		cwd = os.getcwd()
		try:
			return real_index(str(anchored_index))
		finally:
			os.chdir(cwd)

	monkeypatch.setattr(pi, "Index", open_anchored)
	(tmp_path / "PCAS").mkdir()
	(tmp_path / "UMAPS").mkdir()
	monkeypatch.chdir(tmp_path)

	def run(bin_size, eps, neighbors=3, md=0.0, skips=100):
		exec_script("plot_umaps", ["alpha", "chr1", bin_size, eps, neighbors, md, skips], monkeypatch)
		return record

	return run


def test_plot_umaps_prints_bed_like_rows_and_writes_figures(umap_run, figures, capsys, tmp_path, oracle):
	"""One row per bin of chr1: chrom, start, start + bin - 1, cluster, UMAP coordinates; PCA png and UMAP html are written."""
	record = umap_run(500, 15)
	rows = [l.split("\t") for l in capsys.readouterr().out.splitlines()]
	assert len(rows) == 6  # positions 0..2600 in 500 bp bins
	for i, r in enumerate(rows):
		assert r[0] == "chr1"
		assert (int(r[1]), int(r[2])) == (500 * i, 500 * i + 499)
		assert int(r[3]) == 0  # points 10 apart, eps 15: one cluster
		assert float(r[4]) == 10.0 * i and float(r[5]) == 0.0
	assert record["kwargs"]["n_neighbors"] == 3
	assert record["kwargs"]["min_dist"] == 0.0
	# the matrix given to UMAP: bins x genomes, shared k-mers in the bin over the largest count
	m = oracle["chr1"][::100]
	bins = (np.arange(len(m)) * 100) // 500
	expected = np.array([m[bins == b].sum(axis=0) / m[bins == b].sum(axis=0).max() for b in range(6)])
	np.testing.assert_allclose(record["input"], expected, rtol=1e-12)
	assert (tmp_path / "PCAS" / "pca.alpha.chr1.bins500.genome.pos.png").read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
	html = tmp_path / "UMAPS" / "umap.cluster.alpha.chr1.bins500.neighbors3.md0.0.nskips.100.genome.pos.maxbinstest.eps15.0.html"
	assert "plotly" in html.read_text()
	fig, path = figures[-1]
	assert path.endswith("eps15.0.svg")
	np.testing.assert_allclose(np.asarray(fig.data[0].x), 10.0 * np.arange(6))


def test_plot_umaps_reads_the_chromosome_from_zero_to_its_size(umap_run, monkeypatch, capsys, oracle):
	"""The bitmap query covers k-mers 0 .. size of chr1 at the step argument."""
	import panagram.index as pi

	calls = []
	opener = pi.Index  # the fixture's opener of the synthetic index

	def open_spy(path):
		index = opener(path)
		real = index.query_bitmap

		def spy(*args, **kwargs):
			calls.append((args, kwargs))
			return real(*args, **kwargs)

		index.query_bitmap = spy
		return index

	monkeypatch.setattr(pi, "Index", open_spy)
	umap_run(500, 15, skips=100)
	assert calls == [(("alpha", "chr1", 0, len(oracle["chr1"]), 100), {})]


def test_plot_umaps_dbscan_eps_separates_clusters(umap_run, capsys):
	"""With eps 1 and points 10 apart every bin is its own cluster."""
	umap_run(500, 1)
	rows = [l.split("\t") for l in capsys.readouterr().out.splitlines()]
	assert [int(r[3]) for r in rows] == list(range(6))

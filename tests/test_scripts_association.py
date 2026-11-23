"""scripts/association_testing.py: functions and the command line.

Statistics are checked against scipy (pooled t test, Benjamini-Hochberg) or a numpy
generalised least squares written in the test. Haplotype groups are checked against a Ward
dendrogram of the oracle bitmap computed with scipy; windows come from the synthetic
genomes (tests/synthetic.py) and the brute-force k-mer oracle.
"""

import io
import os

os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
import pandas as pd
import pytest
from scipy import stats
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import linkage

if not hasattr(stats, "false_discovery_control"):
	pytest.skip("scipy >= 1.11 is needed for false_discovery_control", allow_module_level=True)
pytest.importorskip("statsmodels")
pytest.importorskip("Bio")
pytest.importorskip("seaborn")

from Bio import Phylo

from helpers_scripts import gene_rows, load_script, oracle_by_chrom, run_script

PHENO = {
	# alpha, beta, gamma, delta, epsilon
	"pheno_a": [0.5, 0.9, 10.2, 9.7, 10.4],
	"pheno_b": [3.1, 7.7, 2.2, 5.9, 4.4],
}


@pytest.fixture(scope="module")
def at():
	return load_script("association_testing")


@pytest.fixture(scope="module")
def oracle(synthetic_genomes):
	return oracle_by_chrom(synthetic_genomes, "alpha")


@pytest.fixture(scope="module")
def pheno_csv(synthetic_genomes, tmp_path_factory):
	path = tmp_path_factory.mktemp("pheno") / "phenotypes.csv"
	df = pd.DataFrame(PHENO, index=synthetic_genomes.names)
	df.index.name = "sample_id"
	df.to_csv(path)
	return path


@pytest.fixture(scope="module")
def tree_file(at, anchored_index, tmp_path_factory):
	path = tmp_path_factory.mktemp("tree") / "bitmap_tree.nwk"
	at.build_bitmap_tree("alpha", str(anchored_index), path)
	return path


# ------------------------------------------------------------ reference computations


def ref_root_split(matrix, names):
	"""The two groups under the root of a Ward dendrogram of the columns of `matrix` (rows x genomes)."""
	z = linkage(np.asarray(matrix, dtype=float).T, method="ward", metric="euclidean")
	n = len(names)

	def leaves(i):
		if i < n:
			return [names[i]]
		return leaves(int(z[i - n, 0])) + leaves(int(z[i - n, 1]))

	return frozenset(
		[frozenset(leaves(int(z[-1, 0]))), frozenset(leaves(int(z[-1, 1])))]
	)


def ref_gls_p(y, h, sigma):
	"""p value of the slope in y = a + b h + e, e ~ N(0, s2 sigma), by weighted least squares and a t test."""
	y, h = np.asarray(y, float), np.asarray(h, float)
	x = np.column_stack([np.ones(len(y)), h])
	w = np.linalg.inv(sigma)
	xtwx = x.T @ w @ x
	beta = np.linalg.solve(xtwx, x.T @ w @ y)
	resid = y - x @ beta
	dof = len(y) - 2
	s2 = resid @ w @ resid / dof
	se = np.sqrt(s2 * np.linalg.inv(xtwx)[1, 1])
	return 2 * stats.t.sf(abs(beta[1] / se), dof)


def ref_kinship(tree_path):
	"""(sample names in tree order, kinship) from root-to-tip depths of the Newick tree."""
	tree = Phylo.read(tree_path, "newick")
	tips = tree.get_terminals()
	depth = tree.depths()
	maxd = 2 * max(depth[t] for t in tips)
	k = np.zeros((len(tips), len(tips)))
	for i, a in enumerate(tips):
		for j, b in enumerate(tips):
			lca = tree.common_ancestor(a, b)
			d = depth[a] + depth[b] - 2 * depth[lca]
			k[i, j] = (maxd - d) / maxd
	k += np.eye(len(tips)) * 1e-6
	return [t.name for t in tips], k / np.mean(np.diag(k))


def ref_feature_rows(features, oracle, names, tree_path, pheno):
	"""(feature, phenotype) -> (partition, sorted [(size, mean)], p) for window features (name, chrom, start, end)."""
	order, kin = ref_kinship(tree_path)
	out = {}
	for name, chrom, start, end in features:
		part = ref_root_split(oracle[chrom][start:end], names)
		for col in pheno.columns:
			groups = sorted(part, key=lambda g: (len(g), sorted(g)))
			y = pheno.loc[order, col].to_numpy(float)
			h = np.array([0.0 if s in groups[0] else 1.0 for s in order])
			gs = [(len(g), pheno.loc[sorted(g), col].mean()) for g in groups]
			out[(name, col)] = (part, sorted(gs), ref_gls_p(y, h, kin))
	return out


def read_assoc(path):
	return pd.read_csv(path, sep="\t", header=None, dtype={0: str, 1: str}).set_axis(
		["feature_id", "phenotype_name", "group1_size", "group2_size", "group1_mean", "group2_mean", "p_value"],
		axis=1,
	)


def assert_rows_match_reference(df, ref):
	"""Every (feature, phenotype) row of `df` has the reference groups, group means and p value."""
	assert len(df) == len(ref)
	for row in df.itertuples():
		part, groups, p = ref[(row.feature_id, row.phenotype_name)]
		got = sorted([(row.group1_size, row.group1_mean), (row.group2_size, row.group2_mean)])
		assert [g[0] for g in got] == [g[0] for g in groups], row.feature_id
		np.testing.assert_allclose([g[1] for g in got], [g[1] for g in groups], rtol=1e-9)
		np.testing.assert_allclose(row.p_value, p, rtol=1e-4, atol=1e-12)


# ------------------------------------------------------------------------ newick, clades


def test_get_newick_reproduces_the_dendrogram(at):
	"""Average linkage of points 0, 1, 10: A and B join at 1, C joins at 9.5; tips are 9.5 from the root."""
	z = linkage(np.array([[0.0], [1.0], [10.0]]), method="average")
	root = hierarchy.to_tree(z, False)
	text = at.get_newick(root, root.dist, ["A", "B", "C"])
	tree = Phylo.read(io.StringIO(text), "newick")
	assert sorted(t.name for t in tree.get_terminals()) == ["A", "B", "C"]
	for tip in tree.get_terminals():
		assert tree.distance(tree.root, tip) == pytest.approx(9.5)
	ab = tree.common_ancestor("A", "B")
	assert sorted(t.name for t in ab.get_terminals()) == ["A", "B"]
	assert tree.distance("A", "B") == pytest.approx(2.0)
	assert text.endswith(");")


def test_get_major_clades_splits_the_root_children(at):
	tree = Phylo.read(io.StringIO("((A:1,B:1):1,C:2,(D:1,E:1,F:1):1);"), "newick")
	assert at.get_major_clades(tree, 1) == {
		"Clade_1": ["A", "B"],
		"Clade_2": ["C"],
		"Clade_3": ["D", "E", "F"],
	}
	# clades below min_size are skipped and the numbering continues without them
	assert at.get_major_clades(tree, 2) == {"Clade_1": ["A", "B"], "Clade_2": ["D", "E", "F"]}
	assert at.get_major_clades(tree, 4) == {}


@pytest.mark.parametrize(
	"attributes, expected",
	[("ID=g1;Name=n", "g1"), ('gene_id "ENS1"; x "y"', "ENS1"), ("Name=zz", "zz")],
)
def test_extract_feature_name(at, attributes, expected):
	assert at.extract_feature_name(attributes) == expected


def test_read_sample_order_skips_the_header(at, tmp_path):
	f = tmp_path / "samples.tsv"
	f.write_text("name\tfasta\tgff\nalpha\ta.fa\t\nbeta\tb.fa\tb.gff\ngamma\tc.fa\t\n")
	assert at.read_sample_order(f) == ["alpha", "beta", "gamma"]


# ----------------------------------------------------------------- window generation


def test_generate_regulatory_windows_by_hand(at):
	"""Gene 1000-2000, flank 500, chromosome 2300: windows start at 500, 600, ..., 2100 (17 of them)."""
	w = at.generate_regulatory_windows("g", 1000, 2000, 2300, flank_size=500)
	assert len(w) == 17
	assert w[0] == ("g_reg_1", 500, 700)
	assert w[-1] == ("g_reg_17", 2100, 2300)
	assert [x[1] for x in w] == list(range(500, 2101, 100))
	assert all(x[2] - x[1] == 200 for x in w)


def test_generate_regulatory_windows_clips_to_the_chromosome(at):
	w = at.generate_regulatory_windows("g", 100, 150, 400, flank_size=500, window_size=100, step_size=50)
	# region 0..400, windows of 100 every 50: starts 0..300
	assert [x[1] for x in w] == [0, 50, 100, 150, 200, 250, 300]
	assert w[-1][2] == 400


def test_generate_regulatory_windows_shorter_than_one_window(at):
	assert at.generate_regulatory_windows("g", 10, 20, 150) == []


def test_generate_chromosome_windows_by_hand(at):
	"""Chromosome of 450 with 200 bp windows every 100: starts 0, 100, 200."""
	assert at.generate_chromosome_windows("chrA", 450) == [
		("chrA_0", "chrA", 0, 200),
		("chrA_100", "chrA", 100, 300),
		("chrA_200", "chrA", 200, 400),
	]
	assert at.generate_chromosome_windows("c", 199) == []
	assert at.generate_chromosome_windows("c", 200) == [("c_0", "c", 0, 200)]
	assert at.generate_chromosome_windows("c", 10, window_size=4, step_size=3) == [
		("c_0", "c", 0, 4),
		("c_3", "c", 3, 7),
		("c_6", "c", 6, 10),
	]


def test_split_into_chunks(at):
	assert at.split_into_chunks([0, 1, 2, 3, 4], 2) == [[0, 1], [2, 3], [4]]
	assert at.split_into_chunks([0, 1], 5) == [[0, 1]]
	assert at.split_into_chunks([], 3) == []
	assert at.split_into_chunks(list(range(6)), 3) == [[0, 1, 2], [3, 4, 5]]


def test_build_feature_windows_all_covers_every_chromosome(at, anchored_index, synthetic_genomes):
	"""alpha chr1 has 2636 15-mers and chr2 1686: 25 and 15 windows of 200 bp every 100."""
	feats = at.build_feature_windows("alpha", str(anchored_index), None, "all")
	chrom_kmers = {
		c: len(s) - synthetic_genomes.k + 1 for c, s in synthetic_genomes.data.seqs["alpha"].items()
	}
	expected = [
		(f"{c}_{s}", c, s, s + 200)
		for c, n in chrom_kmers.items()
		for s in range(0, n - 200 + 1, 100)
	]
	assert feats == expected
	assert (len([f for f in feats if f[1] == "chr1"]), len([f for f in feats if f[1] == "chr2"])) == (25, 15)


def test_build_feature_windows_gene_uses_gff_names_chromosomes_and_ends(at, anchored_index, synthetic_genomes):
	"""Gene features carry the GFF name, chromosome and end. The start is not compared: the script passes the 1-based GFF start on unconverted."""
	feats = at.build_feature_windows("alpha", str(anchored_index), synthetic_genomes.gff["alpha"], "gene")
	assert [(n, c, e) for n, c, s, e in feats] == [(n, c, e) for c, s, e, n in gene_rows(synthetic_genomes, "alpha")]
	assert len(feats) == 6


def test_build_feature_windows_exon_names_are_numbered(at, anchored_index, synthetic_genomes, tmp_path):
	"""Repeated exon names get _1, _2, ... in file order."""
	gff = tmp_path / "x.gff"
	gff.write_text(
		"##gff-version 3\n"
		"chr1\ts\texon\t11\t20\t.\t+\t.\tID=e1\n"
		"chr1\ts\texon\t31\t40\t.\t+\t.\tID=e1\n"
		"chr1\ts\tgene\t1\t100\t.\t+\t.\tID=g\n"
		"chr2\ts\texon\t5\t9\t.\t+\t.\tID=e2\n"
		"chr1\ts\texon\t51\t60\t.\t+\t.\tID=e1\n"
	)
	feats = at.build_feature_windows("alpha", str(anchored_index), gff, "exon")
	# starts are not compared: the 1-based GFF start is passed on unconverted
	assert [(n, c, e) for n, c, s, e in feats] == [
		("e1_1", "chr1", 20),
		("e1_2", "chr1", 40),
		("e2_1", "chr2", 9),
		("e1_3", "chr1", 60),
	]


def test_build_feature_windows_reg_excludes_windows_overlapping_the_gene(at, anchored_index, synthetic_genomes):
	"""Regulatory windows of each gene: every 200 bp window within 500 kb (clipped to the chromosome) that misses the gene body."""
	feats = at.build_feature_windows("alpha", str(anchored_index), synthetic_genomes.gff["alpha"], "reg")
	size = {c: len(s) - synthetic_genomes.k + 1 for c, s in synthetic_genomes.data.seqs["alpha"].items()}
	expected = []
	for chrom, gs1, ge, name in gene_rows(synthetic_genomes, "alpha"):
		gs = gs1 - 1  # the 0-based start; the flank is far larger than the chromosome so lo is 0
		lo, hi = max(gs - 500000, 0), min(ge + 500000, size[chrom])
		for i, ws in enumerate(range(lo, hi - 200 + 1, 100), start=1):
			if not (ws < ge and ws + 200 > gs):
				expected.append((f"{name}_reg_{i}", chrom, ws, ws + 200))
	assert feats == expected
	assert len(feats) > 0
	assert all(not (f[2] < 900 and f[3] > 301) for f in feats if f[0].startswith("alpha_chr1_CORE_1_reg"))


# -------------------------------------------------------- haplotype groups, wide files


def test_create_feature_tree_groups_identical_columns(at):
	bitmap = pd.DataFrame(
		{
			"a": [1, 1, 1, 0, 0, 0],
			"b": [1, 1, 1, 0, 0, 0],
			"c": [0, 0, 0, 1, 1, 1],
			"d": [0, 0, 0, 1, 1, 1],
		}
	)
	clades = at.create_feature_tree(bitmap, ["a", "b", "c", "d"])
	assert sorted(map(sorted, clades.values())) == [["a", "b"], ["c", "d"]]
	assert sorted(clades) == ["Clade_1", "Clade_2"]


def test_create_feature_tree_isolates_an_outlier(at):
	"""Four similar columns and one that differs everywhere: the outlier is its own root clade."""
	rng = np.random.default_rng(3)
	base = rng.integers(0, 2, 40)
	cols = {n: base.copy() for n in "abcd"}
	cols["e"] = 1 - base
	clades = at.create_feature_tree(pd.DataFrame(cols), list("abcde"))
	assert sorted(map(sorted, clades.values())) == [["a", "b", "c", "d"], ["e"]]


def test_create_feature_tree_uses_ward_linkage(at):
	"""Five 8-row columns whose Ward root split ({a, b, e} | {c, d}) differs from the complete, single, average and centroid splits (checked with scipy)."""
	cols = {
		"a": [1, 0, 1, 0, 0, 1, 0, 1],
		"b": [1, 0, 1, 0, 1, 1, 1, 1],
		"c": [0, 0, 1, 1, 0, 1, 0, 0],
		"d": [0, 0, 0, 0, 0, 0, 0, 1],
		"e": [1, 1, 0, 0, 1, 0, 1, 1],
	}
	bitmap = pd.DataFrame(cols)
	for method in ("complete", "single", "average", "centroid"):
		z = linkage(bitmap.to_numpy(float).T, method=method, metric="euclidean")
		root = hierarchy.to_tree(z)
		other = {frozenset("abcde"[i] for i in root.left.pre_order()), frozenset("abcde"[i] for i in root.right.pre_order())}
		assert other != {frozenset("abe"), frozenset("cd")}, method
	clades = at.create_feature_tree(bitmap, list("abcde"))
	assert sorted(map(sorted, clades.values())) == [["a", "b", "e"], ["c", "d"]]


def test_make_haplotype_series_labels_clade_1_zero_and_clade_2_one(at):
	clades = {"Clade_1": ["a", "b"], "Clade_2": ["c"], "Clade_3": ["z"]}
	s = at.make_haplotype_series(clades, ["c", "a", "b", "d", "z"])
	assert s.name == "haplotype"
	assert s.to_dict() == {"c": 1, "a": 0, "b": 0}
	assert list(s.index) == ["c", "a", "b"]


def test_wide_haplotype_header_and_rows(at):
	out = io.StringIO()
	order = ["s1", "s2", "s3"]
	at.write_wide_haplotype_header(out, order)
	at.write_wide_haplotype_row(out, "f1", pd.Series({"s3": 1, "s1": 0, "s2": 1}), order)
	at.write_wide_haplotype_row(out, "f2", pd.Series({"s2": 0}), order)
	assert out.getvalue() == "feature_id\ts1\ts2\ts3\nf1\t0\t1\t1\nf2\tNA\t0\tNA\n"


# ------------------------------------------------------------------------- kinship


def test_build_kinship_matrix_by_hand(at, tmp_path):
	"""Tree ((A:1,B:1):1,C:2): path lengths AB=2, AC=BC=4; max root-to-tip 2 so the scale is 4.

	Similarity (4 - d) / 4 is 1 on the diagonal, 0.5 for A,B and 0 for C; a 1e-6 ridge is added
	and the matrix is divided by the mean diagonal.
	"""
	f = tmp_path / "t.nwk"
	f.write_text("((A:1,B:1):1,C:2);")
	kin, dist, max_samples = at.build_kinship_matrix(f)
	assert list(dist.index) == ["A", "B", "C"] == list(dist.columns)
	np.testing.assert_allclose(dist.to_numpy(), [[1, 0.5, 0], [0.5, 1, 0], [0, 0, 1]])
	s = 1 + 1e-6
	np.testing.assert_allclose(kin, [[1, 0.5 / s, 0], [0.5 / s, 1, 0], [0, 0, 1]], rtol=1e-12)
	assert max_samples == 6


# ----------------------------------------------------------------- one feature (GLS)


def _one_feature(at, hap, pheno_df, sigma, order, phenotypes=None):
	out = io.StringIO()
	dist = pd.DataFrame(np.zeros((len(order), len(order))), index=order, columns=order)
	at.test_one_feature_immediate(
		"feat", hap, sigma, dist, pheno_df, phenotypes or list(pheno_df.columns), out
	)
	rows = [line.split("\t") for line in out.getvalue().splitlines()]
	return rows


def test_one_feature_with_identity_kinship_is_the_pooled_t_test(at):
	"""sigma = I turns GLS into OLS, whose slope test equals the equal-variance two-sample t test."""
	rng = np.random.default_rng(1)
	order = [f"s{i}" for i in range(8)]
	y = rng.normal(size=8)
	pheno = pd.DataFrame({"p1": y}, index=order)
	hap = pd.Series([0, 0, 0, 0, 1, 1, 1, 1], index=order)
	(row,) = _one_feature(at, hap, pheno, np.eye(8), order)
	assert row[:4] == ["feat", "p1", "4", "4"]
	assert float(row[4]) == pytest.approx(y[:4].mean())
	assert float(row[5]) == pytest.approx(y[4:].mean())
	expected = stats.ttest_ind(y[:4], y[4:], equal_var=True).pvalue
	assert float(row[6]) == pytest.approx(expected, rel=1e-8)


def test_one_feature_with_correlated_kinship_matches_weighted_least_squares(at):
	rng = np.random.default_rng(2)
	order = ["a", "b", "c", "d", "e", "f", "g"]
	a = rng.normal(size=(7, 7))
	sigma = a @ a.T + 7 * np.eye(7)
	y = rng.normal(size=7) + np.array([0, 0, 0, 1, 1, 1, 1])
	pheno = pd.DataFrame({"p1": y}, index=order)
	h = np.array([0, 1, 0, 1, 1, 0, 1])
	hap = pd.Series(h, index=order)
	(row,) = _one_feature(at, hap, pheno, sigma, order)
	assert row[2:4] == ["3", "4"]
	assert float(row[4]) == pytest.approx(y[h == 0].mean())
	assert float(row[5]) == pytest.approx(y[h == 1].mean())
	assert float(row[6]) == pytest.approx(ref_gls_p(y, h, sigma), rel=1e-7)


def test_one_feature_aligns_samples_by_name_and_ignores_unshared_ones(at):
	"""Shuffled phenotype rows, an extra phenotype sample and a haplotype-less sample give the same p as sorted input."""
	rng = np.random.default_rng(5)
	order = ["a", "b", "c", "d", "e", "f"]
	a = rng.normal(size=(6, 6))
	sigma = a @ a.T + 6 * np.eye(6)
	y = rng.normal(size=6)
	h = np.array([0, 0, 1, 0, 1, 1])
	pheno = pd.DataFrame({"p1": y}, index=order)
	hap = pd.Series(h, index=order)
	(base,) = _one_feature(at, hap, pheno, sigma, order)

	shuffled = pd.concat([pheno.loc[["e", "a", "f", "c", "b", "d"]], pd.DataFrame({"p1": [99.0]}, index=["extra"])])
	hap2 = pd.concat([hap.loc[["d", "f", "b"]], hap.loc[["a", "c", "e"]]])
	(got,) = _one_feature(at, hap2, shuffled, sigma, order)
	assert got == base

	# a sample of the kinship order with a haplotype but no phenotype row is dropped
	dropped = pheno.drop(index="f")
	(sub,) = _one_feature(at, hap, dropped, sigma, order)
	keep = [0, 1, 2, 3, 4]
	assert sub[2:4] == ["3", "2"]
	assert float(sub[6]) == pytest.approx(
		ref_gls_p(y[keep], h[keep], sigma[np.ix_(keep, keep)]), rel=1e-7
	)


def test_one_feature_leaves_out_missing_phenotypes(at):
	"""A NaN phenotype is dropped from the fit and from the group means; other phenotypes are unaffected."""
	rng = np.random.default_rng(9)
	order = [f"s{i}" for i in range(9)]
	y1 = rng.normal(size=9)
	y2 = rng.normal(size=9)
	y1[3] = np.nan
	pheno = pd.DataFrame({"p1": y1, "p2": y2}, index=order)
	h = np.array([0, 0, 0, 0, 0, 1, 1, 1, 1])
	hap = pd.Series(h, index=order)
	rows = _one_feature(at, hap, pheno, np.eye(9), order)
	assert [r[1] for r in rows] == ["p1", "p2"]
	keep = np.arange(9) != 3
	assert float(rows[0][4]) == pytest.approx(np.nanmean(y1[:5]))
	assert float(rows[0][6]) == pytest.approx(
		stats.ttest_ind(y1[keep & (h == 0)], y1[h == 1]).pvalue, rel=1e-8
	)
	assert float(rows[1][6]) == pytest.approx(stats.ttest_ind(y2[h == 0], y2[h == 1]).pvalue, rel=1e-8)


def test_one_feature_writes_nothing_without_two_groups(at):
	order = ["a", "b", "c", "d"]
	pheno = pd.DataFrame({"p1": [1.0, 2.0, 3.0, 4.0]}, index=order)
	same = pd.Series([0, 0, 0, 0], index=order)
	assert _one_feature(at, same, pheno, np.eye(4), order) == []
	# haplotypes known for samples that have no phenotype row
	other = pd.Series([0, 1], index=["x", "y"])
	assert _one_feature(at, other, pheno, np.eye(4), order) == []
	# a phenotype that is missing for every sample of one group is skipped, the other is written
	pheno2 = pd.DataFrame({"p1": [np.nan, np.nan, 3.0, 4.0], "p2": [1.0, 2.0, 4.0, 8.0]}, index=order)
	hap = pd.Series([0, 0, 1, 1], index=order)
	rows = _one_feature(at, hap, pheno2, np.eye(4), order)
	assert [r[1] for r in rows] == ["p2"]
	allnan = pd.DataFrame({"p1": [np.nan] * 4, "p2": [1.0, 2.0, 4.0, 8.0]}, index=order)
	assert [r[1] for r in _one_feature(at, hap, allnan, np.eye(4), order)] == ["p2"]


# ------------------------------------------------- combining, filtering and BH by file


def _line(feature, pheno, g1, g2, m1, m2, p):
	return f"{feature}\t{pheno}\t{g1}\t{g2}\t{m1}\t{m2}\t{p}\n"


def test_combine_association_files_concatenates_chunks_of_the_feature_type(at, tmp_path):
	(tmp_path / "chunk_000002_gene_associations.tsv").write_text("b1\nb2\n")
	(tmp_path / "chunk_000001_gene_associations.tsv").write_text("a1\n")
	(tmp_path / "chunk_000001_exon_associations.tsv").write_text("exon-row\n")
	(tmp_path / "all_gene_associations.tsv").write_text("stale\n")
	got = at.combine_association_files(tmp_path, "gene")
	assert got == tmp_path / "all_gene_associations.tsv"
	assert got.read_text() == "a1\nb1\nb2\n"


def test_combine_association_files_for_the_all_type_skips_its_own_output(at, tmp_path):
	(tmp_path / "chunk_000001_all_associations.tsv").write_text("w1\n")
	(tmp_path / "chunk_000002_all_associations.tsv").write_text("w2\n")
	at.combine_association_files(tmp_path, "all")
	assert (tmp_path / "all_all_associations.tsv").read_text() == "w1\nw2\n"
	at.combine_association_files(tmp_path, "all")  # a second run must not read the first output
	assert (tmp_path / "all_all_associations.tsv").read_text() == "w1\nw2\n"


def test_filter_by_group_size_keeps_groups_of_at_least_ten_percent(at, tmp_path):
	"""1 of 10 is exactly 10% and stays; 0 of 10 and 9 of 100 go; p = 0 is kept unchanged."""
	src, dst = tmp_path / "in.tsv", tmp_path / "out.tsv"
	keep = [
		_line("f1", "p", 1, 9, 1.5, 2.5, 0.01),
		_line("f2", "p", 50, 50, 1.0, 2.0, 0.0),
		_line("f3", "p", 10, 90, 1.0, 2.0, 0.3),
		_line("f5", "p", 90, 10, 1.0, 2.0, 0.4),
	]
	drop = [
		_line("f4", "p", 0, 10, 1.0, 2.0, 0.2),
		_line("f6", "p", 9, 91, 1.0, 2.0, 0.2),
		_line("f7", "p", 91, 9, 1.0, 2.0, 0.2),
	]
	src.write_text("".join(keep[:2] + drop[:1] + keep[2:] + drop[1:]))
	assert at.filter_associations_by_group_size(src, dst) == 4
	assert dst.read_text() == "".join(keep)


def test_filter_by_group_size_exits_when_nothing_remains(at, tmp_path):
	src = tmp_path / "in.tsv"
	src.write_text(_line("f", "p", 0, 5, 1, 2, 0.5))
	with pytest.raises(SystemExit) as e:
		at.filter_associations_by_group_size(src, tmp_path / "out.tsv")
	assert "No rows remain" in str(e.value)


def test_count_lines(at, tmp_path):
	f = tmp_path / "f"
	f.write_text("a\nb\n\nc")
	assert at.count_lines(f) == 4


def test_run_sort_orders_by_key_and_keeps_all_lines(at, tmp_path, monkeypatch):
	monkeypatch.setenv("LC_ALL", "C")
	src, dst = tmp_path / "in.tsv", tmp_path / "out.tsv"
	src.write_text("x\t10\ny\t9\nz\t1e-5\nw\t0\n")
	at.run_sort(src, dst, ["-k2,2g"])
	assert dst.read_text() == "w\t0\nz\t1e-5\ny\t9\nx\t10\n"
	tmp = tmp_path / "scratch"
	tmp.mkdir()
	at.run_sort(src, dst, ["-k2,2gr"], temp_dir=tmp)
	assert dst.read_text() == "x\t10\ny\t9\nz\t1e-5\nw\t0\n"


def test_add_bh_rank_values_appends_rank_and_raw_bh(at, tmp_path):
	"""Raw BH is p * m / rank: 0.01*4/1, 0.02*4/2, 0.03*4/3, 0.5*4/4; a zero p stays zero."""
	src, dst = tmp_path / "s.tsv", tmp_path / "r.tsv"
	src.write_text(
		_line("f0", "p", 5, 5, 1, 2, 0.0)
		+ _line("f1", "p", 5, 5, 1, 2, 0.01)
		+ _line("f2", "p", 5, 5, 1, 2, 0.02)
		+ _line("f3", "p", 5, 5, 1, 2, 0.03)
		+ _line("f4", "p", 5, 5, 1, 2, 0.5)
	)
	at.add_bh_rank_values(src, dst, 5)
	rows = [l.split("\t") for l in dst.read_text().splitlines()]
	assert [r[7] for r in rows] == ["1", "2", "3", "4", "5"]
	np.testing.assert_allclose(
		[float(r[8]) for r in rows], [0.0, 0.01 * 5 / 2, 0.02 * 5 / 3, 0.03 * 5 / 4, 0.5], rtol=1e-12
	)
	assert [r[:7] for r in rows][1] == _line("f1", "p", 5, 5, 1, 2, 0.01).rstrip("\n").split("\t")


def test_apply_bh_cumulative_min_runs_from_the_largest_rank(at, tmp_path):
	"""Raw values from rank 4 down to 1: 3.0 0.5 0.6 0.04 give 1.0 0.5 0.5 0.04 (capped at 1, then running minimum)."""
	src, dst = tmp_path / "d.tsv", tmp_path / "f.tsv"
	raws = [3.0, 0.5, 0.6, 0.04]
	src.write_text(
		"".join(
			_line(f"f{i}", "p", 5, 5, 1, 2, 0.1 * i).rstrip("\n") + f"\t{4 - i}\t{r}\n"
			for i, r in enumerate(raws)
		)
	)
	at.apply_bh_cumulative_min(src, dst)
	rows = [l.split("\t") for l in dst.read_text().splitlines()]
	assert all(len(r) == 8 for r in rows)
	assert [float(r[7]) for r in rows] == [1.0, 0.5, 0.5, 0.04]
	assert [r[0] for r in rows] == ["f0", "f1", "f2", "f3"]


def test_collect_phenotype_mean_diff_ranges(at, tmp_path):
	f = tmp_path / "f.tsv"
	f.write_text(
		_line("a", "A", 5, 5, 3.0, 1.0, 0.1)  # +2
		+ _line("b", "A", 5, 5, 1.0, 2.0, 0.1)  # -1
		+ _line("c", "A", 5, 5, 4.0, 3.5, 0.1)  # +0.5
		+ _line("a", "B", 5, 5, 5.0, 5.0, 0.1)  # 0
	)
	assert at.collect_phenotype_mean_diff_ranges(f) == {"A": [-1.0, 2.0], "B": [0.0, 0.0]}


def test_collect_top_features_keeps_the_first_rows_per_phenotype(at, tmp_path):
	f = tmp_path / "f.tsv"
	f.write_text(
		_line("a", "A", 5, 5, 1, 2, 0.1)
		+ _line("b", "B", 5, 5, 1, 2, 0.1)
		+ _line("c", "A", 5, 5, 1, 2, 0.1)
		+ _line("d", "A", 5, 5, 1, 2, 0.1)
		+ _line("e", "B", 5, 5, 1, 2, 0.1)
	)
	assert at.collect_top_features_per_phenotype(f, max_points_per_phenotype=2) == {
		("a", "A"),
		("c", "A"),
		("b", "B"),
		("e", "B"),
	}
	assert len(at.collect_top_features_per_phenotype(f)) == 5


def test_add_norm_difference_and_plotted_flag(at, tmp_path):
	"""Difference of group means over the range of differences of that phenotype; 0 when the range is 0."""
	src, dst = tmp_path / "s.tsv", tmp_path / "d.tsv"
	rows = [
		("a", "A", 3.0, 1.0),  # +2
		("b", "A", 1.0, 2.0),  # -1
		("a", "B", 5.0, 5.0),  # 0
	]
	src.write_text("".join(_line(f, p, 5, 5, m1, m2, 0.1).rstrip("\n") + "\t0.25\n" for f, p, m1, m2 in rows))
	at.add_norm_difference_and_plotted_flag(
		src, dst, {"A": [-1.0, 2.0], "B": [0.0, 0.0]}, {("a", "A"), ("a", "B")}
	)
	lines = dst.read_text().splitlines()
	assert lines[0].split("\t") == [
		"feature_id", "phenotype_name", "group1_size", "group2_size", "group1_mean",
		"group2_mean", "p_value", "fdr", "norm_mean_difference", "plotted",
	]
	got = [l.split("\t") for l in lines[1:]]
	assert [float(g[8]) for g in got] == [pytest.approx(2 / 3), pytest.approx(-1 / 3), 0.0]
	assert [g[9] for g in got] == ["1", "0", "1"]


def _random_table(path, n_features=20):
	"""Association rows (two phenotypes) with a zero p value, a tie and some small groups."""
	rng = np.random.default_rng(11)
	lines, keep = [], []
	for i in range(n_features):
		for ph in ("pa", "pb"):
			p = float(rng.uniform())
			if i % 7 == 0:
				p *= 1e-4
			g1, g2 = (6, 4)
			if i == 3:
				g1, g2 = 0, 10  # dropped
			if i == 5:
				g1, g2 = 9, 91  # dropped
			if i == 6:
				g1, g2 = 1, 9  # exactly 10%, kept
			lines.append(_line(f"f{i:02d}", ph, g1, g2, float(i), float(i) * 0.5, repr(p)))
			keep.append((f"f{i:02d}", ph, p, (g1, g2)))
	lines[0] = _line("f00", "pa", 6, 4, 0.0, 0.0, "0.0")  # true zero p value
	keep[0] = ("f00", "pa", 0.0, (6, 4))
	# a tied p value on two rows
	lines[10] = _line("f05", "pa", 6, 4, 5.0, 2.5, "0.3")
	lines[12] = _line("f06", "pa", 6, 4, 6.0, 3.0, "0.3")
	keep[10] = ("f05", "pa", 0.3, (6, 4))
	keep[12] = ("f06", "pa", 0.3, (6, 4))
	path.write_text("".join(lines))
	return keep


def test_run_bh_fdr_disk_sorts_equal_p_values_by_feature_then_phenotype(at, tmp_path, monkeypatch):
	"""The p-value-sorted intermediate file breaks ties by feature id, then phenotype name (numeric p, so 1e-05 < 0.01)."""
	monkeypatch.setenv("LC_ALL", "C")
	combined = tmp_path / "all_gene_associations.tsv"
	combined.write_text(
		_line("fb", "pa", 6, 4, 1.0, 2.0, "0.01")
		+ _line("fa", "pb", 6, 4, 1.0, 2.0, "0.01")
		+ _line("fa", "pa", 6, 4, 1.0, 2.0, "0.01")
		+ _line("fz", "pz", 6, 4, 1.0, 2.0, "1e-05")
		+ _line("fb", "pb", 6, 4, 1.0, 2.0, "0.5")
	)
	at.run_bh_fdr_disk(combined, tmp_path, "gene")
	rows = [l.split("\t")[:2] for l in (tmp_path / "all_gene_associations.sorted_by_p.tmp.tsv").read_text().splitlines()]
	assert rows == [["fz", "pz"], ["fa", "pa"], ["fa", "pb"], ["fb", "pa"], ["fb", "pb"]]


def test_run_bh_fdr_disk_matches_scipy(at, tmp_path, monkeypatch):
	"""FDR over the rows left after the group-size filter equals scipy's Benjamini-Hochberg, zero p stays zero."""
	monkeypatch.setenv("LC_ALL", "C")
	combined = tmp_path / "all_gene_associations.tsv"
	rows = _random_table(combined)
	scratch = tmp_path / "sortdir"
	scratch.mkdir()
	final = at.run_bh_fdr_disk(combined, tmp_path, "gene", sort_temp_dir=scratch, max_points_per_phenotype=3)
	assert final == tmp_path / "all_gene_associations.fdr_sorted.tsv"
	df = pd.read_csv(final, sep="\t")
	assert list(df.columns) == [
		"feature_id", "phenotype_name", "group1_size", "group2_size", "group1_mean",
		"group2_mean", "p_value", "fdr", "norm_mean_difference", "plotted",
	]
	kept = [r for r in rows if min(r[3]) >= 0.1 * sum(r[3])]
	assert len(df) == len(kept) < len(rows)
	expected = dict(
		zip(
			[(r[0], r[1]) for r in kept],
			stats.false_discovery_control([r[2] for r in kept], method="bh"),
		)
	)
	got = {(r.feature_id, r.phenotype_name): r.fdr for r in df.itertuples()}
	assert got.keys() == expected.keys()
	for k in expected:
		assert got[k] == pytest.approx(expected[k], rel=1e-9, abs=0)
	assert df.loc[df["p_value"] == 0, "fdr"].tolist() == [0.0]
	assert (np.diff(df["fdr"].to_numpy()) >= 0).all()
	# top 3 FDR rows per phenotype are marked for plotting
	for ph, sub in df.groupby("phenotype_name"):
		assert sub["plotted"].sum() == 3
		assert (sub["plotted"].to_numpy()[:3] == 1).all()
	# normalised difference: (g1 - g2) over the range of (g1 - g2) in that phenotype
	diff = df["group1_mean"] - df["group2_mean"]
	rng_ = diff.groupby(df["phenotype_name"]).transform(lambda s: s.max() - s.min())
	np.testing.assert_allclose(df["norm_mean_difference"], diff / rng_, rtol=1e-9)


# ---------------------------------------------------------------------- volcano plot


def _final_table(path, rows):
	cols = [
		"feature_id", "phenotype_name", "group1_size", "group2_size", "group1_mean",
		"group2_mean", "p_value", "fdr", "norm_mean_difference", "plotted",
	]
	pd.DataFrame(rows, columns=cols).to_csv(path, sep="\t", index=False)


def test_plot_volcano_draws_plotted_rows_and_keeps_the_table(at, tmp_path, monkeypatch):
	"""Only plotted == 1 rows are drawn; fdr = 0 is drawn at least two orders above the smallest nonzero fdr; the TSV is untouched."""
	import matplotlib.pyplot as plt

	rows = [
		("f1", "ph_a", 5, 5, 2.0, 1.0, 0.0, 0.0, 0.5, 1),
		("f2", "ph_a", 5, 5, 1.0, 1.0, 1e-4, 1e-3, 0.0, 1),
		("f3", "ph_a", 5, 5, 0.5, 1.0, 0.2, 0.4, -0.2, 1),
		("f4", "ph_a", 5, 5, 0.5, 1.0, 0.2, 0.4, -0.2, 0),
		("f1", "ph_b", 5, 5, 2.0, 1.0, 1e-3, 1e-2, 0.5, 1),
		("f2", "ph_b", 5, 5, 3.0, 1.0, 0.5, 0.5, 0.9, 1),
	]
	table = tmp_path / "final.tsv"
	_final_table(table, rows)
	before = table.read_bytes()
	captured = []
	real = plt.savefig

	def spy(path, *a, **k):
		ax = plt.gca()
		pts = [np.asarray(c.get_offsets()) for c in ax.collections if len(c.get_offsets())]
		captured.append((os.path.basename(str(path)), ax.get_xlabel(), np.vstack(pts) if pts else np.zeros((0, 2))))
		return real(path, *a, **k)

	monkeypatch.setattr(plt, "savefig", spy)
	out = tmp_path / "plots" / "all_volcano.png"
	at.plot_volcano(table, out)
	assert table.read_bytes() == before
	names = [c[0] for c in captured]
	assert names[0] == "all_volcano.png"
	assert sorted(names[1:]) == ["ph_a_volcano.png", "ph_b_volcano.png"]
	for n in names:
		data = (out.parent / n).read_bytes()
		assert data[:8] == b"\x89PNG\r\n\x1a\n"
	combined = captured[0]
	assert combined[1] == "Normalized difference in means"
	pts = combined[2]
	assert len(pts) == 5  # f4 is not plotted
	nonzero_y = [-np.log10(f) for f in (1e-3, 0.4, 1e-2, 0.5)]
	star = [p[1] for p in pts if p[0] == 0.5 and p[1] > max(nonzero_y)]
	assert len(star) == 1
	assert star[0] >= -np.log10(1e-3) + 2 - 1e-9
	assert sorted(p[1] for p in pts if p[1] <= max(nonzero_y)) == pytest.approx(sorted(nonzero_y))


def test_plot_volcano_marks_fdr_005_with_a_dashed_line(at, tmp_path, monkeypatch):
	"""The combined plot and each phenotype plot carry a horizontal line at -log10(0.05)."""
	import matplotlib.pyplot as plt

	table = tmp_path / "final.tsv"
	_final_table(
		table,
		[
			("f1", "ph_a", 5, 5, 2.0, 1.0, 1e-4, 1e-3, 0.5, 1),
			("f2", "ph_b", 5, 5, 1.0, 2.0, 0.2, 0.4, -0.5, 1),
		],
	)
	lines = []
	real = plt.axhline

	def spy(*a, **k):
		lines.append((k.get("y", a[0] if a else None), k.get("linestyle")))
		return real(*a, **k)

	monkeypatch.setattr(plt, "axhline", spy)
	at.plot_volcano(table, tmp_path / "v.png")
	assert len(lines) == 3
	assert [y for y, _ in lines] == pytest.approx([-np.log10(0.05)] * 3)
	assert {style for _, style in lines} == {"--"}


def test_plot_volcano_without_plotted_rows_exits(at, tmp_path):
	table = tmp_path / "final.tsv"
	_final_table(table, [("f1", "ph", 5, 5, 2.0, 1.0, 0.1, 0.5, 0.5, 0)])
	with pytest.raises(SystemExit) as e:
		at.plot_volcano(table, tmp_path / "v.png")
	assert "plotted" in str(e.value)


# ----------------------------------------------------------------- tree and chunks


def test_build_bitmap_tree_matches_average_linkage_of_the_oracle(tree_file, oracle, synthetic_genomes):
	"""The tree clusters genomes by UPGMA on the step-100 bitmap rows of every chromosome."""
	names = synthetic_genomes.names
	matrix = np.vstack([oracle[c][::100] for c in oracle]).astype(float)
	z = linkage(matrix.T, method="average")
	root = hierarchy.to_tree(z, False)

	def clades(node):
		if node.is_leaf():
			return {frozenset([names[node.id]])}
		return {frozenset(names[i] for i in node.pre_order())} | clades(node.left) | clades(node.right)

	tree = Phylo.read(tree_file, "newick")
	got = set()
	for clade in tree.find_clades():
		got.add(frozenset(t.name for t in clade.get_terminals()))
	assert got == clades(root)
	# branch lengths are written with two decimals: root-to-tip distance is the root height
	heights = {tree.distance(tree.root, t) for t in tree.get_terminals()}
	assert max(heights) == pytest.approx(root.dist, abs=0.06)


def test_process_feature_chunk_writes_associations_and_haplotypes(
	at, anchored_index, synthetic_genomes, oracle, tree_file, pheno_csv, tmp_path
):
	feats = at.generate_chromosome_windows("chr1", 2636)[:6]
	hap_dir, assoc_dir = tmp_path / "hap", tmp_path / "assoc"
	hap_dir.mkdir()
	assoc_dir.mkdir()
	got = at.process_feature_chunk(
		"000007", feats, "alpha", str(anchored_index), str(synthetic_genomes.samples_tsv),
		str(pheno_csv), "all", str(tree_file), hap_dir, assoc_dir, True,
	)
	assert got == "000007"
	df = read_assoc(assoc_dir / "chunk_000007_all_associations.tsv")
	pheno = pd.read_csv(pheno_csv, index_col=0)
	ref = ref_feature_rows(feats, oracle, synthetic_genomes.names, tree_file, pheno)
	assert_rows_match_reference(df, ref)
	hap = pd.read_csv(hap_dir / "chunk_000007_all_haplotypes.tsv", sep="\t", index_col=0)
	assert list(hap.columns) == synthetic_genomes.names
	for name, row in hap.iterrows():
		groups = frozenset(
			frozenset(row.index[row == v]) for v in (0, 1) if (row == v).any()
		)
		assert groups == ref[(name, "pheno_a")][0]
	assert list(hap.index) == [f[0] for f in feats]


def test_process_feature_chunk_without_haplotypes_writes_no_haplotype_file(
	at, anchored_index, synthetic_genomes, tree_file, pheno_csv, tmp_path
):
	feats = at.generate_chromosome_windows("chr2", 1686)[:2]
	hap_dir, assoc_dir = tmp_path / "hap", tmp_path / "assoc"
	hap_dir.mkdir()
	assoc_dir.mkdir()
	at.process_feature_chunk(
		"1", feats, "alpha", str(anchored_index), str(synthetic_genomes.samples_tsv),
		str(pheno_csv), "all", str(tree_file), hap_dir, assoc_dir, False,
	)
	assert list(hap_dir.iterdir()) == []
	assert (assoc_dir / "chunk_1_all_associations.tsv").exists()


# ------------------------------------------------------------------- command line


def _run(env, tmp_path, out_name, anchored_index, genomes, pheno_csv, *extra, check=True):
	out = tmp_path / out_name
	r = run_script(
		env,
		"association_testing.py",
		"-a", "alpha",
		"-p", anchored_index,
		"-g", genomes.gff["alpha"],
		"-s", genomes.samples_tsv,
		"-m", pheno_csv,
		"-o", out,
		*extra,
		cwd=tmp_path,
		check=check,
	)
	return out, r


def test_cli_help_lists_every_option(panagram_overlay):
	r = run_script(panagram_overlay, "association_testing.py", "--help")
	for flag in (
		"--anchor-name", "--panagram-dir", "--gff", "--sample-file", "--phenotype-file", "--type",
		"--output-dir", "--threads", "--chunk-size", "--keep-haplotypes", "--sort-temp-dir",
		"--max-plot-points-per-phenotype",
	):
		assert flag in r.stdout
	for choice in ("gene", "exon", "reg", "all"):
		assert choice in r.stdout


def test_cli_rejects_an_unknown_type(panagram_overlay, anchored_index, synthetic_genomes, pheno_csv, tmp_path):
	_, r = _run(panagram_overlay, tmp_path, "o", anchored_index, synthetic_genomes, pheno_csv, "-t", "nope", "-n", "1", check=False)
	assert r.returncode == 2


@pytest.fixture(scope="module")
def all_run(panagram_overlay, anchored_index, synthetic_genomes, pheno_csv, tmp_path_factory):
	tmp = tmp_path_factory.mktemp("assoc_all")
	scratch = tmp / "sortscratch"
	scratch.mkdir()
	out, r = _run(
		panagram_overlay, tmp, "out", anchored_index, synthetic_genomes, pheno_csv,
		"-t", "all", "-n", "2", "--chunk-size", "15", "--keep-haplotypes",
		"--sort-temp-dir", scratch, "--max-plot-points-per-phenotype", "10",
	)
	return out, r


def test_cli_all_mode_output_layout(all_run, synthetic_genomes):
	out, r = all_run
	assert sorted(p.name for p in out.iterdir()) == ["associations", "haplotypes", "plots", "tree"]
	assoc = sorted(p.name for p in (out / "associations").glob("chunk_*"))
	# 40 windows in chunks of 15
	assert assoc == [f"chunk_00000{i}_all_associations.tsv" for i in (1, 2, 3)]
	assert sorted(p.name for p in (out / "haplotypes").iterdir()) == [
		f"chunk_00000{i}_all_haplotypes.tsv" for i in (1, 2, 3)
	]
	assert sorted(p.name for p in (out / "plots").iterdir()) == [
		"all_volcano.png", "pheno_a_volcano.png", "pheno_b_volcano.png",
	]
	for p in (out / "plots").iterdir():
		assert p.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
	assert "Total features/windows to test: 40" in r.stderr
	assert "Total chunks: 3" in r.stderr
	tree = Phylo.read(out / "tree" / "bitmap_tree.nwk", "newick")
	assert sorted(t.name for t in tree.get_terminals()) == sorted(synthetic_genomes.names)


def test_cli_all_mode_statistics_match_the_reference(all_run, oracle, synthetic_genomes, pheno_csv):
	out, _ = all_run
	final = pd.read_csv(out / "associations" / "all_all_associations.fdr_sorted.tsv", sep="\t")
	assert len(final) == 80  # 40 windows x 2 phenotypes
	size = {c: len(s) - synthetic_genomes.k + 1 for c, s in synthetic_genomes.data.seqs["alpha"].items()}
	feats = [(f"{c}_{s}", c, s, s + 200) for c, n in size.items() for s in range(0, n - 200 + 1, 100)]
	pheno = pd.read_csv(pheno_csv, index_col=0)
	ref = ref_feature_rows(feats, oracle, synthetic_genomes.names, out / "tree" / "bitmap_tree.nwk", pheno)
	assert_rows_match_reference(final, ref)
	# FDR is Benjamini-Hochberg over all 80 tests
	by = final.set_index(["feature_id", "phenotype_name"])
	expected = stats.false_discovery_control(by["p_value"].to_numpy(), method="bh")
	np.testing.assert_allclose(by["fdr"].to_numpy(), expected, rtol=1e-9)
	assert (np.diff(final["fdr"].to_numpy()) >= 0).all()
	assert final.groupby("phenotype_name")["plotted"].sum().to_dict() == {"pheno_a": 10, "pheno_b": 10}


def test_cli_all_mode_haplotype_files_hold_the_root_split(all_run, oracle, synthetic_genomes):
	out, _ = all_run
	frames = [pd.read_csv(p, sep="\t", index_col=0) for p in sorted((out / "haplotypes").iterdir())]
	hap = pd.concat(frames)
	assert list(hap.columns) == synthetic_genomes.names
	assert len(hap) == 40
	for name, row in hap.iterrows():
		chrom, start = name.rsplit("_", 1)
		start = int(start)
		part = ref_root_split(oracle[chrom][start : start + 200], synthetic_genomes.names)
		got = frozenset(frozenset(row.index[row == v]) for v in (0, 1))
		assert got == part, name
		assert set(row.unique()) == {0, 1}


def test_cli_gene_mode_without_keep_haplotypes(panagram_overlay, anchored_index, synthetic_genomes, pheno_csv, tmp_path):
	out, r = _run(
		panagram_overlay, tmp_path, "out", anchored_index, synthetic_genomes, pheno_csv,
		"-t", "gene", "-n", "1", "--chunk-size", "4",
	)
	assert list((out / "haplotypes").iterdir()) == []
	chunks = sorted(p.name for p in (out / "associations").glob("chunk_*"))
	assert chunks == ["chunk_000001_gene_associations.tsv", "chunk_000002_gene_associations.tsv"]
	final = pd.read_csv(out / "associations" / "all_gene_associations.fdr_sorted.tsv", sep="\t")
	names = {n for _, _, _, n in gene_rows(synthetic_genomes, "alpha")}
	assert set(final["feature_id"]) == names
	assert len(final) == 2 * len(names)
	assert ((final["group1_size"] + final["group2_size"]) == 5).all()
	assert ((final["p_value"] > 0) & (final["p_value"] <= 1)).all()
	expected = stats.false_discovery_control(final["p_value"].to_numpy(), method="bh")
	np.testing.assert_allclose(final["fdr"].to_numpy(), expected, rtol=1e-9)
	# the default of 50000 plotted rows per phenotype marks every row here
	assert (final["plotted"] == 1).all()
	assert (out / "plots" / "gene_volcano.png").exists()


def test_cli_exits_when_the_gff_has_no_features_of_the_type(
	panagram_overlay, anchored_index, synthetic_genomes, pheno_csv, tmp_path
):
	gff = tmp_path / "exons_only.gff"
	gff.write_text("##gff-version 3\nchr1\ts\texon\t1\t50\t.\t+\t.\tID=e\n")
	out = tmp_path / "out"
	r = run_script(
		panagram_overlay, "association_testing.py", "-a", "alpha", "-p", anchored_index, "-g", gff,
		"-s", synthetic_genomes.samples_tsv, "-m", pheno_csv, "-t", "gene", "-o", out, "-n", "1",
		cwd=tmp_path, check=False,
	)
	assert r.returncode == 1
	assert "No features/windows were found" in r.stderr

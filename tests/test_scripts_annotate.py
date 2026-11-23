"""scripts/annotate_significant_windows.py: functions and command line.

Expected values are hand-worked on small gene tables, or taken from the coordinates of the
synthetic GFF3 files (tests/synthetic.py).
"""

import pandas as pd
import pytest

from helpers_scripts import gene_rows, load_script, run_script

GFF_HEADER = "##gff-version 3\n"


@pytest.fixture(scope="module")
def asw():
	return load_script("annotate_significant_windows")


def _genes(rows):
	"""{chrom: DataFrame} the way read_genes_from_gff builds it, from (chrom, start, end, name)."""
	df = pd.DataFrame(rows, columns=["chrom", "start", "end", "gene_name"])
	return {
		c: g.sort_values(["start", "end"]).reset_index(drop=True) for c, g in df.groupby("chrom")
	}


def _write_gff(path, rows, extra=""):
	with open(path, "w") as f:
		f.write(GFF_HEADER)
		for chrom, start, end, name in rows:
			f.write(f"{chrom}\tsrc\tgene\t{start}\t{end}\t.\t+\t.\tID={name};Name=x\n")
		f.write(extra)


# ------------------------------------------------------------ extract_feature_name


@pytest.mark.parametrize(
	"attributes, expected",
	[
		("ID=geneA;Name=other", "geneA"),
		("Name=first;ID=second", "first"),
		('gene_id "G1"; transcript_id "T1"', "G1"),
		("ID=a=b;Name=c", "a"),
	],
)
def test_extract_feature_name_uses_the_first_attribute(asw, attributes, expected):
	"""The value of the first attribute is the name; quotes (GTF style) are stripped."""
	assert asw.extract_feature_name(attributes) == expected


# ---------------------------------------------------------------- read_genes_from_gff


def test_read_genes_from_gff_groups_sorts_and_filters(asw, tmp_path):
	"""Only `gene` lines of 9 or more fields count; chromosomes are separate tables sorted by start then end."""
	gff = tmp_path / "x.gff"
	rows = [
		("chrB", 500, 600, "b2"),
		("chrA", 300, 400, "a2"),
		("chrA", 100, 250, "a1"),
		("chrA", 100, 200, "a0"),
	]
	_write_gff(
		gff,
		rows,
		extra="chrA\tsrc\texon\t1\t50\t.\t+\t.\tID=ex1\nchrA\tsrc\tgene\t1\t2\n#chrA\tsrc\tgene\t7\t8\t.\t+\t.\tID=cm\n",
	)
	got = asw.read_genes_from_gff(gff)
	assert sorted(got) == ["chrA", "chrB"]
	a = got["chrA"]
	assert list(a["gene_name"]) == ["a0", "a1", "a2"]
	assert list(a["start"]) == [100, 100, 300]
	assert list(a["end"]) == [200, 250, 400]
	assert list(got["chrB"]["gene_name"]) == ["b2"]


def test_read_genes_from_gff_sorts_by_start_before_end(asw, tmp_path):
	"""A long gene that starts first precedes a short nested one, though it ends later; equal starts order by end."""
	gff = tmp_path / "nested.gff"
	_write_gff(gff, [("c", 200, 300, "inner"), ("c", 100, 500, "outer"), ("c", 200, 250, "inner_short")])
	got = asw.read_genes_from_gff(gff)["c"]
	assert got["gene_name"].tolist() == ["outer", "inner_short", "inner"]
	assert got["start"].tolist() == [100, 200, 200]


def test_annotate_region_lists_intersecting_genes_in_start_order(asw):
	"""Genes overlapping the region are joined in coordinate order (start, then end)."""
	rows = [("c", 200, 300, "inner"), ("c", 100, 500, "outer")]
	got = asw.annotate_region("c", 250, 260, _genes(rows))
	assert got["intersecting_genes"] == "outer,inner"


def test_read_genes_from_gff_without_genes_raises(asw, tmp_path):
	gff = tmp_path / "none.gff"
	gff.write_text(GFF_HEADER + "chr1\tsrc\texon\t1\t10\t.\t+\t.\tID=e1\n")
	with pytest.raises(ValueError, match="gene"):
		asw.read_genes_from_gff(gff)


def test_read_genes_from_synthetic_gff(asw, synthetic_genomes):
	"""The synthetic alpha GFF3 holds four genes on chr1 and two on chr2, at the layout's coordinates."""
	got = asw.read_genes_from_gff(synthetic_genomes.gff["alpha"])
	expected = gene_rows(synthetic_genomes, "alpha")
	for chrom in ("chr1", "chr2"):
		want = sorted((s, e, n) for c, s, e, n in expected if c == chrom)
		df = got[chrom]
		assert list(zip(df["start"], df["end"], df["gene_name"])) == want
	assert (len(got["chr1"]), len(got["chr2"])) == (4, 2)


# ----------------------------------------------------- parse_all_window_feature_id


def test_parse_window_id_splits_on_the_last_underscore(asw):
	assert asw.parse_all_window_feature_id("chr1_100", 200) == ("chr1", 100, 300)
	assert asw.parse_all_window_feature_id("chr1_RagTag_123400", 200) == ("chr1_RagTag", 123400, 123600)
	assert asw.parse_all_window_feature_id("scaffold_7_0", 50) == ("scaffold_7", 0, 50)


def test_parse_window_id_rejects_non_numeric_start(asw):
	with pytest.raises(ValueError):
		asw.parse_all_window_feature_id("chr1_abc", 200)


# ------------------------------------------------------------------ annotate_region

BLANK = {
	"intersecting_genes": "",
	"upstream_gene": "",
	"upstream_distance": "",
	"downstream_gene": "",
	"downstream_distance": "",
}

# genes A [100, 200] and B [300, 400] on chrX
GENES = [("chrX", 100, 200, "A"), ("chrX", 300, 400, "B")]


def _annot(asw, start, end, rows=GENES, chrom="chrX"):
	return asw.annotate_region(chrom, start, end, _genes(rows))


def test_annotate_region_reports_intersecting_genes_only(asw):
	"""An overlap gives the gene names, comma separated in start order, and blank flank fields."""
	assert _annot(asw, 150, 250) == {**BLANK, "intersecting_genes": "A"}
	assert _annot(asw, 190, 310) == {**BLANK, "intersecting_genes": "A,B"}
	assert _annot(asw, 0, 1000) == {**BLANK, "intersecting_genes": "A,B"}


def test_annotate_region_between_genes_gives_both_flanks(asw):
	"""Region 250-260: upstream A ends at 200 (distance 50); downstream B starts at 300 (distance 40)."""
	assert _annot(asw, 250, 260) == {
		**BLANK,
		"upstream_gene": "A",
		"upstream_distance": "50",
		"downstream_gene": "B",
		"downstream_distance": "40",
	}


def test_annotate_region_touching_a_gene_is_not_intersecting(asw):
	"""Half-open overlap: a region that ends where a gene starts, or starts where one ends, does not intersect."""
	got = _annot(asw, 200, 300)
	assert got["intersecting_genes"] == ""
	assert (got["upstream_gene"], got["upstream_distance"]) == ("A", "0")
	assert (got["downstream_gene"], got["downstream_distance"]) == ("B", "0")


def test_annotate_region_outside_all_genes(asw):
	before = _annot(asw, 0, 50)
	assert before == {**BLANK, "downstream_gene": "A", "downstream_distance": "50"}
	after = _annot(asw, 500, 600)
	assert after == {**BLANK, "upstream_gene": "B", "upstream_distance": "100"}


def test_annotate_region_picks_the_nearest_flank(asw):
	rows = [("c", 10, 20, "far"), ("c", 50, 60, "near"), ("c", 200, 220, "d1"), ("c", 300, 320, "d2")]
	got = _annot(asw, 100, 120, rows=rows, chrom="c")
	assert (got["upstream_gene"], got["upstream_distance"]) == ("near", "40")
	assert (got["downstream_gene"], got["downstream_distance"]) == ("d1", "80")


def test_annotate_region_unknown_chromosome_is_blank(asw):
	assert _annot(asw, 100, 200, chrom="chrY") == BLANK


# ------------------------------------------------- annotate_significant_associations

ASSOC_COLUMNS = [
	"feature_id",
	"phenotype_name",
	"group1_size",
	"group2_size",
	"group1_mean",
	"group2_mean",
	"p_value",
	"fdr",
	"norm_mean_difference",
	"plotted",
]
OUT_COLUMNS = [
	"feature_id",
	"phenotype_name",
	"chrom",
	"start",
	"end",
	"p_value",
	"fdr",
	"intersecting_genes",
	"upstream_gene",
	"upstream_distance",
	"downstream_gene",
	"downstream_distance",
]

# (feature_id, phenotype, p_value, fdr)
ASSOC_ROWS = [
	("chrX_150", "ph1", 0.0, 0.0),
	("chrX_250", "ph1", 1e-6, 0.01),
	("chrX_250", "ph2", 1e-4, 0.05),
	("chrX_900", "ph1", 0.2, 0.0500001),
	("chrX_300", "ph2", 0.5, 0.9),
	("chrY_10", "ph1", 1e-5, 0.02),
]


def _write_assoc(path, rows=ASSOC_ROWS):
	df = pd.DataFrame(
		[(f, p, 3, 2, 1.0, 2.0, pv, fd, 0.1, 1) for f, p, pv, fd in rows], columns=ASSOC_COLUMNS
	)
	df.to_csv(path, sep="\t", index=False)


def _read_out(path):
	return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def test_annotate_significant_associations_end_to_end(asw, tmp_path):
	"""Rows with fdr <= 0.05 are kept (the cutoff itself is kept) and annotated with 200 bp windows."""
	gff = tmp_path / "genes.gff"
	_write_gff(gff, GENES)
	assoc = tmp_path / "assoc.tsv"
	_write_assoc(assoc)
	out = tmp_path / "deep" / "dir" / "out.tsv"
	asw.annotate_significant_associations(assoc, gff, out)
	got = _read_out(out)
	assert list(got.columns) == OUT_COLUMNS
	assert list(zip(got["feature_id"], got["phenotype_name"])) == [
		("chrX_150", "ph1"),
		("chrX_250", "ph1"),
		("chrX_250", "ph2"),
		("chrY_10", "ph1"),
	]
	by = got.set_index(["feature_id", "phenotype_name"])
	# window chrX_150 = [150, 350) intersects A and B
	assert by.loc[("chrX_150", "ph1"), "intersecting_genes"] == "A,B"
	assert (by.loc[("chrX_150", "ph1"), "chrom"], by.loc[("chrX_150", "ph1"), "start"]) == ("chrX", "150")
	assert by.loc[("chrX_150", "ph1"), "end"] == "350"
	# window chrX_250 = [250, 450) intersects B only
	assert by.loc[("chrX_250", "ph2"), "intersecting_genes"] == "B"
	# chrY has no genes: blank annotation
	row = by.loc[("chrY_10", "ph1")]
	assert row["intersecting_genes"] == row["upstream_gene"] == row["downstream_gene"] == ""
	assert float(by.loc[("chrX_250", "ph2"), "p_value"]) == 1e-4
	assert float(by.loc[("chrX_250", "ph2"), "fdr"]) == 0.05


def test_annotate_significant_associations_options(asw, tmp_path):
	"""--fdr-cutoff, --window-size and --chunksize change the selection, window end and chunking."""
	gff = tmp_path / "genes.gff"
	_write_gff(gff, GENES)
	assoc = tmp_path / "assoc.tsv"
	_write_assoc(assoc)
	out = tmp_path / "out.tsv"
	asw.annotate_significant_associations(
		assoc, gff, out, fdr_cutoff=0.015, window_size=40, chunksize=1
	)
	got = _read_out(out)
	assert list(got["feature_id"]) == ["chrX_150", "chrX_250"]
	by = got.set_index("feature_id")
	# window chrX_250 of 40 bp = [250, 290): between A (ends 200) and B (starts 300)
	assert by.loc["chrX_250", "end"] == "290"
	assert (by.loc["chrX_250", "upstream_gene"], by.loc["chrX_250", "upstream_distance"]) == ("A", "50")
	assert (by.loc["chrX_250", "downstream_gene"], by.loc["chrX_250", "downstream_distance"]) == ("B", "10")
	assert by.loc["chrX_150", "intersecting_genes"] == "A"


def test_annotate_significant_associations_nothing_significant_writes_header_only(asw, tmp_path):
	gff = tmp_path / "genes.gff"
	_write_gff(gff, GENES)
	assoc = tmp_path / "assoc.tsv"
	_write_assoc(assoc, [("chrX_1", "ph1", 0.3, 0.8)])
	out = tmp_path / "out.tsv"
	asw.annotate_significant_associations(assoc, gff, out)
	assert out.read_text() == "\t".join(OUT_COLUMNS) + "\n"


def test_annotate_synthetic_windows_against_layout(asw, synthetic_genomes, tmp_path):
	"""Windows of the alpha chr1 layout: inside CORE, and between the CORE and PAIR genes (nearest flanks named)."""
	g = synthetic_genomes
	genes = {n: (s, e) for c, s, e, n in gene_rows(g, "alpha") if c == "chr1"}
	core_s, core_e = genes["alpha_chr1_CORE_1"]
	pair_s, pair_e = genes["alpha_chr1_PAIR_3"]
	assoc = tmp_path / "assoc.tsv"
	_write_assoc(
		assoc,
		[
			("chr1_400", "ph1", 1e-9, 0.001),
			("chr1_950", "ph1", 1e-9, 0.001),
		],
	)
	out = tmp_path / "out.tsv"
	asw.annotate_significant_associations(assoc, g.gff["alpha"], out)
	by = _read_out(out).set_index("feature_id")
	assert by.loc["chr1_400", "intersecting_genes"] == "alpha_chr1_CORE_1"
	row = by.loc["chr1_950"]
	assert row["intersecting_genes"] == ""
	assert (row["upstream_gene"], int(row["upstream_distance"])) == ("alpha_chr1_CORE_1", 950 - core_e)
	# the distance is not compared: the downstream gene start is the 1-based GFF start, unconverted
	assert row["downstream_gene"] == "alpha_chr1_PAIR_3"
	assert int(row["downstream_distance"]) in (pair_s - 1150, pair_s - 1 - 1150)


# -------------------------------------------------------------------- command line


def test_cli_defaults_and_options(panagram_overlay, tmp_path):
	"""The command line: defaults are fdr 0.05 and 200 bp windows; the options override them."""
	gff = tmp_path / "genes.gff"
	_write_gff(gff, GENES)
	assoc = tmp_path / "assoc.tsv"
	_write_assoc(assoc)
	out = tmp_path / "out.tsv"
	run_script(panagram_overlay, "annotate_significant_windows.py", "-a", assoc, "-g", gff, "-o", out)
	default = _read_out(out)
	assert len(default) == 4
	assert set(default["end"].astype(int) - default["start"].astype(int)) == {200}

	out2 = tmp_path / "out2.tsv"
	run_script(
		panagram_overlay,
		"annotate_significant_windows.py",
		"--association-file", assoc,
		"--gff", gff,
		"--output-file", out2,
		"--fdr-cutoff", "0.015",
		"--window-size", "40",
		"--chunksize", "2",
	)
	custom = _read_out(out2)
	assert list(custom["feature_id"]) == ["chrX_150", "chrX_250"]
	assert set(custom["end"].astype(int) - custom["start"].astype(int)) == {40}


def test_cli_requires_its_three_arguments(panagram_overlay, tmp_path):
	r = run_script(panagram_overlay, "annotate_significant_windows.py", "-a", "x", check=False)
	assert r.returncode == 2
	assert "--gff" in r.stderr and "--output-file" in r.stderr

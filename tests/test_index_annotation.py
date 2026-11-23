"""GFF parsing, annotation tables (tabix), anno types, run_annotate and the annotation queries of panagram.index.Genome.

GFF files are written by the tests with hand-known contents. Coordinates of the stored gene/annotation tables are
read back from the tables themselves where the GFF-to-table conversion of `start` is not under test.
"""

import os
import shutil

import numpy as np
import pandas as pd
import pysam
import pytest

import helpers_index as H
import synthetic
from conftest import open_index
from panagram import index as pidx

GFF = """##gff-version 3
# a comment line
chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=g1;Name=GeneOne
chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=t1;Parent=g1
chr1\tsrc\texon\t1\t50\t.\t+\t.\tID=e1;Parent=t1
chr1\tsrc\texon\t1\t50\t.\t+\t.\tID=e1b;Parent=t1
chr1\tsrc\tCDS\t10\t40\t.\t+\t.\tID=c1;Parent=t1
chr1\tsrc\ttranscript\t1\t100\t.\t+\t.\tID=t9;Parent=g1
chr1\tsrc\tgene\t200\t300\t.\t-\t.\tID=g2
chr1\tsrc\texon\t200\t250\t.\t-\t.\tID=e2;Parent=g2
chr1\tsrc\tregion\t1\t1000\t.\t.\t.\tID=r1;Name=Region
chr2\tsrc\tgene\t5\t50\t.\t+\t.\tID=g3;Name=Three
"""


@pytest.fixture
def gff_file(tmp_path):
	p = tmp_path / "x.gff"
	p.write_text(GFF)
	return p


def _gff_genome(tmp_path, gff, **params):
	return H.make_genome(tmp_path / "w", gff=str(gff), **params)


def _anno_rows(g, typ="anno"):
	"""All rows of a written tabix table as tuples of strings."""
	with pysam.TabixFile(g.tabix_fname(typ), parser=pysam.asTuple(), index=g.tabix_idx_fname(typ)) as t:
		return [tuple(r) for r in t.fetch()]


# ---------------------------------------------------------------------- attributes


def test_gffattr_extracts_named_attribute_case_insensitively(tmp_path):
	"""_gffattr pulls NAME=value up to the next semicolon, ignoring case, and fills gaps from another column."""
	g = H.make_genome(tmp_path)
	df = pd.DataFrame(
		{
			"attr": ["ID=a;Name=Alpha", "id=b;name=beta;x=1", "ID=c", "Name=only"],
			"id": ["a", "b", "c", None],
		}
	)
	assert g._gffattr(df, "Name").tolist()[:2] == ["Alpha", "beta"]
	assert pd.isna(g._gffattr(df, "Name")[2])
	assert g._gffattr(df, "ID").tolist()[:3] == ["a", "b", "c"]
	assert pd.isna(g._gffattr(df, "ID")[3])
	filled = g._gffattr(df, "Name", fill="id")
	assert filled.tolist()[:3] == ["Alpha", "beta", "c"]
	assert filled[3] == "only"


def test_iter_gff_skips_comments_and_adds_ids(tmp_path, gff_file):
	"""_iter_gff yields frames with chr, end, type, attr and the parsed ID; comment lines are skipped."""
	g = _gff_genome(tmp_path, gff_file)
	frames = list(g._iter_gff())
	assert len(frames) == 1
	df = frames[0]
	assert set(df.columns) == {"chr", "start", "end", "type", "attr", "id"}
	assert len(df) == 10
	assert df["chr"].tolist() == ["chr1"] * 9 + ["chr2"]
	assert df["type"].tolist() == ["gene", "mRNA", "exon", "exon", "CDS", "transcript", "gene", "exon", "region", "gene"]
	assert df["end"].tolist() == [100, 100, 50, 50, 40, 100, 300, 250, 1000, 50]
	assert df["id"].tolist() == ["g1", "t1", "e1", "e1b", "c1", "t9", "g2", "e2", "r1", "g3"]


def test_iter_gff_reads_in_chunks_of_ten_thousand_lines(tmp_path):
	"""A 25,000-line GFF arrives as chunks of 10,000, 10,000 and 5,000 rows, in file order."""
	p = tmp_path / "big.gff"
	p.write_text("".join(f"chr1\ts\tgene\t{i}\t{i + 5}\t.\t+\t.\tID=g{i}\n" for i in range(1, 25001)))
	g = _gff_genome(tmp_path, p)
	sizes = [len(df) for df in g._iter_gff()]
	assert sizes == [10000, 10000, 5000]
	ids = np.concatenate([df["id"].to_numpy() for df in g._iter_gff()])
	assert ids[0] == "g1" and ids[-1] == "g25000" and len(set(ids)) == 25000


# ---------------------------------------------------------------------- init_gff


def test_init_gff_with_missing_file_returns_none(tmp_path):
	"""A genome whose gff is NA and that is given no filename has nothing to parse."""
	g = H.make_genome(tmp_path)
	assert g.init_gff() is None
	assert g.init_gff(None) is None


def test_init_gff_gene_table(tmp_path, gff_file, capsys):
	"""Genes (type gene) are indexed by (chr, start, end), sorted by chr then start, named by Name, else by ID, with a zero count column per occupancy 0..n."""
	g = _gff_genome(tmp_path, gff_file)
	genes = g.init_gff()
	capsys.readouterr()  # init_gff prints debug tables
	assert genes.index.names == ["chr", "start", "end"]
	assert genes["name"].tolist() == ["GeneOne", "g2", "Three"]
	assert genes["id"].tolist() == ["g1", "g2", "g3"]
	assert genes.index.get_level_values("chr").tolist() == ["chr1", "chr1", "chr2"]
	assert genes.index.get_level_values("end").tolist() == [100, 300, 50]
	for i in range(6):
		assert (genes[i] == 0).all()


def test_init_gff_annotation_table_names_children_after_their_gene(tmp_path, gff_file, capsys):
	"""Non-gene features are written to the anno tabix; children take the name of the gene they descend from.

	exon and CDS under mRNA t1 under gene g1 are named GeneOne; the exon under g2 is named g2 (ID fallback);
	the parentless region keeps its own Name; transcripts are left out; the two identical exons collapse to one row.
	"""
	g = _gff_genome(tmp_path, gff_file)
	g.init_gff()
	capsys.readouterr()
	rows = _anno_rows(g)
	assert sorted((r[0], r[2], r[3], r[4]) for r in rows) == sorted(
		[
			("chr1", "100", "mRNA", "GeneOne"),
			("chr1", "50", "exon", "GeneOne"),
			("chr1", "40", "CDS", "GeneOne"),
			("chr1", "250", "exon", "g2"),
			("chr1", "1000", "region", "Region"),
		]
	)
	assert [int(r[1]) for r in rows] == sorted(int(r[1]) for r in rows)


def test_init_gff_records_annotation_types(tmp_path, gff_file, capsys):
	"""The set of annotation types present is kept on the genome and written to anno_types.txt, one per line."""
	g = _gff_genome(tmp_path, gff_file)
	g.init_gff()
	capsys.readouterr()
	assert g.gff_anno_types == {"mRNA", "exon", "CDS", "region"}
	with open(g.anno_types_fname) as f:
		assert sorted(l.strip() for l in f) == ["CDS", "exon", "mRNA", "region"]


def test_init_gff_gene_types_option_moves_features_to_the_gene_table(tmp_path, gff_file, capsys):
	"""With gff_gene_types=[gene, mRNA] the mRNA is a gene row (named by its ID) and its children are named after it."""
	g = _gff_genome(tmp_path, gff_file, gff_gene_types=["gene", "mRNA"])
	genes = g.init_gff()
	capsys.readouterr()
	assert genes["name"].tolist() == ["GeneOne", "t1", "g2", "Three"]
	assert sorted((r[3], r[4]) for r in _anno_rows(g)) == sorted(
		[("exon", "t1"), ("CDS", "t1"), ("exon", "g2"), ("region", "Region")]
	)


def test_init_gff_name_option_picks_another_attribute(tmp_path, capsys):
	"""gff_name selects which attribute names genes; genes without it fall back to the ID."""
	p = tmp_path / "n.gff"
	p.write_text(
		"chr1\ts\tgene\t1\t10\t.\t+\t.\tID=a;Name=A;gene_id=GA\n"
		"chr1\ts\tgene\t20\t30\t.\t+\t.\tID=b;Name=B\n"
		"chr1\ts\texon\t1\t5\t.\t+\t.\tID=e;Parent=a\n"
	)
	g = _gff_genome(tmp_path, p, gff_name="gene_id")
	genes = g.init_gff()
	capsys.readouterr()
	assert genes["name"].tolist() == ["GA", "b"]
	assert [(r[3], r[4]) for r in _anno_rows(g)] == [("exon", "GA")]


def test_init_gff_resolves_multi_level_parents(tmp_path, capsys):
	"""exon -> mRNA -> gene chains of any depth resolve to the gene's name (here a three-level chain)."""
	p = tmp_path / "d.gff"
	p.write_text(
		"chr1\ts\tgene\t1\t100\t.\t+\t.\tID=g;Name=Top\n"
		"chr1\ts\tmRNA\t1\t100\t.\t+\t.\tID=m;Parent=g\n"
		"chr1\ts\tmiRNA_primary\t1\t90\t.\t+\t.\tID=p;Parent=m\n"
		"chr1\ts\texon\t1\t20\t.\t+\t.\tID=e;Parent=p\n"
	)
	g = _gff_genome(tmp_path, p)
	g.init_gff()
	capsys.readouterr()
	assert {r[3]: r[4] for r in _anno_rows(g)} == {"mRNA": "Top", "miRNA_primary": "Top", "exon": "Top"}


def test_init_gff_sorts_across_chunks(tmp_path, capsys):
	"""A GFF with 12,000 genes in descending order comes out sorted by chromosome and start across chunk boundaries."""
	p = tmp_path / "s.gff"
	lines = [f"chr{c}\ts\tgene\t{i}\t{i + 3}\t.\t+\t.\tID=c{c}g{i}\n" for c in (2, 1) for i in range(6000, 0, -1)]
	p.write_text("".join(lines))
	g = _gff_genome(tmp_path, p)
	genes = g.init_gff()
	capsys.readouterr()
	assert len(genes) == 12000
	assert genes.index.is_monotonic_increasing
	assert genes["id"].iloc[0] == "c1g1" and genes["id"].iloc[-1] == "c2g6000"


# ----------------------------------------------------------------- anno types file


def test_init_anno_types_in_write_mode_is_empty(tmp_path):
	"""While writing an index no annotation types are loaded."""
	g = H.make_genome(tmp_path, write=True)
	open(g.anno_types_fname, "w").write("exon\n")
	g._init_anno_types()
	assert g.gff_anno_types is None and g.anno_type_ids is None


def test_init_anno_types_reads_the_file_and_numbers_the_types(tmp_path):
	"""Read mode: exon, if present, becomes type 0 and the others follow in file order; without exon numbering starts at 1."""
	g = H.make_genome(tmp_path, write=True)
	g.write_mode = False
	with open(g.anno_types_fname, "w") as f:
		f.write("CDS\nmRNA\nexon\n")
	g._init_anno_types()
	assert g.gff_anno_types == {"CDS", "mRNA", "exon"}
	assert g.anno_type_ids.to_dict() == {"exon": 0, "CDS": 1, "mRNA": 2}
	with open(g.anno_types_fname, "w") as f:
		f.write("tRNA\nregion\n")
	g._init_anno_types()
	assert g.anno_type_ids.to_dict() == {"tRNA": 1, "region": 2}
	os.remove(g.anno_types_fname)
	g._init_anno_types()
	assert g.gff_anno_types is None and g.anno_type_ids is None


def test_write_anno_types_writes_one_type_per_line(tmp_path):
	"""_write_anno_types writes the set of types, one per line."""
	g = H.make_genome(tmp_path)
	g.gff_anno_types = {"exon", "CDS", "five_prime_UTR"}
	g._write_anno_types()
	with open(g.anno_types_fname) as f:
		assert sorted(f.read().split("\n")) == ["", "CDS", "exon", "five_prime_UTR"]


# ---------------------------------------------------------------- tabix tables


def test_write_tabix_makes_a_csi_indexed_table_that_load_tabix_reads(tmp_path):
	"""_write_tabix writes bed -> bgzip + csi; _load_tabix opens it, and returns None when the table is missing."""
	g = H.make_genome(tmp_path)
	assert g._load_tabix("gene") is None
	df = pd.DataFrame(
		{"chr": ["c1", "c1", "c2"], "start": [5, 50, 7], "end": [20, 80, 9], "name": ["a", "b", "c"], 1: [3, 4, 5], 5: [1, 0, 2]}
	)
	g._write_tabix(df, "gene")
	assert os.path.isfile(g.tabix_fname("gene")) and os.path.isfile(g.tabix_idx_fname("gene"))
	t = g._load_tabix("gene")
	assert [tuple(r) for r in t.fetch("c1", 0, 30)] == [("c1", "5", "20", "a", "3", "1")]
	assert [r[3] for r in t.fetch("c1", 60, 70)] == ["b"]
	assert [r[3] for r in t.fetch("c2")] == ["c"]
	t.close()


# ----------------------------------------------------------- queries on a real index


def test_query_genes_window_and_columns(anchored_reader, synthetic_genomes):
	"""query_genes returns chr, start, end, name and the unique/universal counts as integers; a window keeps overlapping genes."""
	g = synthetic_genomes
	alpha = anchored_reader["alpha"]
	allg = alpha.query_genes("chr1")
	assert list(allg.columns) == ["chr", "start", "end", "name", 1, 5]
	assert allg[["start", "end", 1, 5]].dtypes.map(lambda t: t.kind).eq("i").all()
	assert len(alpha.query_genes()) == 6  # no chromosome: every gene
	core = [x for x in g.data.genes["alpha"] if x[3] == "alpha_chr1_CORE_1"][0]
	win = alpha.query_genes("chr1", core[1] + 100, core[1] + 200)
	assert win["name"].tolist() == ["alpha_chr1_CORE_1"]
	assert len(alpha.query_genes("chr1", 1300, 1400)) == 1  # PAIR gene
	assert len(alpha.query_genes("chr1", 2400, 2500)) == 0  # past the last gene, inside the chromosome
	assert len(alpha.query_genes("chrZ")) == 0  # unknown chromosome gives an empty table


def test_gene_counts_are_occupancy_counts_over_the_stored_interval(anchored_reader, synthetic_genomes):
	"""Columns 1 and 5 of the gene table count positions in [start, end) shared by one and by all five genomes (brute force)."""
	g = synthetic_genomes
	for genome in ("alpha", "beta", "gamma", "delta"):
		reader = anchored_reader[genome]
		for chrom in g.data.seqs[genome]:
			occ = g.oracle_bits(genome, chrom).sum(axis=1)
			for _, r in reader.query_genes(chrom).iterrows():
				if r["end"] > len(occ):
					continue  # genes reaching the last k-1 bases are not counted (see report)
				window = occ[r["start"] : r["end"]]
				assert r[1] == (window == 1).sum(), (genome, r["name"])
				assert r[5] == (window == 5).sum(), (genome, r["name"])
	core = anchored_reader["alpha"].query_genes("chr1")
	row = core[core["name"] == "alpha_chr1_CORE_1"].iloc[0]
	assert row[5] >= 580 and row[1] <= 20  # a CORE gene is shared by all genomes
	priv = core[core["name"] == "alpha_chr1_private"].iloc[0]
	assert priv[1] >= 290 and priv[5] == 0  # a filler gene is private


def test_bitsum_genes_sum_the_per_gene_histograms(anchored_reader):
	"""bitsum_genes per chromosome is the column-wise sum of the occupancy histograms of its genes; bitfreq_genes are its row fractions."""
	alpha = anchored_reader["alpha"]
	assert list(alpha.bitsum_genes.index) == ["chr1", "chr2"]
	assert list(alpha.bitsum_genes.columns) == ["0", "1", "2", "3", "4", "5"]
	genes1 = alpha.query_genes("chr1")
	assert alpha.bitsum_genes.loc["chr1", "1"] == genes1[1].sum()
	assert alpha.bitsum_genes.loc["chr1", "5"] == genes1[5].sum()
	assert (alpha.bitsum_genes["0"] == 0).all()
	np.testing.assert_allclose(alpha.bitfreq_genes.sum(axis=1).to_numpy(), 1.0)


def test_genome_without_annotation_has_zero_gene_tables(anchored_reader):
	"""An unannotated anchor has all-zero gene tables with the gene-table columns, one row per chromosome."""
	e = anchored_reader[synthetic.NO_GFF]
	assert e.annotated is False
	assert e.bitsum_genes.shape == (1, 6)  # chr1 x gene-table columns [chr, start, end, name, 1, 5]
	assert (e.bitsum_genes == 0).all().all()
	assert e.bitfreq_genes is e.bitsum_genes


def test_query_anno_returns_typed_rows_with_type_ids(anchored_reader, synthetic_genomes):
	"""query_anno gives chr/start/end/type/name plus a type_id, restricted to the window, empty for unknown chromosomes."""
	alpha = anchored_reader["alpha"]
	df = alpha.query_anno("chr1", 0, 10000)
	assert list(df.columns) == ["chr", "start", "end", "type", "name", "type_id"]
	assert len(df) == 4 and set(df["type"]) == {"exon"}
	assert (df["type_id"] == alpha.anno_type_ids["exon"]).all()
	assert df[["start", "end"]].dtypes.map(lambda t: t.kind).eq("i").all()
	assert sorted(df["name"]) == sorted(n for c, _, _, n in synthetic_genomes.data.genes["alpha"] if c == "chr1")
	assert len(alpha.query_anno("chr1", 1300, 1400)) == 1
	assert len(alpha.query_anno("chrZ", 0, 100)) == 0 and list(alpha.query_anno("chrZ", 0, 100).columns)[:5] == ["chr", "start", "end", "type", "name"]


def test_query_anno_without_annotation_is_an_empty_frame(anchored_reader):
	"""A genome without GFF returns an empty table with the annotation columns."""
	df = anchored_reader[synthetic.NO_GFF].query_anno("chr1", 0, 1000)
	assert len(df) == 0 and list(df.columns) == pidx.TABIX_COLS


# ------------------------------------------------------------------ run_annotate


@pytest.fixture(scope="module")
def anchored_copy_factory(anchored_index, tmp_path_factory):
	"""Callable returning a fresh private copy of the anchored index (without KMC working files)."""
	n = [0]

	def make():
		n[0] += 1
		dest = tmp_path_factory.mktemp(f"annot{n[0]}") / "index"
		shutil.copytree(anchored_index, dest, ignore=shutil.ignore_patterns(".snakemake", "kmc", "tmp", "logs"))
		return dest

	return make


def _write_gff(path, genes):
	with open(path, "w") as f:
		for chrom, s, e, name in genes:
			f.write(f"{chrom}\tx\tgene\t{s}\t{e}\t.\t+\t.\tID={name};Name={name}\n")
			f.write(f"{chrom}\tx\texon\t{s}\t{e}\t.\t+\t.\tID={name}.e;Parent={name}\n")


def test_run_annotate_replaces_the_gene_tables(anchored_copy_factory, synthetic_genomes, tmp_path, capsys):
	"""Re-annotating alpha with a new GFF rewrites the gene and annotation tables and the per-chromosome gene histograms.

	Genes: one inside CORE, one inside alpha's private filler, one inside PAIR. Counts must be those of the brute-force
	occupancy over the stored intervals; the type of region decides which column is non-zero.
	"""
	copy = anchored_copy_factory()
	gff = tmp_path / "new.gff"
	_write_gff(gff, [("chr1", 401, 700, "inCore"), ("chr1", 21, 200, "inPrivate"), ("chr1", 1201, 1500, "inPair")])
	idx = open_index(copy)
	os.chdir(copy)
	idx["alpha"].gff = str(gff)  # run_annotate(gff_file) parses the genome's own gff (see report)
	idx["alpha"].run_annotate()
	idx.close()
	capsys.readouterr()

	idx = open_index(copy)
	os.chdir(copy)
	alpha = idx["alpha"]
	genes = alpha.query_genes("chr1")
	assert genes["name"].tolist() == ["inPrivate", "inCore", "inPair"]
	occ = synthetic_genomes.oracle_bits("alpha", "chr1").sum(axis=1)
	for _, r in genes.iterrows():
		window = occ[r["start"] : r["end"]]
		assert r[1] == (window == 1).sum() and r[5] == (window == 5).sum()
	byname = genes.set_index("name")
	assert byname.loc["inCore", 5] >= 290 and byname.loc["inCore", 1] == 0
	assert byname.loc["inPrivate", 1] >= 170 and byname.loc["inPrivate", 5] == 0
	assert byname.loc["inPair", 1] == 0 and byname.loc["inPair", 5] == 0
	hist = alpha.bitsum_genes.loc["chr1"]
	assert hist["1"] == byname[1].sum() and hist["5"] == byname[5].sum()
	assert hist["2"] >= 290  # the PAIR gene is shared by two genomes
	assert hist.sum() == sum(
		(occ[r["start"] : r["end"]]).size for _, r in genes.iterrows()
	)
	assert "chr2" not in alpha.bitsum_genes.index
	annos = alpha.query_anno("chr1", 0, 5000)
	assert sorted(annos["name"]) == ["inCore", "inPair", "inPrivate"]
	idx.close()


def test_run_annotate_nogene_only_writes_annotations(anchored_copy_factory, tmp_path):
	"""nogene=True writes the annotation table and types but leaves the gene table and histograms alone."""
	copy = anchored_copy_factory()
	gff = tmp_path / "new.gff"
	_write_gff(gff, [("chr1", 401, 700, "onlyAnno")])
	before = (copy / "anchor" / "beta" / "gene.bed.gz").read_bytes()
	hist_before = (copy / "anchor" / "beta" / "bitsum.genes.tsv").read_text()
	idx = open_index(copy)
	os.chdir(copy)
	idx["beta"].gff = str(gff)
	idx["beta"].run_annotate(nogene=True)
	idx.close()
	assert (copy / "anchor" / "beta" / "gene.bed.gz").read_bytes() == before
	assert (copy / "anchor" / "beta" / "bitsum.genes.tsv").read_text() == hist_before
	idx = open_index(copy)
	os.chdir(copy)
	assert idx["beta"].query_anno("chr1", 0, 10000)["name"].tolist() == ["onlyAnno"]
	idx.close()


def test_run_annotate_writes_to_a_log_file_when_asked(anchored_copy_factory, tmp_path):
	"""A gene whose end is not after its start is skipped and the warning goes to the log file when logfile is given."""
	import logging

	copy = anchored_copy_factory()
	gff = tmp_path / "bad.gff"
	_write_gff(gff, [("chr1", 401, 700, "ok"), ("chr1", 900, 900, "empty")])
	root = logging.getLogger()
	saved = root.handlers[:]
	root.handlers[:] = []
	idx = open_index(copy)
	os.chdir(copy)
	try:
		idx["alpha"].gff = str(gff)
		idx["alpha"].run_annotate(logfile=str(tmp_path / "ann.log"))
		for h in root.handlers:
			h.flush()
	finally:
		for h in root.handlers[:]:
			h.close()
		root.handlers[:] = saved
		idx.close()
	text = (tmp_path / "ann.log").read_text()
	assert "Skipping gene at chr1:900-900" in text

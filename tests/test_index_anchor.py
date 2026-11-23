"""Genome.run_anchor and the read-side tables it produces, with a fake KMC database (no KMC binaries needed).

The fake database returns counters computed from the position, so every expected table is derived independently
from the formula: bit 0 (the anchor itself) is always set, bit 1 follows (pos // 3) % 2 and bit 2 follows (pos // 5) % 2.
Needs the `bgzip` executable (the same one the real workflow uses); skipped without it.
"""

import os
import random
import shutil
import sys
import types

import numpy as np
import pandas as pd
import pytest

import helpers_index as H
import synthetic
from conftest import open_index
from panagram import index as pidx

K_ = 15
SIZES = {"c1": 400, "c2": 250}


def _bits(npos):
	p = np.arange(npos)
	return np.column_stack([np.ones(npos, int), (p // 3) % 2, (p // 5) % 2])


def _counter_values(npos):
	b = _bits(npos)
	return (b[:, 0] + 2 * b[:, 1] + 4 * b[:, 2]).astype("uint32")


class TableDb:
	"""Fake KMCFile: counters for a sequence are looked up by the sequence text."""

	def __init__(self, table):
		self.table = table

	def GetCountersForRead(self, seq, vec):
		vec.extend(self.table[seq].tolist())


GFF_MAIN = (
	"c1\tsrc\tgene\t11\t200\t.\t+\t.\tID=ga;Name=GA\n"
	"c1\tsrc\texon\t11\t100\t.\t+\t.\tID=ea;Parent=ga\n"
	"c2\tsrc\tgene\t21\t100\t.\t+\t.\tID=gb;Name=GB\n"
)
# c1 has 386 positions: Edge ends exactly at the last one, Over runs past it
GFF_EDGES = (
	"c1\tsrc\tgene\t301\t386\t.\t+\t.\tID=ge;Name=Edge\n"
	"c1\tsrc\tgene\t351\t400\t.\t+\t.\tID=go;Name=Over\n"
)


def _build_fake(root, gff_text):
	"""Write genome g0 (+ gff_text) under root, anchor it with the fake counters, return (idx, root, seqs). Caller closes idx and restores cwd."""
	rng = random.Random(11)
	seqs = {c: "".join(rng.choice("ACGT") for _ in range(n)) for c, n in SIZES.items()}
	synthetic.write_fasta(root / "g0.fa", seqs)
	(root / "g0.gff").write_text(gff_text)
	(root / "samples.tsv").write_text(
		f"name\tfasta\tgff\ng0\t{root / 'g0.fa'}\t{root / 'g0.gff'}\ng1\t{root / 'g1.fa'}\t\ng2\t{root / 'g2.fa'}\t\n"
	)
	idx = pidx.Index(str(root / "samples.tsv"), k=K_, anchor_genomes=["g0"])
	os.makedirs("anchor/g0")
	table = {s: _counter_values(len(s) - K_ + 1) for s in seqs.values()}
	mod = types.SimpleNamespace(KMCFile=lambda: None, CountVec=list)
	import panagram.extra

	sys.modules["panagram.extra.py_kmc_api"] = mod
	panagram.extra.py_kmc_api = mod
	try:
		idx["g0"].run_anchor([TableDb(table)])
	finally:
		del sys.modules["panagram.extra.py_kmc_api"]
		del panagram.extra.py_kmc_api
	return idx, root, seqs


@pytest.fixture(scope="module")
def anchored_fake(tmp_path_factory):
	"""A 3-genome write-mode index where genome g0 (two chromosomes, two genes) was anchored with fake counters."""
	if shutil.which("bgzip") is None:
		pytest.skip("bgzip is not on PATH")
	cwd = os.getcwd()
	try:
		idx, root, seqs = _build_fake(tmp_path_factory.mktemp("fakeanchor"), GFF_MAIN)
		yield idx, root, seqs
		idx.close()
	finally:
		os.chdir(cwd)


@pytest.fixture(scope="module")
def anchored_fake_edges(tmp_path_factory):
	"""Like anchored_fake, with one gene ending exactly at the last position of c1 and one running past it."""
	if shutil.which("bgzip") is None:
		pytest.skip("bgzip is not on PATH")
	cwd = os.getcwd()
	try:
		idx, root, seqs = _build_fake(tmp_path_factory.mktemp("fakeedges"), GFF_EDGES)
		yield idx, root, seqs
		idx.close()
	finally:
		os.chdir(cwd)


def test_run_anchor_counts_a_gene_ending_at_the_last_position_and_skips_one_past_it(anchored_fake_edges):
	"""A gene whose end equals the chromosome's position count is counted; one that ends beyond it is kept with zero counts."""
	idx, root, _ = anchored_fake_edges
	g = idx["g0"]
	occ = _bits(SIZES["c1"] - K_ + 1).sum(axis=1)
	genes = g.query_genes("c1").set_index("name")
	assert sorted(genes.index) == ["Edge", "Over"]
	edge = genes.loc["Edge"]
	window = occ[edge["start"] : edge["end"]]
	assert edge["end"] == 386 and len(window) == 386 - edge["start"] > 0
	assert edge[1] == (window == 1).sum() and edge[3] == (window == 3).sum()
	assert edge[1] + edge[3] > 0
	over = genes.loc["Over"]
	assert over[1] == 0 and over[3] == 0
	assert g.bitsum_genes.loc["c1", ["0", "1", "2", "3"]].sum() == len(window)


def test_run_anchor_writes_every_listed_file(anchored_fake):
	"""All files named by anchor_filenames exist after run_anchor, plus the umap and paircount tables."""
	idx, root, _ = anchored_fake
	g = idx["g0"]
	for f in g.anchor_filenames + [g.chrom_umaps_filename, g.genome_umap_filename, "anchor/g0/total_paircounts.csv"]:
		assert (root / f).is_file(), f
	assert not (root / "anchor" / "g1").exists()  # non-anchor genomes get nothing


def test_run_anchor_bitmap_matches_the_counters(anchored_fake):
	"""Reading the written bitmap back gives the bits the counters encode, at both resolutions."""
	idx, root, _ = anchored_fake
	g = idx["g0"]
	assert list(g.genome_names) == ["g0", "g1", "g2"]
	for chrom, n in SIZES.items():
		npos = n - K_ + 1
		got = g.query(chrom)
		assert list(got.columns) == ["g0", "g1", "g2"]
		np.testing.assert_array_equal(got.to_numpy(), _bits(npos))
		low = g.query(chrom, step=100)
		np.testing.assert_array_equal(low.to_numpy(), _bits(npos)[::100])
		assert list(low.index) == list(range(0, npos, 100))
		np.testing.assert_array_equal(g.query(chrom, 17, 90).to_numpy(), _bits(npos)[17:90])


def test_run_anchor_chromosome_table_counts_genes(anchored_fake):
	"""chrs.tsv records k-mer counts and the number of genes per chromosome."""
	idx, root, _ = anchored_fake
	chrs = pd.read_table(root / "anchor" / "g0" / "chrs.tsv", index_col="name")
	assert chrs["size"].to_dict() == {"c1": 386, "c2": 236}
	assert chrs["gene_count"].to_dict() == {"c1": 1, "c2": 1}
	assert chrs["id"].to_dict() == {"c1": 0, "c2": 1}


def test_run_anchor_bitsum_bins_and_totals(anchored_fake):
	"""bitsum.bins.tsv holds occupancy histograms per bin (386 // 100 = 3 wide), and bitsum_chrs/total sum them."""
	idx, root, _ = anchored_fake
	g = idx["g0"]
	for chrom, n in SIZES.items():
		npos = n - K_ + 1
		occ = _bits(npos).sum(axis=1)
		binlen = npos // 100
		rows = g.bitsum_bins.loc[chrom]
		assert list(rows.index) == list(range(0, npos, binlen))
		for start, row in rows.iterrows():
			np.testing.assert_array_equal(row.to_numpy(), np.bincount(occ[start : start + binlen], minlength=4))
		np.testing.assert_array_equal(g.bitsum_chrs.loc[chrom].to_numpy(), np.bincount(occ, minlength=4))
		np.testing.assert_allclose(g.bitfreq_chrs.loc[chrom].to_numpy(), np.bincount(occ, minlength=4) / npos)
		np.testing.assert_allclose(g.bitfreq_bins.loc[chrom].sum(axis=1).to_numpy(), 1.0)
	all_occ = np.concatenate([_bits(n - K_ + 1).sum(axis=1) for n in SIZES.values()])
	np.testing.assert_array_equal(g.bitsum_total.to_numpy(), np.bincount(all_occ, minlength=4))


def test_read_bitsum_bins_decodes_chromosome_ids(anchored_fake):
	"""_read_bitsum_bins maps the numeric chr column back to chromosome names and gives integer column labels."""
	idx, root, _ = anchored_fake
	g = idx["g0"]
	df = g._read_bitsum_bins()
	assert df.index.names == ["chr", "start"]
	assert sorted(set(df.index.get_level_values("chr"))) == ["c1", "c2"]
	assert list(df.columns) == [0, 1, 2, 3]
	raw = pd.read_table(root / "anchor" / "g0" / "bitsum.bins.tsv")
	assert len(df) == len(raw) and df.to_numpy().sum() == 386 + 236


def test_run_anchor_total_paircounts(anchored_fake):
	"""total_paircounts counts, per genome, the anchor positions whose k-mer that genome carries; frac is relative to the anchor."""
	idx, root, _ = anchored_fake
	g = idx["g0"]
	want = sum(_bits(n - K_ + 1).sum(axis=0) for n in SIZES.values())
	assert list(g.total_paircounts.index) == ["g0", "g1", "g2"]
	np.testing.assert_array_equal(g.total_paircounts["count"].to_numpy(), want)
	np.testing.assert_allclose(g.total_paircounts["frac"].to_numpy(), want / want[0])
	assert g.total_paircounts.loc["g0", "frac"] == 1.0


def test_run_anchor_gene_tables(anchored_fake):
	"""Gene rows carry unique (occupancy 1) and universal (occupancy 3) counts over their stored interval; bitsum_genes sums them."""
	idx, root, _ = anchored_fake
	g = idx["g0"]
	occ = {c: _bits(n - K_ + 1).sum(axis=1) for c, n in SIZES.items()}
	assert g.annotated
	total1 = total3 = 0
	for chrom in SIZES:
		genes = g.query_genes(chrom)
		assert genes["name"].tolist() == [{"c1": "GA", "c2": "GB"}[chrom]]
		r = genes.iloc[0]
		window = occ[chrom][r["start"] : r["end"]]
		assert r[1] == (window == 1).sum() and r[3] == (window == 3).sum()
		assert g.bitsum_genes.loc[chrom, "1"] == r[1] and g.bitsum_genes.loc[chrom, "3"] == r[3]
		assert g.bitsum_genes.loc[chrom].sum() == len(window)
	assert g.gff_anno_types == {"exon"}
	cwd = os.getcwd()
	try:
		annos = open_index(root)["g0"].query_anno("c1", 0, 1000)  # a reader loads the type ids
	finally:
		os.chdir(cwd)
	assert annos["name"].tolist() == ["GA"] and annos["type"].tolist() == ["exon"]
	assert annos["type_id"].tolist() == [0]


def test_run_anchor_loaded_state_matches_files(anchored_fake):
	"""After run_anchor the genome is in read state: the written row counts are recorded and one BGZF block holds each bitmap."""
	idx, root, _ = anchored_fake
	g = idx["g0"]
	assert g.bitmap_lens[1] == 386 + 236
	assert g.bitmap_lens[100] == 4 + 3
	assert g.blocks[1].tolist() == [(0, 0)]


def test_bin_bitsum_with_exactly_min_bin_count_positions():
	"""A chromosome with exactly min_bin_count (100) positions gives 100 one-position bins, each counting one position."""
	g = pidx.Genome.__new__(pidx.Genome)
	g.params = {"max_bin_kbp": 200, "min_bin_count": 100}
	g.ngenomes = 3
	out = g.bin_bitsum(np.ones(100, dtype=int))
	assert out.shape == (100, 4) and (out[1] == 1).all()

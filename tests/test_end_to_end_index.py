"""End to end: synthetic genomes -> panagram index -> bitmap, bitdump, annotation.

Expected values come from how the genomes are built (tests/synthetic.py) and from a
brute-force canonical k-mer oracle, never from a recorded panagram output.
"""

import collections

import numpy as np
import pandas as pd
import pytest

import synthetic
from conftest import K


def _piece(genomes, genome, chrom, kind, nth=0):
	return [p for p in genomes.data.pieces[genome][chrom] if p.kind == kind][nth]


# ---------------------------------------------------------------- the synthetic input


def test_synthetic_layout_matches_its_description(synthetic_genomes):
	"""Segments appear where LAYOUT says, with the right lengths and strand."""
	g = synthetic_genomes
	for genome, chroms in synthetic.LAYOUT.items():
		for chrom, layout in chroms.items():
			assert len(g.data.seqs[genome][chrom]) == sum(
				arg if kind == "rand" else synthetic.SEGMENT_LEN[kind] for kind, arg in layout
			)
	for kind, members in synthetic.MEMBERS.items():
		carriers = {
			genome
			for genome, chroms in synthetic.LAYOUT.items()
			for layout in chroms.values()
			if any(k == kind for k, _ in layout)
		}
		private = {"DIV": {"epsilon"}}.get(kind, set())
		assert carriers - private == members, kind
	core = g.data.seqs["alpha"]["chr1"][300:900]
	assert core in g.data.seqs["gamma"]["chr1"]
	# delta carries CORE reverse-complemented
	assert synthetic.revcomp(core) in g.data.seqs["delta"]["chr1"]
	assert core not in g.data.seqs["delta"]["chr1"]


def test_epsilon_div_differs_from_the_others_in_every_window(synthetic_genomes):
	"""Every k-mer window of epsilon's DIV copy holds a substitution relative to alpha's copy."""
	g = synthetic_genomes
	a = _piece(g, "alpha", "chr1", "DIV")
	e = _piece(g, "epsilon", "chr1", "DIV")
	sa = g.data.seqs["alpha"]["chr1"][a.start : a.end]
	se = g.data.seqs["epsilon"]["chr1"][e.start : e.end]
	diff = np.array([x != y for x, y in zip(sa, se)])
	assert len(sa) == len(se)
	for i in range(len(sa) - K + 1):
		assert diff[i : i + K].any()


# ------------------------------------------------------------------- prepare step


def test_prepare_writes_config_samples_and_snakefile(prepared_index, synthetic_genomes):
	"""`panagram index samples.tsv -k 15 --prepare` writes config, samples (ids, anchors) and a Snakefile."""
	import yaml

	assert (prepared_index / "Snakefile").is_file()
	cfg = yaml.safe_load((prepared_index / "config.yaml").read_text())
	assert cfg["k"] == K
	assert cfg["anchor_genomes"] == synthetic.GENOMES
	samples = pd.read_table(prepared_index / "samples.tsv")
	assert list(samples.columns) == ["name", "fasta", "gff", "id", "anchor"]
	assert list(samples["name"]) == synthetic.GENOMES
	assert list(samples["id"]) == list(range(len(synthetic.GENOMES)))
	assert samples["anchor"].all()
	assert samples["gff"].isna().tolist() == [g == synthetic.NO_GFF for g in synthetic.GENOMES]
	# prepare must not have anchored anything
	assert not (prepared_index / "anchor").exists()


# ---------------------------------------------------------------- anchored outputs


def test_chromosome_sizes_are_kmer_counts(anchored_reader, synthetic_genomes):
	"""Each chromosome's bitmap length is its sequence length minus k plus 1."""
	g = synthetic_genomes
	for genome in g.names:
		chrs = anchored_reader[genome].chrs
		assert list(chrs.index) == list(g.data.seqs[genome])
		for chrom, seq in g.data.seqs[genome].items():
			assert chrs.loc[chrom, "size"] == len(seq) - K + 1


def test_genome_order_gives_bit_order(anchored_reader):
	"""Bit j of a bitmap row is the genome with id j, in samples.tsv order."""
	assert list(anchored_reader.genome_names) == synthetic.GENOMES
	df = anchored_reader.query_bitmap("beta", "chr1", 0, 5)
	assert list(df.columns) == synthetic.GENOMES


def test_bitmap_equals_construction(anchored_reader, synthetic_genomes):
	"""Where the layout determines a position, the bitmap row is exactly what the layout implies.

	Boundary-straddling k-mers are undetermined by the layout; they are covered by the
	brute-force test below. The genome's own bit must be set at every position.
	"""
	g = synthetic_genomes
	for genome in g.names:
		for chrom in g.data.seqs[genome]:
			got = anchored_reader.query_bitmap(genome, chrom).to_numpy()
			expected = g.expected_bits(genome, chrom)
			assert (expected >= 0).mean() > 0.9  # the check covers most of the genome
			g.assert_matches_expected(got, expected)
			assert got[:, g.names.index(genome)].all()


def test_bitmap_equals_bruteforce_canonical_kmer_membership(anchored_reader, synthetic_genomes):
	"""Every position of every anchor chromosome: bit j is set iff the canonical k-mer occurs in genome j."""
	g = synthetic_genomes
	for genome in g.names:
		for chrom in g.data.seqs[genome]:
			got = anchored_reader.query_bitmap(genome, chrom).to_numpy()
			np.testing.assert_array_equal(got, g.oracle_bits(genome, chrom), err_msg=f"{genome} {chrom}")


@pytest.mark.parametrize("genome", synthetic.GENOMES)
def test_core_segment_is_set_for_all_genomes(anchored_reader, synthetic_genomes, genome):
	"""Inside CORE (also when reverse-complemented) every genome's bit is set."""
	g = synthetic_genomes
	for chrom, plist in g.data.pieces[genome].items():
		for p in plist:
			if p.kind != "CORE":
				continue
			bits = anchored_reader.query_bitmap(genome, chrom, p.start, p.end - K + 1).to_numpy()
			assert bits.shape == (synthetic.SEGMENT_LEN["CORE"] - K + 1, 5)
			assert bits.all()


def test_pair_trio_and_private_membership(anchored_reader, synthetic_genomes):
	"""PAIR is set only for alpha and beta, TRIO only for gamma, delta, epsilon, filler only for its own genome."""
	g = synthetic_genomes
	expected_cols = {
		"PAIR": [1, 1, 0, 0, 0],
		"TRIO": [0, 0, 1, 1, 1],
	}
	for genome in g.names:
		for chrom, plist in g.data.pieces[genome].items():
			for p in plist:
				if p.kind in expected_cols:
					bits = anchored_reader.query_bitmap(genome, chrom, p.start, p.end - K + 1)
					assert (bits.to_numpy() == expected_cols[p.kind]).all()
				elif p.kind == "rand":
					# k-mers wholly inside filler belong to this genome alone
					if p.end - p.start < K:
						continue
					bits = anchored_reader.query_bitmap(genome, chrom, p.start, p.end - K + 1)
					want = [int(n == genome) for n in g.names]
					assert (bits.to_numpy() == want).all()


def test_diverged_segment_is_private_in_epsilon_only(anchored_reader, synthetic_genomes):
	"""DIV: alpha..delta share it, epsilon's mutated copy matches nobody else."""
	g = synthetic_genomes
	for genome in g.names:
		p = _piece(g, genome, "chr1", "DIV")
		bits = anchored_reader.query_bitmap(genome, "chr1", p.start, p.end - K + 1).to_numpy()
		want = [1, 1, 1, 1, 0] if genome != "epsilon" else [0, 0, 0, 0, 1]
		assert (bits == want).all(), genome


def test_segment_boundaries(anchored_reader, synthetic_genomes):
	"""The last k-mer wholly inside PAIR is shared; k-mers with 8 or more flank bases are private.

	A k-mer with only 1-2 flank bases may match another genome by chance, so the
	boundary is checked at 8 bases of overlap (a fixed outcome for the fixed seed).
	"""
	g = synthetic_genomes
	p = _piece(g, "alpha", "chr1", "PAIR")
	row = lambda pos: anchored_reader.query_bitmap("alpha", "chr1", pos, pos + 1).to_numpy()[0].tolist()
	assert row(p.end - K) == [1, 1, 0, 0, 0]
	assert row(p.end - K + 8) == [1, 0, 0, 0, 0]  # 8 bases beyond the segment end
	assert row(p.start - 8) == [1, 0, 0, 0, 0]  # window covers 8 bases before the segment... and 7 inside


def test_lowres_bitmap_samples_every_hundredth_position(anchored_reader, synthetic_genomes):
	"""step=100 returns the full-resolution rows at positions 0, 100, 200, ..."""
	g = synthetic_genomes
	for genome in ("alpha", "gamma"):
		for chrom in g.data.seqs[genome]:
			full = g.oracle_bits(genome, chrom)
			got = anchored_reader.query_bitmap(genome, chrom, step=100)
			assert list(got.index) == list(range(0, len(full), 100))
			np.testing.assert_array_equal(got.to_numpy(), full[::100])


def test_bitmap_window_matches_slice(anchored_reader, synthetic_genomes):
	"""A window query returns the same rows as slicing the whole chromosome, with a coordinate index."""
	g = synthetic_genomes
	full = g.oracle_bits("beta", "chr1")
	got = anchored_reader.query_bitmap("beta", "chr1", 777, 1234)
	assert list(got.index) == list(range(777, 1234))
	np.testing.assert_array_equal(got.to_numpy(), full[777:1234])


def test_occupancy_histogram_matches_construction(anchored_reader, synthetic_genomes):
	"""bitsum_chrs[g][chrom][n] counts positions shared by n genomes (counted by brute force)."""
	g = synthetic_genomes
	for genome in g.names:
		for chrom in g.data.seqs[genome]:
			counts = collections.Counter(g.oracle_bits(genome, chrom).sum(axis=1))
			got = anchored_reader[genome].bitsum_chrs.loc[chrom]
			for n in range(1, 6):
				assert got[n] == counts.get(n, 0), (genome, chrom, n)


def test_total_paircounts_count_shared_positions(anchored_reader, synthetic_genomes):
	"""total_paircounts[j] counts the anchor's positions whose k-mer is present in genome j (counted by brute force)."""
	g = synthetic_genomes
	for genome in g.names:
		want = sum(
			(g.oracle_bits(genome, c).sum(axis=0) for c in g.data.seqs[genome]), np.zeros(5, dtype=int)
		)
		tp = anchored_reader[genome].total_paircounts
		assert list(tp.index) == g.names
		np.testing.assert_array_equal(tp["count"].to_numpy(), want)
		np.testing.assert_allclose(tp["frac"].to_numpy(), want / want[g.names.index(genome)])


def test_genome_dist_lists_every_pair(anchored_index):
	"""mash triangle output has one row per unordered genome pair, distances in [0, 1]."""
	rows = pd.read_table(anchored_index / "genome_dist.tsv", header=None)
	pairs = {frozenset(r[:2]) for r in rows.itertuples(index=False)}
	assert len(rows) == 10
	assert pairs == {
		frozenset((a, b)) for i, a in enumerate(synthetic.GENOMES) for b in synthetic.GENOMES[:i]
	}
	assert rows[2].between(0, 1).all()


# ------------------------------------------------------------------ annotation


def test_gene_table_has_the_gff_genes(anchored_reader, synthetic_genomes):
	"""query_genes returns exactly the gene names written to the GFF3, on the right chromosomes."""
	g = synthetic_genomes
	for genome, genes in g.data.genes.items():
		for chrom in g.data.seqs[genome]:
			df = anchored_reader.query_genes(genome, chrom)
			want = sorted(n for c, _, _, n in genes if c == chrom)
			assert sorted(df["name"]) == want
			assert (df["chr"] == chrom).all()
		assert anchored_reader[genome].chrs["gene_count"].sum() == len(genes)


def test_genome_without_gff_has_no_gene_tables(anchored_reader):
	"""A genome with a blank gff column is anchored but not annotated."""
	e = anchored_reader[synthetic.NO_GFF]
	assert e.gene_tabix is None and e.anno_tabix is None
	assert len(e.query_genes("chr1")) == 0
	assert e.chrs["gene_count"].sum() == 0


# ------------------------------------------------------------------------ bitdump


def _bitdump_rows(stdout):
	"""Parse the DataFrame that `panagram bitdump` prints (header of genome names, then index + bits)."""
	lines = stdout.strip().splitlines()
	rows = [l.split() for l in lines[1:]]
	header = [h for h in lines[0].split() if h != "name"]  # pandas prints the columns' name too
	return header, [int(r[0]) for r in rows], np.array([[int(x) for x in r[1:]] for r in rows])


def test_bitdump_window_matches_construction(anchored_index, synthetic_genomes, run_panagram):
	"""`panagram bitdump` over a window spanning the end of CORE then filler prints one 0/1 row per position.

	Windows are kept under pandas' print truncation limit (60 rows).
	"""
	g = synthetic_genomes
	core = _piece(g, "alpha", "chr1", "CORE")
	start, end = core.end - K - 5, core.end - K + 5  # last 5 shared k-mers, then 5 spilling ones
	r = run_panagram("bitdump", anchored_index, "alpha", "chr1", start, end)
	names, index, rows = _bitdump_rows(r.stdout)
	assert names == synthetic.GENOMES
	assert index == list(range(start, end))
	assert (rows[:5] == [1, 1, 1, 1, 1]).all()
	np.testing.assert_array_equal(rows, g.oracle_bits("alpha", "chr1")[start:end])
	assert rows[:, 0].all()  # the anchor's own bit


def test_bitdump_pair_region(anchored_index, synthetic_genomes, run_panagram):
	"""Inside the second copy of PAIR in beta the dump shows the alpha and beta columns set and nothing else."""
	g = synthetic_genomes
	pair = _piece(g, "beta", "chr1", "PAIR", nth=1)
	r = run_panagram("bitdump", anchored_index, "beta", "chr1", pair.start + 10, pair.start + 30)
	names, index, rows = _bitdump_rows(r.stdout)
	assert len(rows) == 20
	assert (rows == [1, 1, 0, 0, 0]).all()


def test_bitdump_step_samples_rows(anchored_index, synthetic_genomes, run_panagram):
	"""bitdump with step 100 over a 5-row-wide stretch prints the rows at multiples of 100 only."""
	g = synthetic_genomes
	full = g.oracle_bits("gamma", "chr1")
	r = run_panagram("bitdump", anchored_index, "gamma", "chr1", 0, 1000, 100)
	_, index, rows = _bitdump_rows(r.stdout)
	assert index == list(range(0, 1000, 100))
	np.testing.assert_array_equal(rows, full[0:1000:100])

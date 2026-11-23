"""panagram.index.Genome without KMC: file names, chromosome tables, bitmap reads, binning, KMC byte packing.

Expected values are worked by hand or come from random truth matrices written by the test itself.
"""

import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from Bio import bgzf

import helpers_index as H
import synthetic
from panagram import index as pidx


# ------------------------------------------------------------------ construction


def test_unanchored_genome_returns_early_with_no_files(tmp_path):
	"""A genome with no fasta is not anchored, has no chrs and no anchor files, and its gff status follows pd.isna."""
	idx = H.stub_index(tmp_path, ["a", "b"])
	g = pidx.Genome(idx, 1, "b")
	assert g.anchored is False
	assert g.chrs is None and g.bitmaps is None
	assert g.annotated is False
	assert g.anchor_filenames == []
	assert g.nbytes == 1 and g.ngenomes == 2
	g.close()  # no bitmaps: must not raise

	annotated = pidx.Genome(idx, 0, "a", gff="x.gff")
	assert annotated.annotated is True and annotated.anchored is False


def test_anchored_flag_defaults_to_having_a_fasta(tmp_path):
	"""anchor=None falls back to 'has a fasta'; an explicit anchor value wins."""
	(tmp_path / "anchor" / "g").mkdir(parents=True)
	fa = tmp_path / "g.fa"
	synthetic.write_fasta(fa, {"c": "ACGTACGTACGTACGTACGT"})
	idx = H.stub_index(tmp_path, ["g"], fastas={"g": str(fa)}, k=5)
	assert pidx.Genome(idx, 0, "g", fasta=str(fa), write=True).anchored is True
	assert pidx.Genome(idx, 0, "g", fasta=str(fa), anchor=False, write=True).anchored is False
	assert pidx.Genome(idx, 0, "g", fasta=None, write=True, anchor=None).anchored is False


def test_init_chrs_sizes_are_kmer_counts_in_fasta_order(tmp_path):
	"""init_chrs gives size = length - k + 1 per chromosome, ids in FASTA order, indexed by name."""
	fa = tmp_path / "g.fa"
	synthetic.write_fasta(fa, {"z": "A" * 30, "a": "C" * 12})
	g = H.make_genome(tmp_path, fasta=str(fa), k=7)
	assert list(g.chrs.index) == ["z", "a"]
	assert list(g.chrs["id"]) == [0, 1]
	assert list(g.chrs["size"]) == [24, 6]
	assert list(g.chrs["gene_count"]) == [0, 0]
	assert g.chr_count == 2
	assert g.seq_len("a") == 6 and g.seq_len("z") == 24


def test_write_and_load_chrs_roundtrip(tmp_path):
	"""write_chrs writes a TSV with a name column that load_chrs reads back unchanged."""
	fa = tmp_path / "g.fa"
	synthetic.write_fasta(fa, {"c1": "ACGT" * 10, "c2": "ACGT" * 5})
	g = H.make_genome(tmp_path, fasta=str(fa), k=5)
	g.chrs["gene_count"] = [3, 1]
	g.write_chrs()
	lines = Path(g.chrs_fname).read_text().splitlines()
	assert lines[0].split("\t") == ["name", "id", "size", "gene_count"]
	g2 = H.make_genome(tmp_path, fasta=None, write=True)  # finds chrs.tsv, loads it
	pd.testing.assert_frame_equal(g2.chrs, g.chrs)
	assert list(g2.sizes) == [36, 16]


def test_set_chrs_offsets_count_rows_per_step(tmp_path):
	"""offsets[step] is the number of rows of earlier chromosomes in that step's bitmap: cumulative ceil(size/step)."""
	g = H.make_genome(tmp_path)
	g.set_chrs(pd.DataFrame({"id": [0, 1, 2], "size": [10, 250, 100]}, index=pd.Index(["a", "b", "c"], name="name")))
	assert list(g.offsets[1]) == [0, 10, 260]
	assert list(g.offsets[100]) == [0, 1, 4]
	assert "gene_count" in g.chrs.columns


def test_steps_are_one_and_lowres_step(tmp_path):
	"""Genome.steps is [1, lowres_step] from the parameters."""
	g = H.make_genome(tmp_path, lowres_step=25)
	g.steps = None
	g._init_steps()
	assert g.steps == [1, 25]


def test_bitsum_index_and_gene_table_columns(tmp_path):
	"""bitsum_index runs 0..ngenomes; the gene table keeps the unique (1) and universal (ngenomes) columns."""
	g = H.make_genome(tmp_path, ngenomes=7)
	assert list(g.bitsum_index) == list(range(8))
	assert g.gene_tabix_cols == ["chr", "start", "end", "name", 1, 7]
	assert g.gene_tabix_types == {"start": int, "end": int, 1: int, 7: int}
	assert g.nbytes == 1
	assert H.make_genome(tmp_path, ngenomes=9, name="x").nbytes == 2


# ---------------------------------------------------------------------- file names


def test_filenames_follow_the_anchor_layout(tmp_path):
	"""Per-genome files live in <prefix>/anchor/<name>/ with the documented names."""
	g = H.make_genome(tmp_path, name="g0")
	base = tmp_path / "anchor" / "g0"
	assert Path(g.prefix) == base
	assert Path(g.chrs_fname) == base / "chrs.tsv"
	assert Path(g.bins_fname) == base / "bitsum.bins.tsv"
	assert Path(g.chr_genes_fname) == base / "bitsum.genes.tsv"
	assert Path(g.anno_types_fname) == base / "anno_types.txt"
	assert Path(g.bitmap_gz_fname(100)) == base / "bitmap.100.gz"
	assert Path(g.bitmap_gzi_fname(1)) == base / "bitmap.1.gzi"
	assert Path(g.bed_tmp_fname("gene")) == base / "gene.bed"
	assert Path(g.tabix_fname("anno")) == base / "anno.bed.gz"
	assert Path(g.tabix_idx_fname("anno")) == base / "anno.bed.gz.csi"
	assert Path(g.chrom_umaps_filename) == base / "chrom_umaps.csv"
	assert Path(g.genome_umap_filename) == base / "genome_umap.csv"


def test_anchor_filenames_list_exactly_the_files_a_complete_anchor_has(anchored_reader, anchored_index):
	"""Annotated anchors list chrs, bins, both bitmaps with indexes, gene and anno tables; all exist on disk.

	A genome without a GFF omits the gene/anno files. Index.anchor_filenames concatenates the genomes' lists.
	"""
	alpha = anchored_reader["alpha"].anchor_filenames
	base = "anchor/alpha/"
	assert alpha == [
		base + "chrs.tsv",
		base + "bitsum.bins.tsv",
		base + "bitmap.1.gz",
		base + "bitmap.1.gzi",
		base + "bitmap.100.gz",
		base + "bitmap.100.gzi",
		base + "bitsum.genes.tsv",
		base + "gene.bed.gz",
		base + "gene.bed.gz.csi",
		base + "anno.bed.gz",
		base + "anno.bed.gz.csi",
	]
	eps = anchored_reader[synthetic.NO_GFF].anchor_filenames
	assert eps == [f.replace("alpha", "epsilon") for f in alpha[:6]]
	everything = anchored_reader.anchor_filenames
	assert everything == sum((anchored_reader[n].anchor_filenames for n in synthetic.GENOMES), [])
	for f in everything:
		assert (anchored_index / f).is_file(), f
	assert anchored_reader.get_anchor_filenames("beta") == anchored_reader["beta"].anchor_filenames


# ---------------------------------------------------------------------- bgzf index


def test_load_bgz_blocks_reads_samtools_gzi(tmp_path):
	"""A .gzi holds n, then n (compressed, uncompressed) offset pairs; the loader prepends the (0, 0) first block."""
	g = H.make_genome(tmp_path)
	f = tmp_path / "x.gzi"
	with open(f, "wb") as out:
		np.array([2], dtype="uint64").tofile(out)
		np.array([[100, 65280], [250, 130560]], dtype="uint64").tofile(out)
	blocks = g.load_bgz_blocks(str(f))
	assert blocks.dtype.names == ("rstart", "dstart")
	assert blocks["rstart"].tolist() == [0, 100, 250]
	assert blocks["dstart"].tolist() == [0, 65280, 130560]
	with open(f, "wb") as out:
		np.array([0], dtype="uint64").tofile(out)
	assert g.load_bgz_blocks(str(f)).tolist() == [(0, 0)]


def test_load_bgz_blocks_matches_the_blocks_of_a_real_bgzf_file(tmp_path):
	"""The block table of a BGZF file written block by block matches the block offsets Biopython reports."""
	g = H.make_genome(tmp_path)
	data = bytes(range(256)) * 4
	blocks = H.write_bgzf_with_gzi(tmp_path / "d.gz", data, 100)
	got = g.load_bgz_blocks(str(tmp_path / "d.gzi"))
	assert [tuple(b) for b in got.tolist()] == blocks
	assert got["dstart"].tolist() == list(range(0, 1024, 100))


# ------------------------------------------------------------------ bitmap reading


def test_bytes_to_bits_is_little_endian_and_trimmed(tmp_path):
	"""Genome j is bit j of the row (byte 0 holds genomes 0-7, lowest bit first); padding bits are dropped."""
	g = H.make_genome(tmp_path, ngenomes=3)
	pac = np.array([[0b101], [0b010], [0b111]], dtype="uint8")
	assert g._bytes_to_bits(pac).tolist() == [[1, 0, 1], [0, 1, 0], [1, 1, 1]]
	g = H.make_genome(tmp_path, ngenomes=11, name="x")
	pac = np.array([[0x01, 0x00], [0x80, 0x02], [0xFF, 0x07]], dtype="uint8")
	got = g._bytes_to_bits(pac)
	assert got.shape == (3, 11)
	assert got[0].tolist() == [1] + [0] * 10
	assert got[1].tolist() == [0] * 7 + [1, 0, 1, 0]
	assert got[2].tolist() == [1] * 11


@pytest.mark.parametrize("ngenomes", [5, 11])
@pytest.mark.parametrize("step", [1, 3, 7, 100, 150, 200])
def test_query_matches_truth_across_bgzf_blocks(tmp_path, ngenomes, step):
	"""Windows that start mid-block and span several 50-byte BGZF blocks return the truth rows, with a coordinate index.

	step 100 and 200 read the low-resolution file; other steps subsample the full bitmap. Starts are multiples of 100
	for the low-resolution cases (see the report on unaligned starts).
	"""
	g, bits = H.build_bitmap_genome(tmp_path, ngenomes, {"c1": 300, "c2": 450})
	assert len(g.blocks[1]) >= 15
	names = list(g.genome_names)
	for chrom, a, b in [("c1", 0, 300), ("c2", 0, 450), ("c2", 100, 400), ("c1", 200, 300), ("c2", 300, 451 - 1)]:
		if step in (100, 200) and a % 100:
			continue
		got = g.query(chrom, a, b, step)
		assert list(got.index) == list(range(a, b, step))
		assert list(got.columns) == names
		np.testing.assert_array_equal(got.to_numpy(), bits[chrom][a:b:step], err_msg=f"{chrom} {a} {b}")


def test_query_unaligned_starts_and_defaults(tmp_path):
	"""At full resolution any start works; with no start and end the whole chromosome is returned."""
	g, bits = H.build_bitmap_genome(tmp_path, 5, {"c1": 300, "c2": 450})
	got = g.query("c2", 123, 321)
	np.testing.assert_array_equal(got.to_numpy(), bits["c2"][123:321])
	np.testing.assert_array_equal(g.query("c1").to_numpy(), bits["c1"])
	np.testing.assert_array_equal(g.query("c2", 17).to_numpy(), bits["c2"][17:])
	np.testing.assert_array_equal(g.query("c2", None, 50).to_numpy(), bits["c2"][:50])
	assert list(g.query("c2", 17).index) == list(range(17, 450))


def test_query_unknown_chromosome_raises_keyerror(tmp_path):
	"""A chromosome that is not in the index raises KeyError."""
	g, _ = H.build_bitmap_genome(tmp_path, 5, {"c1": 300})
	with pytest.raises(KeyError):
		g.query("nope", 0, 10)


def test_query_bytes_returns_packed_rows_and_substeps(tmp_path):
	"""_query_bytes returns raw packed rows: step 1 gives consecutive rows, step 2*bstep keeps every second row."""
	g, bits = H.build_bitmap_genome(tmp_path, 11, {"c1": 300, "c2": 450})
	packed = np.packbits(bits["c2"], axis=1, bitorder="little")
	got = g._query_bytes("c2", 40, 59, 1, 1)
	np.testing.assert_array_equal(got, packed[40:60])
	assert got.shape == (20, 2)
	lowres = np.packbits(bits["c2"][::100], axis=1, bitorder="little")
	got = g._query_bytes("c2", 0, 399, 200, 100)
	np.testing.assert_array_equal(got, lowres[0:4:2])


def test_bitmap_files_of_the_real_index_are_single_block(anchored_reader):
	"""The synthetic index is small enough for one BGZF block per bitmap: the block table is just (0, 0)."""
	g = anchored_reader["alpha"]
	for s in g.steps:
		assert g.blocks[s].tolist() == [(0, 0)]


# --------------------------------------------------------------- bitsum binning


def _hist(values, n):
	return np.bincount(values, minlength=n + 1)


def test_bin_bitsum_uses_min_bin_count_bins_for_short_chromosomes():
	"""Chromosomes shorter than min_bin_count * max_bin_kbp kb are split into len // min_bin_count wide bins.

	1005 positions, 100 minimum bins -> bins of 10: 100 full bins and one 5-wide bin, the end clipped to the length.
	"""
	import types

	g = pidx.Genome.__new__(pidx.Genome)
	g.params = {"max_bin_kbp": 200, "min_bin_count": 100}
	g.ngenomes = 5
	rng = np.random.default_rng(3)
	bitsum = rng.integers(1, 6, 1005)
	out = g.bin_bitsum(bitsum)
	assert out.shape == (101, 6)
	assert list(out.columns) == list(range(6))
	assert [s for s, e in out.index][:3] == [0, 10, 20] and out.index[-1] == (1000, 1005)
	assert out.dtypes.map(lambda t: t.kind).eq("i").all()
	for (s, e), row in out.iterrows():
		np.testing.assert_array_equal(row.to_numpy(), _hist(bitsum[s:e], 5))
	assert out.to_numpy().sum() == 1005


def test_bin_bitsum_floors_the_bin_width():
	"""1099 positions with 100 minimum bins: the width is floor(1099 / 100) = 10, giving 110 bins, the last one 9 wide (1090-1099)."""
	g = pidx.Genome.__new__(pidx.Genome)
	g.params = {"max_bin_kbp": 200, "min_bin_count": 100}
	g.ngenomes = 4
	bitsum = np.arange(1099) % 4 + 1
	out = g.bin_bitsum(bitsum)
	assert len(out) == 110
	assert out.index[1] == (10, 20) and out.index[-1] == (1090, 1099)
	assert out.iloc[-1].tolist() == _hist(bitsum[1090:1099], 4).tolist()
	assert out.to_numpy().sum() == 1099


def test_bin_bitsum_uses_max_bin_width_for_long_chromosomes():
	"""When the chromosome has at least min_bin_count bins of max_bin_kbp, the bin width is max_bin_kbp * 1000."""
	g = pidx.Genome.__new__(pidx.Genome)
	g.params = {"max_bin_kbp": 1, "min_bin_count": 2}
	g.ngenomes = 3
	bitsum = np.array([1, 2, 3] * 833 + [2])  # 2500 positions
	out = g.bin_bitsum(bitsum)
	assert list(out.index) == [(0, 1000), (1000, 2000), (2000, 2500)]
	for (s, e), row in out.iterrows():
		np.testing.assert_array_equal(row.to_numpy(), _hist(bitsum[s:e], 3))


def test_bin_bitsum_matches_bitmap_through_the_real_index(anchored_reader, synthetic_genomes):
	"""The stored bitsum bins of every chromosome equal the occupancy histogram of the brute-force bitmap.

	26-wide bins for chr1 of alpha (2636 positions // 100) and so on; the last bin is clipped.
	"""
	g = synthetic_genomes
	for genome in g.names:
		reader = anchored_reader[genome]
		for chrom in g.data.seqs[genome]:
			occ = g.oracle_bits(genome, chrom).sum(axis=1)
			binlen = len(occ) // 100
			rows = reader.bitsum_bins.loc[chrom]
			assert list(rows.index) == list(range(0, len(occ), binlen))
			for start, row in rows.iterrows():
				np.testing.assert_array_equal(row.to_numpy(), _hist(occ[start : start + binlen], 5))


def test_genome_bitsum_count_is_value_counts():
	"""Genome.bitsum_count counts each occupancy value."""
	g = pidx.Genome.__new__(pidx.Genome)
	got = g.bitsum_count([3, 1, 3, 3, 2])
	assert got.to_dict() == {3: 3, 1: 1, 2: 1}


# ----------------------------------------------------------- KMC-side byte packing


def test_load_kmc_opens_string_databases_and_passes_objects_through(monkeypatch, tmp_path):
	"""_load_kmc opens each path with OpenForRA and keeps already-open databases as they are."""
	import panagram.extra

	fake_dbs = [H.FakeKmcDb([1]), H.FakeKmcDb([2])]
	mod = H.fake_kmc_module(fake_dbs)
	monkeypatch.setitem(sys.modules, "panagram.extra.py_kmc_api", mod)
	monkeypatch.setattr(panagram.extra, "py_kmc_api", mod, raising=False)
	g = H.make_genome(tmp_path)
	already = H.FakeKmcDb([9])
	got = g._load_kmc(["kmc/a.bitvec0", already])
	assert got[0].opened == "kmc/a.bitvec0"
	assert got[1] is already
	assert g.kmc is mod


def test_load_kmc_reports_a_missing_py_kmc_api(monkeypatch, tmp_path):
	"""Without the compiled extension a ModuleNotFoundError naming py_kmc_api is raised."""
	import panagram.extra

	monkeypatch.delattr(panagram.extra, "py_kmc_api", raising=False)
	monkeypatch.setitem(sys.modules, "panagram.extra.py_kmc_api", None)
	g = H.make_genome(tmp_path)
	with pytest.raises(ModuleNotFoundError, match="py_kmc_api"):
		g._load_kmc(["x"])


def _kmc_genome(tmp_path, ngenomes, dbs):
	g = H.make_genome(tmp_path, ngenomes=ngenomes)
	g.kmc = H.fake_kmc_module()
	g.kmc_dbs = dbs
	return g


def test_query_kmc_bytes_takes_little_endian_bytes_of_each_counter(tmp_path):
	"""A counter 0x04030201 contributes bytes 1, 2, 3, 4 (least significant first); up to 4 bytes per database."""
	vals = [0x04030201, 0x08070605]
	g = _kmc_genome(tmp_path, 32, [H.FakeKmcDb(vals)])
	assert g.nbytes == 4
	assert g._query_kmc_bytes(0, "ACGT").tolist() == [[1, 2, 3, 4], [5, 6, 7, 8]]
	g = _kmc_genome(tmp_path, 3, [H.FakeKmcDb(vals)])
	assert g._query_kmc_bytes(0, "ACGT").tolist() == [[1], [5]]
	g = _kmc_genome(tmp_path, 13, [H.FakeKmcDb(vals)])
	assert g._query_kmc_bytes(0, "ACGT").tolist() == [[1, 2], [5, 6]]


def test_query_kmc_bytes_last_database_holds_the_remainder(tmp_path):
	"""With 40 genomes (5 bytes) database 0 gives 4 bytes and the last gives the remaining 1; at 64 genomes both give 4."""
	vals = [0x04030201]
	g = _kmc_genome(tmp_path, 40, [H.FakeKmcDb(vals), H.FakeKmcDb(vals)])
	assert g.nbytes == 5
	assert g._query_kmc_bytes(0, "A").shape == (1, 4)
	assert g._query_kmc_bytes(1, "A").tolist() == [[1]]
	g = _kmc_genome(tmp_path, 64, [H.FakeKmcDb(vals), H.FakeKmcDb(vals)])
	assert g._query_kmc_bytes(1, "A").shape == (1, 4)


def test_write_bitmap_writes_each_step_and_returns_the_bits(tmp_path):
	"""_write_bitmap appends all rows to the full file and every 100th to the low-res file, and returns genome bits."""
	rng = np.random.default_rng(1)
	vals = rng.integers(0, 32, 250).astype("uint32")
	g = _kmc_genome(tmp_path, 5, [H.FakeKmcDb(vals)])
	g.bitmaps = {1: H.Recorder(), 100: H.Recorder()}
	g.seq_lens, g.bitmap_lens = {0: {}}, {1: 0, 100: 0}
	bits = g._write_bitmap("c1", "A" * 300)
	assert g.bitmaps[1].data == vals.astype("uint8").tobytes()
	assert g.bitmaps[100].data == vals[::100].astype("uint8").tobytes()
	assert g.bitmap_lens == {1: 250, 100: 3}
	want = np.array([[(v >> j) & 1 for j in range(5)] for v in vals])
	np.testing.assert_array_equal(bits, want)


def test_write_bitmap_joins_two_databases_side_by_side(tmp_path):
	"""40 genomes need two KMC databases: bytes 0-3 come from the first, byte 4 from the second."""
	rng = np.random.default_rng(2)
	v0 = rng.integers(0, 2**32, 120, dtype="uint64").astype("uint32")
	v1 = rng.integers(0, 256, 120).astype("uint32")
	g = _kmc_genome(tmp_path, 40, [H.FakeKmcDb(v0), H.FakeKmcDb(v1)])
	g.bitmaps = {1: H.Recorder(), 100: H.Recorder()}
	g.seq_lens, g.bitmap_lens = {0: {}}, {1: 0, 100: 0}
	bits = g._write_bitmap("c1", "A")
	rows = np.frombuffer(g.bitmaps[1].data, "uint8").reshape(120, 5)
	for i in (0, 57, 119):
		want = [int(v0[i]) >> j & 1 for j in range(32)] + [int(v1[i]) >> j & 1 for j in range(8)]
		assert bits[i].tolist() == want
		assert rows[i].tolist() == list(int(v0[i]).to_bytes(4, "little")) + [int(v1[i])]
	assert len(g.bitmaps[100].data) == 2 * 5
	assert g.bitmap_lens[1] == 120


def test_run_anchor_skips_non_anchor_genomes(tmp_path, caplog):
	"""run_anchor on a non-anchor genome logs a message and writes nothing."""
	idx = H.stub_index(tmp_path, ["a", "b"])
	g = pidx.Genome(idx, 1, "b")
	with caplog.at_level("INFO", logger="panagram.index"):
		g.run_anchor(["x"])
	assert "Skipping non-anchor genome 'b'" in caplog.text
	assert not (tmp_path / "anchor").exists()


# ---------------------------------------------------------------------- iter_fasta


@pytest.mark.parametrize("ext", ["fa", "fa.gz", "fa.bgz"])
def test_iter_fasta_reads_plain_gzip_and_bgzip(tmp_path, synthetic_genomes, ext):
	"""iter_fasta yields every record with its id and sequence from plain, gzip or bgzip FASTA."""
	chroms = synthetic_genomes.data.seqs["alpha"]
	plain = tmp_path / "a.fa"
	synthetic.write_fasta(plain, chroms)
	path = tmp_path / f"a.{ext}"
	if ext == "fa.gz":
		path.write_bytes(gzip.compress(plain.read_bytes()))
	elif ext == "fa.bgz":
		with bgzf.BgzfWriter(str(path), "wb") as w:
			w.write(plain.read_bytes())
	else:
		path = plain
	g = H.make_genome(tmp_path, write=True)
	g.fasta = str(path)
	recs = list(g.iter_fasta())
	assert [r.id for r in recs] == list(chroms)
	assert [str(r.seq) for r in recs] == list(chroms.values())


def test_close_closes_every_bitmap_handle(tmp_path):
	"""Genome.close closes the open bitmap readers."""
	g, _ = H.build_bitmap_genome(tmp_path, 5, {"c1": 300})
	g.close()
	for s, f in g.bitmaps.items():
		assert f._handle.closed

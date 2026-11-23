"""Helpers for the panagram.index / __main__ / figs tests (not fixtures; plain functions)."""

import os
import types
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import bgzf

import synthetic

DEFAULT_PARAMS = {
	"prefix": "",
	"k": 15,
	"lowres_step": 100,
	"max_bin_kbp": 200,
	"min_bin_count": 100,
	"gff_gene_types": ["gene"],
	"gff_anno_types": None,
	"gff_name": "Name",
}


def stub_index(tmp_path, names, fastas=None, gffs=None, **params):
	"""Minimal stand-in for panagram.index.Index with what Genome.__init__ reads.

	`prefix` is `tmp_path`, so a Genome writes under tmp_path/anchor/<name>/ (create it first).
	"""
	p = dict(DEFAULT_PARAMS, **params)
	p["prefix"] = str(tmp_path)
	samples = pd.DataFrame(
		{
			"fasta": [(fastas or {}).get(n, pd.NA) for n in names],
			"gff": [(gffs or {}).get(n, pd.NA) for n in names],
			"id": range(len(names)),
		},
		index=pd.Index(names, name="name"),
	)
	return types.SimpleNamespace(
		samples=samples, params=p, genome_names=samples.index, lowres_step=p["lowres_step"]
	)


def make_genome(tmp_path, ngenomes=5, name="g0", write=True, fasta=None, gff=None, anchor=True, **params):
	"""A panagram.index.Genome over a stub index with `ngenomes` genomes, anchor dir created."""
	from panagram.index import Genome

	names = [f"g{i}" for i in range(ngenomes)]
	idx = stub_index(tmp_path, names, **params)
	(Path(tmp_path) / "anchor" / name).mkdir(parents=True, exist_ok=True)
	return Genome(idx, 0, name, fasta=fasta, gff=gff, anchor=anchor, write=write)


def write_tsv_fasta(path, chroms):
	synthetic.write_fasta(path, chroms)


def prepared_copy(genomes, dest, **index_kwargs):
	"""Copy the synthetic genomes to `dest` and write config.yaml, samples.tsv, Snakefile in-process.

	Same as `panagram index samples.tsv -k K --prepare` but needs no KMC binaries.
	"""
	import shutil
	from panagram.index import Index

	shutil.copytree(genomes.root, dest)
	cwd = os.getcwd()
	try:
		Index(str(Path(dest) / "samples.tsv"), k=genomes.k, prepare=True, **index_kwargs).run()
	finally:
		os.chdir(cwd)
	return Path(dest)


def open_write(path, **kw):
	"""Open an initialized index directory in write mode (restores the working directory)."""
	from panagram.index import Index

	cwd = os.getcwd()
	try:
		return Index(str(path), mode="w", **kw)
	finally:
		os.chdir(cwd)


def write_bgzf_with_gzi(path, data, chunk):
	"""Write `data` (bytes) as BGZF in blocks of `chunk` bytes each, plus a .gzi index in samtools layout.

	Returns the list of (compressed start, uncompressed start) of every data block.
	"""
	with bgzf.BgzfWriter(str(path), "wb") as w:
		for i in range(0, len(data), chunk):
			w.write(data[i : i + chunk])
			w.flush()
	blocks = []
	with open(path, "rb") as h:
		for start, raw_len, data_start, data_len in bgzf.BgzfBlocks(h):
			if data_len > 0:
				blocks.append((start, data_start))
	with open(str(path)[:-2] + "gzi", "wb") as f:
		np.array([len(blocks) - 1], dtype="uint64").tofile(f)
		np.array(blocks[1:], dtype="uint64").tofile(f)
	return blocks


def build_bitmap_genome(tmp_path, ngenomes, sizes, chunk=50, seed=7, steps=(1, 100)):
	"""A write-mode Genome whose bitmaps are random bits in multi-block BGZF files.

	Returns (genome, bits) with bits[chrom] the (size x ngenomes) uint8 truth matrix. Uses
	the on-disk layout the reader expects: all chromosomes concatenated, one row of
	ceil(ngenomes / 8) bytes per position (little bit order), the low-resolution file
	holding every steps[1]-th row of each chromosome.
	"""
	rng = np.random.default_rng(seed)
	g = make_genome(tmp_path, ngenomes=ngenomes, name="g0", write=True)
	names = list(sizes)
	chrs = pd.DataFrame(
		{"name": names, "id": range(len(names)), "size": [sizes[n] for n in names], "gene_count": 0}
	)
	chrs.to_csv(g.chrs_fname, sep="\t", index=False)
	g.load_chrs()
	bits = {n: rng.integers(0, 2, (sizes[n], ngenomes), dtype="uint8") for n in names}
	g.blocks, g.bitmaps = {}, {}
	for step in steps:
		rows = np.concatenate([np.packbits(bits[n][::step], axis=1, bitorder="little") for n in names])
		write_bgzf_with_gzi(g.bitmap_gz_fname(step), rows.tobytes(), chunk)
		g.blocks[step] = g.load_bgz_blocks(g.bitmap_gzi_fname(step))
		g.bitmaps[step] = bgzf.BgzfReader(g.bitmap_gz_fname(step), "rb")
	return g, bits


class FakeKmcDb:
	"""Stands in for py_kmc_api.KMCFile: GetCountersForRead appends fixed counters."""

	def __init__(self, values):
		self.values = list(values)
		self.opened = None

	def OpenForRA(self, name):
		self.opened = name

	def GetCountersForRead(self, seq, vec):
		vec.extend(self.values)


def fake_kmc_module(dbs=None):
	"""Namespace with KMCFile (hands out `dbs` in order) and CountVec (a list)."""
	dbs = list(dbs or [])
	return types.SimpleNamespace(KMCFile=lambda: dbs.pop(0), CountVec=list)


class Recorder:
	"""Write-only file stand-in that remembers the bytes written."""

	def __init__(self):
		self.data = b""

	def write(self, b):
		self.data += b

	def close(self):
		pass


def occupancy(oracle):
	"""Number of genomes sharing each position, from a brute-force (positions x genomes) matrix."""
	return np.asarray(oracle).sum(axis=1)

"""Helpers for the introgression tests: reference binning, inline executor, file readers.

The reference functions are short, loop-based statements of what the README describes
(kmer similarity of an anchor to each accession per bin), computed from the brute-force
canonical k-mer oracle of the synthetic genomes, never from panagram's own bitmap code.
"""

import concurrent.futures
import copy
import sys

import numpy as np
import pandas as pd

STEP = 100  # panagram's default bitmap step
BIN = 500  # bin size that gives 5 sampled rows per bin on the synthetic chromosomes


def oracle_bitmap(genomes, genome, chrom, step=STEP):
	"""Positions x accessions 0/1 table sampled every `step` k-mer starts (what Genome.query returns)."""
	bits = genomes.oracle_bits(genome, chrom)
	return pd.DataFrame(bits[::step], index=range(0, len(bits), step), columns=genomes.names)


def reference_bins(bitmap, binlen, rmf=False, rmu_keep=None):
	"""Binned similarity table (accessions x bin starts) from a positions x accessions 0/1 table.

	Per bin: count, for every accession, the sampled rows with a 1, divide by the largest count
	in the bin. rmf drops rows that are 1 in every accession first (empty bins become 1).
	rmu_keep: list of columns; rows with no 1 in those columns get a 1 in all of them first.
	"""
	rows = bitmap.copy()
	if rmu_keep is not None:
		for pos in rows.index:
			if rows.loc[pos, rmu_keep].sum() == 0:
				rows.loc[pos, rmu_keep] = 1
	out = {}
	for b in sorted({p // binlen for p in rows.index}):
		sel = rows[[p // binlen == b for p in rows.index]]
		if rmf:
			sel = sel[[not bool((r == 1).all()) for _, r in sel.iterrows()]]
		if len(sel) == 0:
			out[b * binlen] = pd.Series(1.0, index=rows.columns)
		else:
			counts = sel.sum(axis=0).astype(float)
			out[b * binlen] = counts / counts.max()
	return pd.DataFrame(out)


def reference_calls_2way(sims, anchor, ref, thr):
	"""Bin starts where the anchor's similarity to `ref` is below thr."""
	return [int(c) for c in sims.columns if sims.loc[ref, c] < thr]


def runs_to_bed(starts, bin_size, chrom, name):
	"""Merge adjacent bin starts into inclusive BED rows (chrom, start, end, name)."""
	rows = []
	for s in sorted(starts):
		if rows and rows[-1][2] + 1 == s:
			rows[-1][2] = s + bin_size - 1
		else:
			rows.append([chrom, s, s + bin_size - 1, name])
	return [tuple(r) for r in rows]


def read_bed_rows(path):
	"""BED file as a list of (chrom, start, end, name) tuples; [] for an empty file."""
	out = []
	with open(path) as f:
		for line in f:
			if line.strip():
				c, s, e, n = line.rstrip("\n").split("\t")[:4]
				out.append((c, int(s), int(e), n))
	return out


def write_groups(path, mapping):
	"""Write a group.tsv (name, group) from a dict accession -> group."""
	with open(path, "w") as f:
		f.write("name\tgroup\n")
		for name, group in mapping.items():
			f.write(f"{name}\t{group}\n")
	return path


class InlineExecutor:
	"""Stand-in for ProcessPoolExecutor that runs jobs at submit time in this process.

	Arguments are deep-copied first, as pickling would do, because the worker pops keys
	from the preprocessing dictionary it is given.
	"""

	def __init__(self, max_workers=None):
		self.max_workers = max_workers

	def __enter__(self):
		return self

	def __exit__(self, *exc):
		return False

	def submit(self, fn, *args, **kwargs):
		fut = concurrent.futures.Future()
		try:
			fut.set_result(fn(*copy.deepcopy(args), **copy.deepcopy(kwargs)))
		except BaseException as e:
			fut.set_exception(e)
		return fut


def can_render_images(tmp_path):
	"""True when plotly can write an SVG here (kaleido and Chrome present)."""
	try:
		import plotly.express as px

		px.imshow(np.zeros((2, 2))).write_image(str(tmp_path / "probe.svg"))
		return True
	except Exception:
		return False

"""Re-running snakemake on a finished index after an anchor job re-runs.

Rule umaps must be scheduled again so the UMAP csv files are refreshed, and a stale
anchor/<sample>/failed marker must not make the re-run fail.
"""

import os
import shutil
import sys

import pytest

from conftest import ANCHOR_CORES, _run

GENOME = "alpha"
UMAP_FILES = ("chrom_umaps.csv", "genome_umap.csv")


def _copy(anchored_index, dest):
	shutil.copytree(anchored_index, dest)
	return dest


def _forget_anchor(index, genome):
	"""Delete one genome's anchor outputs the way an anchor-only re-run test does."""
	adir = index / "anchor" / genome
	for p in list(adir.glob("bitmap.*")) + [adir / "chrs.tsv"]:
		p.unlink()


def _rerun(index, env):
	return _run([sys.executable, "-m", "snakemake", "--cores", str(ANCHOR_CORES), "all"], index, env)


def _stamps(index):
	return {
		(g.name, f): (g / f).stat().st_mtime_ns
		for g in (index / "anchor").iterdir()
		for f in UMAP_FILES
		if (g / f).exists()
	}


def test_anchor_only_rerun_refreshes_the_umap_files(anchored_index, panagram_overlay, tmp_path):
	"""Deleting one genome's anchor outputs and re-running regenerates every UMAP csv."""
	index = _copy(anchored_index, tmp_path / "index")
	before = _stamps(index)
	assert before
	_forget_anchor(index, GENOME)
	_rerun(index, panagram_overlay)
	after = _stamps(index)
	assert after.keys() == before.keys()
	stale = [k for k in before if after[k] <= before[k]]
	assert not stale, f"UMAP files not rewritten: {stale}"
	assert (index / "anchor" / GENOME / "chrs.tsv").exists()


def test_stale_failed_marker_does_not_fail_the_rerun(anchored_index, panagram_overlay, tmp_path):
	"""A failed marker from an earlier run is removed when the anchor job starts and the re-run succeeds."""
	index = _copy(anchored_index, tmp_path / "index")
	_forget_anchor(index, GENOME)
	failed = index / "anchor" / GENOME / "failed"
	failed.touch()
	_rerun(index, panagram_overlay)
	assert not failed.exists()
	assert (index / "anchor" / GENOME / "chrs.tsv").exists()
	for f in UMAP_FILES:
		assert (index / "anchor" / GENOME / f).stat().st_size > 0


def test_failed_marker_is_written_on_failure_and_cleared_on_success():
	"""The Snakefile text keeps the marker's lifecycle: removed before run_anchor, written on error, removed after."""
	import panagram

	text = (os.path.join(os.path.dirname(panagram.__file__), "workflow", "Snakefile"))
	src = open(text).read()
	start = src.index("rule anchor:")
	block = src[start : src.index("rule umaps:")]
	assert block.count("os.remove(failed)") == 2
	assert block.index("os.remove(failed)") < block.index("run_anchor") < block.rindex("os.remove(failed)")
	assert 'open(failed, "w")' in block
	assert "chrs.tsv" in src[src.index("rule umaps:") : src.index("rule kmc_bitvec:")]


def test_rerun_after_fixing_a_failed_anchor_job(fresh_prepared_index, panagram_overlay):
	"""An anchor job that fails on a bad GFF leaves the marker; after the GFF is fixed the re-run completes."""
	import subprocess

	index = fresh_prepared_index
	gff = index / "GFFS" / f"{GENOME}.gff"
	good = gff.read_text()
	gff.write_text("chr1\tsrc\tgene\tabc\tdef\t.\t+\t.\tID=g1\n")
	cmd = [sys.executable, "-m", "snakemake", "--cores", str(ANCHOR_CORES), "all"]
	r = subprocess.run(cmd, cwd=index, env=panagram_overlay, capture_output=True, text=True, timeout=600)
	assert r.returncode != 0
	assert (index / "anchor" / GENOME / "failed").exists()
	gff.write_text(good)
	_rerun(index, panagram_overlay)
	assert not (index / "anchor" / GENOME / "failed").exists()
	for f in UMAP_FILES:
		assert (index / "anchor" / GENOME / f).stat().st_size > 0

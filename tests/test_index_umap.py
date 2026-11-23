"""UMAP embedding, DBSCAN clustering and the per-genome UMAP files of panagram.index.Genome."""

import os
import shutil

import numpy as np
import pandas as pd
import pytest

import helpers_index as H
from conftest import open_index
from panagram import index as pidx


def _paircounts(n_per_group=20, seed=0):
	"""Two well-separated groups of bins over 5 genomes, indexed by (chrom, start)."""
	rng = np.random.default_rng(seed)
	a = np.abs(rng.normal([1, 0, 0, 0, 0], 0.03, (n_per_group, 5)))
	b = np.abs(rng.normal([0, 0, 0, 0, 1], 0.03, (n_per_group, 5)))
	idx = pd.MultiIndex.from_arrays(
		[["c1"] * (2 * n_per_group), np.arange(2 * n_per_group) * 1000], names=["chrom", "start"]
	)
	return pd.DataFrame(np.vstack([a, b]), index=idx)


def test_run_umap_with_too_few_bins_gives_a_zero_embedding(tmp_path):
	"""At most `neighbors` bins cannot be embedded: every bin gets umap1 = umap2 = cluster = 0, and end = start + bin_size."""
	g = H.make_genome(tmp_path)
	pc = _paircounts(2)  # 4 bins <= neighbors (4)
	out = g.run_umap(pc, pidx.UMAP(neighbors=4, bin_size=777))
	assert list(out.columns) == ["chrom", "start", "end", "umap1", "umap2", "cluster"]
	assert out["chrom"].tolist() == ["c1"] * 4
	assert out["start"].tolist() == [0, 1000, 2000, 3000]
	assert out["end"].tolist() == [777, 1777, 2777, 3777]
	assert (out[["umap1", "umap2", "cluster"]] == 0).all().all()


def test_run_umap_embeds_and_clusters(tmp_path):
	"""With enough bins the output has a 2-D embedding that separates the groups, and clusters are DBSCAN of that embedding."""
	from sklearn.cluster import DBSCAN

	g = H.make_genome(tmp_path)
	pc = _paircounts(20)
	args = pidx.UMAP(neighbors=5, dist=0, eps=1, samples=1, bin_size=1000)
	out = g.run_umap(pc, args)
	assert list(out.columns) == ["chrom", "start", "end", "umap1", "umap2", "cluster"]
	assert len(out) == 40 and out["start"].tolist() == pc.index.get_level_values("start").tolist()
	assert (out["end"] - out["start"] == 1000).all()
	emb = out[["umap1", "umap2"]].to_numpy()
	assert np.isfinite(emb).all() and emb.std(axis=0).min() > 0
	# the two input groups stay apart: every point's 3 nearest neighbours in the embedding come from its own group
	dist = np.linalg.norm(emb[:, None] - emb[None], axis=2)
	np.fill_diagonal(dist, np.inf)
	group = np.arange(40) // 20
	for i in range(40):
		assert (group[np.argsort(dist[i])[:3]] == group[i]).all(), i
	np.testing.assert_array_equal(out["cluster"].to_numpy(), DBSCAN(eps=1, min_samples=1).fit_predict(emb))
	again = g.run_umap(pc, args)
	pd.testing.assert_frame_equal(out, again)  # fixed random_state: repeatable


def test_run_umap_falls_back_to_zeros_when_umap_raises(tmp_path, monkeypatch):
	"""An exception inside UMAP is logged and the chromosome gets the zero embedding."""

	class Boom:
		def __init__(self, **kw):
			pass

		def fit_transform(self, x):
			raise RuntimeError("no embedding")

	monkeypatch.setattr(pidx.umap, "UMAP", Boom)
	g = H.make_genome(tmp_path)
	out = g.run_umap(_paircounts(10), pidx.UMAP(neighbors=3, bin_size=1000))
	assert len(out) == 20 and (out[["umap1", "umap2", "cluster"]] == 0).all().all()
	assert out["end"].tolist() == [s + 1000 for s in out["start"]]


def test_run_umap_passes_parameters_to_umap_and_dbscan(tmp_path, monkeypatch):
	"""neighbors and dist go to UMAP (2 components, one job, random_state 42); eps and samples go to DBSCAN."""
	seen = {}

	class FakeUmap:
		def __init__(self, **kw):
			seen["umap"] = kw

		def fit_transform(self, x):
			seen["shape"] = x.shape
			return np.column_stack([np.arange(len(x)) * 10.0, np.zeros(len(x))])

	class FakeDbscan:
		def __init__(self, **kw):
			seen["dbscan"] = kw

		def fit_predict(self, e):
			return np.arange(len(e)) % 2

	monkeypatch.setattr(pidx.umap, "UMAP", FakeUmap)
	monkeypatch.setattr(pidx, "DBSCAN", FakeDbscan)
	g = H.make_genome(tmp_path)
	out = g.run_umap(_paircounts(5), pidx.UMAP(neighbors=3, dist=0.25, eps=2.5, samples=7, bin_size=50))
	assert seen["umap"] == {"n_neighbors": 3, "min_dist": 0.25, "n_components": 2, "n_jobs": 1, "random_state": 42}
	assert seen["dbscan"] == {"eps": 2.5, "min_samples": 7}
	assert seen["shape"] == (10, 5)
	assert out["umap1"].tolist() == [i * 10.0 for i in range(10)]
	assert out["cluster"].tolist() == [i % 2 for i in range(10)]
	assert out["end"].tolist() == [s + 50 for s in out["start"]]


def test_load_umaps_reads_both_files(tmp_path):
	"""load_umaps reads the chromosome embedding indexed by chrom and the genome embedding as a plain table."""
	g = H.make_genome(tmp_path)
	chrom = pd.DataFrame({"chrom": ["a", "b"], "start": [0, 0], "end": [10, 10], "umap1": [0.5, 1.5], "umap2": [2.0, 3.0], "cluster": [0, 1]})
	chrom.to_csv(g.chrom_umaps_filename, index=False)
	chrom.to_csv(g.genome_umap_filename, index=False)
	g.load_umaps()
	assert g.chrom_umaps.index.name == "chrom" and g.chrom_umaps.index.tolist() == ["a", "b"]
	assert g.chrom_umaps["umap2"].tolist() == [2.0, 3.0]
	assert list(g.genome_umap.columns) == list(chrom.columns) and len(g.genome_umap) == 2


def test_load_umaps_without_files_sets_none(tmp_path):
	"""With neither UMAP file present both attributes are None."""
	g = H.make_genome(tmp_path)
	g.load_umaps()
	assert g.chrom_umaps is None and g.genome_umap is None


def test_umap_files_of_the_real_index(anchored_reader, synthetic_genomes):
	"""Synthetic chromosomes fit in one 100 kb bin, so every embedding is the zero fallback: one row per chromosome."""
	for genome in synthetic_genomes.names:
		g = anchored_reader[genome]
		chroms = list(synthetic_genomes.data.seqs[genome])
		assert g.chrom_umaps.index.tolist() == chroms
		assert g.chrom_umaps["start"].tolist() == [0] * len(chroms)
		assert g.chrom_umaps["end"].tolist() == [100000] * len(chroms)
		assert (g.chrom_umaps[["umap1", "umap2", "cluster"]] == 0).all().all()
		assert g.genome_umap["chrom"].tolist() == chroms


@pytest.fixture
def private_copy(anchored_index, tmp_path):
	dest = tmp_path / "index"
	shutil.copytree(anchored_index, dest, ignore=shutil.ignore_patterns(".snakemake", "kmc", "tmp", "logs"))
	return dest


def test_write_umaps_bins_the_low_resolution_bitmap(private_copy, synthetic_genomes):
	"""With 500-bp bins, write_umaps embeds each chromosome's 100-step bitmap: one row per bin, files written and reloadable.

	alpha chr1 has 2636 positions, so 27 low-resolution rows fall into 6 bins of 500 positions (5 rows each, the last 2).
	"""
	idx = open_index(private_copy)
	os.chdir(private_copy)
	g = idx["alpha"]
	idx.chrom_umap.bin_size = 500
	idx.genome_umap.bin_size = 1000
	g.write_umaps()
	chroms = g.chrom_umaps
	assert chroms.index.unique().tolist() == ["chr1", "chr2"]
	c1 = chroms.loc["chr1"]
	assert c1["start"].tolist() == [0, 500, 1000, 1500, 2000, 2500]
	assert (c1["end"] - c1["start"] == 500).all()
	assert len(chroms.loc["chr2"]) == 4  # 1686 positions -> 17 rows -> bins of 5 rows
	assert list(g.genome_umap.columns) == ["chrom", "start", "end", "umap1", "umap2", "cluster"]
	assert g.genome_umap["chrom"].tolist() == ["chr1"] * 3 + ["chr2"] * 2
	assert (g.genome_umap["end"] - g.genome_umap["start"] == 1000).all()
	idx.close()
	again = open_index(private_copy)
	os.chdir(private_copy)
	pd.testing.assert_frame_equal(again["alpha"].chrom_umaps, chroms)
	assert len(again["alpha"].genome_umap) == 5
	again.close()

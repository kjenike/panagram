"""Four regressions of the single-process umaps rule against the original Snakefile.

An empty samples table, duplicate genome names, invalid UMAP parameters and the head-start
writing UMAP files into a directory snakemake empties when rule umaps starts.
"""

import os
import shutil
import subprocess
import sys

import pandas as pd
import pytest

from conftest import ANCHOR_CORES, K, _run, open_index, panagram_cli


def _prepared(genomes, dest, env, *args, samples=None):
    shutil.copytree(genomes.root, dest)
    if samples is not None:
        (dest / "samples.tsv").write_text(samples)
    _run(panagram_cli("index", "samples.tsv", "-k", K, *args, "--prepare"), dest, env)
    return dest


def _snakemake(index, env):
    cmd = [sys.executable, "-m", "snakemake", "--cores", str(ANCHOR_CORES), "all"]
    return subprocess.run(cmd, cwd=index, env=env, capture_output=True, text=True, timeout=600)


def test_header_only_samples_table_builds_genome_dist(synthetic_genomes, panagram_overlay, tmp_path):
    """A samples table with only a header exits 0 and writes genome_dist.tsv (kmc_count's threads divided by zero)."""
    index = _prepared(synthetic_genomes, tmp_path / "index", panagram_overlay, samples="name\tfasta\tgff\n")
    r = _snakemake(index, panagram_overlay)
    assert r.returncode == 0, r.stderr[-2000:]
    assert (index / "genome_dist.tsv").exists()
    assert "ZeroDivisionError" not in r.stderr


def test_duplicate_anchor_genomes_build_once(synthetic_genomes, panagram_overlay, tmp_path):
    """`--anchor_genomes alpha alpha` exits 0 and writes alpha's UMAP files (duplicate output pattern in rule umaps)."""
    index = _prepared(synthetic_genomes, tmp_path / "index", panagram_overlay, "--anchor_genomes", "alpha", "alpha")
    r = _snakemake(index, panagram_overlay)
    assert r.returncode == 0, r.stderr[-2000:]
    assert "Duplicate output file pattern" not in r.stderr
    for f in ("chrom_umaps.csv", "genome_umap.csv"):
        assert (index / "anchor" / "alpha" / f).stat().st_size > 0


@pytest.mark.parametrize("args", [("--chrom_umap.dist", "-1"), ("--chrom_umap.neighbors", "1")])
def test_invalid_umap_parameters_write_zeros_and_exit_zero(synthetic_genomes, panagram_overlay, tmp_path, args):
    """An invalid UMAP parameter is caught per fit as in the original: zero embeddings, exit 0, not a failed warm-up."""
    index = _prepared(synthetic_genomes, tmp_path / "index", panagram_overlay, *args)
    r = _snakemake(index, panagram_overlay)
    assert r.returncode == 0, r.stderr[-2000:]
    umaps = pd.read_csv(index / "anchor" / "alpha" / "chrom_umaps.csv")
    assert (umaps[["umap1", "umap2"]] == 0).all().all()


def test_head_start_survives_snakemake_emptying_the_umap_outputs(anchored_index, panagram_overlay, tmp_path, monkeypatch):
    """The head-start never writes into anchor/<genome>/: deleting those files while it works loses nothing."""
    from panagram import _umaps_job
    from panagram.index import Genome

    index = tmp_path / "index"
    shutil.copytree(anchored_index, index)
    control = tmp_path / "control"
    control.mkdir()
    for marker in ("go", "anchors"):
        (control / marker).touch()
    monkeypatch.chdir(index)
    original = Genome.write_umaps

    def write_after_snakemake_deleted_the_outputs(self, outdir=None):
        for name in ("chrom_umaps.csv", "genome_umap.csv"):
            path = index / "anchor" / self.name / name
            if path.exists():
                path.unlink()
        return original(self, outdir)

    monkeypatch.setattr(Genome, "write_umaps", write_after_snakemake_deleted_the_outputs)
    assert _umaps_job.main(early=str(control)) == 0
    for genome in (index / "anchor").iterdir():
        for name in ("chrom_umaps.csv", "genome_umap.csv"):
            assert not (genome / name).exists()
            assert (control / "stage" / genome.name / name).stat().st_size > 0


def test_write_umaps_to_another_directory_matches_the_default(anchored_reader, anchored_index, tmp_path):
    """write_umaps(outdir) writes the same two files as the default location, and nothing else changes."""
    genome = anchored_reader["alpha"]
    genome.init_read()
    genome.write_umaps(str(tmp_path))
    for name in ("chrom_umaps.csv", "genome_umap.csv"):
        assert (tmp_path / name).read_bytes() == (anchored_index / "anchor" / "alpha" / name).read_bytes()

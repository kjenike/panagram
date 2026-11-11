"""Entry point of rule umaps: warm up umap's numba kernels, then write every anchor's UMAPs."""
import os
import time
from types import SimpleNamespace

import yaml

from panagram.index import Genome, Index


def main():
    with open("config.yaml") as f:
        cfg = yaml.load(f, yaml.Loader) or {}
    chrom = cfg.get("chrom_umap") or {}
    Genome.warm_umap(SimpleNamespace(chrom_umap=SimpleNamespace(
        neighbors=chrom.get("neighbors", 4), dist=chrom.get("dist", 0))))

    index = Index(".", mode="w")
    for name in index.anchor_genomes:
        done = f"logs/anchor.{name}.benchmark.txt"
        failed = f"anchor/{name}/failed"
        while not os.path.exists(done):
            if os.path.exists(failed):
                raise RuntimeError(f"anchor job for {name} failed; no UMAPs")
            time.sleep(0.05)
        genome = index[name]
        genome.init_read()
        genome.write_umaps()


if __name__ == "__main__":
    main()

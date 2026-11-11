"""Entry point of rule umaps: warm up umap's numba kernels, then write every anchor's UMAPs."""
import os
import shutil
import sys
import time
from types import SimpleNamespace

import yaml

from panagram.index import Genome, Index


def _parent_alive(parent):
    return os.getppid() == parent


def main(early=None):
    """Run the umaps job. With `early` (a control directory), run as the head-start
    process of the Snakefile: warm up, wait for the go file, then write the UMAPs."""
    parent = os.getppid()
    with open("config.yaml") as f:
        cfg = yaml.load(f, yaml.Loader) or {}
    chrom = cfg.get("chrom_umap") or {}
    Genome.warm_umap(SimpleNamespace(chrom_umap=SimpleNamespace(
        neighbors=chrom.get("neighbors", 4), dist=chrom.get("dist", 0))))

    if early:
        while not os.path.exists(os.path.join(early, "go")):
            if not _parent_alive(parent):
                return 1
            time.sleep(0.01)

    index = Index(".", mode="w")
    for name in index.anchor_genomes:
        done = f"logs/anchor.{name}.benchmark.txt"
        failed = f"anchor/{name}/failed"
        while not os.path.exists(done):
            if os.path.exists(failed):
                raise RuntimeError(f"anchor job for {name} failed; no UMAPs")
            if early and not _parent_alive(parent):
                return 1
            time.sleep(0.05)
        genome = index[name]
        genome.init_read()
        genome.write_umaps()
    return 0


def start_early(control):
    """Spawn the head-start process (called from the Snakefile while it is parsed)."""
    import subprocess
    log = open(os.path.join(control, "log.txt"), "w")
    return subprocess.Popen([sys.executable, "-m", "panagram._umaps_job", "--early", control],
        stdout=log, stderr=subprocess.STDOUT)


def attach(control):
    """Rule umaps: release the head-start process and wait for it; True if it finished
    cleanly, False when there is none or it failed (the caller then runs the job itself)."""
    pidfile = os.path.join(control, "pid")
    if not control or not os.path.exists(pidfile):
        return False
    with open(os.path.join(control, "go"), "w"):
        pass
    done = os.path.join(control, "status")
    while not os.path.exists(done):
        time.sleep(0.01)
    with open(done) as f:
        ok = f.read().strip() == "0"
    shutil.rmtree(control, ignore_errors=True)
    return ok


def _early_entry(control):
    status = 1
    try:
        status = main(early=control) or 0
    finally:
        with open(os.path.join(control, "status.tmp"), "w") as f:
            f.write(str(status))
        os.replace(os.path.join(control, "status.tmp"), os.path.join(control, "status"))
    return status


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "--early":
        sys.exit(_early_entry(sys.argv[2]))
    if sys.argv[1:] == ["--attach"]:
        if attach(os.environ.get("PANAGRAM_UMAPS_CONTROL", "")):
            sys.exit(0)
    main()

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


def _stamp(path):
    """Identity of a file as it is now (None if absent), to tell a stale marker from a new one."""
    try:
        st = os.stat(path)
    except OSError:
        return None
    return (st.st_mtime_ns, st.st_ino, st.st_size)


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
    # Markers left by an earlier run are stale until an anchor job rewrites them; once
    # the anchors file exists (rule umaps is scheduled, so every anchor job is done) a
    # benchmark file that never changed belongs to an anchor job that did not re-run.
    stale = {}
    if early:
        for name in index.anchor_genomes:
            stale[name] = (_stamp(f"logs/anchor.{name}.benchmark.txt"), _stamp(f"anchor/{name}/failed"))
    for name in index.anchor_genomes:
        done = f"logs/anchor.{name}.benchmark.txt"
        failed = f"anchor/{name}/failed"
        old_done, old_failed = stale.get(name, (None, None))
        while True:
            now = _stamp(done)
            if now is not None and (now != old_done or (early and os.path.exists(os.path.join(early, "anchors")))):
                break
            now_failed = _stamp(failed)
            if now_failed is not None and now_failed != old_failed:
                raise RuntimeError(f"anchor job for {name} failed; no UMAPs")
            if early and not _parent_alive(parent):
                return 1
            time.sleep(0.05)
        genome = index[name]
        genome.init_read()
        genome.write_umaps()
        if early:
            # Snakemake deletes a rule's existing outputs when the rule starts, and rule
            # umaps starts after the anchor jobs: park the files for attach() to move in.
            stage = os.path.join(early, "stage", name)
            os.makedirs(stage, exist_ok=True)
            for path in (genome.chrom_umaps_filename, genome.genome_umap_filename):
                shutil.move(path, os.path.join(stage, os.path.basename(path)))
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
    with open(os.path.join(control, "anchors"), "w"):
        pass
    done = os.path.join(control, "status")
    while not os.path.exists(done):
        time.sleep(0.01)
    with open(done) as f:
        ok = f.read().strip() == "0"
    stage = os.path.join(control, "stage")
    if ok and os.path.isdir(stage):
        for name in os.listdir(stage):
            for fn in os.listdir(os.path.join(stage, name)):
                shutil.move(os.path.join(stage, name, fn), os.path.join("anchor", name, fn))
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

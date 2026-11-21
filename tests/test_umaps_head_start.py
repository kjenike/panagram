"""Life cycle of the umaps head-start process started while the Snakefile is parsed.

A re-run with nothing to do must leave no process and no panagram-umaps-* directory, a
head-start killed before it wrote its status must fail the build with a clear error
instead of hanging, and attach() removes the control directory on every exit path.
"""

import os
import shutil
import signal
import subprocess
import sys
import time

import pytest

from conftest import ANCHOR_CORES, _run

GENOME = "alpha"


def _copy(anchored_index, dest):
    shutil.copytree(anchored_index, dest)
    return dest


def _env(overlay, tmp):
    """Environment whose temp directory is private, so leftovers are this test's only."""
    tmp.mkdir()
    return dict(overlay, TMPDIR=str(tmp), TEMP=str(tmp), TMP=str(tmp))


def _snakemake():
    return [sys.executable, "-m", "snakemake", "--cores", str(ANCHOR_CORES), "all"]


def _head_starts(tmp):
    """Pids of running head-start processes whose control directory is under tmp."""
    psutil = pytest.importorskip("psutil")
    out = []
    for p in psutil.process_iter(["cmdline", "status"]):
        cmd = p.info["cmdline"] or []
        if any("_umaps_job" in c for c in cmd) and any(str(tmp) in c for c in cmd):
            if p.info["status"] != psutil.STATUS_ZOMBIE:
                out.append(p.pid)
    return out


def _wait_for(cond, seconds):
    end = time.time() + seconds
    while time.time() < end:
        if cond():
            return True
        time.sleep(0.1)
    return cond()


def test_noop_rerun_leaves_no_process_and_no_directory(anchored_index, panagram_overlay, tmp_path):
    """Five re-runs that have nothing to do leave neither a head-start process nor a control directory."""
    index = _copy(anchored_index, tmp_path / "index")
    tmp = tmp_path / "tmp"
    env = _env(panagram_overlay, tmp)
    for _ in range(5):
        r = _run(_snakemake(), index, env)
        assert "Nothing to be done" in r.stdout + r.stderr
    assert _wait_for(lambda: not _head_starts(tmp), 20), _head_starts(tmp)
    assert list(tmp.glob("panagram-umaps-*")) == []


@pytest.mark.skipif(os.name == "nt", reason="the test kills the process with SIGKILL")
def test_killed_head_start_fails_the_build_with_a_clear_error(anchored_index, panagram_overlay, tmp_path):
    """SIGKILL of the head-start after `go` was written: rule umaps fails with a message, within a bound."""
    index = _copy(anchored_index, tmp_path / "index")
    adir = index / "anchor" / GENOME
    for p in list(adir.glob("bitmap.*")) + [adir / "chrs.tsv"]:
        p.unlink()
    tmp = tmp_path / "tmp"
    env = _env(panagram_overlay, tmp)
    proc = subprocess.Popen(_snakemake(), cwd=index, env=env, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, text=True)
    try:
        def go_written():
            return any((d / "go").exists() and (d / "pid").exists() for d in tmp.glob("panagram-umaps-*"))
        assert _wait_for(go_written, 60), "head-start was not started"
        control = next(tmp.glob("panagram-umaps-*"))
        pid = int((control / "pid").read_text())
        assert not (control / "status").exists()
        os.kill(pid, signal.SIGKILL)
        out, _ = proc.communicate(timeout=120)
    except subprocess.TimeoutExpired:
        proc.kill()
        proc.communicate()
        pytest.fail("snakemake hung after the head-start process was killed")
    assert proc.returncode != 0
    assert "head-start process" in out and "died" in out
    assert _wait_for(lambda: not list(tmp.glob("panagram-umaps-*")), 20)


def test_attach_reports_a_dead_process_and_removes_the_directory(tmp_path, monkeypatch):
    """attach() with a pid that is gone and no status file raises and removes the control directory."""
    from panagram import _umaps_job

    dead = subprocess.Popen([sys.executable, "-c", "pass"])
    dead.wait()
    control = tmp_path / "panagram-umaps-x"
    control.mkdir()
    (control / "pid").write_text(str(dead.pid))
    (control / "log.txt").write_text("")
    with pytest.raises(SystemExit) as e:
        _umaps_job.attach(str(control))
    assert "died before" in str(e.value)
    assert not control.exists()


def test_attach_without_a_head_start_returns_false(tmp_path):
    """No pid file: nothing to attach to, the caller runs the job itself."""
    from panagram import _umaps_job

    assert _umaps_job.attach(str(tmp_path)) is False

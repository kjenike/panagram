"""Tests of panagram/introgressions/introgression_runner.py (config parsing, pipeline orchestration,
the `main` entry point and the script's command line).

subprocess.run is replaced by a recorder where only the commands matter; one test runs the real
postprocessing step through the runner on the anchored synthetic index.
"""

import os
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest
import yaml

from helpers_introgression_post import INTROGRESSIONS_DIR, fake_tools, write_bed  # noqa: F401
from panagram.introgressions import introgression_runner as runner

RUNNER_DIR = Path(runner.__file__).resolve().parent


def make_config(tmp_path, **sections):
	"""Write a config. Sections default to everything off; pass dicts to override keys of a section."""
	cfg = {
		"general": {
			"output_dir": str(tmp_path / "out"),
			"index_dir": str(tmp_path / "idx"),
			"tsv": str(tmp_path / "groups.tsv"),
			"bin": 125000,
			"ref": "Reference",
			"threads": 2,
		},
		"calling": {
			"run": False, "rmf": True, "rmu": None, "ogrp": None, "thr": [0.5], "stp": 100, "grp": "OFFSPRING",
			"anc": None, "chr": None, "cmp": ["REF"], "urf": False, "gnm": None, "trm": None, "sft": None,
			"ssz": None, "edg": False, "vis": False,
		},  # fmt: skip
		"postprocessing": {"run": False, "act": ["fgap", "rmbn"], "min": 8, "gap": 4, "map": None, "paf": None},
		"scoring": {
			"run": False, "gdt": None, "act": None, "min": None, "gap": None, "thr": None, "cmp": None, "vis": False,
		},  # fmt: skip
	}
	for name, over in sections.items():
		cfg[name].update(over)
	path = tmp_path / "config.yaml"
	path.write_text(yaml.safe_dump(cfg))
	return path


# --- parse_config ------------------------------------------------------------------------------------


def test_parse_config_with_everything_off_returns_no_flags_and_the_general_settings(tmp_path):
	call, post, score, out, thr, cmp_, threads = runner.parse_config(make_config(tmp_path))
	assert (call, post, score) == (None, None, None)
	assert out == (tmp_path / "out").resolve()
	assert thr == [0.5]
	assert cmp_ == ["REF"]
	assert threads == 2


def test_parse_config_builds_the_calling_flags_in_order(tmp_path):
	cfg = make_config(
		tmp_path,
		calling={
			"run": True, "rmf": True, "grp": "OFFSPRING", "anc": ["A", "B"], "chr": ["chr1", "chr2"],
			"cmp": ["WT", "REF"], "stp": 100, "gnm": 0.9, "trm": 3, "sft": "mean", "ssz": 2, "urf": True,
			"rmu": ["x", "y"], "ogrp": ["WT"], "edg": True, "vis": True,
		},
	)  # fmt: skip
	call, *_ = runner.parse_config(cfg)
	out, idx, tsv = (tmp_path / n for n in ("out", "idx", "groups.tsv"))
	assert call == [
		f"--out {out.resolve()}",
		f"--idx {idx.resolve()}",
		f"--tsv {tsv.resolve()}",
		"--grp OFFSPRING",
		"--anc A B",
		"--chr chr1 chr2",
		"--cmp WT REF",
		"--bin 125000",
		"--stp 100",
		"--gnm 0.9",
		"--trm 3",
		"--sft mean",
		"--ssz 2",
		"--urf",
		"--rmu x y",
		"--ogrp WT",
		"--ref Reference",
		"--rmf",
		"--edg",
		"--vis",
		"--threads 2",
	]


def test_parse_config_omits_unset_calling_options(tmp_path):
	cfg = make_config(tmp_path, calling={"run": True, "grp": None, "rmf": False})
	call, *_ = runner.parse_config(cfg)
	out, idx, tsv = (tmp_path / n for n in ("out", "idx", "groups.tsv"))
	assert call == [
		f"--out {out.resolve()}",
		f"--idx {idx.resolve()}",
		f"--tsv {tsv.resolve()}",
		"--cmp REF",
		"--bin 125000",
		"--stp 100",
		"--ref Reference",
		"--threads 2",
	]


def test_parse_config_rmu_true_becomes_the_word_true(tmp_path):
	call, *_ = runner.parse_config(make_config(tmp_path, calling={"run": True, "rmu": True}))
	assert "--rmu true" in call


def test_parse_config_builds_the_postprocessing_flags(tmp_path):
	cfg = make_config(
		tmp_path,
		postprocessing={"run": True, "act": ["lift", "fgap"], "min": 3, "gap": 0, "map": "-x asm5", "paf": "/p/paf"},
	)
	_, post, *_ = runner.parse_config(cfg)
	assert post == [
		f"--idx {(tmp_path / 'idx').resolve()}",
		"--act lift fgap",
		"--bin 125000",
		"--min 3",
		"--gap 0",
		"--map -x asm5",
		"--paf /p/paf",
		"--ref Reference",
	]


def test_parse_config_postprocessing_drops_unset_options_but_keeps_zero(tmp_path):
	cfg = make_config(tmp_path, postprocessing={"run": True, "act": None, "min": None, "gap": 0})
	_, post, *_ = runner.parse_config(cfg)
	assert post == [f"--idx {(tmp_path / 'idx').resolve()}", "--bin 125000", "--gap 0", "--ref Reference"]


def test_parse_config_builds_the_scoring_flags(tmp_path):
	cfg = make_config(
		tmp_path,
		scoring={"run": True, "gdt": str(tmp_path / "gt"), "act": ["fgap"], "min": 1, "gap": 2, "thr": 0.25,
			"cmp": ["WT"], "vis": True},
	)  # fmt: skip
	_, _, score, *_ = runner.parse_config(cfg)
	assert score == [
		f"--idx {(tmp_path / 'idx').resolve()}",
		f"--gdt {(tmp_path / 'gt').resolve()}",
		"--thr 0.25",
		"--cmp WT",
		"--act fgap",
		"--bin 125000",
		"--min 1",
		"--gap 2",
		"--ref Reference",
		"--vis",
		f"--grp {(tmp_path / 'groups.tsv').resolve()}",
	]


def test_parse_config_scoring_without_vis_has_no_vis_or_grp_flag(tmp_path):
	cfg = make_config(tmp_path, scoring={"run": True, "gdt": str(tmp_path / "gt"), "thr": None, "cmp": None})
	_, _, score, *_ = runner.parse_config(cfg)
	assert score == [
		f"--idx {(tmp_path / 'idx').resolve()}",
		f"--gdt {(tmp_path / 'gt').resolve()}",
		"--bin 125000",
		"--ref Reference",
	]


def test_parse_config_resolves_relative_paths_against_the_working_directory(tmp_path, monkeypatch):
	monkeypatch.chdir(tmp_path)
	cfg = make_config(tmp_path, general={"output_dir": "./rel/out", "index_dir": "./rel/idx", "tsv": "g.tsv"})
	_, _, _, out, *_ = runner.parse_config(cfg)
	assert out == (tmp_path / "rel" / "out").resolve()


def test_parse_config_reads_the_repository_example_configs():
	configs = sorted((INTROGRESSIONS_DIR / "example_configs").glob("*.yaml")) + sorted(
		(INTROGRESSIONS_DIR / "example").glob("*.yaml")
	)
	if not configs:
		pytest.skip("example configs are not part of the installed package")
	for c in configs:
		call, post, score, out, thr, cmp_, threads = runner.parse_config(c)
		raw = yaml.safe_load(c.read_text())
		assert thr == raw["calling"]["thr"]
		assert threads == raw["general"]["threads"]
		assert (call is not None) == bool(raw["calling"]["run"])
		assert (post is not None) == bool(raw["postprocessing"]["run"])
		assert (score is not None) == bool(raw["scoring"]["run"])


# --- run_introgression_pipeline ---------------------------------------------------------------------


class Recorder:
	"""Stand-in for subprocess.run that records commands and returns scripted results."""

	def __init__(self, returncodes=None, stdout="ran\n"):
		self.commands = []
		self.kwargs = []
		self.returncodes = returncodes or {}
		self.stdout = stdout

	def __call__(self, cmd, **kwargs):
		self.commands.append(cmd)
		self.kwargs.append(kwargs)
		rc = 0
		for key, code in self.returncodes.items():
			if key in cmd:
				rc = code
		if kwargs.get("check") and rc != 0:
			raise subprocess.CalledProcessError(rc, cmd)
		return subprocess.CompletedProcess(cmd, rc, stdout=self.stdout, stderr="")


@pytest.fixture
def recorder(monkeypatch):
	rec = Recorder()
	monkeypatch.setattr(runner.subprocess, "run", rec)
	return rec


def _threshold_dirs(out, thrs):
	for t in thrs:
		(out / f"{out.name}_{t}").mkdir(parents=True)


def test_pipeline_runs_calling_once_with_every_threshold_as_floats(tmp_path, recorder):
	out = tmp_path / "res"
	runner.run_introgression_pipeline(["--out x", "--cmp REF"], None, None, out, [1, "0.25"], ["REF"], 1, False)
	assert recorder.commands == [f"python {RUNNER_DIR / 'call_introgressions.py'} --out x --cmp REF --thr 1.0 0.25"]
	assert recorder.kwargs[0]["shell"] is True and recorder.kwargs[0]["check"] is True


def test_pipeline_without_call_flags_skips_calling(tmp_path, recorder):
	runner.run_introgression_pipeline(None, None, None, tmp_path / "res", [0.5], ["REF"], 1, False)
	assert recorder.commands == []


def test_pipeline_sweep_for_two_way_comparison_uses_0_1_to_0_95(tmp_path, recorder, capsys):
	runner.run_introgression_pipeline(["--x"], None, None, tmp_path / "res", [0.3], ["REF"], 1, True)
	expected = [round(0.1 + 0.05 * i, 2) for i in range(18)]
	assert recorder.commands[0].endswith("--thr " + " ".join(str(t) for t in expected))
	assert "2-way" in capsys.readouterr().out


def test_pipeline_sweep_for_three_way_comparison_uses_0_to_0_68(tmp_path, recorder, capsys):
	runner.run_introgression_pipeline(["--x"], None, None, tmp_path / "res", [0.3], ["WT"], 1, True)
	expected = [round(0.04 * i, 2) for i in range(18)]
	assert recorder.commands[0].endswith("--thr " + " ".join(str(t) for t in expected))
	assert "3-way" in capsys.readouterr().out


def test_pipeline_runs_postprocessing_and_scoring_per_threshold_in_the_threshold_folder(tmp_path, recorder):
	out = tmp_path / "res"
	_threshold_dirs(out, [0.2, 0.4])
	runner.run_introgression_pipeline(
		None, ["--idx I", "--act fgap"], ["--idx I", "--gdt G"], out, [0.2, 0.4], ["REF"], 2, False
	)
	assert len(recorder.commands) == 4
	for t in (0.2, 0.4):
		d = out / f"res_{t}"
		post = f"python {RUNNER_DIR / 'postprocess_introgressions.py'} --idx I --act fgap --bed {d / 'raw'} --out {d / 'postprocessed'}"
		score = f"python {RUNNER_DIR / 'score_introgressions.py'} --idx I --gdt G --pre {d / 'postprocessed'} --out {d / 'scored'}"
		assert post in recorder.commands
		assert score in recorder.commands
		assert recorder.commands.index(post) < recorder.commands.index(score)
		assert (d / "caller.log").read_text() == "ran\n\nran\n\n"


@pytest.mark.parametrize("threads,thresholds,workers", [(8, [0.2, 0.4], 2), (1, [0.2, 0.4, 0.6], 1), (3, [0.2, 0.4, 0.6], 3)])
def test_pipeline_uses_at_most_one_worker_per_threshold(tmp_path, recorder, monkeypatch, threads, thresholds, workers):
	seen = []
	real = runner.ThreadPoolExecutor

	def spy(max_workers=None, **kw):
		seen.append(max_workers)
		return real(max_workers=max_workers, **kw)

	monkeypatch.setattr(runner, "ThreadPoolExecutor", spy)
	out = tmp_path / "res"
	_threshold_dirs(out, thresholds)
	runner.run_introgression_pipeline(None, ["--idx I"], None, out, thresholds, ["REF"], threads, False)
	assert seen == [workers]


def test_pipeline_postprocessing_alone_does_not_score(tmp_path, recorder):
	out = tmp_path / "res"
	_threshold_dirs(out, [0.5])
	runner.run_introgression_pipeline(None, ["--idx I"], None, out, [0.5], ["REF"], 1, False)
	assert len(recorder.commands) == 1 and "postprocess_introgressions.py" in recorder.commands[0]


def test_pipeline_requires_the_call_output_folder_of_each_threshold(tmp_path, recorder):
	with pytest.raises(ValueError, match="Expected call output directory"):
		runner.run_introgression_pipeline(None, ["--idx I"], None, tmp_path / "res", [0.5], ["REF"], 1, False)


def test_pipeline_exits_with_status_1_when_a_step_fails_and_names_the_log(tmp_path, monkeypatch, capsys):
	rec = Recorder(returncodes={"postprocess_introgressions.py": 4}, stdout="trace\n")
	monkeypatch.setattr(runner.subprocess, "run", rec)
	out = tmp_path / "res"
	_threshold_dirs(out, [0.5])
	with pytest.raises(SystemExit) as e:
		runner.run_introgression_pipeline(None, ["--idx I"], ["--idx I"], out, [0.5], ["REF"], 1, False)
	assert e.value.code == 1
	assert "caller.log" in capsys.readouterr().out
	assert (out / "res_0.5" / "caller.log").read_text() == "trace\n\n"
	assert not any("score_introgressions.py" in c for c in rec.commands)


def test_pipeline_visualizes_only_for_a_scored_sweep(tmp_path, recorder):
	out = tmp_path / "res"
	thrs = [round(0.1 + 0.05 * i, 2) for i in range(18)]
	_threshold_dirs(out, thrs)
	runner.run_introgression_pipeline(None, None, ["--idx I", "--vis"], out, [0.5], ["REF"], 2, True)
	vis_cmds = [c for c in recorder.commands if "visualize_introgressions.py" in c]
	assert vis_cmds == [
		f"python {RUNNER_DIR / 'visualize_introgressions.py'} -v prc prcc prca mcc shtmp --dir {out} --thresholds "
		+ " ".join(str(t) for t in thrs)
	]


def test_pipeline_does_not_visualize_without_sweep_or_without_vis(tmp_path, recorder):
	out = tmp_path / "res"
	_threshold_dirs(out, [0.5])
	runner.run_introgression_pipeline(None, None, ["--idx I", "--vis"], out, [0.5], ["REF"], 1, False)
	out2 = tmp_path / "res2"
	thrs = [round(0.1 + 0.05 * i, 2) for i in range(18)]
	_threshold_dirs(out2, thrs)
	runner.run_introgression_pipeline(None, None, ["--idx I"], out2, [0.5], ["REF"], 1, True)
	assert not any("visualize_introgressions.py" in c for c in recorder.commands)


def test_pipeline_prints_completion_message(tmp_path, recorder, capsys):
	runner.run_introgression_pipeline(None, None, None, tmp_path / "res", [0.5], ["REF"], 1, False)
	assert "Introgressions analysis complete." in capsys.readouterr().out


# --- main ----------------------------------------------------------------------------------------------


@pytest.fixture
def spy_pipeline(monkeypatch):
	calls = []
	monkeypatch.setattr(runner, "run_introgression_pipeline", lambda *a: calls.append(a))
	return calls


def test_main_copies_the_config_into_the_output_folder_and_runs_the_pipeline(tmp_path, spy_pipeline):
	cfg = make_config(tmp_path, calling={"run": True, "thr": [0.1, 0.2]})
	runner.main(cfg)
	copy = tmp_path / "out" / "intro_config.yaml"
	assert copy.read_text() == cfg.read_text()
	assert len(spy_pipeline) == 1
	call_flags, post, score, out, thr, cmp_, threads, sweep = spy_pipeline[0]
	assert call_flags[0] == f"--out {(tmp_path / 'out').resolve()}"
	assert (post, score, thr, cmp_, threads, sweep) == (None, None, [0.1, 0.2], ["REF"], 2, False)


def test_main_passes_the_sweep_flag(tmp_path, spy_pipeline):
	runner.main(make_config(tmp_path), sweep=True)
	assert spy_pipeline[0][-1] is True


def test_main_exits_when_the_config_is_missing(tmp_path, spy_pipeline, capsys):
	with pytest.raises(SystemExit) as e:
		runner.main(tmp_path / "nope.yaml")
	assert e.value.code == 1
	assert "does not exist" in capsys.readouterr().out
	assert spy_pipeline == []


def _populated_out(tmp_path):
	out = tmp_path / "out"
	out.mkdir()
	(out / "a.txt").write_text("a")
	(out / "b.txt").write_text("b")
	return out


def test_main_deletes_an_existing_output_folder_when_answered_y(tmp_path, spy_pipeline, monkeypatch):
	out = _populated_out(tmp_path)
	monkeypatch.setattr("builtins.input", lambda prompt="": "y")
	runner.main(make_config(tmp_path))
	assert sorted(p.name for p in out.iterdir()) == ["intro_config.yaml"]
	assert len(spy_pipeline) == 1


def test_main_keeps_an_existing_output_folder_when_answered_n(tmp_path, spy_pipeline, monkeypatch):
	out = _populated_out(tmp_path)
	monkeypatch.setattr("builtins.input", lambda prompt="": "n")
	runner.main(make_config(tmp_path))
	assert sorted(p.name for p in out.iterdir()) == ["a.txt", "b.txt", "intro_config.yaml"]
	assert len(spy_pipeline) == 1


def test_main_exits_on_any_other_answer_and_leaves_the_folder(tmp_path, spy_pipeline, monkeypatch):
	out = _populated_out(tmp_path)
	monkeypatch.setattr("builtins.input", lambda prompt="": "maybe")
	with pytest.raises(SystemExit) as e:
		runner.main(make_config(tmp_path))
	assert e.value.code == 1
	assert sorted(p.name for p in out.iterdir()) == ["a.txt", "b.txt"]
	assert spy_pipeline == []


def test_main_does_not_ask_when_the_output_folder_holds_at_most_one_entry(tmp_path, spy_pipeline, monkeypatch):
	out = tmp_path / "out"
	out.mkdir()
	(out / "only.txt").write_text("x")

	def no_input(prompt=""):
		raise AssertionError("asked")

	monkeypatch.setattr("builtins.input", no_input)
	runner.main(make_config(tmp_path))
	assert len(spy_pipeline) == 1


def test_main_requires_a_paf_folder_for_liftover_with_several_thresholds(tmp_path, spy_pipeline):
	cfg = make_config(
		tmp_path,
		calling={"thr": [0.1, 0.2]},
		postprocessing={"run": True, "act": ["lift", "fgap"], "paf": None},
	)
	with pytest.raises(ValueError, match="PAF folder must be provided"):
		runner.main(cfg)
	assert spy_pipeline == []


def test_main_liftover_with_one_threshold_or_a_paf_folder_is_accepted(tmp_path, spy_pipeline):
	runner.main(make_config(tmp_path, postprocessing={"run": True, "act": ["lift"], "paf": None}))
	runner.main(
		make_config(tmp_path, calling={"thr": [0.1, 0.2]}, postprocessing={"run": True, "act": ["lift"], "paf": "/p"})
	)
	assert len(spy_pipeline) == 2


def test_main_liftover_with_sweep_requires_a_paf_folder(tmp_path, spy_pipeline):
	cfg = make_config(tmp_path, postprocessing={"run": True, "act": ["lift"], "paf": None})
	with pytest.raises(ValueError, match="PAF folder"):
		runner.main(cfg, sweep=True)


# --- the script's command line ----------------------------------------------------------------------------


def _script(*args, cwd):
	return subprocess.run(
		[sys.executable, str(RUNNER_DIR / "introgression_runner.py"), *map(str, args)],
		capture_output=True, text=True, cwd=cwd,
	)  # fmt: skip


def test_script_without_arguments_prints_usage_and_exits_zero(tmp_path):
	r = _script(cwd=tmp_path)
	assert r.returncode == 0
	assert "usage: python introgression_runner.py <config.yaml> [--sweep]" in r.stdout


def test_script_with_a_missing_config_exits_with_status_1(tmp_path):
	r = _script(tmp_path / "nope.yaml", cwd=tmp_path)
	assert r.returncode == 1
	assert "does not exist" in r.stdout


def test_script_with_everything_off_creates_the_output_folder_and_copies_the_config(tmp_path):
	cfg = make_config(tmp_path)
	r = _script(cfg, cwd=tmp_path)
	assert r.returncode == 0, r.stderr[-1500:]
	assert (tmp_path / "out" / "intro_config.yaml").is_file()
	assert "Introgressions analysis complete." in r.stdout


def test_script_accepts_the_sweep_flag(tmp_path):
	"""With every step off the sweep flag changes nothing visible; the script must still accept it and finish."""
	r = _script(make_config(tmp_path), "--sweep", cwd=tmp_path)
	assert r.returncode == 0, r.stderr[-1500:]
	assert "Unknown command" not in r.stdout
	assert "Introgressions analysis complete." in r.stdout


def test_script_unknown_second_argument_prints_a_hint(tmp_path):
	cfg = make_config(tmp_path)
	r = _script(cfg, "--bogus", cwd=tmp_path)
	assert "Unknown command. Use '--sweep'" in r.stdout


# --- the real postprocessing step, driven by the runner --------------------------------------------------


def test_runner_drives_postprocessing_for_each_threshold_on_the_real_index(
	tmp_path, anchored_index, fake_tools, monkeypatch
):
	out = tmp_path / "res"
	for thr in ("0.3", "0.6"):
		raw = out / f"res_{thr}" / "raw"
		raw.mkdir(parents=True)
		write_bed(raw / "alpha_chr1_WT.bed", [("chr1", 200, 400, "x"), ("chr1", 600, 800, "x")])
	cfg = make_config(
		tmp_path,
		general={"output_dir": str(out), "index_dir": str(anchored_index), "bin": 100, "threads": 2},
		calling={"thr": [0.3, 0.6]},
		postprocessing={"run": True, "act": ["fgap"], "min": 1, "gap": 2},
	)
	monkeypatch.setattr("builtins.input", lambda prompt="": "n")
	runner.main(cfg)
	for thr in ("0.3", "0.6"):
		done = out / f"res_{thr}" / "postprocessed" / "alpha_chr1_WT.bed"
		assert done.read_text() == "chr1\t200\t799\tWT_intro\n"
		assert "Postprocessing introgressions..." in (out / f"res_{thr}" / "caller.log").read_text()

"""Helpers for the introgression post-processing, visualization, heatmap and runner tests.

Fixtures defined here are imported into the test files that use them.
"""

import os
import shutil
import stat
import sys
from pathlib import Path

import pytest

import panagram

INTROGRESSIONS_DIR = Path(panagram.__file__).resolve().parent / "introgressions"


def add_scripts_to_path():
	"""bed_to_txt.py and create_heatmap.py import their siblings by bare name (they run as scripts)."""
	d = str(INTROGRESSIONS_DIR)
	if d not in sys.path:
		sys.path.insert(0, d)


class FigureLog(list):
	"""(figure, path) pairs recorded by the write_image stub."""

	def paths(self):
		return [Path(p) for _, p in self]

	def by_name(self, name):
		found = [fig for fig, p in self if Path(p).name == name]
		assert len(found) == 1, f"{name}: {len(found)} figures written; wrote {[Path(p).name for _, p in self]}"
		return found[0]


@pytest.fixture
def figure_log(monkeypatch):
	"""Replace Figure.write_image by a stub that records the figure and writes a placeholder file.

	Rendering needs kaleido and Chrome; the tests assert on the figure objects instead.
	"""
	from plotly.basedatatypes import BaseFigure

	log = FigureLog()

	def write_image(self, file, *args, **kwargs):
		log.append((self, str(file)))
		Path(file).write_bytes(b"<svg/>")

	monkeypatch.setattr(BaseFigure, "write_image", write_image)
	return log


_SCRIPT_MINIMAP2 = """#!/bin/sh
# prints its arguments, one per line, so tests can read the command line back
if [ -n "$FAKE_MINIMAP_FAIL" ]; then
	echo "minimap2 exploded" >&2
	exit 3
fi
printf '%s\\n' "$@"
echo "fake minimap2 log" >&2
"""

_SCRIPT_PAFTOOLS = """#!/bin/sh
# records its arguments in $FAKE_PAFTOOLS_ARGS and prints the content of $FAKE_PAFTOOLS_OUT
if [ -n "$FAKE_PAFTOOLS_FAIL" ]; then
	echo "paftools exploded" >&2
	exit 2
fi
printf '%s\\n' "$@" > "$FAKE_PAFTOOLS_ARGS"
cat "$FAKE_PAFTOOLS_OUT"
"""

_SCRIPT_PYTHON = """#!/bin/sh
exec "{python}" "$@"
"""


def _write_script(path, text):
	path.write_text(text)
	path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


@pytest.fixture
def fake_tools(tmp_path, monkeypatch):
	"""Directory on PATH with stand-ins for minimap2, paftools.js and a `python` wrapper.

	minimap2 prints its arguments; paftools.js prints the file named by the
	FAKE_PAFTOOLS_OUT environment variable. Returns a namespace with `args_file`
	and `set_lift(text)` to choose the lifted BED the fake paftools.js prints.
	"""
	bindir = tmp_path / "fakebin"
	bindir.mkdir()
	_write_script(bindir / "minimap2", _SCRIPT_MINIMAP2)
	_write_script(bindir / "paftools.js", _SCRIPT_PAFTOOLS)
	_write_script(bindir / "python", _SCRIPT_PYTHON.format(python=sys.executable))
	monkeypatch.setenv("PATH", str(bindir) + os.pathsep + os.environ.get("PATH", ""))
	args_file = tmp_path / "paftools_args.txt"
	out_file = tmp_path / "paftools_out.bed"
	out_file.write_text("")
	monkeypatch.setenv("FAKE_PAFTOOLS_ARGS", str(args_file))
	monkeypatch.setenv("FAKE_PAFTOOLS_OUT", str(out_file))
	monkeypatch.delenv("FAKE_MINIMAP_FAIL", raising=False)
	monkeypatch.delenv("FAKE_PAFTOOLS_FAIL", raising=False)

	class Tools:
		pass

	t = Tools()
	t.args_file = args_file
	t.set_lift = lambda text: out_file.write_text(text)
	t.bindir = bindir
	return t


@pytest.fixture(scope="module")
def index_copy(anchored_index, tmp_path_factory):
	"""A writable copy of the anchored index (the tools write alignments/ and panagram_visuals/ into it)."""
	dest = tmp_path_factory.mktemp("introgression_index") / "index"
	shutil.copytree(anchored_index, dest)
	return dest


def write_bed(path, rows):
	"""Write BED rows (tuples) tab separated."""
	Path(path).write_text("".join("\t".join(map(str, r)) + "\n" for r in rows))
	return Path(path)

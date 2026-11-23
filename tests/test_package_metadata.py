"""panagram.__about__ and the names re-exported by panagram/__init__.py."""

import re
from pathlib import Path

import pytest

import panagram
from panagram import __about__

REPO = Path(panagram.__file__).resolve().parent.parent
needs_source_tree = pytest.mark.skipif(
	not ((REPO / "setup.py").is_file() and (REPO / "pyproject.toml").is_file()),
	reason="setup.py and pyproject.toml are not next to the installed package",
)


def test_about_fields():
	"""The about module names the package, its authors, repository and license; pins the current output."""
	assert __about__.__title__ == "panagram"
	assert __about__.__uri__ == "https://github.com/kjenike/panagram"
	assert __about__.__license__ == "MIT License"
	assert __about__.__author__ == "Katharine Jenike, Sam Kovaka, Nicole Brown"
	assert __about__.__summary__ == ("Pangenome viewer",)


def test_about_version_is_a_release_number():
	"""The version is dotted integers (major.minor.patch), the form pip orders."""
	assert re.fullmatch(r"\d+\.\d+\.\d+", __about__.__version__)


def test_copyright_names_the_authors_and_a_year():
	"""The copyright line is a four-digit year followed by the author string."""
	m = re.fullmatch(r"(\d{4}) (.+)", __about__.__copyright__)
	assert m and m.group(2) == __about__.__author__
	assert 2020 <= int(m.group(1)) <= 2100


@needs_source_tree
def test_setup_py_reads_the_version_and_first_summary_from_about():
	"""setup.py exec()s __about__.py and uses __version__ and __summary__[0], so both must be reachable that way."""
	ns = {}
	exec((REPO / "panagram" / "__about__.py").read_text(), ns)
	assert ns["__version__"] == panagram.__version__
	assert ns["__summary__"][0] == "Pangenome viewer"
	setup = (REPO / "setup.py").read_text()
	assert 'about["__version__"]' in setup and 'about["__summary__"][0]' in setup


def test_package_reexports_main_and_about_names():
	"""`import panagram` exposes main (the console-script entry point) and every about field."""
	from panagram import __main__

	assert panagram.main is __main__.main
	for name in ("__title__", "__version__", "__summary__", "__author__", "__uri__", "__license__", "__copyright__"):
		assert getattr(panagram, name) == getattr(__about__, name)


@needs_source_tree
def test_console_script_target_exists():
	"""pyproject.toml's console script panagram = panagram.__main__:main points at a callable."""
	text = (REPO / "pyproject.toml").read_text()
	assert 'panagram = "panagram.__main__:main"' in text
	from panagram.__main__ import main

	assert callable(main)

"""Harness that loads panagram.view.view() without serving and exposes its nested objects.

view() builds a Dash app, defines every helper and callback as a closure, and finishes with
app.run(). The harness replaces Dash.run with a function that records the locals of its
caller (view's frame), so tests can call the nested functions directly.
"""

import dataclasses
import sys
import types

import dash


@dataclasses.dataclass
class Params:
	"""Stand-in for panagram.__main__.View with the same fields and defaults."""

	index_dir: str
	genome: str = None
	chrom: str = None
	start: int = None
	end: int = None
	order: list = None
	ndebug: bool = True
	port: str = "8050"
	host: str = "127.0.0.1"
	url_base: str = "/"
	max_chr_bins: int = 350
	bookmarks: str = None


def load_view(index_dir, monkeypatch, params=None, **overrides):
	"""Run panagram.view.view with a recording Dash.run; return a namespace of view's locals.

	params defaults to a Params(index_dir, **overrides); pass a parsed panagram.__main__.View to
	use the command-line object instead. Attributes: app, index, params, every nested function.
	"""
	from panagram import view as view_module

	captured = {}

	def fake_run(self, *args, **kwargs):
		captured["locals"] = dict(sys._getframe(1).f_locals)
		captured["run_kwargs"] = kwargs

	monkeypatch.setattr(dash.Dash, "run", fake_run)
	if params is None:
		params = Params(str(index_dir), **overrides)
	view_module.view(params)
	ns = types.SimpleNamespace(**captured["locals"])
	ns.run_kwargs = captured["run_kwargs"]
	return ns


def set_triggered(prop_id):
	"""Set the Dash callback context so ctx.triggered_id is the id part of prop_id."""
	from dash._callback_context import context_value
	from dash._utils import AttributeDict

	context_value.set(
		AttributeDict(
			**{
				"triggered_inputs": [{"prop_id": prop_id, "value": None}] if prop_id else [],
				"inputs_list": [],
				"state_list": [],
				"outputs_list": [],
				"input_values": {},
				"state_values": {},
				"args_grouping": [],
				"outputs_grouping": [],
			}
		)
	)


def call(func, prop_id, *args):
	"""Call a callback function with ctx.triggered_id set from prop_id."""
	set_triggered(prop_id)
	return func(*args)


# ---------------------------------------------------------------- fixtures shared by the view tests
# Test modules import these names; pytest finds them in the importing module.

import os

import numpy as np
import pandas as pd
import pytest


@pytest.fixture(scope="session")
def view_ns(anchored_index):
	"""Namespace of panagram.view.view() run over the anchored synthetic index (session-scoped).

	Holds app, index, params and every nested helper and callback of view().
	"""
	mp = pytest.MonkeyPatch()
	cwd = os.getcwd()
	try:
		ns = load_view(anchored_index, mp)
	finally:
		mp.undo()
		os.chdir(cwd)
	yield ns
	ns.index.close()


@pytest.fixture
def v(view_ns, anchored_index, monkeypatch):
	"""view_ns with the working directory on the index (tabix files are reopened by relative path)."""
	monkeypatch.chdir(anchored_index)
	return view_ns


UMAP_ROWS = {
	# chrom, start, end, umap1, umap2, cluster: a hand-written table with several windows per chromosome
	"chrom": ["chr1", "chr1", "chr1", "chr2", "chr2"],
	"start": [0, 1000, 2000, 0, 1000],
	"end": [1000, 2000, 3000, 1000, 2000],
	"umap1": [-1.0, 0.5, 2.0, 1.5, -2.5],
	"umap2": [0.0, 1.0, -1.0, 2.0, 0.25],
	"cluster": [0, 1, 1, 0, 2],
}


def patch_umaps(ns, monkeypatch):
	"""Replace alpha's UMAP tables in a view namespace by UMAP_ROWS.

	The synthetic chromosomes are shorter than one 100 kb window, so the index holds a single UMAP
	row per chromosome, which the Chromosome tab cannot plot (see the report).
	"""
	genome = ns.index.genomes["alpha"]
	table = pd.DataFrame(UMAP_ROWS)
	monkeypatch.setattr(genome, "genome_umap", table)
	monkeypatch.setattr(genome, "chrom_umaps", table.set_index("chrom"))
	return table


@pytest.fixture
def umaps(v, monkeypatch):
	"""UMAP_ROWS installed as alpha's UMAP tables of the shared view namespace."""
	return patch_umaps(v, monkeypatch)


def decode_array(value):
	"""Plotly 6+ JSON stores numeric arrays as {"dtype": ..., "bdata": base64}; return a list either way."""
	import base64

	if isinstance(value, dict) and "bdata" in value:
		return np.frombuffer(base64.b64decode(value["bdata"]), dtype=value["dtype"]).tolist()
	return list(value)


def popcount_bins(bits, edges):
	"""Per bin (starts `edges`, last bin open) histogram of row popcounts: array (n_popcounts, n_bins)."""
	bits = np.asarray(bits)
	counts = bits.sum(axis=1)
	ends = list(edges[1:]) + [len(bits)]
	out = np.zeros((bits.shape[1] + 1, len(edges)), dtype=int)
	for b, (s, e) in enumerate(zip(edges, ends)):
		for c in counts[s:e]:
			out[c, b] += 1
	return out


def components(node):
	"""Yield every Dash component below node (depth first), following list and component children."""
	if isinstance(node, (list, tuple)):
		for n in node:
			yield from components(n)
	elif hasattr(node, "_type"):
		yield node
		yield from components(getattr(node, "children", None))

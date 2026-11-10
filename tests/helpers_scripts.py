"""Helpers for the tests of scripts/*.py, scripts/*.sh and cpp/anchor.cpp.

Scripts are loaded by path (they are not part of the installed package). Scripts with a
`__main__` guard are imported with importlib; scripts whose body runs at import time are
executed with `exec_script`, which sets `sys.argv` and returns the namespace.
"""

import importlib.util
import subprocess
import sys
import types
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = REPO_ROOT / "scripts"
CPP = REPO_ROOT / "cpp"

THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS")
_MODULES = {}


def load_script(name):
	"""Import scripts/<name>.py as a module (for scripts guarded by `if __name__ == "__main__"`)."""
	if name not in _MODULES:
		spec = importlib.util.spec_from_file_location(f"panagram_script_{name}", SCRIPTS / f"{name}.py")
		mod = importlib.util.module_from_spec(spec)
		spec.loader.exec_module(mod)
		_MODULES[name] = mod
	return _MODULES[name]


def exec_script(name, argv, monkeypatch, tolerate=(), edits=()):
	"""Run the body of scripts/<name>.py with sys.argv = [name, *argv]; return (namespace, exception).

	Exceptions in `tolerate` raised by the module body are returned instead of propagated, so
	the functions defined before the failure stay usable. `edits` is a list of (old, new) text
	replacements applied to the source first (used only to shrink a size constant to the
	synthetic genomes); each `old` must occur in the script.
	"""
	path = SCRIPTS / f"{name}.py"
	source = path.read_text()
	for old, new in edits:
		assert old in source, old
		source = source.replace(old, new)
	monkeypatch.setattr(sys, "argv", [str(path), *map(str, argv)])
	ns = {"__name__": f"panagram_script_{name}", "__file__": str(path)}
	exc = None
	try:
		exec(compile(source, str(path), "exec"), ns)
	except tolerate as e:
		exc = e
	return ns, exc


def run_script(env, name, *args, cwd=None, timeout=600, check=True):
	"""Run `python scripts/<name> args...` in a subprocess with the overlay environment `env`."""
	cmd = [sys.executable, str(SCRIPTS / name), *map(str, args)]
	env = dict(env, **{v: "1" for v in THREAD_VARS})
	r = subprocess.run(cmd, cwd=cwd, env=env, capture_output=True, text=True, timeout=timeout)
	if check and r.returncode != 0:
		raise RuntimeError(
			f"{' '.join(cmd)} failed ({r.returncode})\nstdout:\n{r.stdout[-3000:]}\nstderr:\n{r.stderr[-3000:]}"
		)
	return r


def install_mycolorpy_shim(monkeypatch):
	"""Provide `mycolorpy.colorlist.gen_color` built on matplotlib.colormaps.

	scripts/query.py and scripts/plot_umaps.py import mycolorpy, which is not a declared
	dependency, and its gen_color calls matplotlib's removed `plt.cm.get_cmap`.
	"""
	import matplotlib
	from matplotlib.colors import to_hex

	def gen_color(cmap, n):
		return [to_hex(matplotlib.colormaps[cmap](x)) for x in np.linspace(0, 1, n)]

	pkg = types.ModuleType("mycolorpy")
	sub = types.ModuleType("mycolorpy.colorlist")
	sub.gen_color = gen_color
	pkg.colorlist = sub
	monkeypatch.setitem(sys.modules, "mycolorpy", pkg)
	monkeypatch.setitem(sys.modules, "mycolorpy.colorlist", sub)


def oracle_by_chrom(genomes, genome):
	"""chrom -> (k-mer starts x genomes) 0/1 matrix, the brute-force canonical k-mer oracle."""
	return {chrom: genomes.oracle_bits(genome, chrom) for chrom in genomes.data.seqs[genome]}


def gene_rows(genomes, genome):
	"""[(chrom, gff_start_1based, gff_end, name)] of the synthetic GFF3 of `genome`."""
	return [(c, s0 + 1, e0, n) for c, s0, e0, n in genomes.data.genes[genome]]

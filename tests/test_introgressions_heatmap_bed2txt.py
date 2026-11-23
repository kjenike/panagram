"""Tests of panagram/introgressions/create_heatmap.py and bed_to_txt.py on the anchored synthetic index.

The k-mer sharing matrix comes from the synthetic ground truth (Genomes.oracle_bits), reduced to
bins by a short reference implementation in this file. Rendering is stubbed (create_heatmap hands
its matrix to call_introgressions.visualize, which is captured).
"""

import os
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

from helpers_introgression_post import (  # noqa: F401
	INTROGRESSIONS_DIR,
	add_scripts_to_path,
	index_copy,
	write_bed,
)

add_scripts_to_path()
import bed_to_txt  # noqa: E402
import create_heatmap  # noqa: E402


# --- reference binning ----------------------------------------------------------------------------


def reference_binned(bits, bin_size, omit_fixed):
	"""bits: positions x genomes 0/1. Per bin: genomes' count of 1s (rows that are all 1 dropped if
	omit_fixed), a bin with no rows left reads 1 everywhere, then each bin is divided by its largest count."""
	n_bins = (len(bits) + bin_size - 1) // bin_size
	cols = []
	for b in range(n_bins):
		block = bits[b * bin_size : (b + 1) * bin_size]
		if omit_fixed:
			block = block[~(block == 1).all(axis=1)]
		if len(block) == 0:
			cols.append(np.ones(bits.shape[1]))
		else:
			cols.append(block.sum(axis=0).astype(float))
	m = np.array(cols).T
	return m / m.max(axis=0)


@pytest.fixture
def captured(monkeypatch):
	"""Replace call_introgressions.visualize inside create_heatmap by a recorder."""
	calls = []

	def visualize(binned_bitmap, output_file, **kwargs):
		calls.append((binned_bitmap, str(output_file), kwargs))

	monkeypatch.setattr(create_heatmap.call, "visualize", visualize)
	return calls


# --- create_heatmap.panagram_heatmap_general -------------------------------------------------------


def test_heatmap_hands_one_binned_matrix_per_chromosome_to_the_plot(index_copy, synthetic_genomes, captured):
	g = synthetic_genomes
	create_heatmap.panagram_heatmap_general(index_copy, "alpha", bitmap_step=1, bin_size=500)
	assert [os.path.basename(c[1]) for c in captured] == ["alpha_chr1_heatmap.svg", "alpha_chr2_heatmap.svg"]
	assert os.path.dirname(captured[0][1]) == str(index_copy / "panagram_visuals")
	assert (index_copy / "panagram_visuals").is_dir()
	for (binned, _, kwargs), chrom in zip(captured, ["chr1", "chr2"]):
		expected = reference_binned(g.oracle_bits("alpha", chrom), 500, True)
		assert list(binned.index) == g.names
		assert list(binned.columns) == [500 * i for i in range(expected.shape[1])]
		np.testing.assert_allclose(binned.values, expected)
		assert kwargs["inverse"] is True
		assert kwargs["groups"] is None
		assert kwargs["title"] == f"alpha - {chrom.capitalize()}"


def test_heatmap_anchor_row_is_the_anchor_genome_column_of_the_matrix(index_copy, synthetic_genomes, captured):
	# every k-mer of the anchor is in the anchor, so its row is 1 only if the anchor is the max in the bin
	create_heatmap.panagram_heatmap_general(index_copy, "gamma", bitmap_step=1, bin_size=500, omit_fixed_kmers=False)
	binned = captured[0][0]
	bits = synthetic_genomes.oracle_bits("gamma", "chr1")
	np.testing.assert_allclose(binned.loc["gamma"].values, reference_binned(bits, 500, False)[2])
	assert (binned.loc["gamma"] == 1.0).all()


def test_heatmap_omit_fixed_kmers_changes_the_matrix(index_copy, synthetic_genomes, captured):
	g = synthetic_genomes
	create_heatmap.panagram_heatmap_general(index_copy, "alpha", bitmap_step=1, bin_size=500, omit_fixed_kmers=False)
	with_fixed = captured[0][0]
	create_heatmap.panagram_heatmap_general(index_copy, "alpha", bitmap_step=1, bin_size=500, omit_fixed_kmers=True)
	without_fixed = captured[2][0]
	bits = g.oracle_bits("alpha", "chr1")
	np.testing.assert_allclose(with_fixed.values, reference_binned(bits, 500, False))
	np.testing.assert_allclose(without_fixed.values, reference_binned(bits, 500, True))
	assert not np.allclose(with_fixed.values, without_fixed.values)


def test_heatmap_bitmap_step_samples_every_step_th_kmer(index_copy, synthetic_genomes, captured):
	create_heatmap.panagram_heatmap_general(index_copy, "alpha", bitmap_step=100, bin_size=500)
	bits = synthetic_genomes.oracle_bits("alpha", "chr1")[::100]
	binned = captured[0][0]
	# sampled positions 0, 100, ... fall five to a 500 bp bin
	assert list(binned.columns) == [500 * i for i in range(binned.shape[1])]
	assert binned.shape[1] == 6
	assert binned.shape[0] == len(synthetic_genomes.names)
	# the first bin holds positions 0, 100, ..., 400: count the 1s per genome by hand from the oracle
	block = bits[:5]
	block = block[~(block == 1).all(axis=1)]
	counts = block.sum(axis=0).astype(float)
	assert counts.max() > 0
	np.testing.assert_allclose(binned[0].values, counts / counts.max())


def test_heatmap_groups_file_is_read_as_a_table_indexed_by_name(index_copy, tmp_path, captured):
	grp = tmp_path / "groups.tsv"
	grp.write_text("name\tgroup\nbeta\tX\nalpha\tX\ngamma\tY\n")
	create_heatmap.panagram_heatmap_general(index_copy, "alpha", groups=grp, bitmap_step=1, bin_size=500)
	groups = captured[0][2]["groups"]
	assert isinstance(groups, pd.DataFrame)
	assert groups.index.tolist() == ["beta", "alpha", "gamma"]
	assert groups["group"].tolist() == ["X", "X", "Y"]


def _trimmed_mean(row, trim_std):
	mean, std = row.mean(), row.std()
	if trim_std == -1:
		return mean
	return row[(row >= mean - trim_std * std) & (row <= mean + trim_std * std)].mean()


def _normalize(binned_list, target):
	"""Reference of genome-similarity normalization: shift bins <= .98 of each genome so its mean reaches target."""
	rounded = [b.round(2) for b in binned_list]
	sims = pd.concat(binned_list, axis=1).apply(_trimmed_mean, trim_std=3, axis=1)
	if target == -1:
		target = sims[sims != 1].max()
	out = []
	for b in rounded:
		b = b.copy()
		for name in b.index:
			row = b.loc[name].copy()
			row[row <= 0.98] += target - sims[name]
			b.loc[name] = row.clip(0, 1)
		out.append(b)
	return out


@pytest.mark.parametrize("gnm", [0.9, -1])
def test_heatmap_gnm_shifts_each_genome_towards_the_target_mean(index_copy, synthetic_genomes, captured, gnm):
	g = synthetic_genomes
	# unnormalized matrices, with the reference binning
	plain = {}
	for chrom in ("chr1", "chr2"):
		m = reference_binned(g.oracle_bits("alpha", chrom), 500, True)
		plain[chrom] = pd.DataFrame(m, index=g.names, columns=[500 * i for i in range(m.shape[1])])
	expected = _normalize([plain["chr1"], plain["chr2"]], gnm)
	create_heatmap.panagram_heatmap_general(index_copy, "alpha", bitmap_step=1, bin_size=500, gnm=gnm)
	for (binned, _, _), want in zip(captured, expected):
		np.testing.assert_allclose(binned.values, want.values, atol=1e-9)


def test_heatmap_command_line_passes_every_option(index_copy, captured, monkeypatch, tmp_path):
	grp = tmp_path / "g.tsv"
	grp.write_text("name\tgroup\nalpha\tX\n")
	monkeypatch.setattr(
		sys,
		"argv",
		["create_heatmap", "--index-dir", str(index_copy), "--anchor", "beta", "--groups", str(grp),
		 "--bitmap-step", "1", "--bin-size", "1000"],
	)  # fmt: skip
	create_heatmap.main()
	assert [os.path.basename(c[1]) for c in captured] == ["beta_chr1_heatmap.svg"]
	assert captured[0][0].shape == (5, 3)  # 2886 - 14 = 2872 positions -> 3 bins of 1000
	assert captured[0][2]["groups"].index.tolist() == ["alpha"]


def test_heatmap_command_line_requires_index_dir_and_anchor(monkeypatch):
	for argv in (["h", "--anchor", "alpha"], ["h", "--index-dir", "x"]):
		monkeypatch.setattr(sys, "argv", argv)
		with pytest.raises(SystemExit) as e:
			create_heatmap.main()
		assert e.value.code == 2


def test_heatmap_command_line_gnm_option_is_a_float(index_copy, captured, monkeypatch):
	seen = {}
	real = create_heatmap.panagram_heatmap_general

	def spy(**kwargs):
		seen.update(kwargs)
		return real(**kwargs)

	monkeypatch.setattr(create_heatmap, "panagram_heatmap_general", spy)
	monkeypatch.setattr(
		sys, "argv", ["h", "--index-dir", str(index_copy), "--anchor", "beta", "--gnm", "-1", "--bitmap-step", "100"]
	)
	create_heatmap.main()
	assert seen["gnm"] == -1.0
	assert seen["bitmap_step"] == 100
	assert seen["bin_size"] == 1000000
	assert seen["omit_fixed_kmers"] is True


# --- bed_to_txt.bed_to_text ---------------------------------------------------------------------------


def _run_bed2txt(monkeypatch, bed, idx, *extra, ref="alpha", wt="beta", group="WT"):
	argv = ["bed_to_txt", "--gt_bed_file", str(bed), "--index_dir", str(idx), "--ref", ref,
		"--wild_type", wt, "--wild_type_group", group, *extra]  # fmt: skip
	monkeypatch.setattr(sys, "argv", argv)
	bed_to_txt.bed_to_text()


def test_bed_to_text_writes_one_table_per_chromosome_with_a_row_per_offspring(index_copy, tmp_path, monkeypatch):
	bed = write_bed(tmp_path / "gt.bed", [("chr1", 200, 400, "x"), ("chr1", 600, 800, "x"), ("chr2", 1500, 1700, "x")])
	_run_bed2txt(monkeypatch, bed, index_copy, "--bin_size", "100")
	assert sorted(p.name for p in tmp_path.glob("*.txt")) == ["chr1_WT.txt", "chr2_WT.txt"]

	t1 = pd.read_csv(tmp_path / "chr1_WT.txt", sep="\t", index_col="Sample")
	# offspring are all genomes except the reference (alpha) and the wild type (beta)
	assert t1.index.tolist() == ["gamma", "delta", "epsilon"]
	assert t1.columns.tolist() == [str(100 * i) for i in range(27)]  # alpha chr1: 2636 positions
	want = np.zeros(27, dtype=int)
	want[[2, 3, 6, 7]] = 1
	for sample in t1.index:
		assert t1.loc[sample].tolist() == want.tolist()

	t2 = pd.read_csv(tmp_path / "chr2_WT.txt", sep="\t", index_col="Sample")
	assert t2.shape == (3, 17)  # alpha chr2: 1686 positions
	assert t2.iloc[0].tolist() == [0] * 15 + [1, 1]


def test_bed_to_text_names_the_file_after_the_wild_type_group_and_writes_next_to_the_bed(
	index_copy, tmp_path, monkeypatch
):
	d = tmp_path / "sub"
	d.mkdir()
	bed = write_bed(d / "gt.bed", [("chr1", 200, 400, "x")])
	_run_bed2txt(monkeypatch, bed, index_copy, "--bin_size", "100", group="GROUP_Z", wt="gamma")
	assert [p.name for p in d.glob("*.txt")] == ["chr1_GROUP_Z.txt"]
	t = pd.read_csv(d / "chr1_GROUP_Z.txt", sep="\t", index_col="Sample")
	assert t.index.tolist() == ["beta", "delta", "epsilon"]  # alpha is the reference, gamma the wild type


def test_bed_to_text_default_bin_size_is_one_million(index_copy, tmp_path, monkeypatch):
	bed = write_bed(tmp_path / "gt.bed", [("chr1", 0, 2000, "x")])
	_run_bed2txt(monkeypatch, bed, index_copy)
	t = pd.read_csv(tmp_path / "chr1_WT.txt", sep="\t", index_col="Sample")
	assert t.columns.tolist() == ["0"]


def test_bed_to_text_requires_all_options(index_copy, tmp_path, monkeypatch):
	bed = write_bed(tmp_path / "gt.bed", [("chr1", 0, 10, "x")])
	full = ["--gt_bed_file", str(bed), "--index_dir", str(index_copy), "--ref", "alpha", "--wild_type", "beta",
		"--wild_type_group", "WT"]  # fmt: skip
	for i in (0, 2, 4, 6, 8):
		argv = ["b"] + full[:i] + full[i + 2 :]
		monkeypatch.setattr(sys, "argv", argv)
		with pytest.raises(SystemExit) as e:
			bed_to_txt.bed_to_text()
		assert e.value.code == 2


def test_bed_to_txt_runs_as_a_script(index_copy, tmp_path):
	bed = write_bed(tmp_path / "gt.bed", [("chr1", 200, 400, "x")])
	r = subprocess.run(
		[sys.executable, str(INTROGRESSIONS_DIR / "bed_to_txt.py"), "--gt_bed_file", str(bed), "--index_dir",
		 str(index_copy), "--ref", "alpha", "--wild_type", "beta", "--wild_type_group", "WT", "--bin_size", "100"],
		capture_output=True, text=True, cwd=tmp_path,
	)  # fmt: skip
	assert r.returncode == 0, r.stderr[-2000:]
	t = pd.read_csv(tmp_path / "chr1_WT.txt", sep="\t", index_col="Sample")
	assert t.iloc[0, 2] == 1 and t.iloc[0, 3] == 1 and t.iloc[0].sum() == 2

"""panagram.introgressions.simulate_introgressions: FASTA/BED helpers, mutation model, introgressions, CLI.

Expected values are hand-worked tiny cases, invariants of the mutation model (the position
mapper must send every surviving base to the same base in the mutated sequence), scipy for the
size distribution, and, for the CLI, the planted introgressions themselves: with a
substitution-only wild relative the offspring must equal the reference outside the BED
intervals and the wild relative inside them.
"""

import gzip
import sys

import numpy as np
import pytest
from scipy.stats import beta

import helpers_introgressions as H
from panagram.introgressions import simulate_introgressions as sim

# ------------------------------------------------------------------ parse_fasta / write_fasta / write_bed


def test_parse_fasta_multiline_blank_lines_and_header_first_token(tmp_path):
	p = tmp_path / "a.fa"
	p.write_text(">chr1 some description\nacgt\nACG\n\n>chr2\nTTTT\n>empty\n>chr3\nGG\n")
	got = sim.parse_fasta(p)
	assert list(got) == ["chr1", "chr2", "empty", "chr3"]
	assert got["chr1"] == "ACGTACG"  # upper-cased, lines joined
	assert got["chr2"] == "TTTT"
	assert got["empty"] == ""
	assert got["chr3"] == "GG"


def test_parse_fasta_reads_gzip_by_suffix(tmp_path):
	p = tmp_path / "a.fasta.gz"
	with gzip.open(p, "wt") as f:
		f.write(">x\nAC\nGT\n")
	assert dict(sim.parse_fasta(p)) == {"x": "ACGT"}
	assert dict(sim.parse_fasta(str(p))) == {"x": "ACGT"}


def test_parse_fasta_empty_file_gives_no_records(tmp_path):
	p = tmp_path / "e.fa"
	p.write_text("")
	assert len(sim.parse_fasta(p)) == 0


def test_write_fasta_wraps_lines(tmp_path):
	p = tmp_path / "o.fa"
	sim.write_fasta({"a": "ACGTACGTAC", "b": "GG", "c": ""}, p, wrap=4)
	assert p.read_text() == ">a\nACGT\nACGT\nAC\n>b\nGG\n>c\n"
	sim.write_fasta({"a": "A" * 130}, p)
	lines = p.read_text().split("\n")
	assert [len(x) for x in lines[1:4]] == [60, 60, 10]


def test_fasta_round_trip(tmp_path):
	rng = np.random.default_rng(3)
	seqs = {f"c{i}": "".join(rng.choice(list("ACGT"), size=n)) for i, n in enumerate([1, 59, 60, 61, 500])}
	sim.write_fasta(seqs, tmp_path / "r.fa")
	assert dict(sim.parse_fasta(tmp_path / "r.fa")) == seqs


def test_write_bed_one_line_per_entry(tmp_path):
	sim.write_bed(["c\t1\t5\tx", "c\t9\t12\ty"], tmp_path / "b.bed")
	assert (tmp_path / "b.bed").read_text() == "c\t1\t5\tx\nc\t9\t12\ty\n"
	sim.write_bed([], tmp_path / "e.bed")
	assert (tmp_path / "e.bed").read_text() == ""


# ------------------------------------------------------------------ position_spacer


def test_position_spacer_keeps_non_overlapping_positions():
	assert sim.position_spacer([10, 50], [5, 5], 100) == [(10, 5), (50, 5)]


def test_position_spacer_pushes_overlapping_position_past_previous_end():
	"""(10, len 5) ends at 15; a second start at 12 moves to 16."""
	assert sim.position_spacer([10, 12], [5, 5], 100) == [(10, 5), (16, 5)]


def test_position_spacer_drops_what_does_not_fit_and_sorts_positions():
	out = sim.position_spacer([50, 10], [5, 5], 52)
	assert out == [(10, 5)]  # sorted to [10, 50]; 50 + 5 > 52
	assert sim.position_spacer([10, 12], [5, 90], 100) == [(10, 5)]  # moved to 16 + 90 > 100


# ------------------------------------------------------------------ generate_skewed_sizes


def test_generate_skewed_sizes_are_ints_in_range():
	sizes = sim.generate_skewed_sizes(2000, 10, 50, np.random.default_rng(0))
	assert len(sizes) == 2000
	assert all(isinstance(x, int) for x in sizes)
	assert min(sizes) >= 10 and max(sizes) <= 50


def test_generate_skewed_sizes_mean_follows_the_beta_distribution():
	"""Beta(0.05, 1) has mean 0.05/1.05, so sizes average size_min + 0.0476 * (max - min); skewed to small sizes."""
	sizes = np.array(sim.generate_skewed_sizes(20000, 1, 1001, np.random.default_rng(1)))
	assert sizes.mean() == pytest.approx(1 + 1000 * 0.05 / 1.05, abs=8)
	assert np.median(sizes) < 50  # median of Beta(0.05, 1) is 0.5**20 ~ 0
	# other shape parameters: Beta(2, 2) is symmetric
	sym = np.array(sim.generate_skewed_sizes(20000, 0, 1000, np.random.default_rng(1), a=2, b=2))
	assert sym.mean() == pytest.approx(500, abs=10)


def test_generate_skewed_sizes_constant_when_min_equals_max():
	assert set(sim.generate_skewed_sizes(50, 7, 7, np.random.default_rng(0))) == {7}


def test_generate_skewed_sizes_reproducible_from_the_seed():
	a = sim.generate_skewed_sizes(30, 1, 100, np.random.default_rng(5))
	b = sim.generate_skewed_sizes(30, 1, 100, np.random.default_rng(5))
	assert a == b
	assert a == [int(round(1 + v * 99)) for v in beta.rvs(0.05, 1, size=30, random_state=np.random.default_rng(5))]


# ------------------------------------------------------------------ create_random_trenches / edge weights


def test_create_random_trenches_no_drop_normalizes_only():
	w = np.array([1.0, 1, 2, 4])
	out = sim.create_random_trenches(w.copy(), np.random.default_rng(0), block_size=2, drop_fraction=0.0)
	np.testing.assert_allclose(out, [0.125, 0.125, 0.25, 0.5])


def test_create_random_trenches_full_drop_scales_whole_blocks_only():
	"""25 weights of 1, block 10, every block dropped: the two full blocks are x0.8, the 5-weight tail is not."""
	out = sim.create_random_trenches(np.ones(25), np.random.default_rng(0), block_size=10, drop_fraction=1.0)
	total = 20 * 0.8 + 5
	np.testing.assert_allclose(out[:20], 0.8 / total)
	np.testing.assert_allclose(out[20:], 1 / total)
	assert out.sum() == pytest.approx(1)


def test_create_random_trenches_without_rng_still_normalizes():
	out = sim.create_random_trenches(np.array([1.0, 3.0]), None, block_size=1, drop_fraction=0.0)
	np.testing.assert_allclose(out, [0.25, 0.75])


def test_generate_edge_tapered_weights_shape():
	"""Sum 1; flat in the middle; highest at both ends where the weight is exactly twice the middle weight."""
	n = 1000
	w = sim.generate_edge_tapered_weights(n, np.random.default_rng(0))
	assert w.sum() == pytest.approx(1)
	assert (w > 0).all()
	assert w[0] == pytest.approx(2 * w[n // 2], rel=1e-9)
	assert w[-1] == pytest.approx(w[0], rel=1e-6)
	np.testing.assert_allclose(w, w[::-1], rtol=1e-6)
	assert w[400:600].max() == pytest.approx(w[400:600].min())
	assert w.max() == w[0]
	assert (w[1 : n // 2] <= w[:-1][: n // 2 - 1] + 1e-15).all()  # non-increasing towards the middle


def test_generate_edge_tapered_weights_extent_of_the_edge_taper():
	"""n = 1001 so x = i / 1000. With edge_fraction 0.3 the taper is 1 - (d / 0.15)**5 for distance d from the nearest
	end (floored at the 0.5 of the flat middle), so relative to the middle: d = 0.1 -> (1 - (2/3)**5) / 0.5, d = 0.14 and
	d = 0.2 -> 1. Both ends match."""
	w = sim.generate_edge_tapered_weights(1001, np.random.default_rng(0))
	mid = w[500]
	assert w[100] / mid == pytest.approx((1 - (2 / 3) ** 5) / 0.5)
	assert w[900] / mid == pytest.approx(w[100] / mid)
	assert w[140] / mid == pytest.approx(1.0)
	assert w[200] / mid == pytest.approx(1.0)
	assert w[50] / mid == pytest.approx((1 - (1 / 3) ** 5) / 0.5)


# ------------------------------------------------------------------ generate_mutation_positions


def test_generate_mutation_positions_deletions_do_not_overlap_and_are_zeroed():
	rng = np.random.default_rng(0)
	w = np.full(1000, 1 / 1000)
	muts, new_w = sim.generate_mutation_positions("deletion", 10, 3, 8, w.copy(), rng)
	assert len(muts) == 10
	spans = sorted((p, p + n) for p, n, t in muts)
	assert all(t == "deletion" for _, _, t in muts)
	assert all(3 <= n <= 8 for _, n, _ in muts)
	assert all(a[1] <= b[0] for a, b in zip(spans, spans[1:]))
	assert spans[-1][1] <= 1000
	zero = np.zeros(1000, dtype=bool)
	for a, b in spans:
		zero[a:b] = True
	assert (new_w[zero] == 0).all()
	assert (new_w[~zero] > 0).all()
	assert new_w.sum() == pytest.approx(1)


@pytest.mark.parametrize("kind", ["insertion", "snp"])
def test_generate_mutation_positions_point_mutations_use_distinct_positions(kind):
	rng = np.random.default_rng(1)
	w = np.full(200, 1 / 200)
	muts, new_w = sim.generate_mutation_positions(kind, 30, 2, 2 if kind == "insertion" else 9, w.copy(), rng)
	positions = [int(p) for p, _, _ in muts]
	assert len(set(positions)) == 30
	assert all(t == kind for _, _, t in muts)
	assert (new_w[positions] == 0).all()
	assert np.count_nonzero(new_w) == 170
	assert new_w.sum() == pytest.approx(1)
	if kind == "insertion":
		assert {n for _, n, _ in muts} == {2}


def test_generate_mutation_positions_caps_at_available_positions():
	"""Only three positions have weight, ten are requested: three are returned, all on the weighted positions."""
	w = np.zeros(50)
	w[[4, 20, 33]] = 1 / 3
	muts, new_w = sim.generate_mutation_positions("snp", 10, 1, 2, w, np.random.default_rng(0))
	assert sorted(int(p) for p, _, _ in muts) == [4, 20, 33]
	assert new_w.sum() == 0  # nothing left; not renormalized


# ------------------------------------------------------------------ mutate_seq_with_indels_and_snps


def _random_seq(n, seed=1):
	return "".join(np.random.default_rng(seed).choice(list("ACGT"), size=n))


def _mutate(seq, sub=0.0, ins=0.0, dele=0.0, ins_size=(3, 9), del_size=(2, 7), seed=0):
	return sim.mutate_seq_with_indels_and_snps(
		seq, sub, ins, dele, ins_size[0], ins_size[1], del_size[0], del_size[1], np.random.default_rng(seed)
	)


def test_mutate_with_zero_rates_is_the_identity():
	seq = _random_seq(1000)
	new, mapper, avail = _mutate(seq)
	assert new == seq
	assert mapper == list(range(1000))
	assert avail == list(range(1000))


def test_mutate_substitutions_only(capsys):
	seq = _random_seq(4000)
	new, mapper, avail = _mutate(seq, sub=0.01)
	assert len(new) == 4000
	assert mapper == list(range(4000))
	changed = [i for i in range(4000) if new[i].upper() != seq[i]]
	assert len(changed) == 40  # int(4000 * 0.01)
	assert len(avail) == 4000 - 40
	assert not set(changed) & set(avail)


def test_mutate_insertions_only_keep_all_original_bases_in_order():
	"""Fixed insertion size 5, rate 0.005: 20 insertions, each after its position; the mapper points at the unchanged bases."""
	seq = _random_seq(4000)
	new, mapper, avail = _mutate(seq, ins=0.005, ins_size=(5, 5))
	assert len(new) == 4000 + 20 * 5
	assert len(mapper) == 4000
	assert all(new[mapper[i]].upper() == seq[i] for i in range(4000))
	assert all(b > a for a, b in zip(mapper, mapper[1:]))
	inserted = set(range(len(new))) - set(mapper)
	assert len(inserted) == 100
	assert len(avail) == 4000 - 20


def test_mutate_deletions_only_remove_runs_marked_minus_one():
	seq = _random_seq(4000)
	new, mapper, avail = _mutate(seq, dele=0.005, del_size=(4, 4))
	deleted = [i for i, m in enumerate(mapper) if m == -1]
	assert len(deleted) % 4 == 0 and 4 * 15 <= len(deleted) <= 80
	assert len(new) == 4000 - len(deleted)
	kept = [i for i in range(4000) if mapper[i] != -1]
	assert "".join(seq[i] for i in kept) == new.upper()
	assert [mapper[i] for i in kept] == list(range(len(new)))
	assert not set(deleted) & set(avail)
	assert len(avail) == 4000 - len(deleted)


def test_mutate_all_kinds_mapper_invariant_over_seeds():
	"""Every surviving original base maps to the same base in the new sequence, except at the substituted positions."""
	seq = _random_seq(4000)
	for seed in range(4):
		new, mapper, avail = _mutate(seq, sub=0.01, ins=0.005, dele=0.005, seed=seed)
		m = np.array(mapper)
		kept = np.nonzero(m != -1)[0]
		mismatches = [i for i in kept if new[m[i]].upper() != seq[i]]
		assert len(mismatches) == 40
		assert (np.diff(m[kept]) > 0).all()
		n_deleted = int((m == -1).sum())
		n_ins_total = len(new) - (4000 - n_deleted)
		assert n_ins_total >= 20 * 3
		assert len(avail) == 4000 - n_deleted - 20 - 40
		assert all(m[p] != -1 for p in avail)


def test_mutate_substituted_base_differs_from_original_and_is_a_nucleotide():
	seq = "A" * 2000
	new, _, _ = _mutate(seq, sub=0.05)
	assert set(new.upper()) <= set("ACGT")
	assert sum(1 for c in new.upper() if c != "A") == 100


def test_mutate_is_reproducible_for_a_seed():
	seq = _random_seq(2000)
	a = _mutate(seq, 0.01, 0.005, 0.005, seed=9)
	b = _mutate(seq, 0.01, 0.005, 0.005, seed=9)
	assert a[0] == b[0] and a[1] == b[1] and a[2] == b[2]
	assert a[0] != _mutate(seq, 0.01, 0.005, 0.005, seed=10)[0]


def test_mutate_reports_lengths_on_stdout(capsys):
	_mutate(_random_seq(500), sub=0.01)
	out = capsys.readouterr().out
	assert "Original length: 500" in out
	assert "New length: 500" in out
	assert "Number of mutations: 5" in out


def test_apply_genome_wide_mutations_handles_each_chromosome():
	seqs = {"c1": _random_seq(2000, 1), "c2": _random_seq(1000, 2)}
	new, mappers, avail = sim.apply_genome_wide_mutations(
		seqs, 0.01, 0.0, 0.0, 1, 1, 1, 1, np.random.default_rng(0)
	)
	assert list(new) == ["c1", "c2"] and list(mappers) == ["c1", "c2"] and list(avail) == ["c1", "c2"]
	assert sum(a != b.upper() for a, b in zip(seqs["c1"], new["c1"])) == 20
	assert sum(a != b.upper() for a, b in zip(seqs["c2"], new["c2"])) == 10
	assert len(mappers["c2"]) == 1000 and len(avail["c2"]) == 990


# ------------------------------------------------------------------ mutate_seq_with_introgressions


REF10 = "0123456789"


def test_mutate_seq_with_introgressions_hand_worked():
	"""Start 2 length 3 with an identity mapping: ref[2:5] is replaced by rel[2:5]; BED end is inclusive (4)."""
	seq, bed = sim.mutate_seq_with_introgressions(
		"c1", REF10, "abcdefghij", list(range(10)), [2], 1, [3], np.random.default_rng(0)
	)
	assert seq == "01cde56789"
	assert bed == ["c1\t2\t4\t3 bp introgression"]


def test_mutate_seq_with_introgressions_two_regions_in_position_order():
	seq, bed = sim.mutate_seq_with_introgressions(
		"c", REF10, "abcdefghij", list(range(10)), [1, 6], 2, [2, 2], np.random.default_rng(0)
	)
	assert seq == "0bc345gh89"
	assert bed == ["c\t1\t2\t2 bp introgression", "c\t6\t7\t2 bp introgression"]


def test_mutate_seq_with_introgressions_shortens_end_that_maps_into_a_deletion():
	"""rel lacks ref base 5. Start 2, length 3: end 5 maps to -1, so the end moves back to 4 -> rel[2:4] = 'cd'."""
	mapper = [0, 1, 2, 3, 4, -1, 5, 6, 7, 8]
	seq, bed = sim.mutate_seq_with_introgressions(
		"c", REF10, "abcdefghi", mapper, [2], 1, [3], np.random.default_rng(0)
	)
	assert seq == "01cd456789"
	assert bed == ["c\t2\t3\t2 bp introgression"]


def test_mutate_seq_with_introgressions_skips_region_that_is_entirely_deleted():
	mapper = [0, 1, 2, -1, -1, -1, 3, 4, 5, 6]
	seq, bed = sim.mutate_seq_with_introgressions(
		"c", REF10, "abcdefg", mapper, [2], 1, [3], np.random.default_rng(0)
	)
	assert seq == REF10
	assert bed == []


def test_mutate_seq_with_introgressions_drops_region_that_runs_off_the_chromosome():
	seq, bed = sim.mutate_seq_with_introgressions(
		"c", REF10, "abcdefghij", list(range(10)), [8], 1, [5], np.random.default_rng(0)
	)
	assert seq == REF10 and bed == []


# ------------------------------------------------------------------ apply_genome_wide_introgressions


def test_apply_genome_wide_introgressions_one_per_chromosome():
	ref = {"c1": "A" * 1000, "c2": "C" * 1000}
	rel = {"c1": "G" * 1000, "c2": "T" * 1000}
	mappers = {c: list(range(1000)) for c in ref}
	avail = {c: list(range(100, 700)) for c in ref}
	new, bed = sim.apply_genome_wide_introgressions(
		ref, rel, mappers, avail, 2, 50, 100, np.random.default_rng(0)
	)
	assert list(new) == ["c1", "c2"]
	assert len(bed) == 2
	for entry, (chrom, base, intro) in zip(bed, [("c1", "A", "G"), ("c2", "C", "T")]):
		c, start, end, label = entry.split("\t")
		start, end = int(start), int(end)
		assert c == chrom and label == f"{end - start + 1} bp introgression"
		assert 50 <= end - start + 1 <= 100
		want = base * start + intro * (end - start + 1) + base * (1000 - end - 1)
		assert new[chrom] == want


def test_apply_genome_wide_introgressions_places_an_equal_share_on_each_chromosome():
	"""4 introgressions over 2 chromosomes: 2 on each, in chromosome order."""
	ref = {"c1": "A" * 2000, "c2": "A" * 2000}
	rel = {"c1": "G" * 2000, "c2": "G" * 2000}
	mappers = {c: list(range(2000)) for c in ref}
	avail = {c: list(range(0, 1000)) for c in ref}
	new, bed = sim.apply_genome_wide_introgressions(
		ref, rel, mappers, avail, 4, 20, 40, np.random.default_rng(0)
	)
	assert [e.split("\t")[0] for e in bed] == ["c1", "c1", "c2", "c2"]
	spans = [(int(e.split("\t")[1]), int(e.split("\t")[2])) for e in bed]
	assert spans[0][1] < spans[1][0] and spans[2][1] < spans[3][0]  # no overlap within a chromosome
	for c in ref:
		assert new[c].count("G") == sum(e - s + 1 for (s, e), b in zip(spans, bed) if b.startswith(c))


# ------------------------------------------------------------------ main()


L_CHROM = 40000


@pytest.fixture(scope="module")
def ref_fasta(tmp_path_factory):
	d = tmp_path_factory.mktemp("simref")
	rng = np.random.default_rng(11)
	seqs = {c: "".join(rng.choice(list("ACGT"), size=L_CHROM)) for c in ("c1", "c2")}
	path = d / "ref.fasta"
	with open(path, "w") as f:
		for c, s in seqs.items():
			f.write(f">{c}\n")
			for i in range(0, len(s), 70):
				f.write(s[i : i + 70] + "\n")
	return path, seqs


def run_sim(monkeypatch, ref, out, *args):
	monkeypatch.setattr(
		sys, "argv", ["simulate_introgressions", "--ref", str(ref), "--out-folder", str(out), *map(str, args)]
	)
	sim.main()
	return out


SNP_ONLY = [
	"--num-introgressions", 2,
	"--introgression-size-min", 50,
	"--introgression-size-max", 100,
	"--rel-sub-rate", 0.2, "--rel-ins-rate", 0, "--rel-del-rate", 0,
	"--mut-rate-start", 0, "--mut-sub-rate", 0.015, "--mut-ins-rate", 0, "--mut-del-rate", 0,
	"--rounds", 3,
	"--seed", 3,
]


@pytest.fixture(scope="module")
def snp_only_run(ref_fasta, tmp_path_factory):
	"""A substitution-only simulation, run once. Rates go linearly from --mut-rate-start 0: 0, 0.0075, 0.015."""
	mp = pytest.MonkeyPatch()
	try:
		out = run_sim(mp, ref_fasta[0], tmp_path_factory.mktemp("sim_snp") / "o", *SNP_ONLY)
	finally:
		mp.undo()
	return out


def _read(out, name):
	return {k: v.upper() for k, v in sim.parse_fasta(out / name).items()}


def test_main_writes_the_documented_files(snp_only_run):
	names = sorted(p.name for p in snp_only_run.iterdir())
	assert names == [
		"ref_0_introgressions.bed",
		"ref_0_offspring.fasta",
		"ref_1_introgressions.bed",
		"ref_1_offspring.fasta",
		"ref_2_introgressions.bed",
		"ref_2_offspring.fasta",
		"ref_3_introgressions.bed",
		"ref_3_offspring.fasta",
		"ref_wildrelative.fasta",
	]


def test_main_wild_relative_has_the_requested_substitution_count(snp_only_run, ref_fasta):
	wild = _read(snp_only_run, "ref_wildrelative.fasta")
	for c, s in ref_fasta[1].items():
		assert len(wild[c]) == L_CHROM
		assert sum(a != b for a, b in zip(s, wild[c])) == int(L_CHROM * 0.2)


def test_main_offspring_is_reference_outside_and_wild_relative_inside_the_planted_introgressions(
	snp_only_run, ref_fasta
):
	"""Truth comes from the BED file: inclusive intervals; sizes within the requested range."""
	wild = _read(snp_only_run, "ref_wildrelative.fasta")
	off = _read(snp_only_run, "ref_0_offspring.fasta")
	rows = H.read_bed_rows(snp_only_run / "ref_0_introgressions.bed")
	assert [r[0] for r in rows] == ["c1", "c2"]
	for chrom, start, end, label in rows:
		assert 50 <= end - start + 1 <= 100
		assert label == f"{end - start + 1} bp introgression"
		expected = ref_fasta[1][chrom][:start] + wild[chrom][start : end + 1] + ref_fasta[1][chrom][end + 1 :]
		assert off[chrom] == expected
	# the wild relative differs from the reference at ~20% of bases, so the introgressions are visible
	assert sum(a != b for c in off for a, b in zip(off[c], ref_fasta[1][c])) > 0


def test_main_later_generations_add_substitutions_at_the_linear_rates(snp_only_run):
	"""Rates 0, 0.0075, 0.015 over three rounds. As the README says, each generation is the F1 offspring (generation 0) plus fresh mutations at that round's rate."""
	g0, g1, g2, g3 = (_read(snp_only_run, f"ref_{i}_offspring.fasta") for i in range(4))
	for c in g0:
		assert g1[c] == g0[c]
		assert sum(a != b for a, b in zip(g0[c], g2[c])) == int(L_CHROM * 0.0075)
		assert sum(a != b for a, b in zip(g0[c], g3[c])) == int(L_CHROM * 0.015)


def test_main_later_generation_beds_keep_coordinates_without_indels(snp_only_run):
	g0 = H.read_bed_rows(snp_only_run / "ref_0_introgressions.bed")
	for i in (1, 2, 3):
		rows = H.read_bed_rows(snp_only_run / f"ref_{i}_introgressions.bed")
		assert [(c, s, e) for c, s, e, _ in rows] == [(c, s, e) for c, s, e, _ in g0]
		assert {n for *_, n in rows} == {"introgression"}


def test_main_same_seed_gives_identical_files_and_other_seed_differs(monkeypatch, ref_fasta, snp_only_run, tmp_path):
	again = run_sim(monkeypatch, ref_fasta[0], tmp_path / "again", *SNP_ONLY)
	for p in snp_only_run.iterdir():
		assert (again / p.name).read_bytes() == p.read_bytes(), p.name
	args = list(SNP_ONLY)
	args[args.index("--seed") + 1] = 4
	other = run_sim(monkeypatch, ref_fasta[0], tmp_path / "other", *args)
	assert (other / "ref_wildrelative.fasta").read_bytes() != (snp_only_run / "ref_wildrelative.fasta").read_bytes()


def test_main_default_seed_is_42(monkeypatch, ref_fasta, tmp_path):
	base = ["--num-introgressions", 2, "--introgression-size-min", 50, "--introgression-size-max", 100,
		"--rel-ins-rate", 0, "--rel-del-rate", 0, "--rounds", 0, "--mut-rate-start", 0]
	a = run_sim(monkeypatch, ref_fasta[0], tmp_path / "a", *base)
	b = run_sim(monkeypatch, ref_fasta[0], tmp_path / "b", *base, "--seed", 42)
	assert (a / "ref_wildrelative.fasta").read_bytes() == (b / "ref_wildrelative.fasta").read_bytes()
	assert (a / "ref_0_introgressions.bed").read_bytes() == (b / "ref_0_introgressions.bed").read_bytes()


def test_main_rounds_zero_writes_only_generation_zero(monkeypatch, ref_fasta, tmp_path):
	out = run_sim(
		monkeypatch, ref_fasta[0], tmp_path / "o",
		"--rounds", 0, "--num-introgressions", 2, "--introgression-size-min", 50,
		"--introgression-size-max", 100, "--rel-ins-rate", 0, "--rel-del-rate", 0,
	)
	assert sorted(p.name for p in out.iterdir()) == [
		"ref_0_introgressions.bed", "ref_0_offspring.fasta", "ref_wildrelative.fasta",
	]


def test_main_reads_gzipped_reference_and_names_files_after_its_stem(monkeypatch, ref_fasta, tmp_path):
	gz = tmp_path / "genome.fa.gz"
	with gzip.open(gz, "wt") as f:
		f.write(">only\n" + ref_fasta[1]["c1"][:5000] + "\n")
	out = run_sim(
		monkeypatch, gz, tmp_path / "nested" / "o",
		"--rounds", 1, "--num-introgressions", 1, "--introgression-size-min", 20,
		"--introgression-size-max", 30, "--rel-sub-rate", 0.01, "--rel-ins-rate", 0, "--rel-del-rate", 0,
		"--mut-rate-start", 0, "--mut-sub-rate", 0, "--mut-ins-rate", 0, "--mut-del-rate", 0,
	)
	assert sorted(p.name for p in out.iterdir()) == [
		"genome_0_introgressions.bed", "genome_0_offspring.fasta",
		"genome_1_introgressions.bed", "genome_1_offspring.fasta", "genome_wildrelative.fasta",
	]
	# rounds 1, all mutation rates 0 (start 0): generation 1 equals generation 0
	assert _read(out, "genome_1_offspring.fasta") == _read(out, "genome_0_offspring.fasta")
	assert list(_read(out, "genome_wildrelative.fasta")) == ["only"]


def test_main_zero_mutation_rate_start_makes_generation_one_equal_to_zero_and_insertions_start_after(
	monkeypatch, ref_fasta, tmp_path
):
	"""--mut-rate-start 0 with --mut-ins-rate 0.001 over 2 rounds: rates [0, 0.001]. Insertions of fixed size 5 appear in generation 2 only."""
	out = run_sim(
		monkeypatch, ref_fasta[0], tmp_path / "o",
		"--rounds", 2, "--num-introgressions", 2, "--introgression-size-min", 50,
		"--introgression-size-max", 100, "--rel-ins-rate", 0, "--rel-del-rate", 0,
		"--mut-rate-start", 0, "--mut-sub-rate", 0, "--mut-del-rate", 0, "--mut-ins-rate", 0.001,
		"--mut-ins-size-min", 5, "--mut-ins-size-max", 5, "--seed", 8,
	)
	g0, g1, g2 = (_read(out, f"ref_{i}_offspring.fasta") for i in range(3))
	assert g1 == g0
	b1 = H.read_bed_rows(out / "ref_1_introgressions.bed")
	b2 = H.read_bed_rows(out / "ref_2_introgressions.bed")
	for c in g0:
		assert len(g2[c]) == len(g1[c]) + 40 * 5
	for (c, s1, e1, _), (c2, s2, e2, _) in zip(b1, b2):
		assert c == c2
		assert s2 >= s1 and e2 >= e1
		# the mapped ends are the same bases as before the insertions
		assert g2[c][s2] == g1[c][s1] and g2[c][e2] == g1[c][e1]


def test_main_later_generation_beds_follow_deletions_in_the_mutated_genome(monkeypatch, ref_fasta, tmp_path):
	"""Rounds are given a hand-made deletion mask built from the generation 0 BED: the first 3 and last 2 bases of every
	introgression, one interior base and 10 bases upstream of its start. The expected BED is the first and last surviving
	base of each interval, renumbered as (old index - deleted bases before it)."""
	out = tmp_path / "o"
	masks = {}

	def fake(seqs, *args, **kwargs):
		if not (out / "ref_0_introgressions.bed").exists():
			# wild relative: unchanged
			return dict(seqs), {c: list(range(len(s))) for c, s in seqs.items()}, {c: list(range(len(s))) for c, s in seqs.items()}
		rows = H.read_bed_rows(out / "ref_0_introgressions.bed")
		new_seqs, mappers = {}, {}
		for c, s in seqs.items():
			dead = np.zeros(len(s), dtype=bool)
			for chrom, st, en, _ in rows:
				if chrom == c:
					dead[st : st + 3] = True
					dead[en - 1 : en + 1] = True
					dead[(st + en) // 2] = True
					dead[st - 10 : st] = True
			masks[c] = dead
			new_idx = np.cumsum(~dead) - 1
			mappers[c] = [-1 if d else int(i) for d, i in zip(dead, new_idx)]
			new_seqs[c] = "".join(b for b, d in zip(s, dead) if not d)
		return new_seqs, mappers, {}

	monkeypatch.setattr(sim, "apply_genome_wide_mutations", fake)
	run_sim(
		monkeypatch, ref_fasta[0], out,
		"--rounds", 1, "--num-introgressions", 4, "--introgression-size-min", 50,
		"--introgression-size-max", 100, "--seed", 6,
	)
	g0 = H.read_bed_rows(out / "ref_0_introgressions.bed")
	g1 = H.read_bed_rows(out / "ref_1_introgressions.bed")
	assert len(g0) == 4 and len(g1) == 4
	for (c, st, en, _), (c1, ns, ne, label) in zip(g0, g1):
		dead = masks[c]
		alive = np.nonzero(~dead)[0]
		first = alive[alive >= st][0]
		last = alive[alive <= en][-1]
		before = lambda i: int(dead[:i].sum())
		assert c1 == c and label == "introgression"
		assert (ns, ne) == (first - before(first), last - before(last))
		assert ns == st + 3 - before(st + 3)


def test_main_relative_indel_options_change_wild_relative_length(monkeypatch, ref_fasta, tmp_path):
	"""--rel-ins-rate 0.001 with fixed size 5: 40 insertions per chromosome, so 200 bases longer. Deletions of fixed size 4 shorten it by a multiple of 4."""
	common = ["--rounds", 0, "--num-introgressions", 2, "--introgression-size-min", 50,
		"--introgression-size-max", 100, "--rel-sub-rate", 0]
	ins = run_sim(monkeypatch, ref_fasta[0], tmp_path / "ins", *common,
		"--rel-ins-rate", 0.001, "--rel-del-rate", 0, "--rel-ins-size-min", 5, "--rel-ins-size-max", 5)
	for c, n in sim.parse_fasta(ins / "ref_wildrelative.fasta").items():
		assert len(n) == L_CHROM + 40 * 5
	dele = run_sim(monkeypatch, ref_fasta[0], tmp_path / "del", *common,
		"--rel-ins-rate", 0, "--rel-del-rate", 0.001, "--rel-del-size-min", 4, "--rel-del-size-max", 4)
	for c, n in sim.parse_fasta(dele / "ref_wildrelative.fasta").items():
		lost = L_CHROM - len(n)
		assert lost % 4 == 0 and 4 * 30 <= lost <= 4 * 40


def test_main_num_introgressions_and_size_options(monkeypatch, ref_fasta, tmp_path):
	out = run_sim(
		monkeypatch, ref_fasta[0], tmp_path / "o",
		"--rounds", 0, "--num-introgressions", 6, "--introgression-size-min", 200,
		"--introgression-size-max", 300, "--rel-ins-rate", 0, "--rel-del-rate", 0, "--seed", 5,
	)
	rows = H.read_bed_rows(out / "ref_0_introgressions.bed")
	assert [r[0] for r in rows] == ["c1"] * 3 + ["c2"] * 3
	assert all(200 <= e - s + 1 <= 300 for _, s, e, _ in rows)


def test_main_rejects_a_reference_without_sequences(monkeypatch, tmp_path):
	empty = tmp_path / "empty.fa"
	empty.write_text("")
	with pytest.raises(ValueError):
		run_sim(monkeypatch, empty, tmp_path / "o")


def test_main_requires_ref_and_out_folder(monkeypatch):
	for argv in (["x", "--out-folder", "o"], ["x", "--ref", "r"]):
		monkeypatch.setattr(sys, "argv", argv)
		with pytest.raises(SystemExit) as e:
			sim.main()
		assert e.value.code == 2

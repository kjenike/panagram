"""Seeded synthetic genomes with shared and diverged segments.

Everything the tests expect about the pan-kmer bitmap follows from this layout,
not from running panagram. Each chromosome is a list of pieces; a piece is either
random filler private to one genome or an instance of a named segment.

Segment membership (which genomes carry the same k-mers in that segment):
	CORE    all five genomes (the copy in delta is reverse-complemented)
	PAIR    alpha and beta (beta carries it twice)
	TRIO    gamma, delta and epsilon
	DIV     all five genomes, but epsilon has a substitution every 10 bp, so
	        with k >= 10 no k-mer of epsilon's copy equals one of the others
"""

import dataclasses
import random

K_MIN = 10  # DIV substitution spacing; every k-mer window needs k >= this

GENOMES = ["alpha", "beta", "gamma", "delta", "epsilon"]
NO_GFF = "epsilon"

SEGMENT_LEN = {"CORE": 600, "PAIR": 400, "TRIO": 300, "DIV": 500}

# Genomes that carry each segment (a non-private instance), hence share its k-mers.
# test_synthetic_layout_matches_its_description checks this against LAYOUT.
MEMBERS = {
	"CORE": {"alpha", "beta", "gamma", "delta", "epsilon"},
	"PAIR": {"alpha", "beta"},
	"TRIO": {"gamma", "delta", "epsilon"},
	"DIV": {"alpha", "beta", "gamma", "delta"},  # epsilon's copy is private
}

# genome -> chrom -> pieces; ("rand", n) private filler, (name, "+"/"-") segment instance.
LAYOUT = {
	"alpha": {
		"chr1": [("rand", 300), ("CORE", "+"), ("rand", 250), ("PAIR", "+"), ("rand", 200),
			("DIV", "+"), ("rand", 400)],
		"chr2": [("rand", 500), ("CORE", "-"), ("rand", 600)],
	},
	"beta": {
		"chr1": [("rand", 200), ("PAIR", "+"), ("rand", 100), ("PAIR", "+"), ("rand", 300),
			("CORE", "+"), ("rand", 150), ("DIV", "+"), ("rand", 250)],
	},
	"gamma": {
		"chr1": [("rand", 250), ("TRIO", "+"), ("rand", 300), ("DIV", "+"), ("rand", 200),
			("CORE", "+"), ("rand", 350)],
		"chr2": [("rand", 700)],
	},
	"delta": {
		"chr1": [("rand", 150), ("CORE", "-"), ("rand", 300), ("TRIO", "+"), ("rand", 200),
			("DIV", "+"), ("rand", 300)],
	},
	"epsilon": {
		"chr1": [("rand", 400), ("TRIO", "+"), ("rand", 200), ("CORE", "+"), ("rand", 300),
			("DIV", "+"), ("rand", 200)],
	},
}

_COMP = str.maketrans("ACGT", "TGCA")


def revcomp(seq):
	return seq.translate(_COMP)[::-1]


def _rand_seq(rng, n):
	return "".join(rng.choice("ACGT") for _ in range(n))


def _mutate(seq, spacing):
	"""Substitute every `spacing`-th base (offset spacing-1) with the next base in ACGT order."""
	out = list(seq)
	for i in range(spacing - 1, len(out), spacing):
		out[i] = "ACGT"["ACGT".index(out[i]) - 3]  # shifts by +1 mod 4
	return "".join(out)


@dataclasses.dataclass
class Piece:
	kind: str  # "rand" or a segment name
	strand: str  # "+" / "-" / ""
	start: int  # 0-based, in the chromosome
	end: int  # exclusive
	private: bool  # no other genome shares its k-mers


@dataclasses.dataclass
class SyntheticSet:
	seqs: dict  # genome -> chrom -> sequence
	pieces: dict  # genome -> chrom -> [Piece]
	genes: dict  # genome -> [(chrom, start0, end0, name)]

	def kmer_owners(self, genome, chrom, k):
		"""Expected bitmap membership per k-mer start, derived from the layout only.

		Entry p is a set of genome names for k-mers wholly inside one segment instance
		(the segment's members) or wholly inside private filler (the genome alone), and
		None for k-mers that straddle a boundary: those share most of their bases with
		a segment and may match another genome by chance (a 1-base spill matches with
		probability 1/4), so the layout does not determine them.
		"""
		n = len(self.seqs[genome][chrom]) - k + 1
		owners = [None] * n
		for p in self.pieces[genome][chrom]:
			if p.private or p.kind == "rand":
				members = {genome}
			elif genome in MEMBERS[p.kind]:
				members = set(MEMBERS[p.kind])
			else:
				members = {genome}
			for pos in range(p.start, p.end - k + 1):
				owners[pos] = members
		return owners


def build(seed=20261109):
	"""Build the genome set. Same seed gives the same sequences."""
	rng = random.Random(seed)
	segs = {name: _rand_seq(rng, n) for name, n in SEGMENT_LEN.items()}
	seqs, pieces = {}, {}
	for genome in GENOMES:
		seqs[genome], pieces[genome] = {}, {}
		for chrom, layout in LAYOUT[genome].items():
			parts, plist, pos = [], [], 0
			for kind, arg in layout:
				if kind == "rand":
					s = _rand_seq(rng, arg)
					plist.append(Piece("rand", "", pos, pos + len(s), True))
				else:
					s = segs[kind]
					private = False
					if kind == "DIV" and genome == "epsilon":
						s, private = _mutate(s, K_MIN), True
					if arg == "-":
						s = revcomp(s)
					plist.append(Piece(kind, arg, pos, pos + len(s), private))
				parts.append(s)
				pos += len(s)
			seqs[genome][chrom] = "".join(parts)
			pieces[genome][chrom] = plist

	# Genes: one covering each shared segment instance on chr1, plus one private filler gene.
	genes = {}
	for genome in GENOMES:
		if genome == NO_GFF:
			continue
		glist = []
		for chrom, plist in pieces[genome].items():
			for i, p in enumerate(plist):
				if p.kind != "rand":
					glist.append((chrom, p.start, p.end, f"{genome}_{chrom}_{p.kind}_{i}"))
			filler = next(p for p in plist if p.kind == "rand")
			glist.append((chrom, filler.start, filler.end, f"{genome}_{chrom}_private"))
		genes[genome] = glist
	return SyntheticSet(seqs, pieces, genes)


def write_fasta(path, chroms, width=60):
	with open(path, "w") as f:
		for chrom, seq in chroms.items():
			f.write(f">{chrom}\n")
			for i in range(0, len(seq), width):
				f.write(seq[i : i + width] + "\n")


def write_gff3(path, genes):
	"""GFF3 with a gene and one exon child per gene (1-based inclusive coordinates)."""
	with open(path, "w") as f:
		f.write("##gff-version 3\n")
		for chrom, start0, end0, name in genes:
			f.write(f"{chrom}\tsynthetic\tgene\t{start0 + 1}\t{end0}\t.\t+\t.\tID={name};Name={name}\n")
			f.write(
				f"{chrom}\tsynthetic\texon\t{start0 + 1}\t{end0}\t.\t+\t.\tID={name}.exon;Parent={name}\n"
			)

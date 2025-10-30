from pathlib import Path
import argparse
from collections import OrderedDict
import numpy as np
from scipy.stats import beta


def parse_fasta(path):
    """Parse FASTA file into an ordered dict of sequences.
    Args:
        path (str or Path): path to FASTA file
    Returns:
        OrderedDict: keys are sequence names, values are sequences (str)
    """

    seqs = OrderedDict()
    name = None
    seq_lines = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(seq_lines).upper()
                name = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        if name is not None:
            seqs[name] = "".join(seq_lines).upper()
    return seqs


def write_fasta(seqs, path):
    """
    Write sequences to a FASTA file.
    Args:
        seqs (dict): keys are sequence names, values are sequences (str)
        path (str or Path): path to output FASTA file
    """

    with open(path, "w") as f:
        for name, seq in seqs.items():
            f.write(f">{name}\n")
            f.write(seq + "\n")


def write_bed(bed_entries, bed_path):
    """Write BED file with introgression entries.
    Args:
        bed_entries (list of str): list of BED entries as strings
        bed_path (str or Path): path to output BED file
    """

    # inclusive coordinates of introgressions
    with open(bed_path, "w") as bed:
        bed.write("#chrom\tstart\tend\tnotes\n")
        for entry in bed_entries:
            bed.write(entry + "\n")
    return


def mutate_seq_with_introgressions(
    chrom,
    ref_seq,
    rel_seq,
    reverse_mapper,
    available_positions,
    n_introgressions,
    introgression_lengths,
    rng=None,
):
    """Apply introgressions by copying segments from rel_seq into ref_seq. Introgression sizes are
    sampled from a skewed distribution between size_min and size_max (inclusive).
    Args:
        chrom (str): chromosome name for bed entries
        ref_seq (str): reference sequence
        rel_seq (str): wild relative sequence
        reverse_mapper (list of int): mapping from ref to rel positions
        available_positions (list of int): available positions for introgressions
        n_introgressions (int): number of introgressions to apply
        introgression_lengths (list of int): lengths of introgressions
        rng (np.random.Generator): random number generator
    Returns:
        new_seq (str): sequence after introgressions applied
        bed_entries (list of str): list of BED entries for introgressions
    """

    bed_entries = []
    length_chr = len(ref_seq)
    introgressed_seq = []

    # choose introgressions on this chromosome from positions that are still valid
    start_positions = rng.choice(available_positions, size=n_introgressions, replace=False)

    # sort them
    start_positions.sort()

    # ensure spacing by adjusting positions
    required_min = start_positions[0] + np.cumsum(np.concatenate(([0], introgression_lengths[:-1])))
    start_positions = np.maximum(start_positions, required_min)
    end_positions = start_positions + introgression_lengths

    current_index = 0
    for i in range(len(start_positions)):

        start_position = int(start_positions[i])
        end_position = int(end_positions[i])

        # check if end position is within chromosome before mapping
        if end_position > length_chr:
            continue  # skip this introgression

        mapped_start = reverse_mapper[start_position]
        mapped_end = reverse_mapper[end_position]

        # try to correct end position if it maps to -1 (deleted region)
        while mapped_end == -1:
            end_position -= 1
            if end_position <= start_position:
                break  # cannot introgress anything
            mapped_end = reverse_mapper[end_position]
        if end_position <= start_position:
            continue  # skip this introgression

        # add unchanged sequence before introgression
        introgressed_seq.extend(ref_seq[current_index:start_position])
        current_index = end_position

        # apply introgression
        introgression = rel_seq[mapped_start:mapped_end]
        introgressed_seq.extend(introgression)

        # add an entry to bed
        bed_entries.append((f"{chrom}\t{start_position}\t{end_position - 1}\tintrogression"))

    # add remaining sequence after last introgression
    introgressed_seq.extend(ref_seq[current_index:])
    return "".join(introgressed_seq), bed_entries


def apply_genome_wide_introgressions(
    ref_seqs,
    rel_seqs,
    reverse_mappers,
    available_positions,
    n_introgressions,
    size_min,
    size_max,
    rng=None,
):
    """Apply introgressions across all chromosomes and save corresponding BED entries.
    Args:
        ref_seqs (dict): reference sequences (chrom -> seq)
        rel_seqs (dict): wild relative sequences (chrom -> seq)
        reverse_mappers (dict): mapping from ref to rel positions (chrom -> list of int)
        available_positions (dict): available positions for introgressions (chrom -> list of int)
        n_introgressions (int): number of introgressions to apply
        size_min (int): minimum introgression size
        size_max (int): maximum introgression size
        rng (np.random.Generator): random number generator
    Returns:
        new_seqs (dict): sequences after introgressions applied (chrom -> seq)
        bed_entries (list of str): list of BED entries for introgressions
    """

    bed_entries = []
    new_seqs = {}
    chroms = list(ref_seqs.keys())
    # generate all introgression lengths ahead of time; gives a better distribution of sizes
    # TODO: make a parameter to control this behavior
    # uniformly sample introgression sizes
    intro_lengths = rng.integers(size_min, size_max + 1, size=n_introgressions).tolist()
    # intro_lengths = generate_skewed_sizes(n_introgressions, size_min, size_max, rng)

    introgressions_completed = 0

    for chrom in chroms:
        print(f"Applying introgressions on chromosome {chrom}...", flush=True)

        n_introgressions_chr = n_introgressions // len(chroms)
        intro_lengths_chr = intro_lengths[
            introgressions_completed : introgressions_completed + n_introgressions_chr
        ]

        introgressed_seq, bed_entries_chr = mutate_seq_with_introgressions(
            chrom=chrom,
            ref_seq=ref_seqs[chrom],
            rel_seq=rel_seqs[chrom],
            reverse_mapper=reverse_mappers[chrom],
            available_positions=available_positions[chrom],
            n_introgressions=n_introgressions_chr,
            introgression_lengths=intro_lengths_chr,
            rng=rng,
        )

        introgressions_completed += n_introgressions_chr
        bed_entries.extend(bed_entries_chr)
        new_seqs[chrom] = introgressed_seq

    return new_seqs, bed_entries


def generate_skewed_sizes(n, size_min, size_max, rng, a=0.05, b=1):
    """Generate n sizes from a beta distribution skewed towards smaller sizes.
    Args:
        n (int): number of sizes to generate
        size_min (int): minimum size
        size_max (int): maximum size
        rng (np.random.Generator): random number generator
        a (float): alpha parameter of beta distribution (default: 0.05)
        b (float): beta parameter of beta distribution (default: 1)
    Returns:
        list of int: generated sizes
    """

    # generate sizes from beta distribution skewed towards smaller sizes
    # used for introgression and indel sizes
    vals = beta.rvs(a, b, size=n, random_state=rng)
    nums = size_min + vals * (size_max - size_min)
    nums = [int(round(x)) for x in nums]
    return nums


def generate_mutation_positions(
    mutation_type, n_positions, size_min, size_max, available_positions, rng
):
    """
    Generate random positions and lengths for mutations of a given type without overlap or replacement.
    Args:
        mutation_type (str): type of mutation ("deletion", "insertion", "snp")
        n_positions (int): number of mutation positions to generate
        size_min (int): minimum size of mutation
        size_max (int): maximum size of mutation
        available_positions (list of int): list of available positions to choose from
        rng (np.random.Generator): random number generator
    Returns:
        mutation_positions (list of tuples): list of (position, length, mutation_type)
        available_positions (list of int): updated list of available positions after mutations
    """

    # create random positions for each mutation type without overlap or replacement
    mutation_positions = []

    # generate lengths from beta distribution skewed towards smaller sizes
    lengths = generate_skewed_sizes(n_positions, size_min, size_max, rng)

    if mutation_type == "deletion":
        # choose n_positions from available_positions with spacing of size_max
        deletion_positions = rng.choice(available_positions, size=n_positions, replace=False)

        # sort them
        deletion_positions.sort()

        # ensure spacing by adjusting positions
        required_min = deletion_positions[0] + np.cumsum(np.concatenate(([0], lengths[:-1])))
        deletion_positions = np.maximum(deletion_positions, required_min)

        # remove any deletions that don’t fully fit in chromosome
        valid = deletion_positions + lengths <= len(available_positions)
        if not np.all(valid):
            deletion_positions = deletion_positions[valid]
            print(f"Warning: could only fit {len(deletion_positions)} deletions.", flush=True)

        # apply deletions and mark positions as unavailable for insertions/SNPs
        deletion_positions = deletion_positions.tolist()
        for i in range(len(deletion_positions)):
            length = lengths[i]
            pos = deletion_positions.pop()
            available_positions[pos:pos + length] = [-1]*length # mark for deletion
            mutation_positions.append((pos, length, mutation_type))

        # clean up available_positions to remove marked positions
        arr = np.asarray(available_positions)
        available_positions = arr[arr != -1].tolist() # have to convert to list to get rid of dtype issues

    else:
        # insertions and SNPs can go anywhere in remaining available positions
        for length in lengths:
            if not available_positions:
                print("Warning: no more available positions for insertions/SNPs.", flush=True)
                break
            pos = available_positions.pop()
            mutation_positions.append((pos, length, mutation_type))

    return mutation_positions, available_positions


def mutate_seq_with_indels_and_snps(
    seq, sub_rate, ins_rate, del_rate, ins_size_min, ins_size_max, del_size_min, del_size_max, rng
):
    """
    Precompute random positions for substitutions, deletions, insertions, then apply them.
    Use apply_genome_wide_mutations() to apply to all sequences in a dict.

    Args:
        seq (str): input sequence
        sub_rate (float): per-base substitution rate
        ins_rate (float): per-base insertion rate
        del_rate (float): per-base deletion rate
        ins_size_min (int): minimum insertion size
        ins_size_max (int): maximum insertion size
        del_size_min (int): minimum deletion size
        del_size_max (int): maximum deletion size
        rng (np.random.Generator): random number generator
    Returns:
        mutated_seq (str): sequence after mutations applied
        reverse_mapper (np.array): mapping from original indices to new indices (-1 for deleted)
        available_positions (list of int): remaining non-mutated positions
    """

    L = len(seq)

    # mutation counts
    n_sub = int(L * sub_rate)
    n_ins = int(L * ins_rate)
    n_del = int(L * del_rate)

    # generate deletions first since they have spacing requirements; shuffling done inside function
    # DO NOT shuffle available positions before deletions - they need to be in order for updating
    available_positions = list(range(L))

    if n_del > 0:
        del_positions, available_positions = generate_mutation_positions(
            "deletion", n_del, del_size_min, del_size_max, available_positions, rng
        )
    else:
        del_positions = []
    # shuffle remaining available positions for insertions and SNPs
    rng.shuffle(available_positions)
    if n_ins > 0:
        ins_positions, available_positions = generate_mutation_positions(
            "insertion", n_ins, ins_size_min, ins_size_max, available_positions, rng
        )
    else:
        ins_positions = []
    if n_sub > 0:
        sub_positions, available_positions = generate_mutation_positions(
            "snp", n_sub, 1, 2, available_positions, rng
        )
    else:
        sub_positions = []
    all_mutations = ins_positions + del_positions + sub_positions
    all_mutations.sort()

    # TODO: change back to upper/lower case if needed
    bases = ["a", "c", "g", "t"]
    current_index = 0
    offset = 0
    new_seq = []
    # index mapper for introgressions/alignment: original_index -> new index
    reverse_mapper = []

    for pos, length, mtype in all_mutations:
        if mtype == "insertion":
            # append everything up to and including pos
            # insertion is AFTER pos
            new_seq.extend(seq[current_index:pos + 1])
            offset_indicies = range(current_index + offset, pos + 1 + offset)
            reverse_mapper.extend(offset_indicies)

            # update current index and offset
            current_index = pos + 1
            offset += length

            # build insertion sequence
            ins_seq = rng.choice(bases, size=length, replace=True)
            new_seq.extend(ins_seq)

        if mtype == "deletion":
            # append everything up to pos
            # deletion starts at and includes pos
            new_seq.extend(seq[current_index:pos])
            offset_indicies = range(current_index + offset, pos + offset)
            reverse_mapper.extend(offset_indicies)

            # update current index and offset to skip deleted region
            current_index = pos + length
            offset -= length

            # pad reverse_mapper for deleted indices
            reverse_mapper.extend([-1]*length)

        if mtype == "snp":
            # append everything up to pos
            new_seq.extend(seq[current_index:pos])
            # include pos index for mapping; the index of a SNP is the same
            offset_indicies = range(current_index + offset, pos + 1 + offset)
            reverse_mapper.extend(offset_indicies)

            # change base at pos
            old_base = seq[pos]
            new_base = rng.choice([b for b in bases if b != old_base])
            new_seq.append(new_base)

            # update current index
            current_index = pos + 1

    # append remaining sequence after last mutation
    new_seq.extend(seq[current_index:])
    offset_indicies = range(current_index + offset, L + offset)
    reverse_mapper.extend(offset_indicies)

    return "".join(new_seq), reverse_mapper, available_positions


def apply_genome_wide_mutations(
    seqs, sub_rate, ins_rate, del_rate, ins_size_min, ins_size_max, del_size_min, del_size_max, rng
):
    """Apply per-base mutation rates (SNPs, indels) to each sequence in seqs dict.
    Args:
        seqs (dict): input sequences (chrom -> seq)
        sub_rate (float): per-base substitution rate
        ins_rate (float): per-base insertion rate
        del_rate (float): per-base deletion rate
        ins_size_min (int): minimum insertion size
        ins_size_max (int): maximum insertion size
        del_size_min (int): minimum deletion size
        del_size_max (int): maximum deletion size
        rng (np.random.Generator): random number generator
    Returns:
        new_seqs (dict): sequences after mutations applied (chrom -> seq)
        mappers (dict): mapping from original indices to new indices for each chrom (chrom -> list of int)
        available_positions (dict): remaining non-mutated positions for each chrom (chrom -> list of int)
    """

    new_seqs = {}
    mappers = {}
    available_positions = {}
    for chrom, seq in seqs.items():
        print(f"Mutating chromosome {chrom}...", flush=True)
        new_seqs[chrom], mappers[chrom], available_positions[chrom] = mutate_seq_with_indels_and_snps(
            seq,
            sub_rate,
            ins_rate,
            del_rate,
            ins_size_min,
            ins_size_max,
            del_size_min,
            del_size_max,
            rng,
        )
    return new_seqs, mappers, available_positions


def main():
    parser = argparse.ArgumentParser(
        description="Simulate introgressions and mutations with indels and SNPs."
    )
    parser.add_argument("--ref", required=True, help="Reference FASTA (input)")
    parser.add_argument("--out-folder", required=True, help="Output folder")
    parser.add_argument(
        "--offspring",
        type=int,
        default=1,
        help="Number of offspring to simulate from created wild relative (default: 1)",
    )
    parser.add_argument(
        "--rounds",
        type=int,
        default=4,
        help="Number of rounds to apply mutations to offspring. Will output offspring for each round (default: 4)",
    )

    # introgressions (counts and size range)
    parser.add_argument(
        "--num-introgressions",
        type=int,
        default=48,
        help="Number of introgressions per generation (default: 48)",
    )
    parser.add_argument(
        "--introgression-size-min", type=int, default=50000, help="Min introgression size (bp)"
    )
    parser.add_argument(
        "--introgression-size-max", type=int, default=10000000, help="Max introgression size (bp)"
    )

    # wild relative mutation rates and size ranges
    # by default, total mutation rate is 1e-1 if you add sub, ins, del rates
    parser.add_argument(
        "--rel-sub-rate",
        type=float,
        default=3.3e-3,
        help="Wild relative per-base substitution rate (default: 3.3e-3)",
    )
    parser.add_argument(
        "--rel-ins-rate",
        type=float,
        default=3.3e-3,
        help="Wild relative per-base insertion rate (default: 3.3e-3)",
    )
    parser.add_argument(
        "--rel-del-rate",
        type=float,
        default=3.3e-3,
        help="Wild relative per-base deletion rate (default: 3.3e-3)",
    )
    parser.add_argument(
        "--rel-ins-size-min", type=int, default=1, help="Wild relative insertion size min (bp)"
    )
    parser.add_argument(
        "--rel-ins-size-max", type=int, default=1000, help="Wild relative insertion size max (bp)"
    )
    parser.add_argument(
        "--rel-del-size-min", type=int, default=1, help="Wild relative deletion size min (bp)"
    )
    parser.add_argument(
        "--rel-del-size-max", type=int, default=1000, help="Wild relative deletion size max (bp)"
    )

    # offspring mutation rates and size ranges
    # by default, the total mutation rate is 1e-2 if you add sub, ins, del rates
    parser.add_argument(
        "--mut-sub-rate",
        type=float,
        default=3.3e-4,
        help="Offspring per-base substitution rate (default: 3.3e-4)",
    )
    parser.add_argument(
        "--mut-ins-rate",
        type=float,
        default=3.3e-4,
        help="Offspring per-base insertion rate (default: 3.3e-4)",
    )
    parser.add_argument(
        "--mut-del-rate",
        type=float,
        default=3.3e-4,
        help="Offspring per-base deletion rate (default: 3.3e-4)",
    )
    parser.add_argument(
        "--mut-ins-size-min", type=int, default=1, help="Offspring insertion size min (bp)"
    )
    parser.add_argument(
        "--mut-ins-size-max", type=int, default=1000, help="Offspring insertion size max (bp)"
    )
    parser.add_argument(
        "--mut-del-size-min", type=int, default=1, help="Offspring deletion size min (bp)"
    )
    parser.add_argument(
        "--mut-del-size-max", type=int, default=1000, help="Offspring deletion size max (bp)"
    )

    parser.add_argument(
        "--seed", type=int, default=42, help="Random seed for reproducibility (default: 42)"
    )

    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    output_folder = Path(args.out_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    # read in reference as an ordered dict with each chr
    reference = Path(args.ref)
    ref_seqs = parse_fasta(reference)
    if not ref_seqs:
        raise ValueError("ERROR: no sequences read from", reference)

    # create wild relative by mutating reference
    rel_seqs, reverse_mappers, available_positions = apply_genome_wide_mutations(
        ref_seqs,
        sub_rate=args.rel_sub_rate,
        ins_rate=args.rel_ins_rate,
        del_rate=args.rel_del_rate,
        ins_size_min=args.rel_ins_size_min,
        ins_size_max=args.rel_ins_size_max,
        del_size_min=args.rel_del_size_min,
        del_size_max=args.rel_del_size_max,
        rng=rng,
    )

    # write wild relative to file
    print("Writing wild relative FASTA...", flush=True)
    write_fasta(rel_seqs, output_folder / f"{reference.stem}_wildrelative.fasta")

    with open(output_folder / f"{reference.stem}_wildrelative.mapper", "w") as mapf:
        for chrom in ref_seqs:
            mapper = reverse_mappers[chrom]
            mapf.write(f">{chrom}_mapper\n")
            mapf.write(" ".join(map(str, mapper)) + "\n")
            mapf.write(f">{chrom}_available_positions\n")
            mapf.write(" ".join(map(str, available_positions[chrom])) + "\n")

    # apply introgressions back to reference to create offspring
    print("Applying introgressions to create generation 1 offspring...", flush=True)
    offspring_seqs, introgressions = apply_genome_wide_introgressions(
        ref_seqs,
        rel_seqs,
        reverse_mappers,
        available_positions,
        args.num_introgressions,
        args.introgression_size_min,
        args.introgression_size_max,
        rng=rng,
    )

    # write offspring fasta
    print("Writing gen. 1 introgressed offspring FASTA...", flush=True)
    write_fasta(offspring_seqs, output_folder / f"{reference.stem}_gen1_offspring.fasta")
    write_bed(introgressions, output_folder / f"{reference.stem}_gen1_introgressions.bed")

    parent_seqs = offspring_seqs
    for gen in range(2, args.rounds + 1):
        print(f"Simulating generation {gen}...", flush=True)

        # apply mutations to offspring from previous generation
        offspring_seqs, _, _ = apply_genome_wide_mutations(
            parent_seqs,
            sub_rate=args.mut_sub_rate,
            ins_rate=args.mut_ins_rate,
            del_rate=args.mut_del_rate,
            ins_size_min=args.mut_ins_size_min,
            ins_size_max=args.mut_ins_size_max,
            del_size_min=args.mut_del_size_min,
            del_size_max=args.mut_del_size_max,
            rng=rng,
        )
        # save per-gen FASTA
        write_fasta(offspring_seqs, output_folder / f"{reference.stem}_gen{gen}_offspring.fasta")
        parent_seqs = offspring_seqs

    print("Simulation finished.")
    return


# # Dummy generate_mutation_positions for testing
# def generate_mutation_positions(mtype, n, size_min, size_max, available_positions, rng):
#     if mtype == "insertion":
#         return [(1, 2, "insertion")], available_positions
#     elif mtype == "deletion":
#         return [(3, 2, "deletion")], available_positions
#     elif mtype == "snp":
#         return [(0, 1, "snp")], available_positions
#     return [], available_positions


def test_apply_mutations():
    # --- Use your function exactly as is, just override the behavior ---
    rng = np.random.default_rng(1)
    seq = "GATTACA"

    mut_seq, reverse_map, remaining_positions = mutate_seq_with_indels_and_snps(
        seq=seq,
        sub_rate=1/7,   # Force 1 SNP
        ins_rate=1/7,   # Force 1 insertion
        del_rate=1/7,   # Force 1 deletion
        ins_size_min=2,
        ins_size_max=2,
        del_size_min=2,
        del_size_max=2,
        rng=rng
    )

    print("Original: ", seq)
    print("Mutated : ", mut_seq)
    print("Reverse map:", reverse_map)
    return


def test_apply_introgressions():
    # test introgressions
    # Reference sequences
    ref_seqs = {
        "chr1": "AAAAAAAAAATTTTTTTTTTCCCCCCCCCC",  # length 30
        "chr2": "GGGGGGGGGGAAAAAAAAAATTTTTTTTTT"
    }

    # Relative sequences (introgressed sequence will come from here)
    rel_seqs = {
        "chr1": "aaaaaaaaaattttttttttcccccccccc",
        "chr2": "ggggggggggaaaaaaaaaatttttttttt"
    }

    # Reverse mapper (identity mapping — assume no deletions applied before)
    reverse_mappers = {
        "chr1": list(range(len(ref_seqs["chr1"]))),
        "chr2": list(range(len(ref_seqs["chr2"])))
    }

    # All positions available (no deletions or insertions earlier)
    available_positions = {
        "chr1": list(range(len(ref_seqs["chr1"]))),
        "chr2": list(range(len(ref_seqs["chr2"])))
    }

    rng = np.random.default_rng(1)

    new_seqs, bed = apply_genome_wide_introgressions(
        ref_seqs=ref_seqs,
        rel_seqs=rel_seqs,
        reverse_mappers=reverse_mappers,
        available_positions=available_positions,
        n_introgressions=3,
        size_min=3,
        size_max=3,
        rng=rng
    )

    print("\nModified Sequences:")
    for chrom in new_seqs:
        print(chrom, ":", new_seqs[chrom], len(new_seqs[chrom]))

    print("\nBED Entries:")
    for b in bed:
        print(b)
    return


def test_apply_introgressions_with_insertion():
    # Reference sequences
    ref_seqs = {
        "chr1": "AAAAAAAAAATTTTTTTTTTCCCCCCCCCC",  # length 30
        "chr2": "GGGGGGGGGGAAAAAAAAAATTTTTTTTTT"
    }

    # Relative sequences with an insertion of length 5 at position 16 on chr1
    insertion_pos = 16  # 0-based
    insertion_seq = "NNNNN"
    chr1 = ref_seqs["chr1"][:insertion_pos] + insertion_seq + ref_seqs["chr1"][insertion_pos:]
    chr1 = chr1.lower()  # make it lowercase to distinguish
    rel_seqs = {
        "chr1": chr1,  # length 35
        "chr2": "ggggggggggaaaaaaaaaatttttttttt"  # unchanged
    }

    # Build reverse mapper: maps ref index → rel index
    reverse_mapper_chr1 = []
    for i in range(len(ref_seqs["chr1"])):
        if i < insertion_pos:
            reverse_mapper_chr1.append(i)           # before insertion: same index
        else:
            reverse_mapper_chr1.append(i + len(insertion_seq))  # after insertion: shift by insertion length

    reverse_mappers = {
        "chr1": reverse_mapper_chr1,
        "chr2": list(range(len(ref_seqs["chr2"])))  # unchanged
    }
    # All positions available (no deletions applied previously)
    available_positions = {
        "chr1": list(range(len(ref_seqs["chr1"]))),
        "chr2": list(range(len(ref_seqs["chr2"])))
    }

    rng = np.random.default_rng(1)

    # Apply introgressions
    new_seqs, bed = apply_genome_wide_introgressions(
        ref_seqs=ref_seqs,
        rel_seqs=rel_seqs,
        reverse_mappers=reverse_mappers,
        available_positions=available_positions,
        n_introgressions=3,  # just one for testing
        size_min=3,
        size_max=3,
        rng=rng
    )

    # Print results
    print("\nModified Sequences:")
    for chrom in new_seqs:
        print(chrom, ":", new_seqs[chrom], len(new_seqs[chrom]))

    print("\nBED Entries:")
    for b in bed:
        print(b)


if __name__ == "__main__":
    # test_apply_mutations()
    # test_apply_introgressions()
    test_apply_introgressions_with_insertion()
    # main()

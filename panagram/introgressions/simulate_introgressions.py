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


def write_fasta(seqs, path, wrap=60):
    """Write sequences to a FASTA file with configurable line wrapping.

    Args:
        seqs (dict): keys are sequence names, values are sequences (str)
        path (str or Path): path to output FASTA file
        wrap (int): number of bases per line (default: 60)
    """

    with open(path, "w") as f:
        for name, seq in seqs.items():
            f.write(f">{name}\n")
            for i in range(0, len(seq), wrap):
                f.write(seq[i : i + wrap] + "\n")
    return


def write_bed(bed_entries, bed_path):
    """Write BED file with introgression entries.

    Args:
        bed_entries (list of str): list of BED entries as strings
        bed_path (str or Path): path to output BED file
    """

    # inclusive coordinates of introgressions
    with open(bed_path, "w") as bed:
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
    rng,
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
        tuple(str, list): sequence after introgressions applied, BED entries for introgressions
    """

    bed_entries = []
    length_chr = len(ref_seq)
    introgressed_seq = []

    # choose introgressions on this chromosome from positions that are still valid
    start_positions = rng.choice(available_positions, size=n_introgressions, replace=False)

    # ensure spacing by adjusting positions
    start_positions = position_spacer(start_positions.tolist(), introgression_lengths, length_chr)
    end_positions = [start + length for start, length in start_positions]
    start_positions = [
        start for start, length in start_positions
    ]  # unzip start positions from tuples

    current_index = 0
    for i in range(len(start_positions)):

        start_position = int(start_positions[i])
        end_position = int(end_positions[i])

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
        actual_length = len(introgression)

        # add an entry to bed
        bed_entries.append(
            (f"{chrom}\t{start_position}\t{end_position - 1}\t{actual_length} bp introgression")
        )

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
    rng,
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
        tuple(dict, list): sequences after introgressions applied (chrom -> seq), BED entries for
        introgressions
    """

    bed_entries = []
    new_seqs = {}
    chroms = list(ref_seqs.keys())
    # generate all introgression lengths ahead of time; gives a better distribution of sizes
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
    """Generate n sizes from a beta distribution skewed towards smaller sizes. Used for sampling
    indel sizes.

    Args:
        n (int): number of sizes to generate
        size_min (int): minimum size
        size_max (int): maximum size
        rng (np.random.Generator): random number generator
        a (float): alpha parameter of beta distribution (default: 0.05)
        b (float): beta parameter of beta distribution (default: 1)

    Returns:
        list: generated sizes
    """

    # generate sizes from beta distribution skewed towards smaller sizes
    vals = beta.rvs(a, b, size=n, random_state=rng)
    nums = size_min + vals * (size_max - size_min)
    nums = [int(round(x)) for x in nums]
    return nums


def create_random_trenches(weights, rng, block_size=1_000_000, drop_fraction=0.3):
    """
    Apply random trenches to the weights by reducing their values in certain blocks.

    Args:
        weights (np.ndarray): array of weights to modify
        block_size (int): size of the blocks to consider for dropping weights
        drop_fraction (float): probability of dropping weights in a block
        rng (np.random.Generator): random number generator

    Returns:
        np.ndarray: modified weights with trenches applied
    """
    if rng is None:
        rng = np.random.default_rng()

    n = len(weights)
    num_blocks = n // block_size

    for i in range(num_blocks):
        if rng.random() < drop_fraction:
            start = i * block_size
            end = start + block_size
            weights[start:end] *= 0.8

    return weights / weights.sum()


def generate_edge_tapered_weights(n, rng, edge_fraction=0.3, edge_power=5):
    """
    Generate n weights with a trench shape: high at edges, low in middle. To be used for sampling
    indel positions with higher probability at edges of the chromosome. Adds random noise to increase
    variability.

    Args:
        n (int): number of weights to generate
        edge_fraction (float): proportion of positions considered "edges"
        edge_power (float): how steep the edge taper is
        rng (np.random.Generator): random number generator

    Returns:
        np.ndarray: array of weights with trench shape
    """
    x = np.linspace(0, 1, n)

    # distance from nearest edge
    d = np.minimum(x, 1 - x)

    # define how far both edges extend
    edge_cut = edge_fraction / 2
    edge_taper = 1 - (d / edge_cut) ** edge_power
    edge_taper = np.clip(edge_taper, 0.5, None)  # limit minimum taper to avoid zero weights

    # edges: high, middle: constant low, smooth transitions controlled by edge_taper
    weights = np.where(d < edge_cut, edge_taper, 0.5)
    weights = weights / weights.sum()

    # add some random trenches to increase variability
    weights = create_random_trenches(weights, rng=rng)

    return weights


def position_spacer(positions, lengths, chr_size):
    """Add space to positions to ensure no overlap given lengths and max size. Removes positions
    that cannot fit within chromosome size.

    Args:
        positions (list of int): candidate positions
        lengths (list of int): corresponding lengths
        chr_size (int): size of chromosome to avoid overflow

    Returns:
        list[tuple]: selected positions and lengths without overlap
    """
    selected_positions = []
    current_end = -1
    positions.sort()
    for position, length in zip(positions, lengths):
        if position + length > chr_size:
            break  # skip other positions; we're at the end of the chromosome
        if position > current_end:
            selected_positions.append((position, length))
            current_end = position + length
        else:
            # adjust position to ensure spacing
            position = current_end + 1
            if position + length > chr_size:
                break  # skip other positions; we're at the end of the chromosome
            else:
                selected_positions.append((position, length))
                current_end = position + length
    if len(selected_positions) < len(positions):
        print(
            f"Warning: could only fit {len(selected_positions)} / {len(positions)} mutations.",
            flush=True,
        )
    return selected_positions


def generate_mutation_positions(
    mutation_type, n_positions, size_min, size_max, position_weights, rng
):
    """Generate random positions and lengths for mutations of a given type without overlap or replacement.

    Args:
        mutation_type (str): type of mutation ("deletion", "insertion", "snp")
        n_positions (int): number of mutation positions to generate
        size_min (int): minimum size of mutation
        size_max (int): maximum size of mutation
        position_weights (array of float): corresponding weights for available positions
        rng (np.random.Generator): random number generator

    Returns:
        tuple(list, np.ndarray): list of (position, length, mutation_type), updated list of
        available weights after mutations
    """

    # create random positions for each mutation type without overlap or replacement
    mutation_positions = []

    # generate lengths from beta distribution skewed towards smaller sizes
    lengths = generate_skewed_sizes(n_positions, size_min, size_max, rng)

    if mutation_type == "deletion":
        # choose n_positions from available_positions with spacing of size_max
        chr_size = len(position_weights)
        deletion_positions = rng.choice(
            chr_size, size=n_positions, replace=False, p=position_weights
        )
        deletion_positions = position_spacer(deletion_positions.tolist(), lengths, chr_size)
        for _ in range(len(deletion_positions)):
            pos, length = deletion_positions.pop()
            position_weights[pos : pos + length] = [0] * length  # ensure weights are zeroed out
            mutation_positions.append((pos, length, mutation_type))

    else:
        # insertions and SNPs can go anywhere in remaining available positions
        num_available_positions = np.count_nonzero(position_weights)
        if n_positions > num_available_positions:
            n_positions = num_available_positions
            print("Warning: no more available positions for insertions/SNPs.", flush=True)
        chr_size = len(position_weights)
        position_indices = rng.choice(chr_size, size=n_positions, replace=False, p=position_weights)

        for pos, length in zip(position_indices, lengths):
            position_weights[pos] = 0  # mark position as used
            mutation_positions.append((pos, length, mutation_type))

    # renormalize weights
    if position_weights.sum() > 0:
        position_weights = position_weights / position_weights.sum()
    return mutation_positions, position_weights


def mutate_seq_with_indels_and_snps(
    seq, sub_rate, ins_rate, del_rate, ins_size_min, ins_size_max, del_size_min, del_size_max, rng
):
    """Precompute random positions for substitutions, deletions, insertions, then apply them.
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
        tuple(str, np.ndarray, list): sequence after mutations applied, mapping from original
        indices to new indices (-1 for deleted), remaining non-mutated positions
    """

    L = len(seq)

    # mutation counts
    n_sub = int(L * sub_rate)
    n_ins = int(L * ins_rate)
    n_del = int(L * del_rate)

    # generate deletions first since they have spacing requirements; shuffling done inside function
    # DO NOT shuffle available positions before deletions - they need to be in order for updating
    position_weights = generate_edge_tapered_weights(L, rng)

    if n_del > 0:
        del_positions, position_weights = generate_mutation_positions(
            "deletion", n_del, del_size_min, del_size_max, position_weights, rng
        )
    else:
        del_positions = []
    if n_ins > 0:
        ins_positions, position_weights = generate_mutation_positions(
            "insertion", n_ins, ins_size_min, ins_size_max, position_weights, rng
        )
    else:
        ins_positions = []
    if n_sub > 0:
        sub_positions, position_weights = generate_mutation_positions(
            "snp", n_sub, 1, 2, position_weights, rng
        )
    else:
        sub_positions = []
    all_mutations = ins_positions + del_positions + sub_positions
    all_mutations.sort()

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
            new_seq.extend(seq[current_index : pos + 1])
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
            reverse_mapper.extend([-1] * length)

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

    # available positions are those that were not mutated (i.e., weights > 0)
    available_positions = np.nonzero(position_weights)[0].tolist()

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
        tuple(dict, dict, dict): sequences after mutations applied (chrom -> seq), mapping from
        original indices to new indices for each chrom (chrom -> list of int), remaining non-mutated
        positions for each chrom (chrom -> list of int)
    """

    new_seqs = {}
    mappers = {}
    available_positions = {}
    for chrom, seq in seqs.items():
        print(f"Mutating chromosome {chrom}...", flush=True)
        new_seqs[chrom], mappers[chrom], available_positions[chrom] = (
            mutate_seq_with_indels_and_snps(
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
        help="Number of offspring to simulate from created wild relative",
    )
    parser.add_argument(
        "--rounds",
        type=int,
        default=5,
        help="Number of rounds to apply mutations to offspring. Will output offspring for each round",
    )

    # introgressions (counts and size range)
    parser.add_argument(
        "--num-introgressions",
        type=int,
        default=48,
        help="Number of introgressions per generation",
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
        default=3e-3,
        help="Wild relative per-base substitution rate",
    )
    parser.add_argument(
        "--rel-ins-rate",
        type=float,
        default=3e-3,
        help="Wild relative per-base insertion rate",
    )
    parser.add_argument(
        "--rel-del-rate",
        type=float,
        default=3e-3,
        help="Wild relative per-base deletion rate",
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
        default=3e-3,
        help="Offspring per-base substitution rate",
    )
    parser.add_argument(
        "--mut-ins-rate",
        type=float,
        default=3e-3,
        help="Offspring per-base insertion rate",
    )
    parser.add_argument(
        "--mut-del-rate",
        type=float,
        default=3e-3,
        help="Offspring per-base deletion rate",
    )
    parser.add_argument(
        "--mut-rate-start",
        type=float,
        default=3e-4,
        help="Starting mutation rate for offspring (will increase linearly to mut_rate over rounds)",
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

    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility")

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

    # write mapper and available positions for debugging
    # with open(output_folder / f"{reference.stem}_wildrelative.mapper", "w") as mapf:
    #     for chrom in ref_seqs:
    #         mapper = reverse_mappers[chrom]
    #         mapf.write(f">{chrom}_mapper\n")
    #         mapf.write(" ".join(map(str, mapper)) + "\n")
    #         mapf.write(f">{chrom}_available_positions\n")
    #         mapf.write(" ".join(map(str, available_positions[chrom])) + "\n")

    # apply introgressions back to reference to create offspring
    print("Applying introgressions to create generation 0 offspring...", flush=True)
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
    print("Writing mutation rate 0, generation 0 introgressed offspring FASTA...", flush=True)
    write_fasta(offspring_seqs, output_folder / f"{reference.stem}_0_offspring.fasta")
    write_bed(introgressions, output_folder / f"{reference.stem}_0_introgressions.bed")

    parent_seqs = offspring_seqs
    sub_rates = np.linspace(args.mut_rate_start, args.mut_sub_rate, args.rounds)
    ins_rates = np.linspace(args.mut_rate_start, args.mut_ins_rate, args.rounds)
    del_rates = np.linspace(args.mut_rate_start, args.mut_del_rate, args.rounds)

    # get intro information for offspring
    introgression_chromosomes = [introgression.split("\t")[0] for introgression in introgressions]
    introgression_starts = [int(introgression.split("\t")[1]) for introgression in introgressions]
    introgression_ends = [int(introgression.split("\t")[2]) for introgression in introgressions]

    for i in range(args.rounds):
        print(
            f"Simulating mutation rate {(sub_rates[i] + ins_rates[i] + del_rates[i]):.2e}, generation {i+1} introgressed offspring...",
            flush=True,
        )

        # apply mutations to offspring from previous generation
        offspring_seqs, reverse_mappers, _ = apply_genome_wide_mutations(
            parent_seqs,
            sub_rate=sub_rates[i],
            ins_rate=ins_rates[i],
            del_rate=del_rates[i],
            ins_size_min=args.mut_ins_size_min,
            ins_size_max=args.mut_ins_size_max,
            del_size_min=args.mut_del_size_min,
            del_size_max=args.mut_del_size_max,
            rng=rng,
        )

        # figure out and save new introgression positions in offspring genome using reverse mappers
        # assumes introgressions are in same order and non-overlapping
        new_introgressions = []
        for j in range(len(introgressions)):
            # get introgression info for this introgression
            intro_chrom = introgression_chromosomes[j]
            intro_start = introgression_starts[j]
            intro_end = introgression_ends[j]

            # get reverse mapper for this chromosome
            reverse_mapper = reverse_mappers[intro_chrom]

            # get new start and end positions in offspring genome
            new_start = reverse_mapper[intro_start]
            # try to correct start/end positions if they map to -1 (deleted region)
            while new_start == -1:
                intro_start += 1
                new_start = reverse_mapper[intro_start]

            # get new end position in offspring genome
            new_end = reverse_mapper[intro_end]
            while new_end == -1:
                intro_end -= 1
                new_end = reverse_mapper[intro_end]

            # save new introgression info
            new_introgressions.append(f"{intro_chrom}\t{new_start}\t{new_end}\tintrogression")

        # save per-gen FASTA
        write_fasta(offspring_seqs, output_folder / f"{reference.stem}_{i+1}_offspring.fasta")
        write_bed(new_introgressions, output_folder / f"{reference.stem}_{i+1}_introgressions.bed")

    print("Simulation finished.")
    return


if __name__ == "__main__":
    main()

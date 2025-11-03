import numpy as np
import panagram.panagram.simulate_introgressions as sim

TEST_REF_FASTA = "./test/Athal_chr4_only.fasta"


def test_generate_positions():
    # set up random number generator
    rng = np.random.default_rng(42)

    # read in test data
    ref_seqs = sim.parse_fasta(TEST_REF_FASTA)
    ref_seq = ref_seqs["chr4"]

    L = len(ref_seq)
    sub_rate=0
    ins_rate=3.3e-3
    del_rate=3.3e-3
    del_size_min=1
    del_size_max=500


    # mutation counts
    n_del = int(L * del_rate)

    available_positions = list(range(L))
    if n_del > 0:
        del_positions, available_positions = sim.generate_mutation_positions(
            "deletion", n_del, del_size_min, del_size_max, available_positions, rng
        )

    # check that deletions do not overlap
    occupied = np.zeros(L, dtype=bool)
    for pos, length, _ in del_positions:
        assert not np.any(occupied[pos:pos + length])
        occupied[pos:pos + length] = True

    # check that available positions do not include deleted positions
    for pos in available_positions:
        assert not occupied[pos]
    return


def test_mapper():
    # set up random number generator
    rng = np.random.default_rng(42)

    # read in test data
    ref_seqs = sim.parse_fasta(TEST_REF_FASTA)

    # run insertions/deletions with no SNPs to generate wild relative
    rel_seqs, mappers, available_positions = sim.apply_genome_wide_mutations(
        ref_seqs,
        sub_rate=0,
        ins_rate=3.3e-3,
        del_rate=3.3e-3,
        ins_size_min=1,
        ins_size_max=1000,
        del_size_min=1,
        del_size_max=500,
        rng=rng,
    )

    # check that for each position in ref, the mapper returns the correct base in the wild relative
    # or -1 if the position was deleted
    for chrom in ref_seqs:
        ref_seq = ref_seqs[chrom]
        rel_seq = rel_seqs[chrom]
        # print(len(ref_seq), len(rel_seq))
        mapper = mappers[chrom]
        replicated_seq = []
        for ref_pos in range(len(ref_seq)):
            # print(ref_pos)
            rel_pos = mapper[ref_pos]
            if rel_pos == -1:
                # position was deleted, check that it is not in the available positions
                # assert ref_pos not in available_positions[chrom]
                replicated_seq.append(ref_seq[ref_pos])
                # print(f"ref_pos: {ref_pos} was deleted, base: {ref_seq[ref_pos]}")
            else:
                # position exists in wild relative, check that the bases match
                # print(f"ref_pos: {ref_pos}, rel_pos: {rel_pos}", ref_seq[ref_pos], rel_seq[rel_pos])
                replicated_seq.append(rel_seq[rel_pos])
                assert ref_seq[ref_pos] == rel_seq[rel_pos]
    assert "".join(replicated_seq) == ref_seq
    return


def test_introgression():
    # set up random number generator
    rng = np.random.default_rng(42)

    # read in test data
    ref_seqs = sim.parse_fasta(TEST_REF_FASTA)

    # run insertions/deletions with no SNPs to generate wild relative
    rel_seqs, mappers, available_positions = sim.apply_genome_wide_mutations(
        ref_seqs,
        sub_rate=0,
        ins_rate=3.3e-3,
        del_rate=3.3e-3,
        ins_size_min=1,
        ins_size_max=1000,
        del_size_min=1,
        del_size_max=500,
        rng=rng,
    )

    # apply introgressions
    introgressed_seqs, bed_entries = sim.apply_genome_wide_introgressions(
        ref_seqs,
        rel_seqs,
        mappers,
        available_positions,
        n_introgressions=2,
        size_min=2000000,
        size_max=4000000,
        rng=rng,
    )

    # get positions and length of introgressions from bed entries
    introgressions = []
    for entry in bed_entries:
        chrom, start, end, length = entry.split("\t")
        length = length.split(" ")[0]  # get length before " bp introgression"
        introgressions.append((int(start), int(end), int(length)))

    # check that introgressed sequences match wild relative at introgressed positions
    for chrom in ref_seqs:
        ref_seq = ref_seqs[chrom]
        introgressed_seq = introgressed_seqs[chrom]
        rel_seqs = rel_seqs[chrom]

        # get current introgression position

        # loop through ref, if position is introgressed, get coorresponding position in wild relative and check that bases match introgressed_seq
        # if not introgressed, check that bases match ref_seq
        # update the position in ref accordingly
        mapper = mappers[chrom]
        current_introgression_index = 0
        ref_pos = 0
        offspring_pos = 0

        try:
            while ref_pos < len(ref_seq):
                current_introgression = introgressions[current_introgression_index] if current_introgression_index < len(introgressions) else None
                current_introgression_start = current_introgression[0] if current_introgression else None
                if not current_introgression or ref_pos < current_introgression_start:
                    # position is not introgressed, should match reference
                    assert introgressed_seq[offspring_pos] == ref_seq[ref_pos]
                    ref_pos += 1
                    offspring_pos += 1
                else:
                    # position is introgressed, should match wild relative
                    introgression_length = current_introgression[2]
                    rel_pos = mapper[ref_pos]
                    for i in range(introgression_length):
                        assert introgressed_seq[offspring_pos] == rel_seqs[rel_pos]
                        rel_pos += 1
                        offspring_pos += 1
                    ref_pos = current_introgression[1] + 1 # move ref_pos to end of introgression + 1 since end is inclusive
                    current_introgression_index += 1
        except Exception as e:
            # print a couple of positions for debugging
            print(f"ref_seq: {ref_seq[ref_pos-5:ref_pos+5]}")
            print(f"introgressed_seq: {introgressed_seq[offspring_pos-5:offspring_pos+5]}")
            print(f"rel_seq: {rel_seqs[rel_pos-5:rel_pos+5]}")
            print(f"ref_pos: {ref_pos}, offspring_pos: {offspring_pos}")
            print(f"Error: {e}")
            raise

    return

# Legacy tests
def old_test_apply_mutations():
    # --- Use your function exactly as is, just override the behavior ---
    rng = np.random.default_rng(1)
    seq = "GATTACA"

    mut_seq, reverse_map, remaining_positions = sim.mutate_seq_with_indels_and_snps(
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


def old_test_apply_introgressions():
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

    new_seqs, bed = sim.apply_genome_wide_introgressions(
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


def old_test_apply_introgressions_with_insertion():
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
    new_seqs, bed = sim.apply_genome_wide_introgressions(
        ref_seqs=ref_seqs,
        rel_seqs=rel_seqs,
        reverse_mappers=reverse_mappers,
        available_positions=available_positions,
        n_introgressions=3,
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
    return


if __name__ == "__main__":
    test_generate_positions()
    test_mapper()
    test_introgression()
    print("All tests passed.")

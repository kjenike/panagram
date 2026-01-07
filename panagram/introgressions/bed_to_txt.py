from pathlib import Path
import argparse
from panagram.index import Index
from postprocess_introgressions import (
    read_bed_file,
    bed_to_bins,
)


def bed_to_text():
    """Convert BED file to text format for each chromosome. Allows usage of ground truth
    introgressions from the simulator in score_introgressions.py.
    """

    parser = argparse.ArgumentParser(
        description="Convert BED file to text format for each chromosome."
    )
    parser.add_argument(
        "--gt_bed_file",
        type=Path,
        help="Path to ground truth introgressions BED file.",
        required=True,
    )
    parser.add_argument(
        "--index_dir",
        type=Path,
        help="Path to Panagram index directory.",
        required=True,
    )
    parser.add_argument(
        "--ref",
        type=str,
        help="Name of the reference genome.",
        required=True,
    )
    parser.add_argument(
        "--wild_type",
        type=str,
        help="Name of the wild type genome.",
        required=True,
    )
    parser.add_argument(
        "--wild_type_group",
        type=str,
        help="Name of the wild type group to use in output text files.",
        required=True,
    )
    parser.add_argument(
        "--bin_size",
        type=int,
        default=1_000_000,
        help="Size of bins to divide chromosomes into.",
    )

    args = parser.parse_args()
    gt_bed_file = args.gt_bed_file
    index_dir = args.index_dir
    ref = args.ref
    wild_type = args.wild_type
    wild_type_group = args.wild_type_group
    bin_size = args.bin_size

    bed_df = read_bed_file(gt_bed_file)
    if bed_df.empty:
        print("No introgressions found in ground truth bed file.")
        return

    # get index and chromosome lengths
    index = Index(index_dir)
    reference_genome = index.genomes[ref]
    offspring_genomes = [
        genome for name, genome in index.genomes.items() if name not in [ref, wild_type]
    ]

    for chr in bed_df["Chromosome"].unique():
        chr_length = reference_genome.sizes[chr]
        bed_chr_df = bed_df[bed_df["Chromosome"] == chr]
        bins_df = bed_to_bins(bed_chr_df, bin_size, chr_length)

        # duplicate introgression column for each offspring genome
        for genome in offspring_genomes:
            bins_df[genome.name] = bins_df["introgression"]

        # remove introgression column and format
        bins_df = bins_df.drop(columns=["introgression"]).T
        bins_df.index.name = "Sample"

        # save to text file
        output_text_file = gt_bed_file.parent / f"{chr}_{wild_type_group}.txt"
        bins_df.to_csv(output_text_file, sep="\t")

    return


if __name__ == "__main__":
    bed_to_text()

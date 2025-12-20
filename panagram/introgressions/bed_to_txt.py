from pathlib import Path
from panagram.index import Index
from postprocess_introgressions import (
    read_bed_file,
    bed_to_bins,
)

GT_BED_FILE = Path(
    "/home/nbrown62/data_mschatz1/nbrown62/panagram_data/arabidopsis_simulated_v3/introgressions_gt/Athal_edited_0_introgressions.bed"
)
INDEX_DIR = Path("/home/nbrown62/data_mschatz1/nbrown62/panagram_data/arabidopsis_simulated_v3")
REF = "Reference"
WILD_TYPE = "WildRelative"
WILD_TYPE_GROUP = "WT"
BIN_SIZE = 1000000


def main():
    bed_df = read_bed_file(GT_BED_FILE)
    if bed_df.empty:
        print("No introgressions found in ground truth bed file.")
        return

    # get index and chromosome lengths
    index = Index(INDEX_DIR)
    reference_genome = index.genomes[REF]
    offspring_genomes = [
        genome for name, genome in index.genomes.items() if name not in [REF, WILD_TYPE]
    ]

    for chr in bed_df["Chromosome"].unique():
        chr_length = reference_genome.sizes[chr]
        bed_chr_df = bed_df[bed_df["Chromosome"] == chr]
        bins_df = bed_to_bins(bed_chr_df, BIN_SIZE, chr_length)
        # duplicate introgression column for each offspring genome
        for genome in offspring_genomes:
            bins_df[genome.name] = bins_df["introgression"]
        # remove introgression column and format
        bins_df = bins_df.drop(columns=["introgression"]).T
        bins_df.index.name = "Sample"

        # save to text file
        output_text_file = GT_BED_FILE.parent / f"{chr}_{WILD_TYPE_GROUP}.txt"
        bins_df.to_csv(output_text_file, sep="\t")

    return


if __name__ == "__main__":
    main()

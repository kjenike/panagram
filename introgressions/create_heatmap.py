import argparse
from pathlib import Path
import pandas as pd
from panagram.index import Index
import call_introgressions as call


def panagram_heatmap_general(
    index_dir,
    anchor,
    groups=None,
    bitmap_step=100,
    bin_size=1000000,
    omit_fixed_kmers=True,
    gnm=None,
):
    """Given a panagram index and an anchor genome, creates heatmaps showing kmer
    similarity across the anchor genome. Each row is a different genome in the index, and each
    column is a different bin across the anchor genome. The color of each cell corresponds to the
    number of shared kmers between that genome and the anchor genome in that bin.

    Args:
        index_dir (str or Path): directory containing the panagram index
        anchor (str): genome name to anchor on for visualization
        groups (str or Path, optional): path to the groups file (default: None)
        bitmap_step (int, optional): Step size for bitmap query (default: 100)
        bin_size (int, optional): Size of bins for visualization (default: 1000000)
        omit_fixed_kmers (bool, optional): Whether to omit fixed kmers (default: True)
        gnm (float, optional): normalize genomes to this mean, or -1 to auto-calc (default: None)
    """

    output_dir = Path(index_dir) / "panagram_visuals"
    output_dir.mkdir(parents=True, exist_ok=True)
    index = Index(index_dir)
    genome = index.genomes[anchor]
    chrs = genome.sizes.keys()
    if groups is not None:
        groups = pd.read_csv(
            groups,
            sep="\t",
            index_col=0,
        )
    if gnm is not None:
        genome_similarities = call.get_genome_similarities(
            genome,
            bitmap_step,
            bin_size,
            omit_fixed_kmers,
            False,
            None,
            None,
            3,
        )

    for chr_name in chrs:
        print(f"Processing {anchor} - {chr_name}...")

        # get an entire chr's bitmap
        chr_size = genome.sizes[chr_name]
        chr_bitmap = genome.query(chr_name, 0, chr_size, step=bitmap_step)
        chr_name = str(chr_name)

        binned_bitmap = call.bitmap_to_bins(chr_bitmap, bin_size, omit_fixed_kmers=omit_fixed_kmers)
        if gnm is not None:
            binned_bitmap = call.preprocess_binned_bitmap(
                binned_bitmap,
                genome_similarities,
                gnm,
                None,
                None,
                False,
            )

        call.visualize(
            binned_bitmap,
            output_dir / f"{anchor}_{chr_name}_heatmap.svg",
            inverse=True,
            groups=groups,
            title=f"{anchor.capitalize()} - {chr_name.capitalize()}",
        )
    return


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--index-dir", help="Directory containing the panagram index", required=True
    )
    parser.add_argument(
        "--anchor", help="Genome name to anchor on for visualization", required=True
    )
    parser.add_argument(
        "--groups", help="group.tsv file to determine row order; not required", required=False
    )
    parser.add_argument(
        "--bitmap-step", help="Step size for bitmap query (default: 100)", type=int, default=100
    )
    parser.add_argument(
        "--bin-size",
        help="Size of bins for visualization (default: 1000000)",
        type=int,
        default=1000000,
    )
    parser.add_argument(
        "--rmf",
        help="Whether to omit fixed kmers (default: True)",
        action="store_true",
        default=True,
    )
    parser.add_argument(
        "--gnm",
        help="Normalize genomes to this mean, or -1 to auto-calc (default: None)",
        type=float,
        default=None,
    )
    args = parser.parse_args()

    panagram_heatmap_general(
        index_dir=args.index_dir,
        anchor=args.anchor,
        groups=args.groups,
        bitmap_step=args.bitmap_step,
        bin_size=args.bin_size,
        omit_fixed_kmers=args.rmf,
        gnm=args.gnm,
    )

    return


if __name__ == "__main__":
    main()

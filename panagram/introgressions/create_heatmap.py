import argparse
from pathlib import Path
import pandas as pd
from panagram.index import Index
import call_introgressions as call


def panagram_heatmap_general(
    index_dir, anchor, groups, bitmap_step=100, bin_size=1000000
):
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

    for chr_name in chrs:
        # get an entire chr's bitmap
        chr_size = genome.sizes[chr_name]
        chr_bitmap = genome.query(chr_name, 0, chr_size, step=bitmap_step)
        chr_name = str(chr_name)

        pair = call.bitmap_to_bins(chr_bitmap, bin_size, omit_fixed_kmers=True)
        call.visualize(
            pair,
            output_dir / f"{anchor}_{chr_name}_heatmap.svg",
            inverse=True,
            groups=groups,
            title=f"{anchor.capitalize()} - {chr_name.capitalize()}",
        )
    return


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--index-dir", help="Directory containing the panagram index", required=True)
    parser.add_argument("--anchor", help="Genome name to anchor on for visualization", required=True)
    parser.add_argument("--groups", help="group.tsv file to determine row order; not required", required=False)
    parser.add_argument("--bitmap-step", help="Step size for bitmap query (default: 100)", type=int, default=100)
    parser.add_argument("--bin-size", help="Size of bins for visualization (default: 1000000)", type=int, default=1000000)
    args = parser.parse_args()

    panagram_heatmap_general(
        index_dir=args.index_dir,
        anchor=args.anchor,
        groups=args.groups,
        bitmap_step=args.bitmap_step,
        bin_size=args.bin_size,
    )

    return

if __name__ == "__main__":
    main()

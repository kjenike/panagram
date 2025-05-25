from pathlib import Path
import argparse
import pandas as pd
from panagram.index import Index
import plotly.express as px


def visualize(pair, output_file, inverse=False):
    # take a look at what pair looks like after manipulation
    if inverse:
        fig = px.imshow(
            pair,
            color_continuous_scale=px.colors.sequential.Greens[
                ::-1
            ],  # px.colors.sequential.Plasma[::-1],
            x=pair.columns,
            y=pair.index,
            aspect="auto",
            zmin=0,
            # zmax=1,
        )
    else:
        fig = px.imshow(
            pair,
            x=pair.columns,
            y=pair.index,
            aspect="auto",
        )
    fig.update_layout(
        xaxis=dict(
            dtick=2000000,
        ),
    )
    fig.write_image(output_file)
    return


def threshold_introgressions(pair, anchor, comp_group, threshold):
    # get the group the anchor belongs to
    anchor_group = pair.loc[anchor, "group"]

    # define df with group info and kmer sims of anchor's group; drop anchor's self-similarity
    pair_anchor_group = (
        pair[pair["group"] == anchor_group].drop(columns=["group"]).drop(anchor, axis=0)
    )

    # define df with group info and kmer sim for comparison group
    pair_comp_group = pair[pair["group"] == comp_group].drop(columns=["group"])

    # throw out all acessions with large global diffs w/ the anchor
    # TODO: maybe make this a parameter
    if len(pair_anchor_group) > 1:
        for acession, row in pair_anchor_group.iterrows():
            # TODO: make this at least 0.85 when omit fixed is not active; its too low otherwise
            if row.mean() < 0.75:
                pair_anchor_group = pair_anchor_group[pair_anchor_group.index != acession]

        if pair_anchor_group.empty:
            print(f"Warning: accession {anchor} is not similar to any of its group members...")
            # run anyway only if comparing to REF
            if comp_group != "REF":
                print("Skipping...")
                return None

    # get mean kmer similarities per window for each group
    group_sims = pair_anchor_group.mean(axis=0).to_frame(name="anchor_sim")
    group_sims["comp_sim"] = pair_comp_group.mean(axis=0)

    # 4 ways to define an introgression:
    # 1. If anchor is in comp_group: introgression is where anchor has low similarity to its group
    # 2. If "REF" is the comp_group: introgression is where anchor has low similarity to REF group
    # 3. If anchor isn't in comp_group: introgression is where anchor is more similar to comp_group
    # 4. Average of all introgression calls for comp_group when comp_group is a list
    # than its own group
    if anchor_group == comp_group:
        group_sims["introgression"] = (group_sims.anchor_sim < threshold).astype(int)
    elif comp_group == "REF":
        # comp_group should only consist of 1 REF
        group_sims["introgression"] = (group_sims.comp_sim < threshold).astype(int)
    # TODO: could refine so that must be at least 10% more similar to comp_sim
    else:
        group_sims["introgression"] = (group_sims.comp_sim >= group_sims.anchor_sim).astype(int)
    return group_sims


def bins_to_bed(bins_df, bin_size, chr_name, comp_group):
    # takes df organized into rows of bins w/ introgression column and converts to bed file format
    # find start/end coordinates
    introgressions = bins_df[bins_df.introgression > 0].copy()
    introgressions["start"] = introgressions.index
    introgressions["end"] = introgressions["start"] + bin_size

    # call adjacent introgressions as the same
    # check if end of the prev col is the same as the start of the current
    introgressions["groups"] = (introgressions.start - introgressions.end.shift(1)).fillna(0)
    # use cumsum to assign adjacent intros to the same group
    introgressions["groups"] = introgressions["groups"].cumsum()
    # get the number of bins in each group and use to calculate new start/end
    bins_in_groups = introgressions.groupby("groups").count()["end"]
    introgressions = introgressions.drop_duplicates(subset="groups", keep="first")
    introgressions = introgressions.drop(columns="end").merge(bins_in_groups, on="groups")
    introgressions["end"] = introgressions.start + (introgressions.end * bin_size)

    # put into bed file format (chr, start, end, name)
    introgressions["chr"] = chr_name
    introgressions["name"] = f"{comp_group}_intro"
    introgressions = introgressions[["chr", "start", "end", "name"]]
    return introgressions


def bitmap_to_bins(bitmap, binlen, omit_fixed_kmers=False):
    # modded version of the same function found in the Index class

    # change index from chr position to the number of the bin the position falls in
    df = bitmap.set_index(bitmap.index // binlen)

    # remove fixed kmers shared by all members of the pangenome (i.e., rows that are all 1)
    if omit_fixed_kmers:
        df = df.loc[~(df == 1).all(axis=1)]

    # get the sum of how many kmers equal 1 in each bin
    paircount_bins = df.groupby(level=0).sum()

    # fix the index's bin number to be the bin's chr position and transpose
    paircount_bins = paircount_bins.set_index(paircount_bins.index * binlen).T

    # divide each kmer sum in each column by the max # of kmers in the bin (some bins are uneven)
    paircount_bins = paircount_bins.div(paircount_bins.max(axis=0), axis=1)
    return paircount_bins


def run_introgression_finder(
    anchor,
    genome,
    chr_name,
    groups,
    comp_groups,
    threshold,
    bitmap_step,
    bin_size,
    omit_fixed_kmers,
    render_vis,
    output_dir,
):

    # Step 1 - choose an anchor and re-create pairwise correlation matrix for it
    # get an entire chr's bitmap
    chr_size = genome.sizes[chr_name]
    chr_bitmap = genome.query(chr_name, 0, chr_size, step=bitmap_step)

    # get correlation matrix
    pair = bitmap_to_bins(chr_bitmap, bin_size, omit_fixed_kmers)

    # get the kmer similarities for the anchor's group and the comparison group
    pair = pair.merge(groups, left_index=True, right_index=True, how="left")

    # Step 2 - identify potential introgressions by thresholding kmer similarities
    merged_sims = None
    for comp_group in comp_groups:
        group_sims = threshold_introgressions(pair, anchor, comp_group, threshold)
        if group_sims is None:
            continue
        elif merged_sims is None:
            merged_sims = group_sims
        else:
            merged_sims += group_sims

        # Step 3 - visualization
        # show user introgressions labeled on the original heatmap from panagram
        # invert values for figure so that introgressions = 0
        pair.loc["Intro. Score"] = (~(group_sims["introgression"].astype(bool))).astype(int)
        if render_vis:
            output_vis = output_dir / f"{anchor}_{chr_name}_{comp_group}_heatmap.png"
            visualize(pair.drop(columns=["group"]), output_vis, inverse=True)

        # Step 4 - save introgressions to bed file
        output_bed = output_dir / f"{anchor}_{chr_name}_{comp_group}.bed"
        introgressions = bins_to_bed(group_sims, bin_size, chr_name, comp_group)
        introgressions.to_csv(output_bed, header=False, index=False, sep="\t")

    # Step 5 - repeat analysis one more time with merged introgressions from all comp_groups
    if merged_sims is not None:
        # divide by the max to get between 0 and 1, do 1 - x to invert
        pair.loc["Intro. Score"] = 1 - (
            merged_sims["introgression"] / merged_sims["introgression"].max()
        )

        if render_vis:
            output_vis = output_dir / f"{anchor}_{chr_name}_merged_heatmap.png"
            visualize(pair.drop(columns=["group"]), output_vis, inverse=True)

        # save introgressions to bed file
        output_bed = output_dir / f"{anchor}_{chr_name}_merged.bed"
        introgressions = bins_to_bed(merged_sims, bin_size, chr_name, "merged")
        introgressions.to_csv(output_bed, header=False, index=False, sep="\t")
    return


def main():
    parser = argparse.ArgumentParser(description="Introgression highlighter tool.")
    parser.add_argument("--stp", type=int, default=100, help="bitmap kmer step size")
    parser.add_argument("--bin", type=int, default=1000000, help="size of bitmap bin in bases")
    parser.add_argument(
        "--thr", type=float, default=0.75, help="introgression kmer similarity lower threshold"
    )
    parser.add_argument("--rmf", action="store_true", help="remove fixed kmers from bitmap")
    parser.add_argument("--vis", action="store_true", help="save pngs of visualized results")
    parser.add_argument("--isc", action="store_true", help="use anchor's group to self compare")
    parser.add_argument(
        "-a", nargs="+", help="name of anchor(s) to mark introgressions for (default: all)"
    )
    parser.add_argument(
        "-g",
        nargs="+",
        help="if -a is not defined, groups in the tsv to mark introgressions for (default: none)",
    )
    parser.add_argument(
        "-c", nargs="+", help="chromosome(s) to mark introgressions for (default: all)"
    )
    parser.add_argument(
        "-p", nargs="+", help="group(s) to compare against anchor(s)", required=True
    )
    parser.add_argument("--idx", type=str, help="path to Panagram index folder", required=True)
    parser.add_argument("--tsv", type=str, help="path to acession group TSV file", required=True)
    parser.add_argument("--out", type=str, help="path to folder to save all outputs", required=True)
    args = parser.parse_args()

    bitmap_step = args.stp
    bin_size = args.bin
    threshold = args.thr
    omit_fixed_kmers = args.rmf
    render_vis = args.vis

    group_tsv = Path(args.tsv)
    if not group_tsv.is_file():
        raise ValueError("TSV file not found. Check --tsv path.")

    index_dir = Path(args.idx)
    if not index_dir.is_dir():
        raise ValueError("Index directory not found. Check --idx path.")
    index = Index(index_dir)
    groups = pd.read_csv(group_tsv, sep="\t", index_col=0)

    # determine if we have a single anchor or multiple to run
    anchors = args.a
    if anchors is None:
        if args.g is None:
            raise ValueError("No anchor selected. Use either -a or -g to specify anchors.")
        # get multiple anchors
        anchors = list(groups[groups.group.isin(args.g)].index)

    comp_groups = args.p
    # ensure every element is only in there once
    comp_groups = list(set(comp_groups))

    output_dir = Path(args.out)
    output_dir.mkdir(parents=True, exist_ok=True)

    for anchor in anchors:
        loop_comp_groups = comp_groups
        anchor_group = groups.loc[anchor, "group"]

        # add or remove self-compare based on user input
        if args.isc:
            if anchor_group not in loop_comp_groups:
                loop_comp_groups.append(anchor_group)
        else:
            if anchor_group in loop_comp_groups:
                loop_comp_groups.remove(anchor_group)

        genome = index.genomes[anchor]
        chromosomes = args.c
        if chromosomes is None:
            chromosomes = genome.sizes.keys()

        for chr_name in chromosomes:
            print("Now running introgression analysis for", anchor, chr_name)
            run_introgression_finder(
                anchor,
                genome,
                chr_name,
                groups,
                comp_groups,
                threshold,
                bitmap_step,
                bin_size,
                omit_fixed_kmers,
                render_vis,
                output_dir,
            )
    return


if __name__ == "__main__":
    main()

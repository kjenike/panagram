from pathlib import Path
import pandas as pd
from panagram.index import Index
import plotly.express as px
import subprocess


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


def fill_gaps(row, rounds=1):
    row = row.values
    for j in range(rounds):
        # Find gaps of 0s surrounded by 1s
        for i in range(1, len(row) - 1):
            if row[i] == 0 and (row[i - 1] >= 1 and row[i + 1] >= 1):
                row[i] = 1
    return row


def threshold_introgressions(pair, anchor, comp_group):
    # get the group the anchor belongs to
    anchor_group = pair.loc[anchor, "group"]

    # make sure anchor's self-similarity gets dropped
    pair_anchor_group = (
        pair[pair["group"] == anchor_group].drop(columns=["group"]).drop(anchor, axis=0)
    )
    pair_comp_group = pair[pair["group"] == comp_group].drop(columns=["group"])

    # get mean kmer similarities per window for each group
    group_sims = pair_anchor_group.mean(axis=0).to_frame(name="anchor_sim")
    group_sims["comp_sim"] = pair_comp_group.mean(axis=0)

    # 4 ways to define an introgression:
    # 1. If anchor is in comp_group: introgression is where anchor has low similarity to its group
    # 2. If "REF" is the comp_group: introgression is where anchor has low similarity to REF group
    # 3. If anchor isn't in comp_group: introgression is where anchor is more similar to comp_group
    # 4. Average of all introgression calls for comp_group when comp_group is a list
    # than its own group
    # Postprocessing: gap fill between nearby introgressions

    if anchor_group == comp_group:
        # print(f"{anchor} is part of group {comp_group}. Comparing using thresholding...")
        group_sims["introgression"] = fill_gaps((group_sims.comp_sim < 0.9).astype(int))
    elif comp_group == "REF":
        group_sims["introgression"] = fill_gaps((group_sims.comp_sim < 0.9).astype(int))
    else:
        group_sims["introgression"] = fill_gaps(
            (group_sims.comp_sim >= group_sims.anchor_sim).astype(int)
        )

    return group_sims


def threshold_to_bed(group_sims, bin_size, chr_name, comp_group, output_file):
    # find start/end coordinates
    introgressions = group_sims[group_sims.introgression > 0].copy()
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

    # write out a bed file of all found introgressions (chr, start, end)
    introgressions["chr"] = chr_name
    introgressions["name"] = f"{comp_group}_intro"
    introgressions = introgressions[["chr", "start", "end", "name"]]
    introgressions.to_csv(
        output_file,
        header=False,
        index=False,
        sep="\t",
    )
    return


def run_introgression_finder(
    index,
    anchor,
    chr_name,
    group_tsv,
    comp_groups,
    bitmap_step,
    bin_size,
    output_dir,
):
    # only generate missing bed files
    # if (output_dir / f"{anchor}_{chr_name}_{comp_group}.bed").exists():
    #     return

    # Step 1 - choose an anchor and re-create pairwise correlation matrix for it
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    genome = index.genomes[anchor]
    groups = pd.read_csv(group_tsv, sep="\t", index_col=0)

    # get an entire chr's bitmap
    chr_size = genome.sizes[chr_name]
    chr_bitmap = genome.query(chr_name, 0, chr_size, step=bitmap_step)

    # get correlation matrix
    _, pair = index.bitmap_to_bins(chr_bitmap, bin_size)

    # get the kmer similarities for the anchor's group and the comparison group
    pair = pair.merge(groups, left_index=True, right_index=True, how="left")

    merged_sims = None
    for comp_group in comp_groups:
        group_sims = threshold_introgressions(pair, anchor, comp_group)
        if merged_sims is None:
            merged_sims = group_sims
        else:
            merged_sims += group_sims

        # show user introgressions labeled on the original heatmap from panagram
        # invert values for figure so that introgressions = 0
        pair.loc["Intro. Score"] = (~(group_sims["introgression"].astype(bool))).astype(int)
        output_vis = output_dir / f"{anchor}_{chr_name}_{comp_group}_heatmap.png"
        visualize(pair.drop(columns=["group"]), output_vis, inverse=True)

        # save introgressions to bed file
        output_bed = output_dir / f"{anchor}_{chr_name}_{comp_group}.bed"
        threshold_to_bed(group_sims, bin_size, chr_name, comp_group, output_bed)

    # repeat analysis one more time with merged introgressions from all comp_groups
    # pair.loc["Intro. Score"] = (~(merged_sims["introgression"].astype(bool))).astype(int)

    # divide by the max to get between 0 and 1, do 1 - x to invert
    pair.loc["Intro. Score"] = 1 - (
        merged_sims["introgression"] / merged_sims["introgression"].max()
    )

    output_vis = output_dir / f"{anchor}_{chr_name}_merged_heatmap.png"
    visualize(pair.drop(columns=["group"]), output_vis, inverse=True)

    # save introgressions to bed file
    output_bed = output_dir / f"{anchor}_{chr_name}_merged.bed"
    threshold_to_bed(merged_sims, bin_size, chr_name, "merged", output_bed)
    return


def main():
    # TODO: allow for a 2-bin gap
    # USER PARAMS
    bitmap_step = 100
    bin_size = 1000000
    index_dir = "/home/nbrown62/data_mschatz1/nbrown62/panagram_data/tomato_sl4"
    group_tsv = "/home/nbrown62/data_mschatz1/nbrown62/panagram_data/tomato_sl4/group.tsv"
    comp_groups = ["SP", "SLC", "SLL", "REF"]
    output_dir = Path(
        "/home/nbrown62/data_mschatz1/nbrown62/panagram_data/tomato_sl4/introgression_analysis_v2/"
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    index = Index(index_dir)

    # For testing with tomato pangenome
    # for anchor in ["M82"]:
    #     genome = index.genomes[anchor]
    #     print(genome.sizes.keys())
    #     for chr_name in ["chr11", "chr4"]:
    #         print("Now running introgression analysis for", anchor, chr_name)
    #         run_introgression_finder(
    #             index,
    #             anchor,
    #             chr_name,
    #             group_tsv,
    #             comp_groups,
    #             bitmap_step,
    #             bin_size,
    #             output_dir,
    #         )

    for anchor in index.genomes.keys():
        genome = index.genomes[anchor]
        for chr_name in genome.sizes.keys():
            print("Now running introgression analysis for", anchor, chr_name)
            run_introgression_finder(
                index,
                anchor,
                chr_name,
                group_tsv,
                comp_groups,
                bitmap_step,
                bin_size,
                output_dir,
            )
    # NOTE: if anchor not in SLL group, we could skip
    return


if __name__ == "__main__":
    main()

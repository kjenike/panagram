from pathlib import Path
import pandas as pd
from panagram.index import Index
import plotly.express as px


def better_dir(item):
    # don't print hidden functions
    methods = dir(item)
    return [method for method in methods if not method.startswith("_")]


def visualize(pair, output_file, inverse=False):
    # take a look at what pair looks like after manipulation
    # pair[pair >= 1] = 10
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
            zmax=1,
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


def merge_adjacent():
    return


def run_introgression_finder(
    index,
    anchor,
    chr_name,
    group_tsv,
    comp_group,
    bitmap_step,
    bin_size,
    output_dir,
):
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

    # show the original heatmap of kmer similarities that panagram shows
    # deeper green = more kmer dissimilarity
    # visualize(pair, output_dir / f"{anchor}_{chr_name}_original_heatmap.png", inverse=True)

    # get the kmer similarities for the anchor's group and the comparison group
    pair = pair.merge(groups, left_index=True, right_index=True, how='left')
    anchor_group = pair.loc[anchor, "group"] # get the group the anchor belongs to
    # make sure anchor's self-similarity gets dropped
    pair_anchor_group = pair[pair["group"] == anchor_group].drop(columns=["group"]).drop(anchor, axis=0)
    pair_comp_group = pair[pair["group"] == comp_group].drop(columns=["group"])

    # get mean similarities per window for each group
    group_sims = pair_anchor_group.mean(axis=0).to_frame(name="anchor_sim")
    group_sims["comp_sim"] = pair_comp_group.mean(axis=0)
    group_sims["introgression"] = (group_sims.comp_sim >= group_sims.anchor_sim)

    # show user introgressions labeled on the original heatmap from panagram
    pair = pair.drop(columns = ["group"])
    pair.loc["Intro?"] = (~group_sims["introgression"]).astype(int)
    visualize(pair, output_dir / f"{anchor}_{chr_name}_{comp_group}_heatmap.png", inverse=True)

    # find start/end coordinates by merging adjacent introgressions
    introgressions = group_sims[group_sims.introgression > 0].copy()
    introgressions['start'] = introgressions.index
    introgressions['end'] = introgressions['start'] + bin_size

    print(introgressions)

    # write out a bed file of all found introgressions (chr, start, end)

    exit()
    return


def main():
    # USER PARAMS
    bitmap_step = 100
    bin_size = 1000000
    index_dir = "/home/nbrown62/data_mschatz1/nbrown62/panagram_data/tomato_sl4"
    group_tsv = "/home/nbrown62/data_mschatz1/nbrown62/panagram_data/tomato_sl4/group.tsv"
    comp_group = "SP"
    output_dir = Path(
        "/home/nbrown62/data_mschatz1/nbrown62/panagram_data/tomato_sl4/introgression_analysis_v2/"
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    index = Index(index_dir)

    # For testing with tomato pangenome
    for anchor in ["SL4"]:
        genome = index.genomes[anchor]
        print(genome.sizes.keys())
        for chr_name in ["chr11", "chr4"]:
            for comp_group in ["SLC", "SP"]:
                print("Now running introgression analysis for", anchor, chr_name, comp_group)
                run_introgression_finder(
                    index,
                    anchor,
                    chr_name,
                    group_tsv,
                    comp_group,
                    bitmap_step,
                    bin_size,
                    output_dir,
                )

    # for anchor in index.genomes.keys():
    #     genome = index.genomes[anchor]
    #     for chr_name in genome.sizes.keys():
    #         print("Now running introgression analysis for", anchor, chr_name)
    #         run_introgression_finder(
    #             index,
    #             anchor,
    #             chr_name,
    #             bitmap_step,
    #             bin_size,
    #             output_dir,
    #         )
    return


if __name__ == "__main__":
    main()

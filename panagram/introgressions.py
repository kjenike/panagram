from pathlib import Path
import math
import numpy as np
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
            color_continuous_scale=px.colors.sequential.Plasma[::-1],
            x=pair.columns,
            y=pair.index,
        )
    else:
        fig = px.imshow(pair, x=pair.columns, y=pair.index)
    fig.write_image(output_file)
    return


def fill_gaps(row, rounds=2):
    row = row.values
    for j in range(rounds):
        # Find gaps of 0s surrounded by 1s
        for i in range(1, len(row) - 1):
            # TODO: if previous row and following row do not share identities, don't add together
            if row[i] == 0 and (row[i - 1] >= 1 and row[i + 1] >= 1):
                row[i] = (row[i - 1] + row[i + 1]) / 2
    return row


def run_introgression_finder(
    index,
    anchor,
    chr_name,
    bitmap_step,
    max_chr_bins,
    k,
    set_difference_threshold,
    output_dir,
):
    # Step 1 - choose an anchor and re-create pairwise correlation matrix for it
    # kmer size is 20-30ish; kmer at position X starts at position X (not centered at position X)
    # there are multiple positions in a bin; there are <bin size> - k + 1 kmers in a bin
    # default: bin_size = ((end_coord - start_coord) // max_chr_bins) + 1; max_chr_bins = 350; step = 100
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # index = Index(index_dir)
    genome = index.genomes[anchor]
    # # print(index.genomes)
    # # print(genome.sizes)

    # get an entire chr's bitmap
    chr_size = genome.sizes[chr_name]
    chr_bitmap = genome.query(chr_name, 0, chr_size, step=bitmap_step)

    # get correlation matrix
    start_coord = 0
    end_coord = chr_size
    bin_size = ((end_coord - start_coord) // max_chr_bins) + 1
    num_kmers_in_bin = bin_size - k + 1

    print("# Positions in a bin", bin_size)
    print("# Kmers in a bin", num_kmers_in_bin)

    pan, pair = index.bitmap_to_bins(chr_bitmap, bin_size)
    visualize(
        pair, output_dir / f"{anchor}_{chr_name}_original_heatmap.png", inverse=True
    )

    # # sanity check
    # print(len(pair.columns))

    # Step 2 - slide through the matrix one genome at a time and calculate average pairwise correlation
    # introgression - area of yellow in panagram where all other accessions/samples are less similar to anchor
    # assumes other accessions are related variants without the introgression
    # sliding window size based on size of the matrix - 5% of the size? alternatively based on bps
    # threshold maybe around <=0.6; if average within threshold, mark all locations within window as part of introgression

    # TODO: decide btwn using the mean/max of outliers to set threshold
    # Or using q3; could increase the threshold and use redundancy as an additional metric for deciding
    values = pair.values.flatten()
    q1 = np.percentile(values, 25)
    q3 = np.percentile(values, 75)
    iqr = q3 - q1
    lower_bound = q1 - 1.5 * iqr
    outliers = values[values < lower_bound]

    if len(outliers) > 0:
        dissimilarity_threshold = outliers.max()
    else:
        dissimilarity_threshold = 0.6
    print("Chosen dissimilarity threshold", dissimilarity_threshold)

    # Histogram of correlations
    fig = px.histogram(values, title="Histogram of Correlation Values")
    fig.add_vline(
        x=dissimilarity_threshold,
        line_width=2,
        line_dash="dash",
        line_color="red",
        annotation_text="Threshold",
        annotation_position="top left",
    )
    fig.write_image(output_dir / f"{anchor}_{chr_name}_corr_hist.png")

    # NOTE: if there are no outliers, or threshold is very high, there are likely no introgressions
    # test for this scenario

    # flip matrix so we can use pandas rolling operation
    transposed_pair = pair.transpose().drop(columns=[anchor])
    pangenome_names = list(transposed_pair.columns)  # save for reference later
    # print(transposed_pair)

    # convert to binary array by applying dissimilarity threshold
    transposed_pair[transposed_pair > dissimilarity_threshold] = 0
    transposed_pair[transposed_pair != 0] = 1

    # Step 3 - create an introgression score for each position based on total dissimilarity
    # for each genome combine nearby locations as the same introgression - column by column
    transposed_pair = transposed_pair.apply(fill_gaps)

    # For each position, figure out which subset of genomes share an introgression at the position
    transposed_pair["genomes"] = transposed_pair.apply(
        lambda row: {col for col in pangenome_names if row[col] > 0}, axis=1
    )
    # Calculate set difference with the next row and store the length
    transposed_pair["genomes_overlap"] = (
        transposed_pair["genomes"]
        .shift(1)
        .combine(
            transposed_pair["genomes"],
            lambda next_row, current_row: (
                len(current_row & (next_row or set()))
                / len(current_row | (next_row or set()))
                if current_row | (next_row or set())
                else 0
            ),
        )
    )

    # introgression score - counter of number of genomes that share the introgression
    # TODO: normalize by number of genomes in the pangenome?
    # transposed_pair["introgression_score"] = transposed_pair[pangenome_names].sum(axis=1)
    # transposed_pair["introgression_score"] = fill_gaps(transposed_pair["introgression_score"].values, transposed_pair["genomes"].values, set_difference_threshold)

    transposed_pair["introgression_score"] = transposed_pair["genomes"].apply(len)
    visualize(
        transposed_pair[pangenome_names + ["introgression_score"]].transpose(),
        output_dir / f"{anchor}_{chr_name}_introgressions_heatmap.png",
    )
    pair.loc["intro_score"] = 1 - (
        transposed_pair["introgression_score"].values
        / max(transposed_pair["introgression_score"].values)
    )
    visualize(
        pair, output_dir / f"{anchor}_{chr_name}_hybrid_heatmap.png", inverse=True
    )

    # Step 4 - report, sorted by size and then by score
    # transposed_pair['introgression_starts'] = (transposed_pair['introgression_score'] != 0) & (transposed_pair['introgression_score'].shift(1, fill_value=0) == 0)
    transposed_pair["introgression_starts"] = (
        (transposed_pair["genomes_overlap"] > 0)
        & (transposed_pair["genomes_overlap"] < 0.9)
    ) | (
        (transposed_pair["genomes_overlap"] == 0)
        & (transposed_pair["introgression_score"] == 1)
    )

    # Create a group identifier for each set of non-zeros
    transposed_pair["introgression_group"] = transposed_pair[
        "introgression_starts"
    ].cumsum()

    # sum of all introgression scores
    introgression_groupby = transposed_pair[
        transposed_pair["introgression_score"] != 0
    ].groupby("introgression_group")
    total_introgression_scores = introgression_groupby["introgression_score"].mean()
    # end index for finding introgression length in bps
    last_indices = introgression_groupby.tail(1).index.values
    # genomes that introgression belongs to (union of all sets in the group)
    introgression_genomes = (
        introgression_groupby["genomes"].agg(lambda sets: set.union(*sets)).values
    )

    introgressions = transposed_pair[
        transposed_pair.introgression_starts == True
    ].copy()

    # print(transposed_pair[transposed_pair['introgression_score'] != 0][transposed_pair.columns[10:20]])
    introgressions["introgression_score"] = total_introgression_scores.values
    introgressions["introgression_end"] = last_indices
    introgressions["introgression_end"] = (
        introgressions["introgression_end"] + bin_size
    )  # account for the fact that the index is the start and not the end of a bin
    introgressions = (
        introgressions[["introgression_end", "introgression_score"]]
        .reset_index()
        .rename(columns={"index": "introgression_start"})
    )
    introgressions["introgression_length"] = (
        introgressions["introgression_end"] - introgressions["introgression_start"]
    )
    introgressions["introgression_genomes"] = introgression_genomes
    print(introgressions.sort_values(by=["introgression_start"], ascending=True))
    introgressions.sort_values(by=["introgression_start"], ascending=True).to_csv(
        output_dir / f"{anchor}_{chr_name}_introgressions.csv", index=False
    )

    return


# USER PARAMS
index_dir = "/home/nbrown62/data_mschatz1/nbrown62/panagram_data/tomato"
# anchor = "BGV006775_MAS2" #"SL5"
# chr_name = "BGV006775_MAS2.0ch11"
bitmap_step = 100
max_chr_bins = 350
size_threshold = 3000000  # NOTE: unused, minimum size in bps of the introgression
k = 31  # TODO: k should be defined somewhere else; don't need from the user
output_dir = "/home/nbrown62/data_mschatz1/nbrown62/panagram_data/tomato/introgression_analysis_v1/"

index = Index(index_dir)
set_difference_threshold = 2  # NOTE: could also make this proportionate to pangenome size int(len(index.genomes) / 10)

# For testing with tomato pangenome
for anchor in ["SL5"]:
    genome = index.genomes[anchor]
    print(genome.sizes.keys())
    for chr_name in [11]:
        print(anchor, chr_name)
        run_introgression_finder(
            index,
            anchor,
            chr_name,
            bitmap_step,
            max_chr_bins,
            k,
            set_difference_threshold,
            output_dir,
        )
        break
    break

# for anchor in index.genomes.keys():
#     genome = index.genomes[anchor]
#     for chr_name in genome.sizes.keys():
#         print("Now running introgression analysis for", anchor, chr_name)
#         run_introgression_finder(index, anchor, chr_name, bitmap_step, max_chr_bins, k, output_dir)

# NOTE: Could look at underlying sequence in found introgressions; compare/cluster? align with annotations? This would be computationally expensive.

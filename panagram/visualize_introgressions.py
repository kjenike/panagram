from pathlib import Path
import numpy as np
import pandas as pd
import plotly.express as px
from postprocess_introgressions import threshold_introgressions


def self_dotplot(seq):
    import matplotlib.pyplot as plt

    # look at a 10kb piece of a sequence
    middle = int(len(seq) / 2)
    seq = seq[middle - 5000 : middle + 5000]
    length = len(seq)

    matrix = np.zeros((length, length))

    for i in range(length):
        for j in range(length):
            if seq[i] == seq[j]:
                matrix[i, j] = 1

    plt.imshow(matrix, cmap="Greys", interpolation="none")
    plt.savefig("./test.png")
    return


def create_heatmap(distances_file):
    distances_file = Path(distances_file)
    output_file = distances_file.parent / (distances_file.stem + ".png")
    distances = pd.read_csv(distances_file, sep="\t", index_col=0).fillna(0)
    # new_index = list(distances.index)

    # for sl4
    # new_index = [i.split(".")[0] for i in new_index]
    # for sl5 - get the name of the acession without number added by Jasmine
    # new_index = [i.split("_")[1] for i in new_index]

    # distances.index = new_index
    distances = distances.sort_index()

    fig = px.imshow(distances, color_continuous_scale="Greens", aspect="auto", zmin=0, zmax=1)
    fig.update_layout(yaxis=dict(tickmode="linear"), title=distances_file.stem)
    fig.write_image(output_file)
    return


def create_scored_heatmap(called_intro_file, gt_intro_file, threshold):
    # read both distances files
    output_file = called_intro_file.parent / (called_intro_file.stem + ".scored.png")

    # threshold and rename the accessions to match
    called_intro_df, gt_intro_df = threshold_introgressions(
        called_intro_file, gt_intro_file, threshold
    )

    # only visualize the accessions that are shared btwn both files
    shared_accessions = list(
        set(called_intro_df.index.values).intersection(set(gt_intro_df.index.values))
    )
    called_intro_df = called_intro_df[called_intro_df.index.isin(shared_accessions)].sort_index()
    gt_intro_df = gt_intro_df[gt_intro_df.index.isin(shared_accessions)].sort_index()

    # label all positions as green for TP, yellow for FP, red for FN, white/blank for TN
    # true pos
    called_intro_df[(called_intro_df == 1) & (gt_intro_df == 1)] = 5

    # true neg
    called_intro_df[(called_intro_df == 0) & (gt_intro_df == 0)] = 4

    # false pos
    called_intro_df[(called_intro_df == 1) & (gt_intro_df == 0)] = 3

    # false neg
    called_intro_df[(called_intro_df == 0) & (gt_intro_df == 1)] = 2

    # visualize
    fig = px.imshow(
        called_intro_df, color_continuous_scale=["red", "yellow", "gray", "green"], aspect="auto"
    )
    fig.update_layout(
        coloraxis_showscale=False, yaxis=dict(tickmode="linear"), title=output_file.stem
    )
    fig.write_image(output_file)
    return


def create_heatmap_runner():
    # NOTE: change folder here
    input_folder = Path(
        "/home/nbrown62/data_mschatz1/nbrown62/panagram_data/tomato_sl4/introgression_analysis_v2/postprocessed"
    )

    distances_files = list(input_folder.glob("chr*.txt"))
    for file in distances_files:
        print(f"Visualizing {file.name}")
        create_heatmap(file)
    return


def create_scored_heatmap_runner():
    # NOTE: change folders here
    called_intros_folder = Path(
        "/home/nbrown62/data_mschatz1/nbrown62/panagram_data/tomato_sl4/introgression_analysis_v3/postprocessed"
    )
    gt_intros_folder = Path(
        "/home/nbrown62/data_mschatz1/nbrown62/CallIntrogressions_data/tomato_sl4_paper"
    )
    introgression_type = "REF"
    threshold = 0.5

    # get gt and called intro files
    called_intros_files = list(called_intros_folder.glob(f"chr*{introgression_type}.txt"))
    gt_intros_files_sp = list(gt_intros_folder.glob(f"chr*SP.txt"))
    gt_intros_files_slc = list(gt_intros_folder.glob(f"chr*SLC.txt"))

    if introgression_type == "SP":
        gt_intros_files = gt_intros_files_sp
    elif introgression_type == "SLC":
        gt_intros_files = gt_intros_files_slc
    elif introgression_type == "merged" or introgression_type == "REF":
        # merge files together
        gt_intros_files_sp.sort()
        gt_intros_files_slc.sort()
        gt_intros_files = list(zip(gt_intros_files_slc, gt_intros_files_sp))

    # sort the files and make sure there are equal numbers of files
    called_intros_files.sort()
    gt_intros_files.sort()

    if len(called_intros_files) != len(gt_intros_files):
        raise ValueError("Unequal numbers of GT and called chromosome files...")

    # pass tuples of files to create_scored_heatmap for each chromosome
    for called_intros_file, gt_intros_file in zip(called_intros_files, gt_intros_files):
        create_scored_heatmap(called_intros_file, gt_intros_file, threshold=0.5)
    return


if __name__ == "__main__":
    # create_heatmap_runner()
    create_scored_heatmap_runner()

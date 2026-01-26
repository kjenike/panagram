import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from panagram.index import Index
from call_introgressions import bins_to_bed
import postprocess_introgressions as post
import plotly.express as px


def merge_bed_files(bed_files, bin_size, chr_length):
    """Convert BED files to binned format and append them together as rows in a DataFrame to match
    the style of the DataFrame returned by read_text_file.

    Args:
        bed_files (list): list of paths to BED files
        bin_size (int): size of the bins
        chr_length (int): length of the chromosome

    Returns:
        pd.DataFrame: merged DataFrame containing all BED files
    """

    merged_df = None
    for bed_file in bed_files:
        parts = bed_file.stem.split("_")
        bed_accession = "_".join(parts[:-2])

        # make sure empty bed files get an empty row
        if post.bed_file_is_empty(bed_file):
            bed_df = post.get_intro_df_template(bin_size, chr_length)
        else:
            bed_df = post.read_bed_file(bed_file)
            bed_df = post.bed_to_bins(bed_df, bin_size, chr_length).T

        bed_df.index = [bed_accession]
        bed_df = bed_df.rename_axis("Sample")

        if merged_df is None:
            merged_df = bed_df
        else:
            # append row to df with other accessions
            merged_df = pd.concat([merged_df, bed_df.iloc[[0]]])

    merged_df = merged_df.sort_index()
    return merged_df


def merge_text_files(text_files):
    """Merge DataFrames from read_text_file together using max operator. The highest value between
    multiple Jaccard similarities is taken as the final value for each bin.

    Args:
        text_files (list): list of paths to text files containing Jaccard similarities

    Returns:
        pd.DataFrame: merged DataFrame containing all text files
    """

    text_dfs = []
    for text_file in text_files:
        text_dfs.append(read_text_file(text_file))
    merged_df = pd.concat(text_dfs).groupby(level=0).max()
    return merged_df


def read_text_file(text_file):
    """Read Jaccard similarities for a chromosome from get_distances.py.

    Args:
        text_file (str): path to the text file containing Jaccard similarities

    Returns:
        pd.DataFrame: DataFrame containing text file information
    """

    intro_df = pd.read_csv(text_file, sep="\t", header=0, index_col=0).fillna(0)
    intro_df.columns = intro_df.columns.astype(int)
    return intro_df


def rescale_introgressions_helper(row, original_bin_size, new_bin_size, chr_length):
    """Rescale DataFrame to a new bin size.

    Args:
        row (pd.Series): Series containing introgression data
        original_bin_size (int): original bin size
        new_bin_size (int): new bin size
        chr_length (int): length of the chromosome

    Returns:
        pd.Series: rescaled introgression data
    """

    bins_df = row.rename("introgression").to_frame()
    bins_df.index = bins_df.index.astype(int)
    bed_df = bins_to_bed(bins_df, original_bin_size, "nan", "nan")
    col_names = ["Chromosome", "Start", "End", "Notes"]
    bed_df.columns = col_names
    bed_df["Sequence"] = None
    if bed_df.empty:
        new_row = post.get_intro_df_template(new_bin_size, chr_length).loc[0]
    else:
        new_row = post.bed_to_bins(bed_df, new_bin_size, chr_length)["introgression"]
    return new_row


def merge_centromere_regions_helper(row, bin_size, chr_length, fasta_df):
    """Uses merge_centromere_regions to fill in centromere regions between introgressions. Helper to
    convert between bins and bed formats.

    Args:
        row (pd.Series): Series containing introgression data
        bin_size (int): bin size
        chr_length (int): length of the chromosome
        fasta_df (pd.DataFrame): DataFrame containing FASTA information from the reference genome

    Returns:
        pd.Series: Merged centromere regions.
    """

    # convert row to bed file - save temporarily
    bins_df = row.rename("introgression").to_frame()
    bins_df.index = bins_df.index.astype(int)

    # match format that merge_centromere_regions expects
    bed_df = bins_to_bed(bins_df, bin_size, "nan", "nan")
    col_names = ["Chromosome", "Start", "End", "Notes"]
    bed_df.columns = col_names
    bed_df["Sequence"] = None

    # merge chrs
    bed_df = post.merge_centromere_regions(bed_df, fasta_df, bin_size)

    # convert back to bins
    new_row = post.bed_to_bins(bed_df, bin_size, chr_length)["introgression"]
    return new_row


def threshold_introgressions_helper(intro_df, threshold):
    """Threshold helper function. Values below the threshold
    are set to 0 (not introgressed), and values equal to or above the threshold are set to 1
    (introgressed).

    Args:
        intro_df (pd.DataFrame): DataFrame containing introgression data
        threshold (float): threshold value for introgression

    Returns:
        pd.DataFrame: thresholded introgression data
    """

    intro_df[intro_df < threshold] = 0
    intro_df[intro_df != 0] = 1
    return intro_df.astype(int)


def threshold_introgressions(pred_df, gt_df, threshold):
    """Threshold Jaccard similarities for scoring and make sure all predicted introgressions
    are labeled as a 1.

    Args:
        pred_df (pd.DataFrame): DataFrame containing predicted introgression data
        gt_df (pd.DataFrame): DataFrame containing ground truth introgression data
        threshold (float): threshold value for introgression

    Returns:
        pd.DataFrame: thresholded predicted and ground truth introgression data
    """

    gt_df = threshold_introgressions_helper(gt_df, threshold)
    pred_df = threshold_introgressions_helper(pred_df, threshold=1)
    return pred_df, gt_df


def score_introgression_overlaps_helper(pred_df_row, gt_df_row, overlap_threshold):
    """Score introgression overlaps between predicted and ground truth data.

    Args:
        pred_df_row (pd.Series): Series containing predicted introgression data
        gt_df_row (pd.Series): Series containing ground truth introgression data
        overlap_threshold (float): min. fraction of overlap required to consider a region a TP

    Returns:
        tuple(dict, pd.Series): dictionary containing overlap scores and gt_df_row with relabeled
        TP/FN/TN/FP regions for visualization
    """

    # label every region in pred as TP, TN, FP, or FN
    # using regions defined by pred means that FPs aren't ignored
    gt_df_row = gt_df_row.to_frame(name="introgression")
    pred_df_row = pred_df_row.to_frame(name="introgression")

    # arcane code to label consecutive strings of 1s and 0s as separate regions
    # finds where there are differences between the previous row and the current
    # where there are differences, labels them as different regions using cumsum
    pred_df_row["region"] = (pred_df_row["introgression"].diff().fillna(1) != 0).cumsum()

    pred_ones = pred_df_row.groupby("region")["introgression"].sum()
    counts = pred_df_row.groupby("region")["introgression"].count()
    gt_ones = gt_df_row.groupby(pred_df_row["region"])["introgression"].sum()

    # Combine the results
    region_df = pd.DataFrame({"gt_ones": gt_ones, "pred_ones": pred_ones, "count": counts})

    # Count TP, TN, FP, and FN regions
    # TP - where the pred has 1s and the gt region has >= overlap_threshold% 1s
    tp_regions = region_df[
        (region_df["pred_ones"] > 0)
        & ((region_df["gt_ones"] / region_df["count"]) >= overlap_threshold)
    ]
    # FN - where the pred has 0s and the gt region has < overlap_threshold% 0s
    fn_regions = region_df[
        (region_df["pred_ones"] == 0)
        & ((region_df["count"] - region_df["gt_ones"]) / region_df["count"] < overlap_threshold)
    ]
    # TN - where the pred has 0s and the gt region has >= overlap_threshold% 0s
    tn_regions = region_df[
        (region_df["pred_ones"] == 0)
        & ((region_df["count"] - region_df["gt_ones"]) / region_df["count"] >= overlap_threshold)
    ]
    # FP - where the pred has 1s and the gt region has < overlap_threshold% 1s
    fp_regions = region_df[
        (region_df["pred_ones"] > 0)
        & ((region_df["gt_ones"] / region_df["count"]) < overlap_threshold)
    ]

    true_pos = len(tp_regions)
    false_pos = len(fp_regions)
    true_neg = len(tn_regions)
    false_neg = len(fn_regions)

    # relabel regions in gt_df for visualization
    if not tp_regions.empty:
        # tps should appear as 1s on gt
        tp_regions = list(tp_regions.index)
        gt_df_row.loc[pred_df_row["region"].isin(tp_regions)] = 1
    if not fn_regions.empty:
        # fns should appear as 1s on gt
        fn_regions = list(fn_regions.index)
        gt_df_row.loc[pred_df_row["region"].isin(fn_regions)] = 1
    if not tn_regions.empty:
        # tns should appear as 0s on gt
        tn_regions = list(tn_regions.index)
        gt_df_row.loc[pred_df_row["region"].isin(tn_regions)] = 0
    if not fp_regions.empty:
        # fps should appear as 0s on gt
        fp_regions = list(fp_regions.index)
        gt_df_row.loc[pred_df_row["region"].isin(fp_regions)] = 0

    # return the metrics for the regions and the modified gt_df_row
    metrics = {
        "True Positive": true_pos,
        "True Negative": true_neg,
        "False Positive": false_pos,
        "False Negative": false_neg,
    }
    gt_df_row = list(gt_df_row["introgression"])
    return metrics, gt_df_row


def score_introgression_overlaps(pred_df, gt_df, overlap_threshold):
    """Get confusion matrix and related metrics for introgressions given ground truth in the same
    coordinate/bin space. Accession names must match between pred_df and gt_df.
    Args:
        pred_df (pd.DataFrame): DataFrame containing predicted introgression data
        gt_df (pd.DataFrame): DataFrame containing ground truth introgression data
        overlap_threshold (float): min. fraction of overlap required to consider a region a TP
    Returns:
        pd.DataFrame: DataFrame containing confusion matrix and related metrics
    """

    # score introgressions per-region rather than per-bin
    shared_cols = list(set(pred_df.index).intersection(set(gt_df.index)))
    pred_df = pred_df.transpose()[shared_cols]
    gt_df = gt_df.transpose()[shared_cols]

    # calculate overlap metrics per column
    metrics_df = gt_df.apply(
        lambda col: score_introgression_overlaps_helper(pred_df[col.name], col, overlap_threshold)[
            0
        ],
        axis=0,
    ).tolist()
    metrics_df = pd.DataFrame(metrics_df, index=gt_df.columns)

    # sum across cols to get totals across all columns
    true_pos = metrics_df["True Positive"].sum()
    true_neg = metrics_df["True Negative"].sum()
    false_pos = metrics_df["False Positive"].sum()
    false_neg = metrics_df["False Negative"].sum()
    total = true_pos + true_neg + false_pos + false_neg

    with np.errstate(invalid="ignore"):
        acc = (true_pos + true_neg) / total
        precision = true_pos / (true_pos + false_pos)
        recall = true_pos / (true_pos + false_neg)
        fpr = false_pos / (false_pos + true_neg)

    # return the metrics as a df
    metrics = {
        "True Positive": true_pos,
        "True Negative": true_neg,
        "False Positive": false_pos,
        "False Negative": false_neg,
        "Accuracy": acc,
        "Precision": precision,
        "Recall": recall,
        "FPR": fpr,
    }
    metrics = pd.DataFrame([metrics])
    return metrics


def score_introgressions(pred_df, gt_df):
    """Get confusion matrix and related metrics for introgressions given ground truth in the same
    coordinate/bin space. Accession names must match between pred_df and gt_df.

    Args:
        pred_df (pd.DataFrame): DataFrame containing predicted introgression data
        gt_df (pd.DataFrame): DataFrame containing ground truth introgression data

    Returns:
        pd.DataFrame: DataFrame containing confusion matrix and related metrics
    """

    # rotate dfs, sort cols, and drop all columns that aren't shared btwn called and gt
    shared_cols = list(set(pred_df.index).intersection(set(gt_df.index)))
    pred_df = pred_df.transpose()[shared_cols]
    gt_df = gt_df.transpose()[shared_cols]

    # calculate confusion matrix
    total = gt_df.size

    # true pos
    true_pos = ((pred_df == 1) & (gt_df == 1)).values.sum()
    # true neg
    true_neg = ((pred_df == 0) & (gt_df == 0)).values.sum()
    # false pos
    false_pos = ((pred_df == 1) & (gt_df == 0)).values.sum()
    # false neg
    false_neg = ((pred_df == 0) & (gt_df == 1)).values.sum()

    with np.errstate(invalid="ignore"):
        acc = (true_pos + true_neg) / total
        precision = true_pos / (true_pos + false_pos)
        recall = true_pos / (true_pos + false_neg)
        fpr = false_pos / (false_pos + true_neg)

    # return the metrics as a df
    metrics = {
        "True Positive": true_pos,
        "True Negative": true_neg,
        "False Positive": false_pos,
        "False Negative": false_neg,
        "Accuracy": acc,
        "Precision": precision,
        "Recall": recall,
        "FPR": fpr,
    }
    metrics = pd.DataFrame([metrics])

    return metrics


def create_scored_heatmap_overlaps(pred_df, gt_df, overlap_threshold, output_file):
    """Helper function to create a scored heatmap for introgression overlaps. Performs modification
    of gt_df before passing to create_scored_heatmap. Only visualizes the accessions that are shared
    between both pred_df and gt_df.

    Args:
        pred_df (pd.DataFrame): DataFrame containing predicted introgression data
        gt_df (pd.DataFrame): DataFrame containing ground truth introgression data
        overlap_threshold (float): minimum overlap threshold for considering a region a true positive
        output_file (str or Path): Path to the output file for the heatmap
    """

    shared_accessions = list(set(pred_df.index.values).intersection(set(gt_df.index.values)))
    pred_df = pred_df[pred_df.index.isin(shared_accessions)].sort_index()
    gt_df = gt_df[gt_df.index.isin(shared_accessions)].sort_index()

    gt_cols = gt_df.columns
    # calculate overlap metrics per column
    gt_df = gt_df.apply(
        lambda row: score_introgression_overlaps_helper(
            pred_df.loc[row.name], row, overlap_threshold
        )[1],
        axis=1,
        result_type="expand",
    )
    gt_df.columns = gt_cols
    create_scored_heatmap(pred_df, gt_df, output_file)
    return


def create_scored_heatmap(pred_df, gt_df, output_file, groups=None, xaxis_dtick=2000000):
    """Create a scored heatmap for introgressions given ground truth in the same coordinate/bin
    space.

    Args:
        pred_df (pd.DataFrame): DataFrame containing predicted introgression data
        gt_df (pd.DataFrame): DataFrame containing ground truth introgression data
        output_file (str or Path): Path to the output file for the heatmap
        groups (pd.Series, optional): accession groups, used to determine ordering of heatmap rows, defaults to None
        xaxis_dtick (int, optional): the x-axis tick interval, defaults to 2000000
    """

    pred_df = pred_df.copy()
    # only visualize the accessions that are shared btwn both files
    shared_accessions = list(set(pred_df.index.values).intersection(set(gt_df.index.values)))
    pred_df = pred_df[pred_df.index.isin(shared_accessions)].sort_index()
    gt_df = gt_df[gt_df.index.isin(shared_accessions)].sort_index()

    # label all positions as green for TP, yellow for FP, red for FN, white/blank for TN
    # this won't work if the column labels are not the same name/type
    # true pos
    pred_df[(pred_df == 1) & (gt_df == 1)] = 5
    # true neg
    pred_df[(pred_df == 0) & (gt_df == 0)] = 4
    # false pos
    pred_df[(pred_df == 1) & (gt_df == 0)] = 3
    # false neg
    pred_df[(pred_df == 0) & (gt_df == 1)] = 2

    if groups is not None:
        # reindex pred_df to match groups
        ordered_names = groups.index.tolist()
        ordered_names = [name for name in ordered_names if name in pred_df.index]
        pred_df = pred_df.reindex(index=ordered_names)

    # visualize
    fig = px.imshow(
        pred_df,
        color_continuous_scale=["#EE6677", "#CCBB44", "#BBBBBB", "#228833"],
        range_color=[2, 5],
        aspect="auto",
    )
    fig.update_layout(
        font=dict(family="Arial", color="black"),
        coloraxis_showscale=False,
        yaxis=dict(tickmode="linear", title=""),
        title=dict(
            text=output_file.stem,
            font=dict(size=20),
        ),
        xaxis=dict(
            dtick=xaxis_dtick,
            title=dict(
                text="Genomic Position",
                font=dict(size=16),
            ),
        ),
    )
    fig.write_image(output_file)

    # save out pred df to make more visuals with if needed
    pred_df.to_csv(output_file.with_suffix(".csv"), sep="\t")
    return


def main():
    # scores introgressions using ground truth derived from SV calling and get_distances
    parser = argparse.ArgumentParser(description="Introgression scoring.")
    parser.add_argument(
        "--pre",
        type=str,
        help="path to file/folder of predicted introgression bed files",
        required=True,
    )
    parser.add_argument(
        "--gdt",
        type=str,
        help="path to file/folder of ground truth introgression text files",
        required=True,
    )
    parser.add_argument(
        "--how",
        type=str,
        help="whether to use bins or overlaps in scoring calculations",
        required=True,
    )
    parser.add_argument("--idx", type=str, help="path to Panagram index folder", required=True)
    parser.add_argument("--ref", type=str, help="name of reference in Panagram", required=True)
    parser.add_argument("--out", type=str, help="path to folder to save all outputs", required=True)
    parser.add_argument("--vis", action="store_true", help="save visualized results")
    parser.add_argument(
        "--grp", type=str, help="path to groups file for heatmap visualization", default=None
    )
    parser.add_argument(
        "--bin",
        type=int,
        help="size of bitmap bin used during calling - gt is rescaled to match this",
        default=1000000,
    )
    parser.add_argument(
        "--min",
        type=int,
        help="minimum number of bins an introgression must be; all smaller are clipped by rmbn",
        default=4,
    )
    parser.add_argument(
        "--gap",
        type=int,
        help="maximum number of bins to fill between introgressions for fgap",
        default=1,
    )
    parser.add_argument(
        "--thr",
        type=float,
        default=0.5,
        help="ground truth Jaccard similarity lower threshold to be considered an introgression",
    )
    parser.add_argument(
        "--othr",
        type=float,
        default=0.8,
        help="min. overlap threshold for an overlap to be considered a TP/TN",
    )
    parser.add_argument(
        "--cmp",
        nargs="+",
        help="for REF/merged introgression types, list all comp groups",
    )
    parser.add_argument(
        "--act",
        nargs="+",
        help="action(s) to perform on introgression file: fgap, fcen, rmbn",
    )
    print("Scoring introgressions...", flush=True)
    args = parser.parse_args()

    # input checking
    gt_intros = Path(args.gdt)
    if not gt_intros.is_dir() and not gt_intros.is_file():
        raise ValueError(f"Ground truth file/directory not found. Check --gdt path.")

    ref_accession = args.ref
    bin_size = args.bin

    index_dir = Path(args.idx)
    if not index_dir.is_dir():
        raise ValueError("Index directory not found. Check --idx path.")
    index = Index(index_dir)
    ref_genome = index.genomes[ref_accession]

    actions = args.act
    if actions is not None:
        for action in actions:
            if action not in ["fgap", "fcen", "rmbn"]:
                raise ValueError(
                    f"Unrecognized action {action}. Check --act flag for valid actions."
                )

    # figure out if user provided a file or folder
    bed_files = Path(args.pre)
    if bed_files.is_file():
        bed_files = [bed_files]
    elif bed_files.is_dir():
        bed_files = list(bed_files.glob(f"*.bed"))
    else:
        raise ValueError("Bed file/folder not found. Check --bed path.")

    # figure out names of chrs and intros to loop through
    chrs = set()
    intro_types = set()
    for bed_file in bed_files:
        parts = bed_file.stem.split("_")
        bed_intro_type = parts[-1]
        bed_chr = parts[-2]
        chrs.add(bed_chr)
        intro_types.add(bed_intro_type)
    chrs = list(chrs)
    intro_types = list(intro_types)
    chrs.sort()
    intro_types.sort()

    intro_groups = args.cmp
    if ("REF" in intro_types or "merged" in intro_types) and (intro_groups is None):
        raise ValueError("--cmp is required to process REF/merged bed files.")

    output_dir = Path(args.out)
    output_dir.mkdir(parents=True, exist_ok=True)
    pred_dir = output_dir / "pred"
    pred_dir.mkdir(parents=True, exist_ok=True)
    gt_dir = output_dir / "gt_postprocessed"
    gt_dir.mkdir(parents=True, exist_ok=True)

    threshold = args.thr
    min_size = args.min
    gap_size = args.gap
    how_to_score = args.how
    if how_to_score not in ["bins", "overlaps"]:
        raise ValueError("--how must be either 'bins' or 'overlaps'.")
    if how_to_score == "overlaps":
        overlap_threshold = args.othr

    render_vis = args.vis
    if render_vis:
        vis_dir = output_dir / "heatmaps"
        vis_dir.mkdir(parents=True, exist_ok=True)

    all_metrics = {}
    for chr in chrs:
        for intro_type in intro_types:
            # find the correct gt_intro_file if it is a dir
            # read in corresponding GT introgressions file
            if gt_intros.is_file():
                gt_intro_file = gt_intros
            else:
                if intro_type in ["REF", "merged"]:
                    gt_intro_files = []
                    for intro_group in intro_groups:
                        gt_intro_files += list(gt_intros.glob(f"{chr}_{intro_group}.txt"))
                    if len(gt_intro_files) != len(intro_groups):
                        raise ValueError(
                            f"--gt directory specified does not contain all files needed for REF/merged scoring."
                        )
                    gt_df = merge_text_files(gt_intro_files)
                else:
                    gt_intro_file = list(gt_intros.glob(f"{chr}_{intro_type}.txt"))
                    if not gt_intro_file:
                        raise ValueError(
                            f"--gt directory specified does not contain {chr}_{intro_type}.txt"
                        )
                    gt_df = read_text_file(gt_intro_file[0])

            # get chr length from ref
            chr_length = ref_genome.sizes[chr]

            # find and merge predicted bed files; this rescales automatically if needed
            accession_bed_files = [
                path for path in bed_files if path.name.endswith(f"_{chr}_{intro_type}.bed")
            ]
            pred_df = merge_bed_files(accession_bed_files, bin_size, chr_length)

            # save off merged file as text
            pred_output_file = pred_dir / f"{chr}_{intro_type}.txt"
            pred_df.to_csv(pred_output_file, sep="\t")

            # perform actions on gt
            pred_df, gt_df = threshold_introgressions(pred_df, gt_df, threshold)

            # rescale gt to match pred if needed
            # original_bin_size = int(gt_df.columns[1])
            # if bin_size != original_bin_size:
            #     gt_df = gt_df.apply(
            #         rescale_introgressions_helper,
            #         original_bin_size=original_bin_size,
            #         new_bin_size=bin_size,
            #         chr_length=chr_length,
            #         axis=1,
            #     )

            # rescale pred to match gt if needed
            gt_bin_size = int(gt_df.columns[1])
            if bin_size != gt_bin_size:
                pred_df = pred_df.apply(
                    rescale_introgressions_helper,
                    original_bin_size=bin_size,
                    new_bin_size=gt_bin_size,
                    chr_length=chr_length,
                    axis=1,
                )
                bin_size = gt_bin_size

            if actions:
                for action in actions:
                    if action in "fgap":
                        # fgap - merge regions >=1 bin size apart
                        # save col names - they can be broken after an apply
                        col_names = gt_df.columns
                        gt_df = gt_df.apply(
                            post.fill_gaps,
                            gap_size=gap_size,
                            axis=1,
                            result_type="expand",
                        )
                        gt_df.columns = col_names

                    elif action in "rmbn":
                        # rmbn - perform singleton bin removal
                        # save col names - they can be broken after an apply
                        col_names = gt_df.columns
                        gt_df = gt_df.apply(
                            post.remove_small_regions,
                            min_size=min_size,
                            axis=1,
                            result_type="expand",
                        )
                        gt_df.columns = col_names

                    elif action == "fcen":
                        # fcen - fill in centromeres between introgression regions - requires ref
                        # will have to do this per row
                        fasta_file = index_dir / ref_genome.fasta
                        fasta_df = post.read_fasta(fasta_file)
                        gt_df = gt_df.apply(
                            merge_centromere_regions_helper,
                            bin_size=bin_size,
                            fasta_df=fasta_df,
                            chr_length=chr_length,
                            axis=1,
                        )

                # save
                gt_output_file = gt_dir / f"{chr}_{intro_type}.txt"
                gt_df.to_csv(gt_output_file, sep="\t")

            # score introgressions
            if how_to_score == "bins":
                metrics = score_introgressions(pred_df, gt_df)
                # visualize introgressions
                if render_vis:
                    vis_output_file = vis_dir / f"{chr}_{intro_type}.png"
                    if args.grp is not None:
                        groups = pd.read_csv(args.grp, sep="\t", index_col=0)
                        create_scored_heatmap(pred_df, gt_df, vis_output_file, groups=groups)
                    else:
                        create_scored_heatmap(pred_df, gt_df, vis_output_file)
            else:
                metrics = score_introgression_overlaps(pred_df, gt_df, overlap_threshold)
                # visualize introgressions
                if render_vis:
                    vis_output_file = vis_dir / f"{chr}_{intro_type}.png"
                    create_scored_heatmap_overlaps(
                        pred_df, gt_df, overlap_threshold, vis_output_file
                    )
            metrics.index = [chr]

            # aggregate metrics
            if intro_type in all_metrics:
                all_metrics[intro_type] = pd.concat([all_metrics[intro_type], metrics.iloc[[0]]])
            else:
                all_metrics[intro_type] = metrics

    # save aggregated metrics
    for intro_type, all_metrics_df in all_metrics.items():
        all_metrics_file = output_dir / f"metrics_{intro_type}.tsv"
        all_metrics_df.to_csv(all_metrics_file, sep="\t")

    return


if __name__ == "__main__":
    main()

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import plotly.express as px
from PIL import Image
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.utils import ImageReader
import io


def create_heatmap(text_file, groups=None):
    """Create a heatmap of Jaccard similarities from get_distances.py.

    Args:
        text_file (str or Path): path to distances text file
        groups (pd.Series, optional): accession groups, used to determine ordering of heatmap rows, defaults to None
    """

    text_file = Path(text_file)
    output_file = text_file.parent / (text_file.stem + ".svg")
    distances = pd.read_csv(text_file, sep="\t", index_col=0).fillna(0)
    distances = distances.sort_index()
    distances.columns = distances.columns.astype(int)

    if groups is not None:
        ordered_names = groups.index.tolist()
        # drop any names not in the distances index
        ordered_names = [name for name in ordered_names if name in distances.index]
        distances = distances.reindex(index=ordered_names)

    fig = px.imshow(distances, color_continuous_scale="Greens", aspect="auto", zmin=0, zmax=1)
    fig.update_layout(
        font=dict(family="Helvetica Bold", color="black"),
        coloraxis_colorbar=dict(
            title=dict(
                text="Jaccard Similarity",
                side="top",  # Title above the colorbar
            ),
            orientation="h",  # horizontal colorbar
            x=0.3,  # center the colorbar
            xanchor="left",
            y=1.08,  # move above the plot
            len=0.7,  # shorten the colorbar
            thickness=20,  # make it thicker
        ),
        title=dict(
            text=text_file.stem,
            # x=0.5,  # Center the title
            # xanchor="right",
            font=dict(size=20),
        ),
        xaxis=dict(
            dtick=2000000,
            title=dict(
                text="Genomic Position",
                font=dict(size=16),
            ),
        ),
        width=700,
        height=500,
        yaxis=dict(tickmode="linear", title=""),
    )
    fig.write_image(output_file)
    return


def create_heatmap_runner(input_dir, groups=None):
    """Create heatmaps for all distance files in the input directory.

    Args:
        input_dir (Path): path to input directory containing distance files
        groups (pd.Series, optional): accession groups, used to determine ordering of heatmap rows, defaults to None.
    """

    distances_files = list(input_dir.glob("chr*.txt"))
    distances_files = [file for file in distances_files if "max_species" not in file.name]
    for file in distances_files:
        print(f"Visualizing {file.name}")
        create_heatmap(file, groups=groups)
    return


def calculate_pr_auc(input_dir, intro_type, how_to_score, thresholds):
    """Calculate precision-recall AUC for introgressions.

    Args:
        input_dir (Path): path to input directory containing distance files
        intro_type (str): introgression type (could be a group name, 'REF', or 'merged')
        how_to_score (str): method for scoring introgressions ('overlaps' or 'bins')
        thresholds (list): list of thresholds to evaluate

    Returns:
        pd.DataFrame: DataFrame containing PR and related metrics for each threshold
    """

    scored_dir_name = f"scored_{how_to_score}"

    # find the file we need
    results = []
    for threshold in thresholds:
        all_metrics_file = (
            f"{input_dir}/{input_dir.name}_{threshold}/{scored_dir_name}/metrics_{intro_type}.tsv"
        )

        # create df of precisions and recalls across thresholds
        metrics = pd.read_csv(all_metrics_file, sep="\t", index_col=0)
        tp = metrics["True Positive"].sum()
        fp = metrics["False Positive"].sum()
        fn = metrics["False Negative"].sum()
        tn = metrics["True Negative"].sum()
        results.append({"threshold": threshold, "TP": tp, "FP": fp, "FN": fn, "TN": tn})

    results_df = pd.DataFrame(results).fillna(0)
    results_df["precision"] = (results_df["TP"] / (results_df["TP"] + results_df["FP"])).fillna(0)
    results_df["recall"] = (results_df["TP"] / (results_df["TP"] + results_df["FN"])).fillna(0)
    results_df["MCC"] = (
        (results_df["TP"] * results_df["TN"] - results_df["FP"] * results_df["FN"])
        / np.sqrt(
            (results_df["TP"] + results_df["FP"])
            * (results_df["TP"] + results_df["FN"])
            * (results_df["TN"] + results_df["FP"])
            * (results_df["TN"] + results_df["FN"])
        )
    ).fillna(0)

    # Anchor NaN points at (1,0) for plotting
    # for i in range(len(results_df)):
    #     if results_df["recall"][i] == 0 and results_df["precision"][i] == 0:
    #         results_df.at[i, "recall"] = 0
    #         results_df.at[i, "precision"] = 1

    auc_pr = np.trapz(results_df["precision"].values, results_df["recall"].values)
    print(f"PR AUC {intro_type}:", auc_pr)
    print(f"Max MCC {intro_type}:", results_df["MCC"].max())
    return results_df


def create_mcc_curve(input_dir, intro_type, how_to_score, thresholds):
    """Create a Matthews Correlation Coefficient (MCC) curve.

    Args:
        input_dir (Path): path to input directory containing distance files
        intro_type (str): introgression type (could be a group name, 'REF', or 'merged')
        how_to_score (str): method for scoring introgressions ('overlaps' or 'bins')
        thresholds (list): list of thresholds to evaluate
    """

    results_df = calculate_pr_auc(input_dir, intro_type, how_to_score, thresholds)
    fig = px.line(
        results_df,
        x="threshold",
        y="MCC",
    )

    fig.update_traces(textposition="top center")
    fig.update_layout(
        font=dict(family="Helvetica Bold", color="black"),
        template="plotly_white",
        xaxis_title=dict(text="Threshold"),
        yaxis_title=dict(text="Matthews Correlation Coefficient (MCC)"),
        xaxis=dict(range=[0, 1.01], ticks="outside", linecolor="black"),
        yaxis=dict(range=[0, 1.01], ticks="outside", linecolor="black"),
        width=500,
        height=500,
        margin=dict(l=10, r=10, t=10, b=10),  # tight layout
    )
    output_file = f"{input_dir}/{how_to_score}_{intro_type}_mcc.svg"
    fig.write_image(output_file)
    return


def create_pr_curve(input_dir, intro_type, how_to_score, thresholds):
    """Create a precision-recall curve.

    Args:
        input_dir (Path): path to input directory containing distance files
        intro_type (str): introgression type (could be a group name, 'REF', or 'merged')
        how_to_score (str): method for scoring introgressions ('overlaps' or 'bins')
        thresholds (list): list of thresholds to evaluate
    """

    results_df = calculate_pr_auc(input_dir, intro_type, how_to_score, thresholds)
    results_df = results_df[~((results_df["recall"] == 0) & (results_df["precision"] == 1))]
    results_df = results_df[~((results_df["recall"] == 0) & (results_df["precision"] == 0))]

    fig = px.line(
        results_df,
        x="recall",
        y="precision",
    )
    fig.update_traces(textposition="top center")
    fig.update_layout(
        font=dict(family="Helvetica Bold", color="black"),
        template="plotly_white",
        xaxis_title=dict(text="Recall"),
        yaxis_title=dict(text="Precision"),
        xaxis=dict(range=[0, 1.01], ticks="outside", linecolor="black"),
        yaxis=dict(range=[0, 1.01], ticks="outside", linecolor="black"),
        width=500,
        height=500,
        margin=dict(l=10, r=10, t=10, b=10),  # tight layout
    )
    output_file = f"{input_dir}/{how_to_score}_{intro_type}_prc.svg"
    fig.write_image(output_file)
    return


def create_pr_curve_accessions(input_dir, intro_type, how_to_score, thresholds):
    """Create a precision-recall curve with individual lines for each accession.

    Args:
        input_dir (Path): path to input directory containing distance files
        intro_type (str): introgression type (could be a group name, 'REF', or 'merged')
        how_to_score (str): method for scoring introgressions ('overlaps' or 'bins')
        thresholds (list): list of thresholds to evaluate
    """

    scored_dir_name = f"scored_{how_to_score}"
    output_file = f"{input_dir}/{how_to_score}_{intro_type}_prca.svg"

    # find the file we need
    results = []
    for threshold in thresholds:
        threshold_path = input_dir / f"{input_dir.name}_{threshold}" / scored_dir_name / "heatmaps"
        metrics_files = list(threshold_path.glob(f"*_{intro_type}.csv"))
        counts_across_all_chrs = None

        # add up TP, TN, FP, FN across all chrs for each accession
        for metrics_file in metrics_files:
            # Count labels per sample and rename columns for clarity
            accession_metrics = pd.read_csv(metrics_file, sep="\t", index_col=0)
            counts = accession_metrics.apply(lambda row: row.value_counts(), axis=1).fillna(0)
            counts = counts.rename(columns={5: "TP", 4: "TN", 3: "FP", 2: "FN"})

            # Ensure all expected columns are present
            for col in ["TP", "TN", "FP", "FN"]:
                if col not in counts:
                    counts[col] = 0

            if counts_across_all_chrs is None:
                counts_across_all_chrs = counts
            else:
                counts_across_all_chrs += counts

        # Compute precision and recall
        counts_across_all_chrs["Precision"] = counts_across_all_chrs["TP"] / (
            counts_across_all_chrs["TP"] + counts_across_all_chrs["FP"]
        )
        counts_across_all_chrs["Recall"] = counts_across_all_chrs["TP"] / (
            counts_across_all_chrs["TP"] + counts_across_all_chrs["FN"]
        )

        # Handle division by zero (e.g., no TP or FP)
        counts_across_all_chrs["Precision"] = counts_across_all_chrs["Precision"].fillna(0)
        counts_across_all_chrs["Recall"] = counts_across_all_chrs["Recall"].fillna(0)
        # if both precision and recall are 0, set precision to 1 for plotting
        counts_across_all_chrs.loc[
            (counts_across_all_chrs["Precision"] == 0) & (counts_across_all_chrs["Recall"] == 0),
            "Precision",
        ] = 1
        counts_across_all_chrs["Threshold"] = threshold
        counts_across_all_chrs["Sample"] = counts_across_all_chrs.index
        results.append(counts_across_all_chrs[["Sample", "Threshold", "Precision", "Recall"]])

    results = pd.concat(results, ignore_index=True)
    # TODO: change back after testing
    # results["Sample Type"] = [
    #     "HiFi",
    #     "HiFi",
    #     "HiFi",
    #     "HiFi",
    #     "HiFi",
    #     "HiFi",
    #     "HiFi",
    #     "ONT",
    #     "ONT",
    #     "ONT",
    #     "ONT",
    #     "ONT",
    #     "ONT",
    #     "ONT",
    #     "ONT",
    #     "HiFi",
    #     "HiFi",
    #     "HiFi",
    # ]*len(thresholds)

    # results = results[~((results["Recall"] == 0) & (results["Precision"] == 1))]

    # Plot PR curves: one line per sample
    fig = px.line(
        results,
        x="Recall",
        y="Precision",
        color="Sample",
        # color_discrete_sequence=results["Sample Type"].map(
        #     {"HiFi": "#000000", "ONT": "#636EFA"}
        # ),
        color_discrete_sequence=px.colors.qualitative.Light24,
        # markers=True,
        # text="Threshold",
        # title="Precision-Recall Curves per Sample",
    )
    fig.update_traces(textposition="top center")
    fig.update_layout(
        font=dict(family="Helvetica Bold", color="black"),
        template="plotly_white",
        xaxis_title=dict(text="Recall"),
        yaxis_title=dict(text="Precision"),
        xaxis=dict(range=[0, 1.01], ticks="outside", linecolor="black"),
        yaxis=dict(range=[0, 1.01], ticks="outside", linecolor="black"),
        width=700,
        height=500,
        margin=dict(l=10, r=10, t=10, b=10),  # tight layout
    )
    fig.write_image(output_file)
    return


def create_pr_curve_chromosomes(input_dir, intro_type, how_to_score, thresholds):
    """Create a precision-recall curve with individual lines for each chromosome.

    Args:
        input_dir (Path): path to input directory containing distance files
        intro_type (str): introgression type (could be a group name, 'REF', or 'merged')
        how_to_score (str): method for scoring introgressions ('overlaps' or 'bins')
        thresholds (list): list of thresholds to evaluate
    """

    scored_dir_name = f"scored_{how_to_score}"
    output_file = f"{input_dir}/{how_to_score}_{intro_type}_prcc.svg"

    # find the file we need
    results = []
    for threshold in thresholds:
        all_metrics_file = (
            f"{input_dir}/{input_dir.name}_{threshold}/{scored_dir_name}/metrics_{intro_type}.tsv"
        )

        # create df of precisions and recalls across thresholds
        metrics = pd.read_csv(all_metrics_file, sep="\t", index_col=0).fillna(0)
        metrics["Chromosome"] = metrics.index
        metrics["Threshold"] = threshold
        results.append(metrics[["Chromosome", "Threshold", "Precision", "Recall"]])

    results = pd.concat(results, ignore_index=True)

    # Plot PR curves: one line per chromosome
    fig = px.line(
        results,
        x="Recall",
        y="Precision",
        color="Chromosome",
        markers=True,
        # text="Threshold",
        title="Precision-Recall Curves per Chromosome",
    )
    fig.update_traces(textposition="top center")
    fig.update_layout(xaxis_title="Recall", yaxis_title="Precision")
    fig.write_image(output_file)
    return


def create_scored_heatmap_collage(input_dir, intro_type, how_to_score, thresholds):
    """Create a collage of heatmaps for each chromosome through each threshold.

    Args:
        input_dir (Path): path to input directory containing distance files
        intro_type (str): introgression type (could be a group name, 'REF', or 'merged')
        how_to_score (str): method for scoring introgressions ('overlaps' or 'bins')
        thresholds (list): list of thresholds to evaluate
    """

    scored_dir_name = f"scored_{how_to_score}"
    output_file = f"{input_dir}/{how_to_score}_{intro_type}_scored_heatmaps.pdf"

    # image_batches: list of lists, each with 9 image paths
    image_batches = []
    for threshold in thresholds:
        threshold_path = input_dir / f"{input_dir.name}_{threshold}" / scored_dir_name / "heatmaps"
        threshold_heatmaps = list(threshold_path.glob(f"*_{intro_type}.png"))
        threshold_heatmaps.sort()
        image_batches.append(threshold_heatmaps)

    # get heatmaps organized by chr instead of by threshold
    # assumes every chr's heatmap exists
    image_batches = list(zip(*image_batches))
    final_batches = []
    for lst in image_batches:
        final_batches.append(lst[:9])
        final_batches.append(lst[9:])

    # Setup
    pagesize = landscape(A4)
    page_w, page_h = pagesize
    grid_size = (3, 3)

    # Set high-res tiling canvas (e.g., 1000x1000 per cell)
    hires_cell_w, hires_cell_h = 500, 500
    tiled_w = hires_cell_w * grid_size[0]
    tiled_h = hires_cell_h * grid_size[1]

    c = canvas.Canvas(output_file, pagesize=pagesize)

    for batch in final_batches:
        tiled = Image.new("RGB", (tiled_w, tiled_h), color="white")

        for idx, path in enumerate(batch):
            img = Image.open(path).resize((hires_cell_w, hires_cell_h))
            x = (idx % grid_size[0]) * hires_cell_w
            y = (idx // grid_size[0]) * hires_cell_h
            tiled.paste(img, (x, y))

        buf = io.BytesIO()
        tiled.save(buf, format="PNG")
        buf.seek(0)

        # Draw scaled-down image to PDF (full page)
        c.drawImage(ImageReader(buf), 0, 0, width=page_w, height=page_h)
        c.showPage()

    c.save()
    return


def main():
    parser = argparse.ArgumentParser(description="Introgression visualization.")
    parser.add_argument(
        "-v",
        nargs="+",
        help="visual type(s) to create: htmp, shtmp, prc, prcc, prca",
        required=True,
    )
    parser.add_argument(
        "--dir",
        type=str,
        help="path to introgression results folder for visualization",
        required=True,
    )
    parser.add_argument(
        "--grp",
        type=str,
        help="path to groups file for heatmap visualization",
        default=None,
    )
    parser.add_argument(
        "--how",
        type=str,
        help="whether bins or overlaps were used during scoring",
    )
    print("Visualizing results...")
    args = parser.parse_args()

    input_dir = Path(args.dir)
    if not input_dir.is_dir():
        raise ValueError("Input directory not found. Check --dir path.")

    vis_functions = args.v
    for func in vis_functions:
        if func not in ["htmp", "shtmp", "prc", "prcc", "prca", "mcc"]:
            raise ValueError("Unknown visualization function specified. Check args.v.")

    if any(func in vis_functions for func in ("shtmp", "prc", "prcc", "prca")):
        how_to_score = args.how
        if how_to_score not in ["bins", "overlaps"]:
            raise ValueError("--how must be either 'bins' or 'overlaps'.")

        thresholds = [
            0.1,
            0.15,
            0.2,
            0.25,
            0.3,
            0.35,
            0.4,
            0.45,
            0.5,
            0.55,
            0.6,
            0.65,
            0.7,
            0.75,
            0.8,
            0.85,
            0.9,
            0.95,
        ]

        # take the first threshold path and check for introgression types
        scored_dir_name = f"scored_{how_to_score}"
        scored_dir_path = input_dir / f"{input_dir.name}_{thresholds[0]}" / scored_dir_name
        score_tsvs = scored_dir_path.glob("*.tsv")

        intro_types = set()
        for score_tsv in score_tsvs:
            intro_type = score_tsv.stem.split("_")[1]
            intro_types.add(intro_type)

    for func in vis_functions:
        if func == "prc":
            for intro_type in intro_types:
                create_pr_curve(input_dir, intro_type, how_to_score, thresholds)
        elif func == "prcc":
            for intro_type in intro_types:
                create_pr_curve_chromosomes(input_dir, intro_type, how_to_score, thresholds)
        elif func == "prca":
            for intro_type in intro_types:
                create_pr_curve_accessions(input_dir, intro_type, how_to_score, thresholds)
        elif func == "shtmp":
            for intro_type in intro_types:
                create_scored_heatmap_collage(input_dir, intro_type, how_to_score, thresholds)
        elif func == "mcc":
            for intro_type in intro_types:
                create_mcc_curve(input_dir, intro_type, how_to_score, thresholds)
        elif func == "htmp":
            if args.grp is not None:
                groups = pd.read_csv(args.grp, sep="\t", index_col=0)
                create_heatmap_runner(input_dir, groups=groups)
            else:
                create_heatmap_runner(input_dir)
    print("Done.")
    return


if __name__ == "__main__":
    main()

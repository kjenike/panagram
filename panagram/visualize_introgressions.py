"""
Extra functions for visualization of introgressions.
"""

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
    # visualize ground truth bins from get_distances.py by themselves
    distances_file = Path(distances_file)
    output_file = distances_file.parent / (distances_file.stem + ".png")
    distances = pd.read_csv(distances_file, sep="\t", index_col=0).fillna(0)
    distances = distances.sort_index()

    fig = px.imshow(distances, color_continuous_scale="Greens", aspect="auto", zmin=0, zmax=1)
    fig.update_layout(yaxis=dict(tickmode="linear"), title=distances_file.stem)
    fig.write_image(output_file)
    return


def create_heatmap_runner(input_dir):
    distances_files = list(input_dir.glob("chr*.txt"))
    distances_files = [file for file in distances_files if "max_species" not in file.name]
    for file in distances_files:
        print(f"Visualizing {file.name}")
        create_heatmap(file)
    return


def create_pr_curve(input_dir, intro_type, how_to_score, thresholds):
    scored_dir_name = f"scored_{how_to_score}"
    output_file = f"{input_dir}/{how_to_score}_{intro_type}_prc.png"

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

    # plot results
    results_df = pd.DataFrame(results).fillna(0)
    results_df["precision"] = results_df["TP"] / (results_df["TP"] + results_df["FP"]).fillna(0)
    results_df["recall"] = results_df["TP"] / (results_df["TP"] + results_df["FN"]).fillna(0)

    fig = px.line(
        results_df,
        x="recall",
        y="precision",
        markers=True,
        text="threshold",
        title="Precision-Recall Curve",
    )
    fig.update_traces(textposition="top center")
    fig.update_layout(xaxis_title="Recall", yaxis_title="Precision")
    fig.write_image(output_file)
    return


def create_pr_curve_accessions(input_dir, intro_type, how_to_score, thresholds):
    scored_dir_name = f"scored_{how_to_score}"
    output_file = f"{input_dir}/{how_to_score}_{intro_type}_prca.png"

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
        counts_across_all_chrs["Threshold"] = threshold
        counts_across_all_chrs["Sample"] = counts_across_all_chrs.index
        results.append(counts_across_all_chrs[["Sample", "Threshold", "Precision", "Recall"]])

    results = pd.concat(results, ignore_index=True)

    # Plot PR curves: one line per sample
    fig = px.line(
        results,
        x="Recall",
        y="Precision",
        color="Sample",
        markers=True,
        # text="Threshold",
        title="Precision-Recall Curves per Sample",
    )
    fig.update_traces(textposition="top center")
    fig.update_layout(xaxis_title="Recall", yaxis_title="Precision")
    fig.write_image(output_file)
    return


def create_pr_curve_chromosomes(input_dir, intro_type, how_to_score, thresholds):
    scored_dir_name = f"scored_{how_to_score}"
    output_file = f"{input_dir}/{how_to_score}_{intro_type}_prcc.png"

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
    # create a collage of heatmaps; each page should have a chromosome through each threshold
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
    print("Visualizing results...")
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
        "--how",
        type=str,
        help="whether bins or overlaps were used during scoring",
    )
    args = parser.parse_args()

    input_dir = Path(args.dir)
    if not input_dir.is_dir():
        raise ValueError("Input directory not found. Check --dir path.")

    vis_functions = args.v
    for func in vis_functions:
        if func not in ["htmp", "shtmp", "prc", "prcc", "prca"]:
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
        if func == "prca":
            for intro_type in intro_types:
                create_pr_curve_accessions(input_dir, intro_type, how_to_score, thresholds)
        elif func == "shtmp":
            for intro_type in intro_types:
                create_scored_heatmap_collage(input_dir, intro_type, how_to_score, thresholds)
        elif func == "htmp":
            create_heatmap_runner(input_dir)
    print("Done.")
    return


if __name__ == "__main__":
    main()

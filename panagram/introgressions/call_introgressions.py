from pathlib import Path
import argparse
import numpy as np
from scipy.ndimage import uniform_filter1d, median_filter
import pandas as pd
from panagram.index import Index
import plotly.express as px
from concurrent.futures import ProcessPoolExecutor, as_completed


def bitmap_to_bins(
    bitmap,
    binlen,
    omit_fixed_kmers=False,
    omit_unique_kmers=False,
    ref_genome_name=None,
    outgroup_accessions=None,
):
    """Modded version of the same function found in the Index class. Converts a bitmap to binned
    kmer similarities on a scale of 0 to 1.

    Args:
        bitmap (pd.DataFrame): the bitmap to convert
        binlen (int): the length of the bins to create
        omit_fixed_kmers (bool, optional): whether to omit fixed kmers, defaults to False
        omit_unique_kmers (bool, optional): whether to omit unique kmers (those not present in neither a reference nor outgroup), defaults to False
        ref_genome_name (str, optional): the name of the reference genome, defaults to None
    outgroup_accessions (list, optional): the accessions of the introgression donor genomes, defaults to None

    Returns:
        pd.DataFrame: the binned kmer similarities
    """

    # change index from chr position to the number of the bin the position falls in
    # every column is the name of the accession
    # every row is an anchor kmer's presence/absence in other accessions
    bitmap_with_bin_idx = bitmap.set_index(bitmap.index // binlen)

    if omit_unique_kmers:
        # omit kmers that aren't present in the reference or the out groups
        kmers_to_keep_columns = outgroup_accessions
        kmers_to_keep_columns.append(ref_genome_name)

        # find kmers that aren't present in the reference or the wild accessions and set their reference column to 1
        bitmap_with_bin_idx.loc[
            bitmap_with_bin_idx[kmers_to_keep_columns].sum(axis=1) == 0, ref_genome_name
        ] = 1

    all_bins = bitmap_with_bin_idx.index.get_level_values(0).unique()

    # remove fixed kmers shared by all members of the pangenome (i.e., rows that are all 1)
    if omit_fixed_kmers:
        bitmap_with_bin_idx = bitmap_with_bin_idx.loc[~(bitmap_with_bin_idx == 1).all(axis=1)]

    # get the sum of how many kmers equal 1 in each bin
    binned_bitmap = bitmap_with_bin_idx.groupby(level=0).sum()

    # reindex to include bins that may now be empty
    binned_bitmap = binned_bitmap.reindex(all_bins, fill_value=1)

    # fix the index's bin number to be the bin's chr position and transpose
    binned_bitmap = binned_bitmap.set_index(binned_bitmap.index * binlen).T

    # divide each kmer sum in each column by the max # of kmers in the bin (some bins are uneven)
    binned_bitmap = binned_bitmap.div(binned_bitmap.max(axis=0), axis=1)
    return binned_bitmap


def row_trimmed_mean(row, trim_std):
    """Calculate the trimmed mean of a row. Note that setting the trim_std to anything above 2 has
    very little effect on the final results, as most values are within 2 standard deviations of
    the mean.

    Args:
        row (pd.Series): the row to calculate the trimmed mean for
        trim_std (float): the number of standard deviations to trim; if -1, return the untrimmed mean

    Returns:
        float: the trimmed mean of the row
    """

    mean = row.mean()
    std = row.std()
    if trim_std == -1:
        # If trim_std is -1, return the untrimmed mean
        return mean
    # Trim the row to values within trim_std standard deviations from the mean
    trimmed = row[(row >= mean - trim_std * std) & (row <= mean + trim_std * std)]
    return trimmed.mean()


def get_genome_similarities(
    genome,
    bitmap_step,
    bin_size,
    omit_fixed_kmers,
    omit_unique_for,
    ref_genome_name,
    outgroup_accessions,
    trim_std,
):
    """For a given anchor, get its mean similarity to all other accessions across the genome. Note
    that this function calculates a trimmed mean for each accession. Setting trim_std to -1 will
    return the untrimmed mean.

    Args:
        genome (panagram.index.Genome): the genome to anchor on and calculate similarities for
        bitmap_step (int): step size for sampling kmers in the bitmap (e.g., every 100th kmer)
        bin_size (int): size of the bins to create
        omit_fixed_kmers (bool): whether to omit fixed kmers
        omit_unique_for (list): list of accessions to omit unique kmers for
        ref_genome_name (str): name of the reference genome
        outgroup_accessions (list): list of outgroup accessions
        trim_std (float): number of standard deviations to trim; if -1, return the untrimmed mean

    Returns:
        pd.Series: mean similarity to all other accessions
    """

    chromosomes = genome.sizes.keys()

    # Collect all bin values for each accession across all chromosomes
    all_bins = []

    for chr_name in chromosomes:
        chr_size = genome.sizes[chr_name]
        chr_bitmap = genome.query(chr_name, 0, chr_size, step=bitmap_step)
        binned_bitmap = bitmap_to_bins(
            chr_bitmap,
            bin_size,
            omit_fixed_kmers,
            omit_unique_for,
            ref_genome_name,
            outgroup_accessions,
        )
        all_bins.append(binned_bitmap)

    # Concatenate all bins along columns (axis=1)
    all_bins_df = pd.concat(all_bins, axis=1)

    # For each accession (row), calculate the trimmed mean across all bins
    trimmed_means = all_bins_df.apply(row_trimmed_mean, trim_std=trim_std, axis=1)
    return trimmed_means


def smooth_row(row, filter_type, filter_size):
    """Smooth the row using the specified filter.

    Args:
        row (pd.Series): the row to smooth
        filter_type (str): the type of filter to use ('mean' or 'median')
        filter_size (int): the size of the filter (number of bins)

    Returns:
        pd.Series: the smoothed row
    """

    if filter_type == "mean":
        smoothed = uniform_filter1d(row.values, size=filter_size)
    elif filter_type == "median":
        smoothed = median_filter(row.values, size=filter_size)
    return pd.Series(smoothed, index=row.index)


def edge_tapered_row_normalization(df, intensity=0.1):
    """
    Normalize each row of the DataFrame more heavily at the beginning and end,
    decreasing intensity towards the center using a Gaussian window.

    Args:
        df (pd.DataFrame): DataFrame to normalize
        intensity (float, optional): intensity of the edge tapering effect, defaults to 0.1

    Returns:
        pd.DataFrame: normalized DataFrame with edge tapering
    """

    n_cols = df.shape[1]
    x = np.linspace(-1, 1, n_cols)
    window = np.exp(-4 * x**2)

    center_boost = intensity * (window / window.max())  # high in center

    norm_df = df.copy()
    for idx, row in norm_df.iterrows():
        norm_df.loc[idx] = row * (1 + center_boost)  # boost center more
    norm_df = norm_df.clip(0, 1)
    # increase intensity of all values to compensate for losses
    norm_df = norm_df.where(norm_df == 1, norm_df - 0.2)
    norm_df = norm_df.clip(0, 1)  # ensure values are between 0 and 1
    return norm_df


def preprocess_binned_bitmap(
    binned_bitmap,
    genome_similarities,
    similarity_normalization_mean,
    smoothing_filter,
    smoothing_filter_size,
    edge_normalization,
):
    """Preprocess the binned bitmapto make introgression calling easier/more consistent between
    accessions in a pangenome.

    Args:
        binned_bitmap (pd.DataFrame): binned k-mer similarities (rows are accessions, columns are genomic bins)
        genome_similarities (pd.Series): Series containing genome similarities
        similarity_normalization_mean (float): target mean for similarity normalization
        smoothing_filter (str): type of smoothing filter to apply
        smoothing_filter_size (int): size of the smoothing filter
        edge_normalization (bool): whether to apply edge normalization
    Returns:
        pd.DataFrame: Preprocessed DataFrame
    """

    binned_bitmap = binned_bitmap.copy()
    binned_bitmap = binned_bitmap.round(2)
    if genome_similarities is not None:
        if similarity_normalization_mean == -1:
            # if no target mean is provided, use the max average kmer sim. as target mean
            similarity_normalization_mean = genome_similarities[genome_similarities != 1].max()
        delta = similarity_normalization_mean - genome_similarities

        # Add delta to each row only where binned bitmap < 1
        for idx in binned_bitmap.index:
            row = binned_bitmap.loc[idx]
            mask = row <= 0.98  # don't change bins that are already 1
            row[mask] += delta[idx]
            binned_bitmap.loc[idx] = row.clip(0, 1)

    # TODO: allow edge normalization to be parameterized
    if edge_normalization:
        binned_bitmap = edge_tapered_row_normalization(binned_bitmap)

    if smoothing_filter:
        binned_bitmap = binned_bitmap.apply(
            smooth_row, axis=1, filter_type=smoothing_filter, filter_size=smoothing_filter_size
        )

    return binned_bitmap


def threshold_introgressions(binned_bitmap, anchor, comp_group, threshold):
    """Threshold binned bitmap based on k-mer similarity. Used for threshold logic that requires
    multiple rows from the bitmap (e.g., comparing similarity between the anchor's group and a
    comparison group). All bins are labeled as introgressed (1) or not introgressed (0) based on the
    following criteria:
    2-way comparison: If "REF" is the comp_group, introgression is where anchor has similarity lower
    than threshold to REF group
    3-way comparison: For all other comp_groups, introgression is where anchor is more similar to
    comp_group than REF by at least the given threshold

    Args:
        binned_bitmap (pd.DataFrame): the binned k-mer similarities
        anchor (str): the anchor accession
        comp_group (str): the comparison group
        threshold (float): the similarity threshold for introgression

    Returns:
        pd.DataFrame: DataFrame with introgression calls
    """

    # get the group the anchor belongs to
    anchor_group = binned_bitmap.loc[anchor, "group"]

    # define df with group info and kmer sims of anchor's group; drop anchor's self-similarity
    binned_bitmap_only_anchor_group = (
        binned_bitmap[binned_bitmap["group"] == anchor_group]
        .drop(columns=["group"])
        .drop(anchor, axis=0)
    )

    # define df with group info and kmer sim for comparison group
    binned_bitmap_only_comp_group = binned_bitmap[binned_bitmap["group"] == comp_group].drop(
        columns=["group"]
    )

    # get mean/max kmer similarities per window for each group
    group_sims = binned_bitmap_only_anchor_group.mean(axis=0).to_frame(name="anchor_sim")
    group_sims["comp_sim"] = binned_bitmap_only_comp_group.max(axis=0)

    # apply thresholding logic
    if comp_group == "REF":  # comparing to reference group
        group_sims["introgression"] = (group_sims.comp_sim < threshold).astype(int)
    else:  # comparing to different/wild group along with reference group
        binned_bitmap_only_ref = binned_bitmap[binned_bitmap["group"] == "REF"].drop(
            columns=["group"]
        )
        group_sims["ref_sim"] = binned_bitmap_only_ref.mean(axis=0)
        group_sims["introgression"] = (
            (group_sims.ref_sim < 0.95) & (group_sims.comp_sim >= group_sims.ref_sim + threshold)
        ).astype(int)
    return group_sims


def threshold_introgressions_simple(binned_bitmap, anchor, threshold):
    """Apply a simple threshold to introgressions. All bins are labeled as introgressed (1) or not
    introgressed (0) based on whether the anchor's similarity is below the threshold.

    Args:
        binned_bitmap (pd.DataFrame): the binned k-mer similarities
        anchor (str): the anchor accession
        threshold (float): the similarity threshold for introgression

    Returns:
        pd.DataFrame: DataFrame with introgression calls
    """

    group_sims = binned_bitmap.drop(columns=["group"]).loc[anchor].to_frame(name="anchor_sim")
    group_sims["comp_sim"] = pd.NA
    group_sims["introgression"] = (group_sims.anchor_sim < threshold).astype(int)
    return group_sims


def bins_to_bed(bins_df, bin_size, chr_name, comp_group):
    """Takes df organized into rows of bins with an introgression column and converts to bed file
    format. Adjacent introgressed bins are merged into a single introgression.

    Args:
        bins_df (pd.DataFrame): DataFrame containing binned introgression data
        bin_size (int): Size of each bin
        chr_name (str): Name of the chromosome
        comp_group (str): Name of the comparison group

    Returns:
        pd.DataFrame: DataFrame in BED format
    """

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
    introgressions["end"] = introgressions.start + (introgressions.end * bin_size) - 1

    # put into bed file format (chr, start, end, name)
    introgressions["chr"] = chr_name
    introgressions["name"] = f"{comp_group}_intro"
    introgressions = introgressions[["chr", "start", "end", "name"]]
    return introgressions


def visualize(
    binned_bitmap, output_file, inverse=False, title=None, groups=None, xaxis_dtick=2000000
):
    """Visualize the binned k-mer similarities as a heatmap.
    Args:
        binned_bitmap (pd.DataFrame): the binned k-mer similarities
        output_file (str): the output file path
        inverse (bool, optional): whether to invert the color scale, defaults to False
        title (str, optional): the title of the plot, defaults to None
        groups (pd.Series, optional): accession groups, used to determine ordering of heatmap rows, defaults to None
        xaxis_dtick (int, optional): the x-axis tick interval, defaults to 2000000
    """

    # if groups is not None, sort by group; groups should be in the order that you want them to appear
    if groups is not None:
        ordered_names = groups.index.tolist()
        intros = None
        if "Introgressions" in binned_bitmap.index:
            intros = binned_bitmap.loc["Introgressions"].copy()
            binned_bitmap = binned_bitmap.drop("Introgressions")
        binned_bitmap = binned_bitmap.reindex(index=ordered_names)
        if intros is not None:
            binned_bitmap.loc["Introgressions"] = intros

    if inverse:
        fig = px.imshow(
            binned_bitmap,
            color_continuous_scale=px.colors.sequential.Plasma[::-1],
            x=binned_bitmap.columns,
            y=binned_bitmap.index,
            aspect="auto",
            zmin=0,
            zmax=1,
        )
    else:
        fig = px.imshow(
            binned_bitmap,
            x=binned_bitmap.columns,
            y=binned_bitmap.index,
            aspect="auto",
            title=title,
            zmin=0,
            zmax=1,
        )

    genome_size = binned_bitmap.columns[-1]
    height = 12 * binned_bitmap.shape[0] + 200
    width = int(700 * (genome_size / 100_000_000))
    height = max(height, 500)
    width = max(width, 500)

    fig.update_layout(
        plot_bgcolor="white",
        margin=dict(l=0, r=0, t=100, b=80),
        autosize=True,
        font=dict(family="Arial", color="black"),
        xaxis=dict(
            title=dict(text="Position (Bp)", font=dict(size=12)),
            dtick=xaxis_dtick,
            tickangle=270,
            tickfont=dict(size=12),
        ),
        yaxis=dict(
            title="",  # Omit y-axis label
            tickfont=dict(size=12),
            tickmode="linear",
        ),
        title=dict(
            text=title,
            font=dict(size=14),
            x=0.5,
            xanchor="center",
        ),
        coloraxis_colorbar=dict(
            title=dict(
                text="Kmer Similarity",
                side="right",
                font=dict(size=12),
            ),
            tickfont=dict(size=12),
            outlinecolor="black",
            outlinewidth=1,
            thickness=15,
            # len=0.3,
            ticks="outside",
            tickcolor="black",
            dtick=0.2,
        ),
        width=width,
        height=height,
    )
    fig.write_image(output_file)
    return


def run_introgression_finder(
    anchor,
    genome,
    ref_genome,
    chr_name,
    groups,
    comp_groups,
    thresholds,
    bitmap_step,
    bin_size,
    using_ref_space,
    preprocessing_args,
    genome_similarities,
    ref_genome_similarities,
    render_vis,
    output_dir,
):
    """Run the introgression finder pipeline for a given anchor genome.

    Args:
        anchor (str): the anchor genome accession name
        genome (panagram.index.Genome): the anchor genome
        ref_genome (panagram.index.Genome): the reference genome (None if not using ref space)
        chr_name (str): the chromosome name
        groups (pd.Series): the groups for each accession
        comp_groups (list): the groups of accessions to compare against
        thresholds (list): the similarity thresholds for introgression detection
        bitmap_step (int): step size for sampling kmers in the bitmap (e.g., every 100th kmer)
        bin_size (int): the size of the bins for bitmap creation
        using_ref_space (bool): whether to use the reference genome space instead of accession space
        preprocessing_args (dict): arguments for preprocessing the bitmap
        genome_similarities (pd.DataFrame): the genome similarities (None if not used for preprocessing)
        ref_genome_similarities (pd.DataFrame): the reference genome similarities (None if not used for preprocessing)
        render_vis (bool): whether to render visualizations
        output_dir (Path): the output directory
    """

    # Step 1 - choose an anchor and re-create pairwise correlation matrix for it
    chr_size = genome.sizes[chr_name]
    chr_bitmap = genome.query(chr_name, 0, chr_size, step=bitmap_step)
    omit_fixed_kmers = preprocessing_args.pop("omit_fixed_kmers")
    omit_unique_kmers = preprocessing_args.pop("omit_unique_kmers")
    ref_genome_name = preprocessing_args.pop("ref_genome_name")
    outgroup_accessions = preprocessing_args.pop("outgroup_accessions")

    if using_ref_space:
        # get reference bitmap
        ref_chr_size = ref_genome.sizes[chr_name]
        ref_chr_bitmap = ref_genome.query(chr_name, 0, ref_chr_size, step=bitmap_step)
        # Note omit unique kmers doesn't make sense here since we are using REF's view
        ref_binned_bitmap = bitmap_to_bins(ref_chr_bitmap, bin_size, omit_fixed_kmers)

        # preprocess
        ref_binned_bitmap = preprocess_binned_bitmap(
            ref_binned_bitmap, ref_genome_similarities, **preprocessing_args
        )

        ref_binned_bitmap = ref_binned_bitmap.merge(
            groups, left_index=True, right_index=True, how="left"
        )
        bitmap_for_visualization = ref_binned_bitmap

    else:
        # get bitmap for anchor
        binned_bitmap = bitmap_to_bins(
            chr_bitmap,
            bin_size,
            omit_fixed_kmers,
            omit_unique_kmers,
            ref_genome_name,
            outgroup_accessions,
        )

        # preprocess
        binned_bitmap = preprocess_binned_bitmap(
            binned_bitmap, genome_similarities, **preprocessing_args
        )

        # get group information for the anchor's group and the comparison group
        binned_bitmap = binned_bitmap.merge(groups, left_index=True, right_index=True, how="left")
        bitmap_for_visualization = binned_bitmap

    # Step 2 - identify potential introgressions by thresholding kmer similarities
    for threshold in thresholds:
        merged_introgressions = None
        threshold_dir = output_dir / f"{output_dir.name}_{threshold}"
        threshold_dir.mkdir(parents=True, exist_ok=True)

        threshold_dir_raw = threshold_dir / "raw"
        threshold_dir_raw.mkdir(parents=True, exist_ok=True)

        if render_vis:
            output_dir_vis = threshold_dir / "heatmaps"
            output_dir_vis.mkdir(parents=True, exist_ok=True)

        for comp_group in comp_groups:
            if using_ref_space:
                # get introgressions by looking at differences with REF
                introgressions = threshold_introgressions_simple(
                    ref_binned_bitmap, anchor, threshold
                )
            else:
                # get introgressions by looking at REF and comparison groups if applicable
                introgressions = threshold_introgressions(
                    binned_bitmap, anchor, comp_group, threshold
                )
                # change comp_group to REFA (REF in accession space) for vis and for liftover purposes
                if comp_group == "REF":
                    comp_group = "REFA"

            # merge introgressions from different comparison groups when not using REF
            if len(comp_groups) > 1:
                if merged_introgressions is None:
                    merged_introgressions = introgressions
                else:
                    merged_introgressions += introgressions

            # Step 3 - visualization
            # show user introgressions labeled on the original heatmap from panagram
            if render_vis:
                # invert values for figure so that introgressions = 0
                bitmap_for_visualization.loc["Introgressions"] = (
                    ~(introgressions["introgression"].astype(bool))
                ).astype(int)
                output_vis = (
                    threshold_dir / "heatmaps" / f"{anchor}_{chr_name}_{comp_group}_heatmap.svg"
                )
                title = f"{anchor} {chr_name} Introgressions Called with {comp_group}"
                visualize(
                    bitmap_for_visualization.drop(columns=["group"]),
                    output_vis,
                    inverse=True,
                    title=title,
                    groups=groups,
                )

            # Step 4 - save introgressions to bed file
            output_bed = threshold_dir / "raw" / f"{anchor}_{chr_name}_{comp_group}.bed"
            introgressions = bins_to_bed(introgressions, bin_size, chr_name, comp_group)
            introgressions.to_csv(output_bed, header=False, index=False, sep="\t")

        # Step 5 - repeat analysis one more time with merged introgressions from all comp_groups
        if merged_introgressions is not None:
            if render_vis:
                # divide by the max to get between 0 and 1, do 1 - x to invert
                bitmap_for_visualization.loc["Introgressions"] = 1 - (
                    merged_introgressions["introgression"]
                    / merged_introgressions["introgression"].max()
                )
                output_vis = threshold_dir / "heatmaps" / f"{anchor}_{chr_name}_merged_heatmap.svg"
                title = f"{anchor} {chr_name} Merged Introgressions"
                visualize(
                    bitmap_for_visualization.drop(columns=["group"]),
                    output_vis,
                    inverse=True,
                    title=title,
                    groups=groups,
                )

            # save introgressions to bed file
            output_bed = threshold_dir / "raw" / f"{anchor}_{chr_name}_merged.bed"
            introgressions = bins_to_bed(merged_introgressions, bin_size, chr_name, "merged")
            introgressions.to_csv(output_bed, header=False, index=False, sep="\t")
    return


def run_introgression_finder_worker(
    anchor,
    ref_genome_name,
    chr_name,
    comp_groups,
    thresholds,
    bitmap_step,
    bin_size,
    using_ref_space,
    preprocessing_args,
    genome_similarities,
    ref_genome_similarities,
    render_vis,
    output_dir,
    index_dir,
    group_tsv_path,
):
    """Worker function for parallel chromosome processing. Creates large objects inside the process
    to avoid pickling issues when directly using run_introgression_finder with multiprocessing.

    Args:
        anchor (str): the anchor genome accession name
        ref_genome (panagram.index.Genome): the reference genome (None if not using ref space)
        chr_name (str): the chromosome name
        comp_groups (list): the groups of accessions to compare against
        thresholds (list): the similarity thresholds for introgression detection
        bitmap_step (int): step size for sampling kmers in the bitmap (e.g., every 100th kmer)
        bin_size (int): size of bitmap bin in bases
        using_ref_space (bool): whether to use the reference genome space instead of accession space
        preprocessing_args (dict): arguments for preprocessing the bitmap
        genome_similarities (pd.DataFrame): the genome similarities (None if not used for preprocessing)
        ref_genome_similarities (pd.DataFrame): the reference genome similarities (None if not used for preprocessing)
        render_vis (bool): whether to render visualizations
        output_dir (Path): the output directory for results
        index_dir (Path): the directory containing the genome index
        group_tsv_path (Path): the path to the TSV file containing group information
    """

    # Reconstruct index, genome, ref_genome, and groups inside the worker
    index = Index(index_dir)
    genome = index.genomes[anchor]
    ref_genome = index.genomes[ref_genome_name] if ref_genome_name else None
    groups = pd.read_csv(group_tsv_path, sep="\t", index_col=0)

    # Call the main finder function
    run_introgression_finder(
        anchor,
        genome,
        ref_genome,
        chr_name,
        groups,
        comp_groups,
        thresholds,
        bitmap_step,
        bin_size,
        using_ref_space,
        preprocessing_args,
        genome_similarities,
        ref_genome_similarities,
        render_vis,
        output_dir,
    )
    return


def main():
    parser = argparse.ArgumentParser(description="Introgression highlighter tool.")
    parser.add_argument("--threads", type=int, default=1, help="number of threads to use")
    parser.add_argument("--stp", type=int, default=100, help="bitmap kmer step size")
    parser.add_argument("--bin", type=int, default=1000000, help="size of bitmap bin in bases")
    parser.add_argument(
        "--gnm",
        type=float,
        help="target mean kmer similarity to normalize each genome to (set to -1 to use max average kmer sim. as target mean)",
    )
    parser.add_argument(
        "--trm",
        type=float,
        default=3.0,
        help="Number of standard deviations for trimmed mean normalization (default: 3.0)",
    )
    parser.add_argument("--sft", type=str, help="filter type for smoothing (mean or median)")
    parser.add_argument(
        "--ssz",
        type=int,
        default=5,
        help="filter size for smoothing",
    )
    parser.add_argument(
        "--edg",
        action="store_true",
        help="perform edge tapered normalization on binned bitmap",
    )
    parser.add_argument("--rmf", action="store_true", help="remove fixed kmers from bitmap")
    parser.add_argument("--rmu", nargs="+", help="remove unique kmers from given bitmaps")
    parser.add_argument("--ogrp", nargs="+", help="group(s) to use as outgroup when using --rmu")
    parser.add_argument("--vis", action="store_true", help="save svgs of visualized results")
    parser.add_argument(
        "--urf", action="store_true", help="when using REF as comp group, use the reference's view"
    )
    parser.add_argument("--ref", type=str, help="name of reference genome if using --rmu or --urf")
    parser.add_argument(
        "--anc", nargs="+", help="name of anchor(s) to mark introgressions for (default: all)"
    )
    parser.add_argument(
        "--grp",
        nargs="+",
        help="if --anc is not defined, groups in the tsv to mark introgressions for (default: none)",
    )
    parser.add_argument(
        "--chr", nargs="+", help="chromosome(s) to mark introgressions for (default: all)"
    )
    parser.add_argument(
        "--cmp", nargs="+", help="group(s) to compare against anchor(s)", required=True
    )
    parser.add_argument(
        "--thr",
        type=float,
        nargs="+",
        help="threshold for 2-way or 3-way introgression calling (see threshold_introgressions function for logic)",
    )
    parser.add_argument("--idx", type=str, help="path to Panagram index folder", required=True)
    parser.add_argument("--tsv", type=str, help="path to accession group TSV file", required=True)
    parser.add_argument("--out", type=str, help="path to folder to save all outputs", required=True)
    args = parser.parse_args()

    # input checking
    bitmap_step = args.stp
    bin_size = args.bin
    thresholds = args.thr
    omit_fixed_kmers = args.rmf
    omit_unique_for = args.rmu
    outgroups = args.ogrp
    render_vis = args.vis
    n_threads = args.threads
    using_ref_space = args.urf

    if not thresholds:
        raise ValueError("At least one threshold must be provided with --thr.")

    group_tsv = Path(args.tsv)
    if not group_tsv.is_file():
        raise ValueError("TSV file not found. Check --tsv path.")

    index_dir = Path(args.idx)
    if not index_dir.is_dir():
        raise ValueError("Index directory not found. Check --idx path.")
    index = Index(index_dir)
    groups = pd.read_csv(group_tsv, sep="\t", index_col=0)

    outgroup_accessions = []
    if omit_unique_for is not None:
        if args.ref is None:
            raise ValueError("Reference genome must be provided using --ref when using --rmu.")
        # get list of outgroup accessions
        outgroup_accessions = groups[groups.group.isin(outgroups)].index.tolist()
        if "REF" in outgroup_accessions:
            raise ValueError(
                "REF cannot be used as an outgroup accession. Please remove REF from the groups specified in --ogrp."
            )
        if any(acc in outgroup_accessions for acc in omit_unique_for):
            raise ValueError(
                "Accessions specified in --rmu cannot be in the outgroup. Please remove any accessions specified in --rmu from the groups specified in --ogrp."
            )

    # get preprocessing arguments
    preprocessing_args_base = {}
    preprocessing_args_base["similarity_normalization_mean"] = args.gnm
    if args.sft not in (None, "mean", "median"):
        raise ValueError("Invalid smoothing filter selected. Can be mean, median, or None.")
    preprocessing_args_base["smoothing_filter"] = args.sft
    preprocessing_args_base["smoothing_filter_size"] = args.ssz
    preprocessing_args_base["edge_normalization"] = args.edg
    preprocessing_args_base["omit_fixed_kmers"] = omit_fixed_kmers
    trim_std = args.trm

    # determine if we have a single anchor or multiple to run
    anchors = args.anc
    if anchors is None:
        if args.grp is None:
            raise ValueError("No anchor selected. Use either --anc or --grp to specify anchors.")
        # get multiple anchors
        anchors = list(groups[groups.group.isin(args.grp)].index)
    else:
        # make sure anchors and groups are not both defined
        if args.grp is not None:
            raise ValueError("Cannot use both --anc and --grp. Use one or the other.")

    comp_groups = args.cmp
    # ensure every element is only in there once
    comp_groups = list(set(comp_groups))

    if "REF" in comp_groups and comp_groups != ["REF"]:
        raise ValueError(
            "Error: REF must be the only comparison group specified so that a 2-way comparison can be run."
        )

    # set up ref genome if using its view for intro calling
    ref_genome = None
    reference = args.ref

    ref_genome_similarities = None
    if using_ref_space:
        if comp_groups != ["REF"]:
            raise ValueError(
                "REF must be the only comparison group specified with --cmp if using --urf"
            )
        ref_genome = index.genomes[reference]
        if args.gnm:
            ref_genome_similarities = get_genome_similarities(
                ref_genome,
                bitmap_step,
                bin_size,
                omit_fixed_kmers,
                omit_unique_for=None,
                ref_genome_name=None,
                outgroup_accessions=None,
                trim_std=trim_std,
            )

    output_dir = Path(args.out)

    futures = []
    with ProcessPoolExecutor(max_workers=n_threads) as executor:
        for anchor in anchors:
            print("Now running introgression analysis for", anchor)
            anchor_group = groups.loc[anchor, "group"]

            # don't compare anchor against its own group
            loop_comp_groups = [group for group in comp_groups if group != anchor_group]
            if not loop_comp_groups:
                print(
                    f"Skipping {anchor}: no comparison groups left after removing {anchor_group}."
                )
                continue

            # determine whether or not omit_unique_kmers overrides using the reference view for this anchor
            preprocessing_args = dict(preprocessing_args_base)
            if omit_unique_for and (anchor in omit_unique_for):
                print("Note that this accession will output REFA files to allow rmu to run.")
                loop_using_ref_space = False
                preprocessing_args["omit_unique_kmers"] = True
                preprocessing_args["ref_genome_name"] = args.ref
                preprocessing_args["outgroup_accessions"] = outgroup_accessions
            else:
                loop_using_ref_space = using_ref_space
                preprocessing_args["omit_unique_kmers"] = False
                preprocessing_args["ref_genome_name"] = None
                preprocessing_args["outgroup_accessions"] = None

            genome = index.genomes[anchor]

            genome_similarities = None
            if args.gnm and not loop_using_ref_space:
                genome_similarities = get_genome_similarities(
                    genome,
                    bitmap_step,
                    bin_size,
                    omit_fixed_kmers,
                    preprocessing_args["omit_unique_kmers"],
                    preprocessing_args["ref_genome_name"],
                    preprocessing_args["outgroup_accessions"],
                    trim_std,
                )

            chromosomes = args.chr
            if chromosomes is None:
                chromosomes = list(genome.sizes.keys())

            for chr_name in chromosomes:
                futures.append(
                    executor.submit(
                        run_introgression_finder_worker,
                        anchor,
                        reference,
                        chr_name,
                        loop_comp_groups,
                        thresholds,
                        bitmap_step,
                        bin_size,
                        loop_using_ref_space,
                        preprocessing_args,
                        genome_similarities,
                        ref_genome_similarities,
                        render_vis,
                        output_dir,
                        index_dir,
                        group_tsv,
                    )
                )

        # Wait for all jobs across all anchors/chromosomes and raise exceptions if any
        for future in as_completed(futures):
            future.result()
    print("Done.")
    return


if __name__ == "__main__":
    main()

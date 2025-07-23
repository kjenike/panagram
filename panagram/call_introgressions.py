from pathlib import Path
import argparse
import numpy as np
from scipy.ndimage import uniform_filter1d, median_filter
import pandas as pd
from panagram.index import Index
import plotly.express as px


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


def row_trimmed_mean(row, trim_std):
    mean = row.mean()
    std = row.std()
    if trim_std == -1:
        # If trim_std is -1, return the untrimmed mean
        return mean
    # Trim the row to values within trim_std standard deviations from the mean
    trimmed = row[(row >= mean - trim_std * std) & (row <= mean + trim_std * std)]
    return trimmed.mean()


def get_genome_similarities(genome, bitmap_step, bin_size, omit_fixed_kmers, trim_std):
    # For a given anchor, get its mean similarity to all other accessions across the genome
    # Note that this function calculates a trimmed mean for each accession
    # Setting trim_std to -1 will return the untrimmed mean
    # Note that setting the trim_std to anything above 2 has very little effect
    # on the final results, as most values are within 2 standard deviations of the mean
    chromosomes = genome.sizes.keys()

    # Collect all bin values for each accession across all chromosomes
    all_bins = []

    for chr_name in chromosomes:
        chr_size = genome.sizes[chr_name]
        chr_bitmap = genome.query(chr_name, 0, chr_size, step=bitmap_step)
        pair = bitmap_to_bins(chr_bitmap, bin_size, omit_fixed_kmers)
        all_bins.append(pair)

    # Concatenate all bins along columns (axis=1)
    all_bins_df = pd.concat(all_bins, axis=1)

    # For each accession (row), calculate the trimmed mean across all bins
    trimmed_means = all_bins_df.apply(row_trimmed_mean, trim_std=trim_std, axis=1)
    return trimmed_means


def smooth_row(row, filter_type, filter_size):
    if filter_type == "mean":
        smoothed = uniform_filter1d(row.values, size=filter_size)
    elif filter_type == "median":
        smoothed = median_filter(row.values, size=filter_size)
    return pd.Series(smoothed, index=row.index)


def percentile_contrast_stretch(row, lower=2, upper=98):
    p_low = np.percentile(row, lower)
    p_high = np.percentile(row, upper)
    if p_high == p_low:
        return pd.Series(1, index=row.index)
    stretched = (row - p_low) / (p_high - p_low)
    return stretched.clip(0, 1)


def edge_tapered_row_normalization(pair, intensity=0.1):
    """
    Normalize each row of the DataFrame 'pair' more heavily at the beginning and end,
    decreasing intensity towards the center using a Gaussian window.
    """
    n_cols = pair.shape[1]
    x = np.linspace(-1, 1, n_cols)
    window = np.exp(-4 * x**2)
    center_boost = intensity * (window / window.max())  # high in center

    norm_pair = pair.copy()
    for idx, row in norm_pair.iterrows():
        norm_pair.loc[idx] = row * (1 + center_boost)  # boost center more
    norm_pair = norm_pair.clip(0, 1)
    # increase intensity of all values to compensate for losses
    norm_pair = norm_pair.where(norm_pair == 1, norm_pair - 0.2)
    norm_pair = norm_pair.clip(0, 1)  # ensure values are between 0 and 1
    return norm_pair


def preprocess_pair(
    pair,
    genome_similarities,
    similarity_normalization_mean,
    smoothing_filter,
    smoothing_filter_size,
    contrast_stretching,
):
    pair = pair.copy()
    if genome_similarities is not None:
        if similarity_normalization_mean == -1:
            # if no target mean is provided, use the max average kmer sim. as target mean
            similarity_normalization_mean = genome_similarities[genome_similarities != 1].max()
        delta = similarity_normalization_mean - genome_similarities
        # add delta to data and clip values that are out of bounds
        pair = pair.add(delta, axis=0).clip(0, 1)

    if contrast_stretching:
        pair = edge_tapered_row_normalization(pair)
        # pair = pair.apply(percentile_contrast_stretch, axis=1)

    if smoothing_filter:
        pair = pair.apply(
            smooth_row, axis=1, filter_type=smoothing_filter, filter_size=smoothing_filter_size
        )

    return pair


def threshold_introgressions(pair, anchor, comp_group, threshold, row_mean_threshold=0.75):
    # a version of threshold introgressions used when pair come's from the anchor's view
    # used for threshold logic that requires multiple rows from pair (e.g., comparing similarity between the anchor's group and a comparison group)
    # get the group the anchor belongs to
    anchor_group = pair.loc[anchor, "group"]

    # define df with group info and kmer sims of anchor's group; drop anchor's self-similarity
    pair_anchor_group = (
        pair[pair["group"] == anchor_group].drop(columns=["group"]).drop(anchor, axis=0)
    )

    # define df with group info and kmer sim for comparison group
    pair_comp_group = pair[pair["group"] == comp_group].drop(columns=["group"])

    # throw out all acessions with large global diffs w/ the anchor
    if len(pair_anchor_group) > 1 and anchor_group == comp_group:
        for acession, row in pair_anchor_group.iterrows():
            if row.mean() < row_mean_threshold:
                pair_anchor_group = pair_anchor_group[pair_anchor_group.index != acession]

        if pair_anchor_group.empty:
            print(
                f"Warning: accession {anchor} is not similar to any of its group members. Cannot use --isc."
            )
            group_sims = pair_anchor_group.mean(axis=0).to_frame(name="anchor_sim")
            group_sims["comp_sim"] = pd.NA
            group_sims["introgression"] = 0
            return group_sims

    # get mean/max kmer similarities per window for each group
    group_sims = pair_anchor_group.mean(axis=0).to_frame(name="anchor_sim")
    group_sims["comp_sim"] = pair_comp_group.max(axis=0)

    # 4 ways to define an introgression:
    # 1. If anchor is in comp_group: introgression is where anchor has similarity lower than threshold to its group
    # 2. If "REF" is the comp_group: introgression is where anchor has similarity lower than threshold to REF group
    # 3. If anchor isn't in comp_group: introgression is where anchor is more similar to comp_group and >= threshold
    # 4. Average of all introgression calls for comp_group when comp_group is a list
    # than its own group
    if anchor_group == comp_group:
        group_sims["introgression"] = (group_sims.anchor_sim < threshold).astype(int)
    elif comp_group == "REF":
        # comp_group should only consist of 1 REF
        group_sims["introgression"] = (group_sims.comp_sim < threshold).astype(int)
    else:
        group_sims["introgression"] = (group_sims.comp_sim >= (threshold)).astype(int)
        pair_ref_group = pair[pair["group"] == "REF"].drop(columns=["group"])
        group_sims["ref_sim"] = pair_ref_group.mean(axis=0)
        group_sims["introgression"] = (
            (group_sims.ref_sim < threshold) & (group_sims.comp_sim > group_sims.comp_sim.min())
        ).astype(int)
    return group_sims


def threshold_introgressions_simple(pair, anchor, threshold):
    # applys a threshold to anchor's row; matches format of
    # used when pair comes from REF's view
    group_sims = pair.drop(columns=["group"]).loc[anchor].to_frame(name="anchor_sim")
    group_sims["comp_sim"] = pd.NA
    group_sims["introgression"] = (group_sims.anchor_sim < threshold).astype(int)
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
    introgressions["end"] = introgressions.start + (introgressions.end * bin_size) - 1

    # put into bed file format (chr, start, end, name)
    introgressions["chr"] = chr_name
    introgressions["name"] = f"{comp_group}_intro"
    introgressions = introgressions[["chr", "start", "end", "name"]]
    return introgressions


def visualize(pair, output_file, inverse=False, title=None, groups=None):
    # take a look at what pair looks like after manipulation
    # if groups is not None, sort by group; groups should be in the order that you want them to appear
    if groups is not None:
        ordered_names = groups.index.tolist()
        intros = None
        if "Introgressions" in pair.index:
            intros = pair.loc["Introgressions"].copy()
            pair = pair.drop("Introgressions")
        pair = pair.reindex(index=ordered_names)
        if intros is not None:
            pair.loc["Introgressions"] = intros

    if inverse:
        fig = px.imshow(
            pair,
            color_continuous_scale=px.colors.sequential.Plasma[::-1],
            x=pair.columns,
            y=pair.index,
            aspect="auto",
        )
    else:
        fig = px.imshow(
            pair, x=pair.columns, y=pair.index, aspect="auto", title=title, zmin=0, zmax=1
        )
    fig.update_layout(
        font=dict(family="Helvetica Bold", color="black"),
        xaxis=dict(
            dtick=2000000,
            title=dict(
                text="Genomic Position",
                font=dict(size=16),
            ),
        ),
        yaxis=dict(
            title="",  # Omit y-axis label
        ),
        title=dict(
            text=title,
            x=0.5,  # Center the title
            xanchor="center",
            font=dict(size=20),
        ),
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
    threshold,
    bitmap_step,
    bin_size,
    omit_fixed_kmers,
    using_ref_space,
    preprocessing_args,
    genome_similarities,
    ref_genome_similarities,
    render_vis,
    output_dir,
):

    # Step 1 - choose an anchor and re-create pairwise correlation matrix for it
    # get an entire chr's bitmap
    chr_size = genome.sizes[chr_name]
    chr_bitmap = genome.query(chr_name, 0, chr_size, step=bitmap_step)
    pair = bitmap_to_bins(chr_bitmap, bin_size, omit_fixed_kmers)

    # preprocess bitmap
    pair = preprocess_pair(pair, genome_similarities, **preprocessing_args)

    # get group information for the anchor's group and the comparison group
    pair = pair.merge(groups, left_index=True, right_index=True, how="left")

    # Step 2 - identify potential introgressions by thresholding kmer similarities
    merged_sims = None
    for comp_group in comp_groups:
        ref_with_ref_space = (comp_group == "REF") and using_ref_space
        if ref_with_ref_space:
            # get reference bitmap
            ref_chr_size = ref_genome.sizes[chr_name]
            ref_chr_bitmap = ref_genome.query(chr_name, 0, ref_chr_size, step=bitmap_step)
            ref_pair = bitmap_to_bins(ref_chr_bitmap, bin_size, omit_fixed_kmers)

            # preprocess ref_pair in the same way as pair
            ref_pair = preprocess_pair(ref_pair, ref_genome_similarities, **preprocessing_args)

            ref_pair = ref_pair.merge(groups, left_index=True, right_index=True, how="left")
            pair_for_visualization = ref_pair

            # get introgressions by looking at differences with REF
            group_sims = threshold_introgressions_simple(ref_pair, anchor, threshold)

        else:
            # TODO: revisit this - hard-coding a threshold might not work well post-normalization
            # make this at least 0.85 when omit fixed is not active; its too low otherwise
            row_mean_threshold = 0.75
            if not omit_fixed_kmers:
                row_mean_threshold = 0.85
            group_sims = threshold_introgressions(
                pair, anchor, comp_group, threshold, row_mean_threshold
            )
            pair_for_visualization = pair

        # don't allow merged_sims to include results from REF when using_ref_space
        if len(comp_groups) > 1 and not ref_with_ref_space:
            if merged_sims is None:
                merged_sims = group_sims
            else:
                merged_sims += group_sims

        # Step 3 - visualization
        # show user introgressions labeled on the original heatmap from panagram
        # invert values for figure so that introgressions = 0
        pair_for_visualization.loc["Introgressions"] = (
            ~(group_sims["introgression"].astype(bool))
        ).astype(int)
        if render_vis:
            output_vis = output_dir / "heatmaps"
            output_vis.mkdir(parents=True, exist_ok=True)
            output_vis = output_vis / f"{anchor}_{chr_name}_{comp_group}_heatmap.svg"
            title = f"{anchor} {chr_name} Introgressions Called with {comp_group}"
            visualize(
                pair_for_visualization.drop(columns=["group"]),
                output_vis,
                inverse=True,
                title=title,
                groups=groups,
            )

        # Step 4 - save introgressions to bed file
        output_bed = output_dir / f"{anchor}_{chr_name}_{comp_group}.bed"
        introgressions = bins_to_bed(group_sims, bin_size, chr_name, comp_group)
        introgressions.to_csv(output_bed, header=False, index=False, sep="\t")

    # Step 5 - repeat analysis one more time with merged introgressions from all comp_groups
    if merged_sims is not None:
        # divide by the max to get between 0 and 1, do 1 - x to invert
        pair_for_visualization.loc["Introgressions"] = 1 - (
            merged_sims["introgression"] / merged_sims["introgression"].max()
        )

        if render_vis:
            output_vis = output_dir / "heatmaps" / f"{anchor}_{chr_name}_merged_heatmap.svg"
            title = f"{anchor} {chr_name} Merged Introgressions"
            visualize(
                pair_for_visualization.drop(columns=["group"]),
                output_vis,
                inverse=True,
                title=title,
                groups=groups,
            )

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
        "--thr",
        type=float,
        default=0.75,
        help="introgressions marked as regions with kmer similarity lower than this threshold",
    )
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
        "--cst",
        action="store_true",
        help="perform percentile-based contrast stretching; better highlights introgressions and normalizes kmer similarities",
    )
    parser.add_argument("--rmf", action="store_true", help="remove fixed kmers from bitmap")
    parser.add_argument("--vis", action="store_true", help="save svgs of visualized results")
    parser.add_argument("--isc", action="store_true", help="use anchor's group to self compare")
    parser.add_argument(
        "--urf", type=str, help="when calling in REF mode, use this accession's view"
    )
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
    parser.add_argument("--idx", type=str, help="path to Panagram index folder", required=True)
    parser.add_argument("--tsv", type=str, help="path to acession group TSV file", required=True)
    parser.add_argument("--out", type=str, help="path to folder to save all outputs", required=True)
    args = parser.parse_args()

    # input checking
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

    # get preprocessing arguments
    preprocessing_args = {}
    preprocessing_args["similarity_normalization_mean"] = args.gnm
    if args.sft not in (None, "mean", "median"):
        raise ValueError("Invalid smoothing filter selected. Can be mean, median, or None.")
    preprocessing_args["smoothing_filter"] = args.sft
    preprocessing_args["smoothing_filter_size"] = args.ssz
    preprocessing_args["contrast_stretching"] = args.cst
    trim_std = args.trm

    # determine if we have a single anchor or multiple to run
    anchors = args.anc
    if anchors is None:
        if args.grp is None:
            raise ValueError("No anchor selected. Use either -anc or -grp to specify anchors.")
        # get multiple anchors
        anchors = list(groups[groups.group.isin(args.grp)].index)

    comp_groups = args.cmp
    # ensure every element is only in there once
    comp_groups = list(set(comp_groups))

    # set up ref genome if using its view for intro calling
    ref_genome = None
    using_ref_space = False
    reference = args.urf

    ref_genome_similarities = None
    if reference is not None:
        if "REF" not in comp_groups:
            raise ValueError("REF must be included as a comparison group using -p if using --urf")
        ref_genome = index.genomes[reference]
        using_ref_space = True
        if args.gnm:
            ref_genome_similarities = get_genome_similarities(
                ref_genome, bitmap_step, bin_size, omit_fixed_kmers, trim_std
            )

    output_dir = Path(args.out)
    output_dir.mkdir(parents=True, exist_ok=True)

    for anchor in anchors:
        print("Now running introgression analysis for", anchor)
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

        genome_similarities = None
        if args.gnm:
            genome_similarities = get_genome_similarities(
                genome, bitmap_step, bin_size, omit_fixed_kmers, trim_std
            )

        chromosomes = args.chr
        if chromosomes is None:
            chromosomes = genome.sizes.keys()

        for chr_name in chromosomes:
            run_introgression_finder(
                anchor,
                genome,
                ref_genome,
                chr_name,
                groups,
                comp_groups,
                threshold,
                bitmap_step,
                bin_size,
                omit_fixed_kmers,
                using_ref_space,
                preprocessing_args,
                genome_similarities,
                ref_genome_similarities,
                render_vis,
                output_dir,
            )

    print("Done.")
    return


if __name__ == "__main__":
    main()

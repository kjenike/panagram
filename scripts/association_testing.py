import sys
import argparse
import subprocess
from pathlib import Path
from io import StringIO
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from Bio import Phylo
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import linkage
import statsmodels.api as sm
import seaborn as sns
from panagram.index import Index


def get_newick(node, parent_dist, leaf_names, newick=""):
    """
    Convert scipy.cluster.hierarchy.to_tree() output to Newick format.
    """

    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parent_dist - node.dist, newick)

    if len(newick) > 0:
        newick = "):%.2f%s" % (parent_dist - node.dist, newick)
    else:
        newick = ");"

    newick = get_newick(node.get_left(), node.dist, leaf_names, newick=newick)
    newick = get_newick(node.get_right(), node.dist, leaf_names, newick=",%s" % newick)
    newick = "(%s" % newick

    return newick


def get_major_clades(tree, min_size=2):
    """
    Divide root-level tips of a tree into major clades.
    """

    clades = {}
    counter = 1

    for clade in tree.root.clades:
        tips = clade.get_terminals()

        if len(tips) >= min_size:
            clade_name = f"Clade_{counter}"
            clades[clade_name] = [tip.name for tip in tips]
            counter += 1

    return clades


def extract_feature_name(attributes):
    """
    Extract the feature name from the first field of the GFF attributes column.

    NOTE: GFF names are not standardized, so this may need editing for some GFFs.
    """

    feature_name = attributes.split(";")[0]

    if "=" in feature_name:
        feature_name = feature_name.split("=")[1]
    else:
        feature_name = feature_name.split()[1]

    feature_name = feature_name.replace('"', "")

    return feature_name


def read_sample_order(sample_f):
    """
    Read sample names from the first column of a tab-delimited sample file.

    The first line is treated as a header and skipped.
    """

    sample_order = []

    with open(sample_f, "r") as f:
        line = f.readline()
        line = f.readline()

        while line:
            sample_order.append(line.strip().split("\t")[0])
            line = f.readline()

    return sample_order


def generate_regulatory_windows(
    gene_name,
    gene_start,
    gene_end,
    chrom_size,
    flank_size=500000,
    window_size=200,
    step_size=100,
):
    """
    Generate overlapping regulatory windows around one gene.
    """

    reg_start = max(gene_start - flank_size, 0)
    reg_end = min(gene_end + flank_size, chrom_size)

    windows = []
    window_counter = 1

    for window_start in range(reg_start, reg_end - window_size + 1, step_size):
        window_end = window_start + window_size
        window_name = f"{gene_name}_reg_{window_counter}"

        windows.append((window_name, window_start, window_end))
        window_counter += 1

    return windows


def generate_chromosome_windows(
    chromosome,
    chrom_size,
    window_size=200,
    step_size=100,
):
    """
    Generate overlapping windows across an entire chromosome/contig.
    """

    windows = []

    for window_start in range(0, chrom_size - window_size + 1, step_size):
        window_end = window_start + window_size
        window_name = f"{chromosome}_{window_start}"

        windows.append((window_name, chromosome, window_start, window_end))

    return windows


def build_bitmap_tree(anchor_name, panagram_dir, tree_file):
    """
    Build a genome-wide sample tree from Panagram bitmap data.
    """

    n_skips = 100
    max_sampled_rows = 2000000

    index = Index(panagram_dir)
    chroms = list(index[anchor_name].chrs.index)

    bitmap = []

    for chrom in chroms:
        start_coord = 0
        end_coord = index[anchor_name].chrs.loc[chrom, "size"]

        bittmp = index.query_bitmap(anchor_name, chrom, start_coord, end_coord, n_skips)

        print(f"Tree bitmap query complete: {chrom}", file=sys.stderr, flush=True)
        bitmap.append(bittmp)

    result = pd.concat(bitmap)
    result_sample = result.sample(n=min(len(result), max_sampled_rows), random_state=42)

    matrix = linkage(result_sample.transpose(), method="average")

    tree_tmp = hierarchy.to_tree(matrix, False)
    treedata = get_newick(tree_tmp, tree_tmp.dist, index.genome_names)

    with open(tree_file, "w") as f:
        f.write(treedata)

    print(f"Tree written to: {tree_file}", file=sys.stderr, flush=True)


def create_feature_tree(bitmap, genome_names):
    """
    Cluster samples by their bitmap profiles for one feature/window.
    """

    matrix = linkage(bitmap.transpose(), method="ward", metric="euclidean")

    tree_tmp = hierarchy.to_tree(matrix, False)
    treedata = get_newick(tree_tmp, tree_tmp.dist, genome_names)

    tree = Phylo.read(StringIO(treedata), "newick")

    clades = get_major_clades(tree, 1)

    return clades


def make_haplotype_series(clades, sample_order):
    """
    Convert clades into a per-sample haplotype series.
    """

    haplotypes = {}

    for sample in sample_order:
        if sample in clades.get("Clade_1", []):
            haplotypes[sample] = 0
        elif sample in clades.get("Clade_2", []):
            haplotypes[sample] = 1

    return pd.Series(haplotypes, name="haplotype")


def write_wide_haplotype_header(out_f, sample_order):
    """
    Write the wide-format haplotype TSV header.
    """

    out_f.write("feature_id\t" + "\t".join(sample_order) + "\n")


def write_wide_haplotype_row(out_f, feature_name, haplotype_series, sample_order):
    """
    Write one wide-format haplotype row.
    """

    values = []

    for sample in sample_order:
        if sample in haplotype_series.index:
            values.append(str(int(haplotype_series.loc[sample])))
        else:
            values.append("NA")

    out_f.write(feature_name + "\t" + "\t".join(values) + "\n")


def build_kinship_matrix(tree_file):
    """
    Build a tree-derived kinship-like covariance matrix.
    """

    tree = Phylo.read(tree_file, "newick")

    tips = tree.get_terminals()
    sample_names = [tip.name for tip in tips]
    max_samples = len(sample_names) + 3

    max_distance = max(tree.distance(tree.root, tip) for tip in tips)
    max_distance = max_distance * 2

    distance_matrix = pd.DataFrame(index=sample_names, columns=sample_names, dtype=float)

    for tip_a in tips:
        for tip_b in tips:
            distance_matrix.loc[tip_a.name, tip_b.name] = (
                max_distance - tree.distance(tip_a, tip_b)
            ) / max_distance

    kinship_matrix = distance_matrix.values.copy()
    kinship_matrix += np.eye(kinship_matrix.shape[0]) * 1e-6
    kinship_matrix = kinship_matrix / np.mean(np.diag(kinship_matrix))

    return kinship_matrix, distance_matrix, max_samples


def test_one_feature_immediate(
    feature_name,
    haplotype_series,
    kinship_matrix,
    distance_matrix,
    phenotype_df,
    phenotype_list,
    out_f,
):
    """
    Test one feature/window against all phenotypes using GLS.
    """

    sample_order = distance_matrix.index.tolist()

    shared_samples = [
        sample
        for sample in sample_order
        if sample in haplotype_series.index and sample in phenotype_df.index
    ]

    if len(shared_samples) == 0:
        return

    haplotypes = haplotype_series.loc[shared_samples].astype(float)
    phenotypes = phenotype_df.loc[shared_samples]

    group1_samples = haplotypes[haplotypes == 0].index
    group2_samples = haplotypes[haplotypes == 1].index

    group1_size = len(group1_samples)
    group2_size = len(group2_samples)

    if min(group1_size, group2_size) == 0:
        return

    shared_positions = [sample_order.index(sample) for sample in shared_samples]

    for phenotype_name in phenotype_list:
        phenotype_values_all = phenotypes[phenotype_name].values.astype(float)
        haplotype_values_all = haplotypes.values.astype(float)

        valid_indices = ~np.isnan(phenotype_values_all)

        if valid_indices.sum() == 0:
            continue

        phenotype_values = phenotype_values_all[valid_indices]
        haplotype_values_raw = haplotype_values_all[valid_indices]

        if len(np.unique(haplotype_values_raw)) < 2:
            continue

        group1_phenotype_values = phenotypes.loc[group1_samples, phenotype_name].dropna()
        group2_phenotype_values = phenotypes.loc[group2_samples, phenotype_name].dropna()

        if len(group1_phenotype_values) == 0 or len(group2_phenotype_values) == 0:
            continue

        haplotype_values = haplotype_values_raw.reshape(-1, 1)
        haplotype_values = sm.add_constant(haplotype_values)

        valid_positions = np.array(shared_positions)[valid_indices]
        kinship_matrix_subset = kinship_matrix[np.ix_(valid_positions, valid_positions)]

        model = sm.GLS(
            phenotype_values,
            haplotype_values,
            sigma=kinship_matrix_subset,
        )

        results = model.fit()

        group1_mean = group1_phenotype_values.mean()
        group2_mean = group2_phenotype_values.mean()

        out_f.write(
            feature_name
            + "\t"
            + phenotype_name
            + "\t"
            + str(group1_size)
            + "\t"
            + str(group2_size)
            + "\t"
            + str(group1_mean)
            + "\t"
            + str(group2_mean)
            + "\t"
            + str(results.pvalues[1])
            + "\n"
        )


def build_feature_windows(
    anchor_name,
    panagram_dir,
    gff,
    feature_type,
):
    """
    Build the full list of features/windows to test.
    """

    index = Index(panagram_dir)
    features = []

    if feature_type == "all":
        chroms = list(index[anchor_name].chrs.index)

        for chrom in chroms:
            chrom_size = index[anchor_name].chrs.loc[chrom, "size"]
            chrom_windows = generate_chromosome_windows(chrom, chrom_size)
            features.extend(chrom_windows)

        return features

    feature_counts = {}

    with open(gff, "r") as gff_f:
        line = gff_f.readline()

        while line:
            if line.startswith("#"):
                line = gff_f.readline()
                continue

            fields = line.strip().split("\t")

            if len(fields) < 9:
                line = gff_f.readline()
                continue

            chrom = fields[0]

            if feature_type == "reg":
                if fields[2] != "gene":
                    line = gff_f.readline()
                    continue
            else:
                if fields[2] != feature_type:
                    line = gff_f.readline()
                    continue

            start_coord = int(fields[3])
            end_coord = int(fields[4])
            feature_name = extract_feature_name(fields[8])

            if feature_type == "reg":
                chrom_size = index[anchor_name].chrs.loc[chrom, "size"]

                windows = generate_regulatory_windows(
                    feature_name,
                    start_coord,
                    end_coord,
                    chrom_size,
                )

                for window_name, window_start, window_end in windows:
                    overlaps_gene = window_start < end_coord and window_end > start_coord

                    if overlaps_gene:
                        continue

                    features.append((window_name, chrom, window_start, window_end))

            else:
                if feature_type == "exon":
                    feature_counts[feature_name] = feature_counts.get(feature_name, 0) + 1
                    feature_name = f"{feature_name}_{feature_counts[feature_name]}"

                features.append((feature_name, chrom, start_coord, end_coord))

            line = gff_f.readline()

    return features


def split_into_chunks(items, chunk_size):
    """
    Split a list into fixed-size chunks.
    """

    chunks = []

    for start in range(0, len(items), chunk_size):
        end = start + chunk_size
        chunks.append(items[start:end])

    return chunks


def process_feature_chunk(
    chunk_id,
    feature_chunk,
    anchor_name,
    panagram_dir,
    sample_f,
    phenotype_file,
    feature_type,
    tree_file,
    haplotype_dir,
    association_dir,
    keep_haplotypes,
):
    """
    Process one chunk of features/windows.
    """

    print(f"Starting chunk: {chunk_id}", file=sys.stderr, flush=True)

    index = Index(panagram_dir)
    sample_order = read_sample_order(sample_f)
    phenotype_df = pd.read_csv(phenotype_file, index_col=0)
    phenotype_list = phenotype_df.columns.values.tolist()

    kinship_matrix, distance_matrix, max_samples = build_kinship_matrix(tree_file)

    haplotype_tsv = haplotype_dir / f"chunk_{chunk_id}_{feature_type}_haplotypes.tsv"
    association_tsv = association_dir / f"chunk_{chunk_id}_{feature_type}_associations.tsv"

    cntr = 0

    if keep_haplotypes:
        hap_f = open(haplotype_tsv, "w")
        write_wide_haplotype_header(hap_f, sample_order)
    else:
        hap_f = None

    with open(association_tsv, "w") as assoc_f:
        for test_feature_name, chrom, test_start, test_end in feature_chunk:
            bitmap = index.query_bitmap(
                anchor_name,
                chrom,
                test_start,
                test_end,
                1,
            )

            clades = create_feature_tree(bitmap, index.genome_names)
            haplotype_series = make_haplotype_series(clades, sample_order)

            if keep_haplotypes:
                write_wide_haplotype_row(
                    hap_f,
                    test_feature_name,
                    haplotype_series,
                    sample_order,
                )

            test_one_feature_immediate(
                test_feature_name,
                haplotype_series,
                kinship_matrix,
                distance_matrix,
                phenotype_df,
                phenotype_list,
                assoc_f,
            )

            if cntr % 10 == 0:
                assoc_f.flush()

                if keep_haplotypes:
                    hap_f.flush()

            cntr += 1

    if keep_haplotypes:
        hap_f.close()

    print(f"Finished chunk: {chunk_id}", file=sys.stderr, flush=True)

    return chunk_id


def combine_association_files(association_dir, feature_type):
    """
    Combine chunk-level association files into one association TSV.
    """

    association_dir = Path(association_dir)
    combined_file = association_dir / f"all_{feature_type}_associations.tsv"

    association_files = sorted(association_dir.glob(f"*_{feature_type}_associations.tsv"))
    association_files = [f for f in association_files if not f.name.startswith("all_")]

    with open(combined_file, "w") as out_f:
        for association_file in association_files:
            print(f"Combining: {association_file}", file=sys.stderr, flush=True)

            with open(association_file, "r") as in_f:
                for line in in_f:
                    out_f.write(line)

    print(f"Combined associations: {combined_file}", file=sys.stderr, flush=True)

    return combined_file


def filter_associations_by_group_size(combined_file, filtered_file):
    """
    Apply the 10% haplotype group-size filter before FDR.

    This function does NOT replace p_value == 0.

    The filtered file keeps original unmodified association values, including:
        p_value = 0

    FDR correction is then run directly on these original p-values.
    """

    retained_rows = 0
    zero_pvalue_rows = 0

    with open(combined_file, "r") as in_f, open(filtered_file, "w") as out_f:
        for line in in_f:
            fields = line.rstrip("\n").split("\t")

            group1_size = float(fields[2])
            group2_size = float(fields[3])
            p_value = float(fields[6])

            total_size = group1_size + group2_size

            if group1_size < 0.1 * total_size or group2_size < 0.1 * total_size:
                continue

            retained_rows += 1

            if p_value == 0:
                zero_pvalue_rows += 1

            out_f.write("\t".join(fields) + "\n")

    if retained_rows == 0:
        sys.exit("ERROR: No rows remain after filtering by haplotype group size.")

    print(f"Rows after 10% group-size filter: {retained_rows}", file=sys.stderr, flush=True)
    print(
        f"Rows with p_value == 0 retained for FDR: {zero_pvalue_rows}", file=sys.stderr, flush=True
    )

    return retained_rows


def count_lines(input_file):
    """
    Count lines in a file.
    """

    count = 0

    with open(input_file, "r") as f:
        for _ in f:
            count += 1

    return count


def run_sort(input_file, output_file, sort_key, temp_dir=None):
    """
    Run Unix sort.
    """

    cmd = ["sort", "-t", "\t"]

    if temp_dir is not None:
        cmd.extend(["-T", str(temp_dir)])

    cmd.extend(sort_key)
    cmd.extend([str(input_file), "-o", str(output_file)])

    print("Running:", " ".join(cmd), file=sys.stderr, flush=True)
    subprocess.run(cmd, check=True)


def add_bh_rank_values(sorted_by_p_file, ranked_file, total_tests):
    """
    Add rank and raw BH values to p-value-sorted associations.

    If p_value == 0, raw_bh is also 0.
    """

    with open(sorted_by_p_file, "r") as in_f, open(ranked_file, "w") as out_f:
        rank = 1

        for line in in_f:
            line = line.rstrip("\n")
            fields = line.split("\t")

            p_value = float(fields[6])
            raw_bh = p_value * total_tests / rank

            out_f.write(line + "\t" + str(rank) + "\t" + str(raw_bh) + "\n")

            rank += 1


def apply_bh_cumulative_min(ranked_desc_file, fdr_unsorted_file):
    """
    Apply BH cumulative minimum while streaming from largest rank to smallest rank.
    """

    running_min = 1.0

    with open(ranked_desc_file, "r") as in_f, open(fdr_unsorted_file, "w") as out_f:
        for line in in_f:
            line = line.rstrip("\n")
            fields = line.split("\t")

            raw_bh = float(fields[8])
            running_min = min(running_min, raw_bh)
            fdr = min(running_min, 1.0)

            original_fields = fields[:7]
            out_f.write("\t".join(original_fields) + "\t" + str(fdr) + "\n")


def collect_phenotype_mean_diff_ranges(fdr_sorted_file):
    """
    Collect phenotype-specific min/max mean differences from corrected associations.
    """

    ranges = {}

    with open(fdr_sorted_file, "r") as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")

            phenotype_name = fields[1]
            group1_mean = float(fields[4])
            group2_mean = float(fields[5])
            mean_difference = group1_mean - group2_mean

            if phenotype_name not in ranges:
                ranges[phenotype_name] = [mean_difference, mean_difference]
            else:
                ranges[phenotype_name][0] = min(
                    ranges[phenotype_name][0],
                    mean_difference,
                )
                ranges[phenotype_name][1] = max(
                    ranges[phenotype_name][1],
                    mean_difference,
                )

    return ranges


def collect_top_features_per_phenotype(fdr_sorted_file, max_points_per_phenotype=50000):
    """
    Collect feature/phenotype pairs to plot.
    """

    counts = defaultdict(int)
    plotted_features = set()

    with open(fdr_sorted_file, "r") as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")

            feature_id = fields[0]
            phenotype_name = fields[1]

            if counts[phenotype_name] < max_points_per_phenotype:
                plotted_features.add((feature_id, phenotype_name))
                counts[phenotype_name] += 1

    return plotted_features


def add_norm_difference_and_plotted_flag(
    fdr_sorted_file,
    final_file,
    phenotype_ranges,
    plotted_features,
):
    """
    Add norm_mean_difference and plotted columns to the final corrected file.
    """

    header = [
        "feature_id",
        "phenotype_name",
        "group1_size",
        "group2_size",
        "group1_mean",
        "group2_mean",
        "p_value",
        "fdr",
        "norm_mean_difference",
        "plotted",
    ]

    with open(fdr_sorted_file, "r") as in_f, open(final_file, "w") as out_f:
        out_f.write("\t".join(header) + "\n")

        for line in in_f:
            fields = line.rstrip("\n").split("\t")

            feature_id = fields[0]
            phenotype_name = fields[1]
            group1_mean = float(fields[4])
            group2_mean = float(fields[5])

            mean_difference = group1_mean - group2_mean

            pheno_min, pheno_max = phenotype_ranges[phenotype_name]

            if pheno_max == pheno_min:
                norm_mean_difference = 0.0
            else:
                norm_mean_difference = mean_difference / (pheno_max - pheno_min)

            if (feature_id, phenotype_name) in plotted_features:
                plotted = "1"
            else:
                plotted = "0"

            out_f.write(
                "\t".join(fields) + "\t" + str(norm_mean_difference) + "\t" + plotted + "\n"
            )


def run_bh_fdr_disk(
    combined_association_file,
    association_dir,
    feature_type,
    sort_temp_dir=None,
    max_points_per_phenotype=50000,
):
    """
    Run exact global BH/FDR correction using disk-based sorting.

    Important:
        p_value == 0 values are NOT replaced before FDR.
        True zero p-values are used directly in BH/FDR.
        The final TSV preserves p_value == 0 and fdr == 0 where applicable.
    """

    association_dir = Path(association_dir)

    filtered_file = association_dir / f"all_{feature_type}_associations.filtered.tmp.tsv"
    sorted_by_p_file = association_dir / f"all_{feature_type}_associations.sorted_by_p.tmp.tsv"
    ranked_file = association_dir / f"all_{feature_type}_associations.ranked.tmp.tsv"
    ranked_desc_file = association_dir / f"all_{feature_type}_associations.ranked_desc.tmp.tsv"
    fdr_unsorted_file = association_dir / f"all_{feature_type}_associations.fdr_unsorted.tmp.tsv"
    fdr_sorted_no_header_file = (
        association_dir / f"all_{feature_type}_associations.fdr_sorted.no_header.tmp.tsv"
    )

    final_file = association_dir / f"all_{feature_type}_associations.fdr_sorted.tsv"

    total_tests = filter_associations_by_group_size(
        combined_association_file,
        filtered_file,
    )

    print(f"Total association tests used for FDR: {total_tests}", file=sys.stderr, flush=True)

    # Sort filtered associations by true p-value, then feature_id, then phenotype_name.
    # True p_value == 0 rows are kept and sorted first.
    run_sort(
        filtered_file,
        sorted_by_p_file,
        ["-k7,7g", "-k1,1", "-k2,2"],
        temp_dir=sort_temp_dir,
    )

    add_bh_rank_values(sorted_by_p_file, ranked_file, total_tests)

    run_sort(
        ranked_file,
        ranked_desc_file,
        ["-k8,8nr"],
        temp_dir=sort_temp_dir,
    )

    apply_bh_cumulative_min(ranked_desc_file, fdr_unsorted_file)

    run_sort(
        fdr_unsorted_file,
        fdr_sorted_no_header_file,
        ["-k8,8g"],
        temp_dir=sort_temp_dir,
    )

    phenotype_ranges = collect_phenotype_mean_diff_ranges(fdr_sorted_no_header_file)

    plotted_features = collect_top_features_per_phenotype(
        fdr_sorted_no_header_file,
        max_points_per_phenotype=max_points_per_phenotype,
    )

    add_norm_difference_and_plotted_flag(
        fdr_sorted_no_header_file,
        final_file,
        phenotype_ranges,
        plotted_features,
    )

    print(f"Final FDR-corrected associations: {final_file}", file=sys.stderr, flush=True)

    return final_file


def plot_volcano(corrected_association_file, output_plot):
    """
    Make volcano plots using only rows marked plotted == 1.

    True p_value == 0 rows are shown with star markers.

    Zero FDR values are replaced only in memory for plotting because:
        -log10(0)
    is not finite.

    The final TSV is not modified.
    """

    output_plot = Path(output_plot)
    output_dir = output_plot.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    chunks = []

    for chunk in pd.read_csv(corrected_association_file, sep="\t", chunksize=500000):
        chunk = chunk[chunk["plotted"] == 1]

        if len(chunk) > 0:
            chunks.append(chunk)

    if len(chunks) == 0:
        sys.exit("ERROR: No rows marked plotted == 1.")

    df = pd.concat(chunks, ignore_index=True)

    df["group1_mean"] = pd.to_numeric(df["group1_mean"])
    df["group2_mean"] = pd.to_numeric(df["group2_mean"])
    df["p_value"] = pd.to_numeric(df["p_value"])
    df["fdr"] = pd.to_numeric(df["fdr"])
    df["norm_mean_difference"] = pd.to_numeric(df["norm_mean_difference"])

    # Track original true zero p-values before any plotting-only replacement.
    df["original_p_value_zero"] = df["p_value"] == 0

    # Replace fdr == 0 only in memory for plotting.
    nonzero_fdr = df[df["fdr"] > 0]["fdr"].min()

    if pd.isna(nonzero_fdr):
        replacement_fdr = 1e-300
    else:
        exponent = int(np.floor(np.log10(nonzero_fdr))) - 2
        replacement_fdr = 10**exponent

    df.loc[df["fdr"] == 0, "fdr"] = replacement_fdr

    df["mean_difference"] = df["group1_mean"] - df["group2_mean"]
    df["log_fdr"] = -np.log10(df["fdr"])

    unique_phenotypes = df["phenotype_name"].unique()
    palette = sns.color_palette("husl", n_colors=len(unique_phenotypes))
    color_map = dict(zip(unique_phenotypes, palette))

    sig_threshold_fdr = 0.05
    sig_threshold_log_fdr = -np.log10(sig_threshold_fdr)

    groups = df.groupby("phenotype_name")

    # Combined volcano plot.
    plt.rcParams["figure.figsize"] = (4.5, 5)
    fig, ax = plt.subplots()

    for name, group in groups:
        nonzero_group = group[~group["original_p_value_zero"]]
        zero_group = group[group["original_p_value_zero"]]

        if len(nonzero_group) > 0:
            ax.scatter(
                nonzero_group["norm_mean_difference"],
                nonzero_group["log_fdr"],
                alpha=0.5,
                c=[color_map[name]],
                marker="o",
                label=name.replace("_", " "),
            )

        if len(zero_group) > 0:
            ax.scatter(
                zero_group.norm_mean_difference,
                zero_group.log_fdr,
                alpha=0.5,
                c=[color_map[name]],
                marker="*",
                s=200,
            )

    plt.axhline(y=sig_threshold_log_fdr, color="r", linestyle="--", label="FDR 0.05")
    plt.xlabel("Normalized difference in means")
    plt.ylabel("-log10(FDR corrected p-value)")
    plt.legend(loc="upper right")
    plt.savefig(output_plot, dpi=300, bbox_inches="tight")
    plt.close()

    # Individual phenotype volcano plots.
    for name, group in groups:
        nonzero_group = group[~group["original_p_value_zero"]]
        zero_group = group[group["original_p_value_zero"]]

        if len(nonzero_group) > 0:
            plt.scatter(
                nonzero_group["mean_difference"],
                nonzero_group["log_fdr"],
                alpha=0.9,
                label=name.replace("_", " "),
                c=[color_map[name]],
                marker="o",
            )

        if len(zero_group) > 0:
            plt.scatter(
                zero_group.mean_difference,
                zero_group.log_fdr,
                alpha=0.9,
                c=[color_map[name]],
                marker="*",
                s=200,
            )

        plt.axhline(y=sig_threshold_log_fdr, color="r", linestyle="--", label="FDR 0.05")
        plt.legend(loc="upper right")
        plt.xlabel("Difference in means")
        plt.ylabel("-log10(FDR corrected p-value)")

        individual_plot_path = output_dir / f"{name}_volcano.png"
        plt.savefig(individual_plot_path, dpi=300, bbox_inches="tight")
        plt.close()

    print(f"Volcano plots saved to {output_dir}.", file=sys.stderr, flush=True)


def main():
    """
    Run the full Panagram association testing pipeline.
    """

    parser = argparse.ArgumentParser(
        description="Run Panagram gene/exon/regulatory-window association testing."
    )

    parser.add_argument(
        "-a",
        "--anchor-name",
        required=True,
        help="Anchor genome name in the Panagram index.",
    )

    parser.add_argument(
        "-p",
        "--panagram-dir",
        required=True,
        help="Path to the Panagram index directory.",
    )

    parser.add_argument(
        "-g",
        "--gff",
        required=True,
        help="Path to the GFF annotation file. Ignored when --type all is used.",
    )

    parser.add_argument(
        "-s",
        "--sample-file",
        required=True,
        help="Path to the tab-delimited sample file. First column should contain sample IDs.",
    )

    parser.add_argument(
        "-m",
        "--phenotype-file",
        required=True,
        help="Path to the phenotype CSV file. Sample IDs should be in the first column.",
    )

    parser.add_argument(
        "-t",
        "--type",
        required=True,
        choices=["gene", "exon", "reg", "all"],
        dest="feature_type",
        help="Feature type to test: gene, exon, reg, or all.",
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        required=True,
        type=Path,
        help="Output directory.",
    )

    parser.add_argument(
        "-n",
        "--threads",
        required=True,
        type=int,
        dest="n_threads",
        help="Number of worker processes to use.",
    )

    parser.add_argument(
        "--chunk-size",
        default=10000,
        type=int,
        help="Number of features/windows per worker chunk.",
    )

    parser.add_argument(
        "--keep-haplotypes",
        action="store_true",
        help="Write wide-format haplotype TSV files. Default is not to save haplotypes.",
    )

    parser.add_argument(
        "--sort-temp-dir",
        default=None,
        type=Path,
        help="Optional temporary directory for Unix sort.",
    )

    parser.add_argument(
        "--max-plot-points-per-phenotype",
        default=50000,
        type=int,
        help="Number of top FDR rows per phenotype to mark for plotting.",
    )

    args = parser.parse_args()

    anchor_name = args.anchor_name
    panagram_dir = args.panagram_dir
    gff = args.gff
    sample_f = args.sample_file
    phenotype_file = args.phenotype_file
    feature_type = args.feature_type
    output_dir = args.output_dir
    n_threads = args.n_threads
    chunk_size = args.chunk_size
    keep_haplotypes = args.keep_haplotypes
    sort_temp_dir = args.sort_temp_dir
    max_plot_points_per_phenotype = args.max_plot_points_per_phenotype

    tree_dir = output_dir / "tree"
    haplotype_dir = output_dir / "haplotypes"
    association_dir = output_dir / "associations"
    plot_dir = output_dir / "plots"

    tree_dir.mkdir(parents=True, exist_ok=True)
    haplotype_dir.mkdir(parents=True, exist_ok=True)
    association_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    tree_file = tree_dir / "bitmap_tree.nwk"
    volcano_plot = plot_dir / f"{feature_type}_volcano.png"

    index = Index(panagram_dir)
    chroms = list(index[anchor_name].chrs.index)

    print("Chromosomes/contigs in Panagram index:", file=sys.stderr, flush=True)

    for chrom in chroms:
        print(f"  {chrom}", file=sys.stderr, flush=True)

    print("Building genome-wide bitmap tree...", file=sys.stderr, flush=True)
    build_bitmap_tree(anchor_name, panagram_dir, tree_file)

    print("Building feature/window list...", file=sys.stderr, flush=True)

    features = build_feature_windows(
        anchor_name,
        panagram_dir,
        gff,
        feature_type,
    )

    print(f"Total features/windows to test: {len(features)}", file=sys.stderr, flush=True)

    if len(features) == 0:
        sys.exit("ERROR: No features/windows were found to test.")

    feature_chunks = split_into_chunks(features, chunk_size)

    print(f"Total chunks: {len(feature_chunks)}", file=sys.stderr, flush=True)
    print(f"Chunk size: {chunk_size}", file=sys.stderr, flush=True)

    print("Running chunk-level analyses in parallel...", file=sys.stderr, flush=True)

    with ProcessPoolExecutor(max_workers=n_threads) as executor:
        futures = []

        for chunk_index, feature_chunk in enumerate(feature_chunks, start=1):
            chunk_id = str(chunk_index).zfill(6)

            future = executor.submit(
                process_feature_chunk,
                chunk_id,
                feature_chunk,
                anchor_name,
                panagram_dir,
                sample_f,
                phenotype_file,
                feature_type,
                tree_file,
                haplotype_dir,
                association_dir,
                keep_haplotypes,
            )

            futures.append(future)

        for future in as_completed(futures):
            chunk_id = future.result()
            print(f"Completed chunk: {chunk_id}", file=sys.stderr, flush=True)

    print("Combining association files...", file=sys.stderr, flush=True)
    combined_association_file = combine_association_files(association_dir, feature_type)

    print("Running disk-based BH/FDR correction...", file=sys.stderr, flush=True)
    corrected_association_file = run_bh_fdr_disk(
        combined_association_file,
        association_dir,
        feature_type,
        sort_temp_dir=sort_temp_dir,
        max_points_per_phenotype=max_plot_points_per_phenotype,
    )

    print("Making volcano plots...", file=sys.stderr, flush=True)
    plot_volcano(corrected_association_file, volcano_plot)

    print("Pipeline complete.", file=sys.stderr, flush=True)
    print(f"Tree: {tree_file}", file=sys.stderr, flush=True)
    print(f"Haplotypes: {haplotype_dir}", file=sys.stderr, flush=True)
    print(f"Associations: {association_dir}", file=sys.stderr, flush=True)
    print(f"Corrected associations: {corrected_association_file}", file=sys.stderr, flush=True)
    print(f"Plots: {plot_dir}", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()

import sys
import argparse
from pathlib import Path
from io import StringIO
from concurrent.futures import ProcessPoolExecutor, as_completed

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from Bio import Phylo
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import linkage
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm
import seaborn as sns
from panagram.index import Index


def get_newick(node, parent_dist, leaf_names, newick=""):
    """
    Convert scipy.cluster.hierarchy.to_tree() output to Newick format.

    Parameters:
        node: output of scipy.cluster.hierarchy.to_tree()
        parent_dist: distance of the parent node
        leaf_names: list of leaf names
        newick: used internally during recursion

    Returns:
        Newick-formatted tree string
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
    Divide the tips of a tree into major clades.

    Parameters:
        tree: Bio.Phylo tree object
        min_size: minimum number of tips for a clade to be considered major

    Returns:
        dict: {clade_name: [tip1, tip2, ...]}
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
    NOTE: GFF names are not standardized, so this may need editing to parse certain GFFs correctly
    """

    feature_name = attributes.split(";")[0]

    # Split on "=" or space to extract the feature name.
    if "=" in feature_name:
        feature_name = feature_name.split("=")[1]
    else:
        feature_name = feature_name.split()[1]

    # Remove quotes.
    feature_name = feature_name.replace('"', "")

    return feature_name


def read_sample_order(sample_f):
    """
    Read sample names from the first column of a tab-delimited sample file.

    The first line is treated as a header and skipped.
    """

    sample_order = []

    with open(sample_f, "r") as f:
        line = f.readline()  # Skip header.
        line = f.readline()

        while line:
            sample_order.append(line.strip().split("\t")[0])
            line = f.readline()

    return sample_order


def build_bitmap_tree(anchor_name, panagram_dir, tree_file):
    """
    Build a genome-wide sample tree from Panagram bitmap data.
    """

    # Hard-coded sampling parameters.
    n_skips = 100
    max_sampled_rows = 2000000

    # Load the Panagram index.
    index = Index(panagram_dir)

    # Get all chromosomes/contigs for the anchor genome.
    chroms = list(index[anchor_name].chrs.index)

    bitmap = []

    # Query bitmap data for each chromosome/contig.
    for chrom in chroms:
        start_coord = 0
        end_coord = index[anchor_name].chrs.loc[chrom, "size"]

        bittmp = index.query_bitmap(anchor_name, chrom, start_coord, end_coord, n_skips)

        print(f"Tree bitmap query complete: {chrom}", file=sys.stderr, flush=True)
        bitmap.append(bittmp)

    # Combine chromosome-level bitmaps into one genome-wide bitmap.
    result = pd.concat(bitmap)

    # Randomly sample rows for clustering, using a fixed seed for reproducibility.
    result_sample = result.sample(n=min(len(result), max_sampled_rows), random_state=42)

    # Cluster samples by their genome-wide bitmap profiles.
    matrix = linkage(result_sample.transpose(), method="average")

    # Convert the clustering result to Newick format.
    tree_tmp = hierarchy.to_tree(matrix, False)
    treedata = get_newick(tree_tmp, tree_tmp.dist, index.genome_names)

    # Write the Newick tree to the output file.
    with open(tree_file, "w") as f:
        f.write(treedata)

    print(f"Tree written to: {tree_file}", file=sys.stderr, flush=True)


def create_feature_tree(bitmap, genome_names):
    """
    Cluster samples by their bitmap profiles for one feature and return top-level clades.
    """

    # Cluster samples using their bitmap profiles for this interval.
    matrix = linkage(bitmap.transpose(), method="ward", metric="euclidean")

    # Convert the clustering result to a Newick tree string.
    tree_tmp = hierarchy.to_tree(matrix, False)
    treedata = get_newick(tree_tmp, tree_tmp.dist, genome_names)

    # Read the Newick tree directly from memory.
    tree = Phylo.read(StringIO(treedata), "newick")

    # Use the top-level clades as haplotype groups.
    clades = get_major_clades(tree, 1)

    return clades


def write_haplotype_output(clades, sample_order, feature_name, out_f):
    """
    Write haplotype group assignments for each sample.

    Samples in Clade_1 are labeled 0.
    Samples in Clade_2 are labeled 1.
    """

    for sample in sample_order:
        if sample in clades["Clade_1"]:
            hap = "0"
        elif sample in clades["Clade_2"]:
            hap = "1"

        out_f.write(feature_name + "\t" + sample + "\t" + hap + "\n")


def query_feature_haplotypes(
    anchor_name, chromosome, panagram_dir, gff, sample_f, feature_type, output_tsv
):
    """
    Query Panagram bitmaps for genes or exons on one chromosome and assign haplotype groups.
    NOTE: this may need to be adjusted for different GFFs.
    """

    # Load the Panagram index inside this process.
    index = Index(panagram_dir)

    # Read the sample file and store sample order.
    sample_order = read_sample_order(sample_f)

    # Track how many times each feature name has been seen.
    # This is used to number exons as feature_1, feature_2, etc.
    feature_counts = {}

    cntr = 0

    # Read in the GFF file and write haplotype assignments to the output TSV.
    with open(output_tsv, "w") as out_f, open(gff, "r") as f:
        line = f.readline()

        while line:
            # Skip comments.
            if line.startswith("#"):
                line = f.readline()
                continue

            tmp = line.strip().split("\t")

            # Skip malformed lines.
            if len(tmp) < 9:
                line = f.readline()
                continue

            # Skip features that do not match the requested type.
            if tmp[2] != feature_type:
                line = f.readline()
                continue

            # Extract chromosome and feature coordinates.
            chrom = tmp[0]
            start_coord = int(tmp[3])
            end_coord = int(tmp[4])

            # Extract the feature name from the GFF attributes column.
            feature_name = extract_feature_name(tmp[8])

            # For exons, add _1, _2, etc. based on how many times we have seen this feature.
            # For genes, keep the original feature name.
            if feature_type == "exon":
                feature_counts[feature_name] = feature_counts.get(feature_name, 0) + 1
                feature_name = f"{feature_name}_{feature_counts[feature_name]}"

            # Only process features on the requested chromosome/contig.
            if chrom == chromosome:
                # Query the Panagram bitmap for this interval.
                bitmap = index.query_bitmap(anchor_name, chrom, start_coord, end_coord, 1)

                # Cluster samples and assign haplotype groups.
                clades = create_feature_tree(bitmap, index.genome_names)
                write_haplotype_output(clades, sample_order, feature_name, out_f)

                # Periodically flush output for long runs.
                if cntr % 10 == 0:
                    out_f.flush()

                cntr += 1

            line = f.readline()


def build_kinship_matrix(tree_file):
    """
    Build a tree-derived kinship-like covariance matrix.
    """

    # Load tree.
    tree = Phylo.read(tree_file, "newick")

    # Get terminal taxa/tips.
    tips = tree.get_terminals()
    sample_names = [tip.name for tip in tips]
    max_samples = len(sample_names) + 3

    # Total root-to-tip distance, used for normalization.
    max_distance = max(tree.distance(tree.root, tip) for tip in tips)
    max_distance = max_distance * 2

    # Compute pairwise tree-based similarity values.
    distance_matrix = pd.DataFrame(index=sample_names, columns=sample_names, dtype=float)

    for tip_a in tips:
        for tip_b in tips:
            distance_matrix.loc[tip_a.name, tip_b.name] = (
                max_distance - tree.distance(tip_a, tip_b)
            ) / max_distance

    # Convert the distance matrix into a kinship-like covariance matrix.
    kinship_matrix = distance_matrix.values.copy()
    kinship_matrix += np.eye(kinship_matrix.shape[0]) * 1e-6
    kinship_matrix = kinship_matrix / np.mean(np.diag(kinship_matrix))

    return kinship_matrix, distance_matrix, max_samples


def test_one_feature(
    feature_name,
    kinship_matrix,
    feature_haplotype_df,
    group1_size,
    group2_size,
    group1_samples,
    group2_samples,
    phenotype_df,
    phenotype_list,
    distance_matrix,
    max_samples,
    out_f,
):
    """
    Test one gene/exon against all phenotypes using GLS.
    """

    # Convert phenotype sample IDs from index to a column for merging.
    phenotype_for_merge = phenotype_df.copy()
    phenotype_for_merge.index.name = "sample_id"
    phenotype_for_merge = phenotype_for_merge.reset_index()

    # Merge phenotype values with haplotype calls for this feature.
    merged_df = phenotype_for_merge.merge(feature_haplotype_df, on="sample_id", validate="1:1")

    # Put sample IDs back as the index and order samples to match the kinship matrix.
    merged_df = merged_df.set_index("sample_id")
    sample_order = distance_matrix.index.tolist()
    merged_df = merged_df.loc[sample_order]

    def test_each_phenotype(
        phenotype_name,
        kinship_matrix,
        feature_name,
        group1_size,
        group2_size,
        group1_samples,
        group2_samples,
    ):
        """
        Test one phenotype for one gene/exon.
        """

        phenotype_values = merged_df[phenotype_name].values
        haplotype_values = merged_df["haplotype"].values.reshape(-1, 1)
        haplotype_values = sm.add_constant(haplotype_values)

        valid_indices = ~np.isnan(phenotype_values)

        phenotype_values = phenotype_values[valid_indices]
        haplotype_values = haplotype_values[valid_indices]

        group1_phenotype_values = list(
            merged_df[merged_df.index.isin(group1_samples)][phenotype_name].dropna()
        )
        group2_phenotype_values = list(
            merged_df[merged_df.index.isin(group2_samples)][phenotype_name].dropna()
        )

        if (
            len(valid_indices) < max_samples
            and len(haplotype_values) < max_samples
            and len(phenotype_values) < max_samples
            and len(haplotype_values) > 0
            and len(phenotype_values) > 0
            and len(group1_phenotype_values) > 0
            and len(group2_phenotype_values) > 0
        ):
            kinship_matrix_subset = kinship_matrix[np.ix_(valid_indices, valid_indices)]

            model = sm.GLS(phenotype_values, haplotype_values, sigma=kinship_matrix_subset)

            results = model.fit()

            group1_mean = sum(group1_phenotype_values) / len(group1_phenotype_values)
            group2_mean = sum(group2_phenotype_values) / len(group2_phenotype_values)

            out_f.write(
                feature_name
                + "\t"
                + phenotype_name
                + "\t"
                + group1_size
                + "\t"
                + group2_size
                + "\t"
                + str(group1_mean)
                + "\t"
                + str(group2_mean)
                + "\t"
                + str(results.pvalues[1])
                + "\n"
            )

    for phenotype_column in phenotype_list:
        test_each_phenotype(
            phenotype_column,
            kinship_matrix,
            feature_name,
            group1_size,
            group2_size,
            list(group1_samples.index),
            list(group2_samples.index),
        )


def run_association_tests(
    tree_file, chromosome, haplotype_file, gff, phenotype_file, feature_type, output_tsv
):
    """
    Test associations between gene/exon haplotype groups and phenotypes for one chromosome.
    """

    # Build kinship matrix from tree.
    kinship_matrix, distance_matrix, max_samples = build_kinship_matrix(tree_file)

    # Read phenotype and haplotype files.
    # Phenotype file is expected to be CSV with sample IDs in the first column.
    phenotype_df = pd.read_csv(phenotype_file, index_col=0)

    column_names = ["gene_id", "sample_id", "haplotype"]

    haplotype_df = pd.read_csv(haplotype_file, sep="\t", index_col=1, names=column_names)

    phenotype_list = phenotype_df.columns.values.tolist()

    counter = 0
    feature_counts = {}

    # Read the GFF and write association results to the output TSV.
    with open(output_tsv, "w") as out_f, open(gff, "r") as file:
        line = file.readline()

        while line:
            counter += 1

            # Skip comments.
            if line.startswith("#"):
                line = file.readline()
                continue

            fields = line.strip().split("\t")

            # Skip malformed lines.
            if len(fields) < 9:
                line = file.readline()
                continue

            line_chromosome = fields[0]

            if line_chromosome != chromosome:
                line = file.readline()
                continue

            # Skip features that do not match the requested type.
            if fields[2] != feature_type:
                line = file.readline()
                continue

            # Extract feature name from GFF attributes.
            feature_name = extract_feature_name(fields[8])

            # For exons, add _1, _2, etc. to match the haplotype output.
            # For genes, keep the original feature name.
            if feature_type == "exon":
                feature_counts[feature_name] = feature_counts.get(feature_name, 0) + 1
                feature_name = f"{feature_name}_{feature_counts[feature_name]}"

            feature_haplotype_df = haplotype_df[haplotype_df["gene_id"] == feature_name]

            group1_count = len(feature_haplotype_df[feature_haplotype_df["haplotype"] == 0])
            group2_count = len(feature_haplotype_df[feature_haplotype_df["haplotype"] == 1])

            if min(group1_count, group2_count) > 0:
                test_one_feature(
                    feature_name,
                    kinship_matrix,
                    feature_haplotype_df,
                    str(group1_count),
                    str(group2_count),
                    feature_haplotype_df[feature_haplotype_df["haplotype"] == 0],
                    feature_haplotype_df[feature_haplotype_df["haplotype"] == 1],
                    phenotype_df,
                    phenotype_list,
                    distance_matrix,
                    max_samples,
                    out_f,
                )

            if counter % 10 == 0:
                out_f.flush()

            line = file.readline()


def run_one_chromosome(
    chromosome,
    anchor_name,
    panagram_dir,
    gff,
    sample_f,
    phenotype_file,
    feature_type,
    tree_file,
    haplotype_dir,
    association_dir,
):
    """
    Run haplotype querying and association testing for one chromosome.
    """

    haplotype_tsv = haplotype_dir / f"{chromosome}_{feature_type}_haplotypes.tsv"
    association_tsv = association_dir / f"{chromosome}_{feature_type}_associations.tsv"

    print(f"Starting chromosome: {chromosome}", file=sys.stderr, flush=True)

    # Query haplotypes for this chromosome.
    query_feature_haplotypes(
        anchor_name, chromosome, panagram_dir, gff, sample_f, feature_type, haplotype_tsv
    )

    print(f"Haplotypes complete: {chromosome}", file=sys.stderr, flush=True)

    # Run association tests for this chromosome.
    run_association_tests(
        tree_file, chromosome, haplotype_tsv, gff, phenotype_file, feature_type, association_tsv
    )

    print(f"Associations complete: {chromosome}", file=sys.stderr, flush=True)

    return chromosome


def plot_volcano(input_folder, output_plot):
    """
    Read association result TSVs, apply FDR correction, and make volcano plots.
    """

    output_plot = Path(output_plot)
    output_dir = output_plot.parent

    # Create output directory if it does not already exist.
    output_dir.mkdir(parents=True, exist_ok=True)

    # Find all TSV files in the input folder.
    tsv_files = list(Path(input_folder).glob("*.tsv"))

    # Stop early if no TSV files are found.
    if len(tsv_files) == 0:
        sys.exit(f"ERROR: No .tsv files found in {input_folder}")

    # Expected columns in the association result TSVs.
    column_names = [
        "gene_name",
        "phenotype_name",
        "group1_size",
        "group2_size",
        "group1_mean",
        "group2_mean",
        "p_value",
    ]

    # Read and concatenate all TSV files.
    dataframes = []

    for tsv_file in tsv_files:
        df_temp = pd.read_csv(tsv_file, sep="\t", names=column_names, dtype={"p_value": float})
        dataframes.append(df_temp)

    df = pd.concat(dataframes, ignore_index=True)

    # Drop rows where either haplotype group is less than 10% of the total group size.
    total_size = round(df["group1_size"] + df["group2_size"])
    df = df[(df["group1_size"] >= 0.1 * total_size) & (df["group2_size"] >= 0.1 * total_size)]

    # Stop if filtering removed all rows.
    if len(df) == 0:
        sys.exit("ERROR: No rows remain after filtering by haplotype group size.")

    # Mark original p-values of 0 before replacing them.
    # These will be plotted with star markers.
    df["is_zero_p"] = df["p_value"] == 0

    # Replace p-values of 0 so -log10(p) can be calculated.
    min_nonzero_p = df[df["p_value"] > 0]["p_value"].min()

    if pd.isna(min_nonzero_p):
        # If all p-values are zero, use a fixed very small value.
        replacement_p = 1e-300
    else:
        # If p = 0, set it to 1e<exponent - 2>.
        exponent = int(np.floor(np.log10(min_nonzero_p))) - 2
        replacement_p = 10**exponent

    df.loc[df["p_value"] == 0, "p_value"] = replacement_p

    print(f"Minimum non-zero p-value: {min_nonzero_p}", file=sys.stderr, flush=True)
    print(f"Setting zero p-values to: {replacement_p}", file=sys.stderr, flush=True)

    # Calculate derived columns.
    df["mean_difference"] = df["group1_mean"] - df["group2_mean"]
    df["log_pval"] = -np.log10(df["p_value"])
    df["fdr_corrected"] = multipletests(df["p_value"], method="fdr_bh")[1]
    df["log_fdr"] = -np.log10(df["fdr_corrected"])

    # Normalize mean differences per phenotype.
    # This keeps the original normalization approach.
    groups = df.groupby("phenotype_name")
    norm_diffs = []

    for name, group in groups:
        g_max = group["mean_difference"].max()
        g_min = group["mean_difference"].min()

        # Avoid division by zero if all mean differences are identical.
        if g_max == g_min:
            group_norm = group["mean_difference"] * 0
        else:
            group_norm = group["mean_difference"] / (g_max - g_min)

        norm_diffs.append(group_norm)

    df["norm_mean_difference"] = pd.concat(norm_diffs).sort_index()

    # Save sorted results.
    df_sorted = df.sort_values("log_fdr", ascending=False)
    sorted_output = output_dir / (output_plot.stem + "_sorted.tsv")
    df_sorted.to_csv(sorted_output, sep="\t", index=False)

    print(f"Dataframe saved to: {sorted_output}", file=sys.stderr, flush=True)

    # Assign visually distinct colors to all phenotypes.
    unique_phenotypes = df["phenotype_name"].unique()
    palette = sns.color_palette("husl", n_colors=len(unique_phenotypes))
    color_map = dict(zip(unique_phenotypes, palette))

    # FDR threshold.
    sig_threshold_fdr = 0.05
    sig_threshold_log_fdr = -np.log10(sig_threshold_fdr)

    # Recreate groups after adding derived columns.
    groups = df.groupby("phenotype_name")

    # Combined plot with normalized differences.
    plt.rcParams["figure.figsize"] = (4.5, 5)
    fig, ax = plt.subplots()

    for name, group in groups:
        non_zero = group[~group["is_zero_p"]]
        zero_p = group[group["is_zero_p"]]

        ax.scatter(
            non_zero.norm_mean_difference,
            non_zero.log_fdr,
            alpha=0.5,
            c=[color_map[name]],
            label=name.replace("_", " "),
        )

        ax.scatter(
            zero_p.norm_mean_difference,
            zero_p.log_fdr,
            alpha=0.5,
            c=[color_map[name]],
            marker="*",
            s=200,
        )

    plt.axhline(y=sig_threshold_log_fdr, color="r", linestyle="--", label="FDR")

    plt.xlabel("Difference in means")
    plt.ylabel("-log10(FDR corrected pvalue)")
    plt.legend(loc="upper right")
    plt.savefig(output_plot, dpi=300, bbox_inches="tight")
    plt.close()

    # Individual plots per phenotype with raw differences.
    for name, group in groups:
        non_zero = group[~group["is_zero_p"]]
        zero_p = group[group["is_zero_p"]]

        plt.scatter(
            non_zero.mean_difference,
            non_zero.log_fdr,
            alpha=0.9,
            label=name.replace("_", " "),
            c=[color_map[name]],
        )

        plt.scatter(
            zero_p.mean_difference,
            zero_p.log_fdr,
            alpha=0.9,
            c=[color_map[name]],
            marker="*",
            s=200,
        )

        plt.axhline(y=sig_threshold_log_fdr, color="r", linestyle="--", label="FDR")

        plt.legend(loc="upper right")
        plt.xlabel("Difference in means")
        plt.ylabel("-log10(FDR corrected pvalue)")

        individual_plot_path = output_dir / f"{name}_volcano.png"
        plt.savefig(individual_plot_path, dpi=300, bbox_inches="tight")
        plt.close()

    print(f"Volcano plots saved to {output_dir}.", file=sys.stderr, flush=True)


def main():
    """
    Run the full Panagram association testing pipeline.
    """

    parser = argparse.ArgumentParser(
        description="Run the full Panagram gene/exon association testing pipeline."
    )

    parser.add_argument(
        "-a", "--anchor-name", required=True, help="Anchor genome name in the Panagram index."
    )

    parser.add_argument(
        "-p", "--panagram-dir", required=True, help="Path to the Panagram index directory."
    )

    parser.add_argument("-g", "--gff", required=True, help="Path to the GFF annotation file.")

    parser.add_argument(
        "-s",
        "--sample-file",
        required=True,
        help="Path to the tab-delimited sample file. The first column should contain sample IDs.",
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
        choices=["gene", "exon"],
        dest="feature_type",
        help="Feature type to test. Must be either 'gene' or 'exon'.",
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        required=True,
        type=Path,
        help="Output directory for tree, haplotypes, association results, and plots.",
    )

    parser.add_argument(
        "-n",
        "--threads",
        required=True,
        type=int,
        dest="n_threads",
        help="Number of chromosomes/contigs to process in parallel.",
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

    # Make output directories.
    tree_dir = output_dir / "tree"
    haplotype_dir = output_dir / "haplotypes"
    association_dir = output_dir / "associations"
    plot_dir = output_dir / "plots"

    tree_dir.mkdir(parents=True, exist_ok=True)
    haplotype_dir.mkdir(parents=True, exist_ok=True)
    association_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    # Define output files.
    tree_file = tree_dir / "bitmap_tree.nwk"
    volcano_plot = plot_dir / f"{feature_type}_volcano.png"

    # Load the Panagram index once to get chromosome names.
    index = Index(panagram_dir)
    chroms = list(index[anchor_name].chrs.index)

    print("Chromosomes/contigs to process:", file=sys.stderr, flush=True)
    for chrom in chroms:
        print(f"  {chrom}", file=sys.stderr, flush=True)

    # Build the genome-wide bitmap tree once.
    print("Building genome-wide bitmap tree...", file=sys.stderr, flush=True)
    build_bitmap_tree(anchor_name, panagram_dir, tree_file)

    # Run haplotype querying and association testing by chromosome in parallel.
    print("Running chromosome-level analyses in parallel...", file=sys.stderr, flush=True)

    with ProcessPoolExecutor(max_workers=n_threads) as executor:
        futures = []

        for chrom in chroms:
            future = executor.submit(
                run_one_chromosome,
                chrom,
                anchor_name,
                panagram_dir,
                gff,
                sample_f,
                phenotype_file,
                feature_type,
                tree_file,
                haplotype_dir,
                association_dir,
            )

            futures.append(future)

        for future in as_completed(futures):
            chrom = future.result()
            print(f"Finished chromosome: {chrom}", file=sys.stderr, flush=True)

    # Make FDR-corrected results and volcano plots.
    print("Making volcano plots...", file=sys.stderr, flush=True)
    plot_volcano(association_dir, volcano_plot)

    print("Pipeline complete.", file=sys.stderr, flush=True)
    print(f"Tree: {tree_file}", file=sys.stderr, flush=True)
    print(f"Haplotypes: {haplotype_dir}", file=sys.stderr, flush=True)
    print(f"Associations: {association_dir}", file=sys.stderr, flush=True)
    print(f"Plots: {plot_dir}", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()

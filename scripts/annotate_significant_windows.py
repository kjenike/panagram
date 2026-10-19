#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


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


def read_genes_from_gff(gff_file):
    """
    Read gene coordinates from a GFF file.

    Returns a dictionary:
        genes_by_chrom[chrom] = DataFrame with columns:
            gene_name, start, end

    Coordinates are kept in the same coordinate style as your pipeline.
    """

    genes = []

    with open(gff_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 9:
                continue

            chrom = fields[0]
            feature_type = fields[2]

            if feature_type != "gene":
                continue

            start = int(fields[3])
            end = int(fields[4])
            gene_name = extract_feature_name(fields[8])

            genes.append(
                {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "gene_name": gene_name,
                }
            )

    gene_df = pd.DataFrame(genes)

    if gene_df.empty:
        raise ValueError("No gene features were found in the GFF.")

    genes_by_chrom = {}

    for chrom, chrom_gene_df in gene_df.groupby("chrom"):
        chrom_gene_df = chrom_gene_df.sort_values(["start", "end"]).reset_index(drop=True)
        genes_by_chrom[chrom] = chrom_gene_df

    return genes_by_chrom


def parse_all_window_feature_id(feature_id, window_size):
    """
    Parse an all-mode window feature ID.

    Expected format:
        chrname_startcoordinate

    This uses rsplit so chromosome names may contain underscores.

    Example:
        chr1_RagTag_123400
    """

    chrom, start_string = feature_id.rsplit("_", 1)

    start = int(start_string)
    end = start + window_size

    return chrom, start, end


def annotate_region(chrom, start, end, genes_by_chrom):
    """
    Annotate one region with intersecting genes.

    If the region intersects a gene:
        return intersecting genes and leave upstream/downstream blank.

    If the region does not intersect a gene:
        return nearest upstream and downstream genes.
    """

    blank_annotation = {
        "intersecting_genes": "",
        "upstream_gene": "",
        "upstream_distance": "",
        "downstream_gene": "",
        "downstream_distance": "",
    }

    if chrom not in genes_by_chrom:
        return blank_annotation

    gene_df = genes_by_chrom[chrom]

    # Intersecting genes:
    # gene_start < region_end and gene_end > region_start
    intersecting = gene_df[(gene_df["start"] < end) & (gene_df["end"] > start)]

    if len(intersecting) > 0:
        annotation = blank_annotation.copy()
        annotation["intersecting_genes"] = ",".join(intersecting["gene_name"].astype(str))
        return annotation

    # If the region is not inside/intersecting a gene, find nearest neighbors.
    upstream = gene_df[gene_df["end"] <= start]
    downstream = gene_df[gene_df["start"] >= end]

    annotation = blank_annotation.copy()

    if len(upstream) > 0:
        nearest_upstream = upstream.iloc[upstream["end"].argmax()]
        annotation["upstream_gene"] = str(nearest_upstream["gene_name"])
        annotation["upstream_distance"] = str(start - int(nearest_upstream["end"]))

    if len(downstream) > 0:
        nearest_downstream = downstream.iloc[downstream["start"].argmin()]
        annotation["downstream_gene"] = str(nearest_downstream["gene_name"])
        annotation["downstream_distance"] = str(int(nearest_downstream["start"]) - end)

    return annotation


def annotate_significant_associations(
    association_file,
    gff_file,
    output_file,
    fdr_cutoff=0.05,
    window_size=200,
    chunksize=500000,
):
    """
    Annotate significant all-mode window associations with nearby genes.

    Significant rows are:
        fdr <= fdr_cutoff
    """

    genes_by_chrom = read_genes_from_gff(gff_file)

    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    output_columns = [
        "feature_id",
        "phenotype_name",
        "chrom",
        "start",
        "end",
        "p_value",
        "fdr",
        "intersecting_genes",
        "upstream_gene",
        "upstream_distance",
        "downstream_gene",
        "downstream_distance",
    ]

    # Cache annotations by feature_id so repeated regions across phenotypes do not need to be re-annotated.
    annotation_cache = {}

    with open(output_file, "w") as out_f:
        out_f.write("\t".join(output_columns) + "\n")

        for chunk in pd.read_csv(association_file, sep="\t", chunksize=chunksize):
            # Keep only significant associations.
            chunk["fdr"] = pd.to_numeric(chunk["fdr"])
            sig_chunk = chunk[chunk["fdr"] <= fdr_cutoff]

            if len(sig_chunk) == 0:
                continue

            for _, row in sig_chunk.iterrows():
                feature_id = row["feature_id"]

                if feature_id in annotation_cache:
                    chrom, start, end, annotation = annotation_cache[feature_id]
                else:
                    chrom, start, end = parse_all_window_feature_id(
                        feature_id,
                        window_size,
                    )

                    annotation = annotate_region(
                        chrom,
                        start,
                        end,
                        genes_by_chrom,
                    )

                    annotation_cache[feature_id] = (
                        chrom,
                        start,
                        end,
                        annotation,
                    )

                out_fields = [
                    str(feature_id),
                    str(row["phenotype_name"]),
                    str(chrom),
                    str(start),
                    str(end),
                    str(row["p_value"]),
                    str(row["fdr"]),
                    annotation["intersecting_genes"],
                    annotation["upstream_gene"],
                    annotation["upstream_distance"],
                    annotation["downstream_gene"],
                    annotation["downstream_distance"],
                ]

                out_f.write("\t".join(out_fields) + "\n")


def main():
    """
    Annotate significant all-mode association windows with nearby genes.
    """

    parser = argparse.ArgumentParser(
        description="Annotate significant all-mode association windows with nearby genes."
    )

    parser.add_argument(
        "-a",
        "--association-file",
        required=True,
        help="FDR-sorted association TSV from the main pipeline.",
    )

    parser.add_argument(
        "-g",
        "--gff",
        required=True,
        help="GFF annotation file.",
    )

    parser.add_argument(
        "-o",
        "--output-file",
        required=True,
        help="Output annotated significant association TSV.",
    )

    parser.add_argument(
        "--fdr-cutoff",
        default=0.05,
        type=float,
        help="FDR cutoff for significant associations. Default: 0.05.",
    )

    parser.add_argument(
        "--window-size",
        default=200,
        type=int,
        help="Window size used in all-mode association testing. Default: 200.",
    )

    parser.add_argument(
        "--chunksize",
        default=500000,
        type=int,
        help="Number of association rows to read at a time. Default: 500000.",
    )

    args = parser.parse_args()

    annotate_significant_associations(
        association_file=args.association_file,
        gff_file=args.gff,
        output_file=args.output_file,
        fdr_cutoff=args.fdr_cutoff,
        window_size=args.window_size,
        chunksize=args.chunksize,
    )


if __name__ == "__main__":
    main()

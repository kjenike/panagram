import subprocess
import argparse
from pathlib import Path
import pandas as pd


def parse_gff(gff_file, omit_unknown=True):
    """Parse GFF file to extract only genes.

    Args:
        gff_file (Path): path to the GFF file
        omit_unknown (bool): whether to omit genes with unknown function

    Returns:
        pd.DataFrame: gene annotations from GFF file
    """

    gff_df = pd.read_csv(gff_file, sep="\t", comment="#", header=None)
    gff_df.columns = [
        "seqid",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes",
    ]
    gff_df = gff_df[gff_df["type"] == "gene"]

    if omit_unknown:
        gff_df = gff_df[
            gff_df["attributes"].str.contains("unknown|hypothetical", case=False, na=False) == False
        ]
    return gff_df


def annotate_bed_with_genes(bed_file, gff_file, output_file):
    """Call bedtools intersect to annotate bed with genes from GFF file.

    Args:
        bed_file (Path): path to the bed file
        gff_file (Path): path to the GFF file
        output_file (Path): path to the output file
    """

    command = f"bedtools intersect -a {bed_file} -b {gff_file} -wa -wb > {output_file}"
    subprocess.run(command, shell=True, check=True)
    return


def extract_attribute_from_bed(bed_file, attribute):
    """Extract specified attribute from annotated bed file.

    Args:
        bed_file (Path): path to the bed file
        attribute (str): attribute to extract

    Returns:
        str: extracted attribute values
    """

    bed_df = pd.read_csv(bed_file, sep="\t", header=None)
    bed_df.columns = [
        "chrom_a",
        "start_a",
        "end_a",
        "notes",
        "seqid",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes",
    ]

    # extract specified attribute (everything after '=' up to ';')
    bed_df["attributes"] = bed_df["attributes"].str.extract(f"{attribute}=([^;]+)")

    attributes = bed_df["attributes"].dropna().tolist()
    attributes = "\n".join(attributes).replace(",", "\n").replace(" ", "")

    return attributes


def main():
    """Annotate introgression bed files with gene information from GFF files. Creates annotated bed
    files and a text file with GO terms."""

    parser = argparse.ArgumentParser(description="Annotate introgression bed files.")
    parser.add_argument("--gff", type=Path, required=True, help="Path to the GFF file.")
    parser.add_argument("--bed", type=Path, required=True, help="Path to the bed folder.")
    parser.add_argument("--accession", type=str, required=True, help="Accession name.")
    parser.add_argument(
        "--term", type=str, default="sollyc4.0_ID", help="Attribute to extract from GFF."
    )
    args = parser.parse_args()

    gff_file = args.gff
    bed_folder = args.bed
    accession = args.accession
    term = args.term

    # create output folders
    gff_output_file = gff_file.parent / f"{gff_file.stem}_genes_only.gff3"
    bed_output_folder = bed_folder.parent / "annotated_beds"
    go_output_folder = bed_folder.parent / "go_terms"
    bed_output_folder.mkdir(exist_ok=True)
    go_output_folder.mkdir(exist_ok=True)

    gff_df = parse_gff(gff_file)
    gff_df.to_csv(gff_output_file, sep="\t", index=False, header=False)

    bed_files = list(bed_folder.glob(f"{accession}*.bed"))
    bed_files.sort()

    # process bed files separately
    all_go_terms = ""
    for bed_file in bed_files:
        print(f"Processing {bed_file.name}...")
        bed_output_file = bed_output_folder / bed_file.name
        annotate_bed_with_genes(bed_file, gff_output_file, bed_output_file)

        # extract terms from annotated bed files
        go_attributes = extract_attribute_from_bed(bed_output_file, attribute=term)
        go_output_file = go_output_folder / f"{bed_file.stem}_{term}.txt"
        with open(go_output_file, "w") as go_out:
            go_out.write(go_attributes)
        all_go_terms += go_attributes + "\n"

    # write all GO terms to a single file
    all_go_output_file = go_output_folder / f"{accession}_all_{term}.txt"
    with open(all_go_output_file, "w") as all_go_out:
        all_go_out.write(all_go_terms)

    return


if __name__ == "__main__":
    main()

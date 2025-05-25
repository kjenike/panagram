import argparse
from pathlib import Path
import subprocess
import math
import numpy as np
import pandas as pd
from panagram.index import Index
from call_introgressions import bins_to_bed


def read_fasta_generator(fp):
    """Generator for fasta file - allows reading of one record at a time to save memory.
    Or, all records can be read at once using read_fasta().

    :param fp: file pointer
    :type fp: _io.TextIOWrapper
    :yield: name and seq in the file
    :rtype: tuple(name, seq)
    """

    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                yield (name, "".join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name:
        yield (name, "".join(seq))


def read_fasta(input_file):
    """Read FASTA into a pandas dataframe.

    :param input_file: FASTA file path
    :type input_file: str
    :return: dataframe of FASTA records with "Name" and "Sequence" columns
    :rtype: pd.DataFrame
    """

    input_file = Path(input_file)
    # print(f"Loading {input_file.name}...")
    names = []
    seqs = []
    with open(input_file) as f:
        for name, seq in read_fasta_generator(f):
            names.append(name[1:])  # remove ">" from name
            seqs.append(seq)
    fasta_df = pd.DataFrame({"Name": names, "Sequence": seqs})
    return fasta_df


def extract_bed_region(bed_row, fasta_df):
    """Grab the sequence corresponding to the bed file coordinates listed.

    :param bed_row: row of a bed file's dataframe
    :type bed_row: pd.DataFrame
    :param fasta_df: fasta information containing sequences that the bed file references
    :type fasta_df: pd.DataFrame
    :return: sequence corresponding to bed coordinates
    :rtype: str
    """
    # extract a region from a FASTA that a BED file points to
    # follows bedtools getfasta logic (0-index, inclusive start coord, exclusive end coord)
    chromosome_sequence = fasta_df.loc[fasta_df.Name == bed_row["Chromosome"], "Sequence"].iloc[0]
    start = int(bed_row["Start"])
    end = int(bed_row["End"])
    introgression_sequence = chromosome_sequence[start:end]
    return introgression_sequence


def bed_file_is_empty(bed_file):
    try:
        bed_df = pd.read_csv(
            bed_file,
            sep="\t",
            header=None,
        )
    except pd.errors.EmptyDataError as e:
        return True
    return False


def read_bed_file(
    bed_file,
    accession=None,
    samples_file=None,
    index_dir=None,
    fasta_file=None,
):
    """Read in a given bed file. If given a corresponding FASTA file, copy the sequence data into
    the dataframe. Finds the FASTA file location given the acession, samples_file and index_dir if
    its path is not provided.

    :param bed_file: path to bed file with introgression locations
    :type bed_file: str or Path
    :param accession: id used in panagram for the genome, defaults to None
    :type accession: str, optional
    :param samples_file: samples.tsv file provided to panagram, defaults to None
    :type samples_file: str or Path, optional
    :param index_dir: panagram index folder; required if samples_file contains relative paths, defaults to None
    :type index_dir: str or Path, optional
    :param fasta_file: path to FASTA file; if provided, samples_file is not required, defaults to None
    :type fasta_file: str or Path, optional
    :return: dataframes representing bed file and FASTA file
    :rtype: tuple(pd.DataFrame, pd.DataFrame)
    """

    if bed_file_is_empty(bed_file):
        return None, None

    bed_df = pd.read_csv(
        bed_file,
        sep="\t",
        header=None,
    )

    bed_df = bed_df.iloc[:, 0:4]
    col_names = ["Chromosome", "Start", "End", "Notes"]
    bed_df.columns = col_names
    bed_df["Sequence"] = None

    if fasta_file:
        accession_fasta = fasta_file
    elif samples_file:
        # get the fasta file associated with the given accession
        samples_df = pd.read_csv(samples_file, sep="\t", header=0, index_col=0)
        accession_fasta = samples_df.loc[accession, "fasta"]
        # convert relative path to absolute path if a relative path was given in the samples.tsv
        if index_dir:
            accession_fasta = Path(index_dir) / accession_fasta
    # return bed file without Sequence data if no FASTA is provided
    else:
        return bed_df, None

    # read fasta and bed file
    fasta_df = read_fasta(accession_fasta)

    # extract the sequence from the FASTA and put it into the bed df
    bed_df["Sequence"] = bed_df.apply(extract_bed_region, axis=1, fasta_df=fasta_df)
    return bed_df, fasta_df


def merge_centromere_regions(bed_df, fasta_df, bin_size):
    # detect potential centromere regions (regions that are 2 bins apart)
    bed_df["End_Shifted"] = bed_df.End.shift(1)
    bed_df["New_Start"] = bed_df.Start.shift(1)
    bed_df["New_End"] = bed_df.End
    bed_df["Centromere"] = ((bed_df.Start - bed_df.End_Shifted) / bin_size) == 2

    # where two intros are 2 bins apart, get the in-between region
    centromere_df = bed_df[bed_df.Centromere == True].copy()
    # if there are no such regions, return with no changes
    if centromere_df.empty:
        return bed_df

    # replace End with Start and Start with End_Shifted to get centromere coordinates
    centromere_df["End"] = centromere_df["Start"]
    centromere_df["Start"] = centromere_df["End_Shifted"]

    # get centromere's sequence
    centromere_df["Sequence"] = centromere_df.apply(extract_bed_region, axis=1, fasta_df=fasta_df)

    # drop all regions not containing a string of at least 50 N's in a row
    centromere_df = centromere_df[centromere_df.Sequence.str.contains("N" * 50)]

    # merge adjacent regions that are confirmed to be centromeres, using New_Start and New_End
    centromere_df["End"] = centromere_df["New_End"]
    centromere_df["Start"] = centromere_df["New_Start"]

    # add centromere_df rows to the bottom of bed_df
    bed_df = pd.concat([bed_df, centromere_df], ignore_index=True)

    # drop original unmerged rows up at the top
    bed_df = bed_df.drop_duplicates(subset="Start", keep="last").drop_duplicates(
        subset="End", keep="last"
    )
    bed_df = bed_df[["Chromosome", "Start", "End", "Notes", "Sequence"]].reset_index(drop=True)

    # fix sequences now that centromere has been added
    bed_df["Sequence"] = bed_df.apply(extract_bed_region, axis=1, fasta_df=fasta_df)

    # fix Start/End since they may have been converted to floats
    bed_df["Start"] = bed_df["Start"].astype(int)
    bed_df["End"] = bed_df["End"].astype(int)
    return bed_df


def align_to_reference(
    reference_file,
    query_file,
    minimap_flags,
    output_file,
):
    # Call minimap2 command - requires minimap2 to be in path
    minimap_command = ["minimap2"] + minimap_flags + [reference_file, query_file]

    # Run the command
    try:
        result = subprocess.run(minimap_command, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        print("Error:", e.stderr)
        exit(1)

    # Save the result as a bed file
    log_file = output_file.parent / f"{output_file.stem}.log"

    with open(log_file, "w") as f:
        f.write(result.stderr)
    with open(output_file, "w") as f:
        f.write(result.stdout)
    return


def liftover_to_reference(
    bed_file,
    paf_file,
    output_file,
):
    # Call liftover script - requires paftools.js to be in path
    liftover_command = [
        "paftools.js",
        "liftover",
        paf_file,
        bed_file,
        "output.liftover",
    ]

    # Run the command
    try:
        result = subprocess.run(liftover_command, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        print("Error:", e.stderr)
        exit(1)

    # Save the result as a bed file
    with open(output_file, "w") as f:
        f.write(result.stdout)
    return


def get_intro_df_template(bin_size, chr_length):
    # create pandas df full of 0s
    n_bins = math.ceil(chr_length / bin_size)
    bin_names = [i * bin_size for i in range(n_bins)]
    zeros = np.zeros((1, n_bins), dtype=int)
    intro_df = pd.DataFrame(zeros, columns=bin_names)
    return intro_df


def bed_to_bins(bed_file, bin_size, chr_length):
    # convert bedfile coordinates back into an index of bins for a chromosome
    # match format of txt files from original introgressions script
    # with a per-bin 0/1 introgression score

    # read in bed file
    bed_df, _ = read_bed_file(bed_file)
    intro_df = get_intro_df_template(bin_size, chr_length)

    # convert coordinates to bin labels
    bed_df["Start Bin"] = (bed_df["Start"] // bin_size) * bin_size  # round down to nearest bin
    bed_df["End Bin"] = ((bed_df["End"] // bin_size) * bin_size) + bin_size  # round up
    # get all column labels between start and end bin
    bed_df["Bin Labels"] = bed_df.apply(
        lambda row: list(range(row["Start Bin"], row["End Bin"], bin_size)), axis=1
    )

    # for each bedfile entry, add 1 to corresponding bins in df
    bin_labels = [pos for sublist in bed_df["Bin Labels"] for pos in sublist]
    for bin_label in bin_labels:
        # the last few chr coords were clipped in get_distances.py if they didn't form a full bin
        if bin_label > intro_df.columns[-1]:
            continue
        intro_df.loc[:, bin_label] += 1

    # match the format of the bins_df from call_introgressions
    intro_df = intro_df.T
    intro_df.columns = ["introgression"]
    return intro_df


def fill_gaps(row, rounds=1):
    # Fill small gaps of 0s surrounded by 1s
    row = row.values
    for j in range(rounds):
        for i in range(1, len(row) - 1):
            if row[i] == 0 and (row[i - 1] >= 1 and row[i + 1] >= 1):
                row[i] = 1
    return row


def remove_single_bins(row):
    # omit introgressions that are only a single bin wide
    # recommend running fill gaps prior to this func after liftover, which can fragment intros
    for i in range(1, len(row) - 1):
        if row[i] >= 1 and (row[i - 1] == 0 and row[i + 1] == 0):
            row[i] = 0
    return row


def postprocess_introgressions():
    # does postprocessing on specified file or all files in a folder
    # allows for liftover, gap filling, centromere filling, and lone bin removal in any order
    parser = argparse.ArgumentParser(description="Introgression postprocessing.")
    parser.add_argument(
        "--bed", type=str, help="path to introgression bed file or folder", required=True
    )
    parser.add_argument("--idx", type=str, help="path to Panagram index folder", required=True)
    parser.add_argument(
        "-a",
        nargs="+",
        help="action(s) to perform on introgression file: lift, fgap, fcen, rmbn",
        required=True,
    )
    parser.add_argument("--ref", type=str, help="name of reference in Panagram")
    parser.add_argument(
        "--map",
        type=str,
        help="minimap flags to use",
        default="-x asm20 -c -t 3",
    )
    parser.add_argument(
        "--paf",
        type=str,
        help="path to minimap paf file or folder for liftover to ref space",
    )
    parser.add_argument(
        "--bin",
        type=int,
        help="size of bitmap bin used during calling",
        default=1000000,
    )
    args = parser.parse_args()

    actions = args.a
    for action in actions:
        if action not in ["lift", "fgap", "fcen", "rmbn"]:
            raise ValueError(f"Unrecognized action {action}. Check -a flag for valid actions.")

    index_dir = Path(args.idx)
    if not index_dir.is_dir():
        raise ValueError("Index directory not found. Check --idx path.")
    index = Index(index_dir)

    # figure out if user provided a file or folder
    bed_files = Path(args.bed)
    if bed_files.is_file():
        bed_files = [bed_files]
    elif bed_files.is_dir():
        bed_files = bed_files.glob("*.bed")
    else:
        raise ValueError("Bed file/folder not found. Check --bed path.")

    for bed_file in bed_files:
        # TODO: handle empty bed files by warning and outputting empty processed bed file and skipping
        # accession is inferred using the bed file name
        bed_accession = bed_file.name.split("_")[0]
        bed_chr = bed_file.name.split("_")[1]
        bed_intro_type = bed_file.stem.split("_")[2]
        bed_genome = index.genomes[bed_accession]
        bed_output = bed_file.parent / (bed_file.stem + f"_postprocessed.bed")

        for action in actions:
            if action == "lift":
                # lift - perform liftover - requires ref; paf optional
                if args.ref is None:
                    raise ValueError("--ref must be specified for liftover.")
                reference_accession = args.ref
                reference_genome = index.genomes[reference_accession]
                reference_file = index_dir / reference_genome.fasta
                if not reference_file.is_file():
                    raise ValueError("Reference file not found. Check --ref path.")

                # run minimap if paf files are not already provided
                if args.paf is None:
                    # find FASTA file associated with bed file
                    query_file = index_dir / bed_genome.fasta
                    minimap_flags = args.map.split(" ")
                    output_dir = index_dir / "alignments"
                    output_dir.mkdir(parents=True, exist_ok=True)
                    paf_file = output_dir / f"{bed_accession}_{reference_accession}.paf"
                    align_to_reference(reference_file, query_file, minimap_flags, paf_file)
                else:
                    paf_file = Path(args.paf)
                    # find paf file with the same name of the bed file
                    if paf_file.is_dir():
                        paf_file = paf_file / f"{bed_accession}_{reference_accession}.paf"
                    if not paf_file.is_file():
                        raise ValueError("Cannot find paf file. Check --paf path.")

                # save bed file lifted over to reference coordinate space
                liftover_to_reference(bed_file, paf_file, bed_output)
                bed_file = bed_output
                # bed genome is now the reference since we are in reference space
                bed_genome = reference_genome

            elif action in ["fgap", "rmbn"]:
                # fgap - merge regions >=1 bin size apart - requires bin size
                # rmbn - perform singleton bin removal - requires bin size
                # similar actions - just differ in function applied to the row of bins
                apply_func = fill_gaps
                if action == "rmbn":
                    apply_func = remove_single_bins

                bin_size = args.bin
                chr_length = bed_genome.sizes[bed_chr]
                bins_df = bed_to_bins(bed_file, bin_size, chr_length)
                bins_df["introgression"] = bins_df.apply(apply_func, axis=0)

                # save
                bed_df = bins_to_bed(bins_df, bin_size, bed_chr, bed_intro_type)
                bed_df.to_csv(bed_output, header=False, index=False, sep="\t")
                bed_file = bed_output

            elif action == "fcen":
                # fcen - fill in centromeres between introgression regions - requires bin size
                fasta_file = index_dir / bed_genome.fasta
                bed_df, fasta_df = read_bed_file(bed_file, fasta_file=fasta_file)
                bed_df = merge_centromere_regions(bed_df, fasta_df, bin_size)

                # save
                bed_df.to_csv(bed_output, header=False, index=False, sep="\t")
                bed_file = bed_output
    return


if __name__ == "__main__":
    postprocess_introgressions()

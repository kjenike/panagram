import argparse
from pathlib import Path
import subprocess
import math
import numpy as np
import pandas as pd
from panagram.index import Index
from call_introgressions import bins_to_bed
from concurrent.futures import ProcessPoolExecutor, as_completed


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


def read_bed_file(bed_file):
    """Read in a given bed file.

    :param bed_file: path to bed file with introgression locations
    :type bed_file: str or Path
    :return: dataframe representing bed file
    :rtype: pd.DataFrame
    """

    if bed_file_is_empty(bed_file):
        return None

    bed_df = pd.read_csv(
        bed_file,
        sep="\t",
        header=None,
    )

    bed_df = bed_df.iloc[:, 0:4]
    col_names = ["Chromosome", "Start", "End", "Notes"]
    bed_df.columns = col_names
    bed_df["Sequence"] = None

    return bed_df


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


def bed_to_bins(bed_df, bin_size, chr_length):
    # convert bedfile coordinates back into an index of bins for a chromosome
    # match format of txt files from original introgressions script
    # with a per-bin 0/1 introgression score

    # read in bed file
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
        # clip accidental extra bins
        if bin_label > intro_df.columns[-1]:
            continue
        intro_df.loc[:, bin_label] = 1

    # match the format of the bins_df from call_introgressions
    intro_df = intro_df.T
    intro_df.columns = ["introgression"]
    return intro_df


def fill_gaps(row, gap_size):
    # take gap size argument; count 0s between intros; if the gap is <= gap_size, fill it
    # recommend running fill gaps prior to remove_small_regions after liftover, which can fragment intros
    row = row.values.copy()
    i = 0
    while i < len(row):
        # find an introgression
        if row[i] == 1:
            # run past the introgression into the 0 space
            while i < len(row) and row[i] == 1:
                i += 1
            # define start of 0 space between 2 intros
            region_start = i
            # measure the length of the 0 spaces
            while i < len(row) and row[i] == 0:
                i += 1
            # 0 space ends when we hit another 1 (start of a new intro)
            region_end = i

            # fill the gap if it is small enough
            region_length = region_end - region_start
            if region_length <= gap_size:
                row[region_start:region_end] = 1
        else:
            i += 1
    return row


def remove_small_regions(row, min_size):
    # omit introgressions that are smaller than the minimum size
    # recommend running fill gaps prior to this func after liftover, which can fragment intros
    row = row.values.copy()
    i = 0
    while i < len(row):
        if row[i] == 1:
            region_start = i
            # figure out how many 1s in a row there are
            while i < len(row) and row[i] == 1:
                i += 1
            region_end = i
            region_length = region_end - region_start

            if region_length < min_size:
                row[region_start:region_end] = 0
        else:
            i += 1

    return row


def postprocess_introgressions():
    # does postprocessing on specified file or all files in a folder
    # allows for liftover, gap filling, centromere filling, and lone bin removal in any order
    parser = argparse.ArgumentParser(description="Introgression postprocessing.")
    parser.add_argument(
        "--bed", type=str, help="path to introgression bed file or folder", required=True
    )
    parser.add_argument("--idx", type=str, help="path to Panagram index folder", required=True)
    parser.add_argument("--out", type=str, help="path to folder to save all outputs", required=True)
    parser.add_argument(
        "--act",
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
    print("Postprocessing introgressions...", flush=True)
    args = parser.parse_args()

    actions = args.act
    for action in actions:
        if action not in ["lift", "fgap", "fcen", "rmbn"]:
            raise ValueError(f"Unrecognized action {action}. Check --act flag for valid actions.")

    index_dir = Path(args.idx)
    if not index_dir.is_dir():
        raise ValueError("Index directory not found. Check --idx path.")
    index = Index(index_dir)

    # figure out if user provided a file or folder
    bed_files = Path(args.bed)
    if bed_files.is_file():
        bed_files = [bed_files]
    elif bed_files.is_dir():
        bed_files = list(bed_files.glob("*.bed"))
    else:
        raise ValueError("Bed file/folder not found. Check --bed path.")

    # perform liftover per-accession, so we can use the same results for all files
    if "lift" in actions:
        if args.ref is None:
            raise ValueError("--ref must be specified for liftover.")
        reference_accession = args.ref
        reference_genome = index.genomes[reference_accession]
        reference_file = index_dir / reference_genome.fasta
        if not reference_file.is_file():
            raise ValueError("Reference file not found. Check --ref path.")

        # run minimap in parallel for each accession if paf files are not already provided
        if args.paf is None:
            # Get unique bed_accessions from bed_files
            unique_accessions = set()
            for bed_file in bed_files:
                bed_accession = bed_file.name.split("_")[0]
                unique_accessions.add(bed_accession)

            paf_dir = index_dir / "alignments"
            paf_dir.mkdir(parents=True, exist_ok=True)
            minimap_flags = args.map.split(" ")

            with ProcessPoolExecutor() as executor:
                futures = []
                for accession in unique_accessions:
                    bed_genome = index.genomes[accession]
                    query_file = index_dir / bed_genome.fasta
                    paf_file = paf_dir / f"{accession}_{reference_accession}.paf"
                    futures.append(
                        executor.submit(
                            align_to_reference, reference_file, query_file, minimap_flags, paf_file
                        )
                    )

                for future in as_completed(futures):
                    future.result()
        else:
            paf_dir = Path(args.paf)

    output_dir = Path(args.out)
    output_dir.mkdir(parents=True, exist_ok=True)

    for bed_file in bed_files:
        # information is inferred using the bed file name
        bed_accession = bed_file.name.split("_")[0]
        bed_chr = bed_file.name.split("_")[1]
        bed_intro_type = bed_file.stem.split("_")[2]
        bed_genome = index.genomes[bed_accession]
        bed_output = output_dir / bed_file.name

        for action in actions:
            # handle empty bed files by warning and outputting empty processed bed file and skipping
            # note that liftover and rmbn can produce empty bed files, so we have to check after each action
            if bed_file_is_empty(bed_file):
                # print(f"Warning: {bed_file.name} is empty. Outputting empty postprocessed bed file.")
                bed_output.touch()
                break

            if action == "lift":
                # lift - perform liftover - requires ref; paf optional
                # find paf file with the same name of the bed file
                if paf_dir.is_dir():
                    paf_file = paf_dir / f"{bed_accession}_{reference_accession}.paf"
                elif paf_dir.is_file():
                    paf_file = paf_dir
                else:
                    raise ValueError("Cannot find paf file/folder. Check --paf path.")

                # save bed file lifted over to reference coordinate space
                liftover_to_reference(bed_file, paf_file, bed_output)
                bed_file = bed_output
                # bed genome is now the reference since we are in reference space
                bed_genome = reference_genome

            elif action == "fgap":
                # fgap - merge regions >=1 bin size apart - requires bin size
                bin_size = args.bin
                chr_length = bed_genome.sizes[bed_chr]
                gap_size = args.gap
                bed_df = read_bed_file(bed_file)
                bins_df = bed_to_bins(bed_df, bin_size, chr_length)
                bins_df["introgression"] = bins_df.apply(fill_gaps, gap_size=gap_size, axis=0)

                # save
                bed_df = bins_to_bed(bins_df, bin_size, bed_chr, bed_intro_type)
                bed_df.to_csv(bed_output, header=False, index=False, sep="\t")
                bed_file = bed_output

            elif action == "rmbn":
                # rmbn - perform singleton bin removal - requires bin size, min. introgression size
                bin_size = args.bin
                chr_length = bed_genome.sizes[bed_chr]
                min_size = args.min
                bed_df = read_bed_file(bed_file)
                bins_df = bed_to_bins(bed_df, bin_size, chr_length)
                bins_df["introgression"] = bins_df.apply(
                    remove_small_regions, min_size=min_size, axis=0
                )

                # save
                bed_df = bins_to_bed(bins_df, bin_size, bed_chr, bed_intro_type)
                bed_df.to_csv(bed_output, header=False, index=False, sep="\t")
                bed_file = bed_output

            elif action == "fcen":
                # fcen - fill in centromeres between introgression regions - requires bin size
                fasta_file = index_dir / bed_genome.fasta
                fasta_df = read_fasta(fasta_file)
                bed_df = read_bed_file(bed_file)
                bed_df = merge_centromere_regions(bed_df, fasta_df, bin_size)

                # save
                bed_df.to_csv(bed_output, header=False, index=False, sep="\t")
                bed_file = bed_output
    return


if __name__ == "__main__":
    postprocess_introgressions()

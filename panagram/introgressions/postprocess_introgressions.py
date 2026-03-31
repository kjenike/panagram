import argparse
from pathlib import Path
import shutil
from io import StringIO
import subprocess
import math
import numpy as np
import pandas as pd
from panagram.index import Index
from panagram.introgressions.call_introgressions import bins_to_bed
from concurrent.futures import ProcessPoolExecutor, as_completed


def read_fasta_generator(fp):
    """Generator for fasta file - allows reading of one record at a time to save memory.
    Or, all records can be read at once using read_fasta().

    Args:
        fp (_io.TextIOWrapper): file pointer

    Yields:
        tuple: name and seq in the file
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

    Args:
        input_file (str): FASTA file path

    Returns:
        pd.DataFrame: dataframe of FASTA records with "Name" and "Sequence" columns
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
    """Grab the sequence corresponding to the bed file coordinates listed. Follows bedtools getfasta
    logic (0-index, inclusive start coord, exclusive end coord).

    Args:
        bed_row (pd.DataFrame): row of a bed file's dataframe
        fasta_df (pd.DataFrame): fasta information containing sequences that the bed file references

    Returns:
        str: sequence corresponding to bed coordinates
    """

    chromosome_sequence = fasta_df.loc[fasta_df.Name == bed_row["Chromosome"], "Sequence"].iloc[0]
    start = int(bed_row["Start"])
    end = int(bed_row["End"])
    introgression_sequence = chromosome_sequence[start:end]
    return introgression_sequence


def bed_file_is_empty(bed_file):
    """Check if a bed file is empty.

    Args:
        bed_file (str or Path): path to bed file

    Returns:
        bool: True if bed file is empty, False otherwise
    """

    try:
        _ = pd.read_csv(
            bed_file,
            sep="\t",
            header=None,
        )
    except pd.errors.EmptyDataError:
        return True
    return False


def read_bed_file(bed_file):
    """Read in a given bed file.

    Args:
        bed_file (str or Path): path to bed file with introgression locations

    Returns:
        pd.DataFrame: dataframe representing bed file
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
    """Merge centromere regions in a bed dataframe.

    Args:
        bed_df (pd.DataFrame): DataFrame containing bed file information (Chromosome, Start, End, Notes, Sequence)
        fasta_df (pd.DataFrame): DataFrame containing fasta file information
        bin_size (int): size of the bins used for merging

    Returns:
        pd.DataFrame: DataFrame with centromere regions merged
    """

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


def align_assemblies(
    reference_file,
    query_file,
    minimap_flags,
    output_file,
    direction="query_to_reference",
):
    """Align query assembly to reference using minimap2.

    Args:
        reference_file (str): path to the reference genome file (FASTA format)
        query_file (str): path to the query genome file (FASTA format)
        minimap_flags (list): list of flags to pass to minimap2
        output_file (str): path to the output file (PAF format)
        direction (str, optional): direction of the alignment. Defaults to "query_to_reference"

    Raises:
        ValueError: If direction is not "query_to_reference" or "reference_to_query".
    """

    # Call minimap2 command - requires minimap2 to be in path
    if direction == "query_to_reference":
        minimap_command = ["minimap2"] + minimap_flags + [reference_file, query_file]
    elif direction == "reference_to_query":
        minimap_command = ["minimap2"] + minimap_flags + [query_file, reference_file]
    else:
        raise ValueError("Direction must be 'query_to_reference' or 'reference_to_query'")

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


def liftover_coordinates(
    bed_file,
    paf_file,
):
    """Liftover coordinates from one genome assembly to another. Requires paftools.js to be in path.

    Args:
        bed_file (str): Path to the BED file to liftover
        paf_file (str): Path to the PAF file indicating mapping between assemblies

    Returns:
        pd.DataFrame: DataFrame containing lifted coordinates in BED format (Chromosome, Start, End, Notes)
    """

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

    new_coordinates_df = pd.read_csv(
        StringIO(result.stdout),
        sep="\t",
        header=None,
        usecols=[0, 1, 2, 3],
        names=["Chromosome", "Start", "End", "Notes"],
    )

    return new_coordinates_df


def run_alignments(
    unique_accessions,
    paf,
    index_dir,
    index,
    reference_accession,
    reference_file,
    map_args,
    n_threads,
):
    """Run minimap2 alignments in parallel for each accession if paf files are not already provided.
    If paf files are provided, skip this step and return the paf directory.

    Args:
        unique_accessions (list): List of unique accessions for which to run alignments
        paf (str): Path to the PAF directory containing alignment information
        index_dir (Path): Path to the index directory
        index (Index): Index object containing genome information
        reference_accession (str): Accession name of the reference genome
        reference_file (Path): Path to the reference genome file
        map_args (str): String of arguments for minimap2
        n_threads (int): Number of threads to use for parallel processing

    Returns:
        Path: Path to the PAF directory containing alignment files
    """

    # run minimap in parallel for each accession if paf files are not already provided
    if paf is None:
        paf_dir = index_dir / "alignments"
        paf_dir.mkdir(parents=True, exist_ok=True)
        minimap_flags = map_args.split(" ")

        with ProcessPoolExecutor(max_workers=n_threads) as executor:
            futures = []
            for accession in unique_accessions:
                bed_genome = index.genomes[accession]
                query_file = index_dir / bed_genome.fasta
                paf_file = paf_dir / f"{accession}_{reference_accession}.paf"
                futures.append(
                    executor.submit(
                        align_assemblies, reference_file, query_file, minimap_flags, paf_file
                    )
                )

            for future in as_completed(futures):
                future.result()
    else:
        paf_dir = Path(paf)

    return paf_dir


def run_liftovers(
    bed_files,
    unique_accessions,
    paf_dir,
    reference_accession,
    output_dir,
):
    """Run liftover for each bed file to convert coordinates to reference space. If bed file is
    empty, save an empty liftover file. After liftover, concat all bed files together, split by
    chromosome and introgression type, and save.

    Args:
        bed_files (list): List of paths to BED files to liftover
        unique_accessions (list): List of unique accession names
        paf_dir (Path): Path to the PAF directory containing alignment files
        reference_accession (str): Accession name of the reference genome
        output_dir (Path): Path to the output directory

    Raises:
        ValueError: If paf file corresponding to a bed file is not found in the paf directory.

    Returns:
        list: List of paths to the lifted BED files
    """

    # NOTE: we might end up with more bed files than we started if 1 chr splits into multiple
    # could mess with scoring if you only intended to score a subset of chromosomes
    bed_df = pd.DataFrame(
        columns=["Chromosome", "Start", "End", "Notes", "Accession", "Introgression Type"]
    )
    for bed_file in bed_files:
        bed_chr, bed_accession, bed_intro_type = get_bed_pieces(bed_file)

        # find paf file with the same name of the bed file
        if paf_dir.is_dir():
            paf_file = paf_dir / f"{bed_accession}_{reference_accession}.paf"
        else:
            raise ValueError(
                f"Cannot find paf file {paf_dir / f'{bed_accession}_{reference_accession}.paf'}. Check --paf path."
            )

        # check for empty bed files
        if bed_file_is_empty(bed_file):
            # create placeholder row for new coordinates
            new_coordinates_df = pd.DataFrame(
                {"Chromosome": [bed_chr], "Start": [0], "End": [0], "Notes": ["placeholder"]}
            )
        # liftover if needed
        elif bed_intro_type != "REF":
            new_coordinates_df = liftover_coordinates(bed_file, paf_file)
        # if no liftover needed, just read in the bed file
        else:
            new_coordinates_df = read_bed_file(bed_file)
            new_coordinates_df = new_coordinates_df[["Chromosome", "Start", "End", "Notes"]]

        # change intro type to REF since now introgressions are in reference space
        if bed_intro_type == "REFA":
            bed_intro_type = "REF"

        new_coordinates_df["Accession"] = bed_accession
        new_coordinates_df["Introgression_Type"] = bed_intro_type
        bed_df = pd.concat([bed_df, new_coordinates_df], ignore_index=True)

    new_bed_files = []
    # concat liftover results together - at minimum, there should be a file for each file that we started with, even if empty
    for accession in unique_accessions:
        # get all results for each chr, remove any placeholders, sort, and save
        accession_df = bed_df[bed_df.Accession == accession]

        chromosomes = accession_df.Chromosome.unique()
        for chromosome in chromosomes:
            introgression_types = accession_df.Introgression_Type.unique()
            for introgression_type in introgression_types:
                chromosome_df = accession_df[
                    (accession_df.Chromosome == chromosome)
                    & (accession_df.Introgression_Type == introgression_type)
                ]
                chromosome_df = chromosome_df[chromosome_df.Notes != "placeholder"]
                chromosome_df = chromosome_df.sort_values(by="Start")

                # save
                bed_output = output_dir / f"{accession}_{chromosome}_{introgression_type}.bed"
                chromosome_df["Introgression_Type"] = chromosome_df["Introgression_Type"] + "_intro"
                chromosome_df[["Chromosome", "Start", "End", "Introgression_Type"]].to_csv(
                    bed_output, header=False, index=False, sep="\t"
                )
                # update bed_files list to point to lifted files for downstream actions
                new_bed_files.append(bed_output)

    return new_bed_files


def get_intro_df_template(bin_size, chr_length):
    """Create a template DataFrame with the correct number of bins for introgressions.

    Args:
        bin_size (int): size of the bins used for introgressions
        chr_length (int): length of the chromosome

    Returns:
        pd.DataFrame: template DataFrame with a row for introgressions
    """

    # create pandas df full of 0s
    n_bins = math.ceil(chr_length / bin_size)
    bin_names = [i * bin_size for i in range(n_bins)]
    zeros = np.zeros((1, n_bins), dtype=int)
    intro_df = pd.DataFrame(zeros, columns=bin_names)
    return intro_df


def bed_to_bins(bed_df, bin_size, chr_length):
    """Convert BED file coordinates to bin indices. Bins with introgressions are marked with a 1.
    Rounds coordinates to nearest bin to reduce number of extra bins marked as introgressions.

    Args:
        bed_df (pd.DataFrame): DataFrame containing BED file coordinates
        bin_size (int): size of the bins used for introgressions
        chr_length (int): length of the chromosome

    Returns:
        pd.DataFrame: DataFrame with an "introgression" column and bin indices as rows.
    """

    # read in bed file
    intro_df = get_intro_df_template(bin_size, chr_length)

    # convert coordinates to bin labels
    # round to nearest bin - reduces number of extra bins marked as introgressions
    bed_df["Start Bin"] = ((bed_df["Start"] / bin_size).round() * bin_size).astype(int)
    bed_df["End Bin"] = ((bed_df["End"] / bin_size).round() * bin_size).astype(int)

    # alternative rounding methods - tend to create extra bins:
    # bed_df["Start Bin"] = (bed_df["Start"] // bin_size) * bin_size  # round down to nearest bin
    # bed_df["End Bin"] = ((bed_df["End"] // bin_size) * bin_size) + bin_size  # round up

    # get all column labels between start and end bin
    bed_df["Bin Labels"] = bed_df.apply(
        lambda row: list(range(row["Start Bin"], row["End Bin"], bin_size)), axis=1
    )

    # if no bin labels (e.g. if start and end are in the same bin),
    # check if the start/end are at least 25% of the bin size apart, and if so, mark the start bin as an introgression
    bed_df.loc[
        (bed_df["Bin Labels"].str.len() == 0)
        & ((bed_df["End"] - bed_df["Start"]) >= (bin_size / 4)),
        "Bin Labels",
    ] = bed_df["Start Bin"].apply(lambda x: [x])
    # alternatively, add the start bin as the bin label every time
    # bed_df.loc[bed_df["Bin Labels"].str.len() == 0, "Bin Labels"] = bed_df["Start Bin"].apply(lambda x: [x])

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
    """Fill gaps between introgressions in a row. Running fill gaps is recommended prior to running
    remove_small_regions after liftover, which can fragment introgressions.

    Args:
        row (pd.Series): a row of the DataFrame containing introgression data
        gap_size (int): the maximum gap size to fill

    Returns:
        list: the modified row with filled gaps
    """

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
    """Omit introgressions that are smaller than the minimum size. Run fill gaps prior to using this
    function after liftover, which can fragment introgressions.

    Args:
        row (pd.Series): a row of the DataFrame containing introgression data
        min_size (int): the minimum size an introgression must be to avoid removal

    Returns:
        list: the modified row with small introgressions removed
    """

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


def get_bed_pieces(bed_file):
    """Get the chromosome, accession, and introgression type from a bed file name. Assumes bed file
    names are in the format "accession_chr_introgressiontype.bed".

    Args:
        bed_file (str or Path): path to bed file

    Returns:
        tuple: chromosome, accession, and introgression type
    """

    parts = bed_file.stem.split("_")
    bed_intro_type = parts[-1]
    bed_chr = parts[-2]
    bed_accession = "_".join(parts[:-2])
    return bed_chr, bed_accession, bed_intro_type


def postprocess_introgressions():
    """Postprocess introgression bed files with various actions: liftover, fill gaps, fill
    centromeres, remove small regions.
    """

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
        default="-x asm20 -c -t 1",  # preset for cross-species full-genome alignment
    )
    parser.add_argument(
        "--threads",
        type=int,
        help="number of alignments to run in parallel for liftover",
        default=1,
    )
    parser.add_argument(
        "--paf",
        type=str,
        help="path to folder with minimap paf files for liftover to ref space",
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
    if actions:
        for action in actions:
            if action not in ["lift", "fgap", "fcen", "rmbn"]:
                raise ValueError(
                    f"Unrecognized action {action}. Check --act flag for valid actions."
                )

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

    output_dir = Path(args.out)
    output_dir.mkdir(parents=True, exist_ok=True)

    if "lift" in actions:
        # check for reference file
        if args.ref is None:
            raise ValueError("--ref must be specified for liftover.")
        reference_accession = args.ref
        reference_genome = index.genomes[reference_accession]
        reference_file = index_dir / reference_genome.fasta
        if not reference_file.is_file():
            raise ValueError("Reference file not found. Check --ref path.")

        # Get unique bed_accessions from bed_files
        unique_accessions = set()
        for bed_file in bed_files:
            bed_accession = bed_file.name.split("_")[0]
            unique_accessions.add(bed_accession)

        # perform alignment if needed
        paf_dir = run_alignments(
            unique_accessions,
            args.paf,
            index_dir,
            index,
            reference_accession,
            reference_file,
            args.map,
            args.threads,
        )

        # perform liftover
        bed_files = run_liftovers(
            bed_files,
            unique_accessions,
            paf_dir,
            reference_accession,
            output_dir,
        )

    for bed_file in bed_files:
        bed_chr, bed_accession, bed_intro_type = get_bed_pieces(bed_file)
        bed_genome = index.genomes[bed_accession]
        if "lift" in actions:
            bed_genome = index.genomes[args.ref]
        if bed_intro_type == "REF":
            if args.ref is None:
                raise ValueError("--ref must be specified if REF files are present.")
            bed_genome = index.genomes[args.ref]
        bed_output = output_dir / bed_file.name

        if actions:
            for action in actions:
                # handle empty bed files by warning and outputting empty processed bed file and skipping
                # note that rmbn can produce empty bed files, so we have to check after each action
                if bed_file_is_empty(bed_file):
                    bed_output.touch()
                    break

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
        else:
            # copy file to output directory if no actions specified
            shutil.copy(bed_file, bed_output)
    return


if __name__ == "__main__":
    postprocess_introgressions()

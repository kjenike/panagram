# take in:
# samples.tsv
# a folder where bed files live
# a chromosome

# extract introgressions using bed files
# convert to kmer dictionary
# while doing this, append kmers to a set
# convert all dicts to vectors of kmer counts
# cluster vectors using umap/tsne/dbscan/etc.
# color by acession

from pathlib import Path
import subprocess
import numpy as np
import pandas as pd
from panagram.index import Index


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

    try:
        bed_df = pd.read_csv(
            bed_file,
            sep="\t",
            header=None,
        )
    except pd.errors.EmptyDataError as e:
        print("Bed file is empty. Returning None.")
        return None, None
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


def liftover_to_reference(
    bed_file,
    paf_file,
    output_file,
    paftools_script_path="/home/nbrown62/data_mschatz1/nbrown62/minimap2/misc/paftools.js",
):
    # Call liftover script - requires k8 to be in path
    liftover_command = [
        "k8",
        paftools_script_path,
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

    return intro_df


def get_intro_df_template(bin_size, chr_length):
    # create pandas df full of 0s
    n_bins = chr_length // bin_size
    bin_names = [i * bin_size for i in range(n_bins)]
    zeros = np.zeros((1, n_bins), dtype=int)
    intro_df = pd.DataFrame(zeros, columns=bin_names)
    return intro_df


def read_introgressions(introgression_file, fix_names=False):
    # read introgressions for a chromosome
    intro_df = pd.read_csv(introgression_file, sep="\t", header=0, index_col=0).fillna(0)

    if fix_names:
        new_index = list(intro_df.index)
        # changing a specific name
        new_index = ["Fla8924" if i == "Fla.8924.ont.s" else i.split(".")[0] for i in new_index]
        intro_df.index = new_index
    return intro_df


def threshold_introgressions_helper(intro_df, threshold):
    intro_df[intro_df < threshold] = 0
    intro_df[intro_df != 0] = 1
    return intro_df


def threshold_introgressions(called_intro_file, gt_intro_file, threshold):
    # if given multiple gt files, merge together
    if type(gt_intro_file) == list or type(gt_intro_file) == tuple:
        gt_intro_df1 = read_introgressions(gt_intro_file[0], fix_names=True)
        gt_intro_df2 = read_introgressions(gt_intro_file[1], fix_names=True)
        # threshold and merge multiple intro types
        gt_intro_df1 = threshold_introgressions_helper(gt_intro_df1, threshold)
        gt_intro_df2 = threshold_introgressions_helper(gt_intro_df2, threshold)
        gt_intro_df = gt_intro_df1 + gt_intro_df2
        gt_intro_df[gt_intro_df > 0] = 1
    else:
        gt_intro_df = read_introgressions(gt_intro_file, fix_names=True)
        threshold_introgressions_helper(gt_intro_df, threshold)

    # make sure all called introgressions are scored as a 1
    called_intro_df = read_introgressions(called_intro_file)
    called_intro_df = threshold_introgressions_helper(called_intro_df, threshold=1)
    return called_intro_df, gt_intro_df


def score_introgressions(called_intro_file, gt_intro_file, threshold):
    # get confusion matrix for introgressions given ground truth in the same coordinate/bin space
    called_intro_df, gt_intro_df = threshold_introgressions(
        called_intro_file, gt_intro_file, threshold
    )

    # rotate dfs, sort cols, and drop all columns that aren't shared btwn called and gt
    shared_cols = list(set(called_intro_df.index).intersection(set(gt_intro_df.index)))
    called_intro_df = called_intro_df.transpose()[shared_cols]
    gt_intro_df = gt_intro_df.transpose()[shared_cols]

    # calculate confusion matrix
    total = gt_intro_df.size

    # true pos
    true_pos = ((called_intro_df == 1) & (gt_intro_df == 1)).values.sum()

    # true neg
    true_neg = ((called_intro_df == 0) & (gt_intro_df == 0)).values.sum()

    # false pos
    false_pos = ((called_intro_df == 1) & (gt_intro_df == 0)).values.sum()

    # false neg
    false_neg = ((called_intro_df == 0) & (gt_intro_df == 1)).values.sum()

    # accuracy
    acc = (true_pos + true_neg) / total
    precision = true_pos / (true_pos + false_pos)
    recall = true_pos / (true_pos + false_neg)

    # return the metrics as a df
    metrics = {
        "True Positive": true_pos,
        "True Negative": true_neg,
        "False Positive": false_pos,
        "False Negative": false_neg,
        "Accuracy": acc,
        "Precision": precision,
        "Recall": recall,
    }
    metrics = pd.DataFrame([metrics])

    return metrics


def score_all_introgressions():
    # NOTE: change parameters here
    index_dir = Path("/home/nbrown62/data_mschatz1/nbrown62/panagram_data/tomato_sl4")
    index = Index(index_dir)
    reference = "SL4"
    genome = index.genomes[reference]
    introgression_type = "REF"
    bin_size = 1000000
    threshold = 0.5
    samples_file = Path(
        "/home/nbrown62/data_mschatz1/nbrown62/panagram_data/tomato_sl4/samples.tsv"
    )
    output_dir = Path(
        "/home/nbrown62/data_mschatz1/nbrown62/panagram_data/tomato_sl4/introgression_analysis_v2/postprocessed"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    all_metrics_df = None
    for chr_name in genome.sizes.keys():
        chr_intro_df = None
        chr_length = genome.sizes[chr_name]

        # NOTE: change paths to point to correct folders
        gt_introgression_sp = Path(
            f"/home/nbrown62/data_mschatz1/nbrown62/CallIntrogressions_data/tomato_sl4_paper/{chr_name}.SP.txt"
        )
        gt_introgression_slc = Path(
            f"/home/nbrown62/data_mschatz1/nbrown62/CallIntrogressions_data/tomato_sl4_paper/{chr_name}.SLC.txt"
        )

        # figure out which ground truth to use
        if introgression_type == "SP":
            gt_intro_file = gt_introgression_sp
        elif introgression_type == "SLC":
            gt_intro_file = gt_introgression_slc
        elif introgression_type in ["REF", "merged"]:
            gt_intro_file = [gt_introgression_slc, gt_introgression_sp]
        else:
            raise ValueError("Invalid introgression type selected.")

        for anchor in index.genomes.keys():
            print(f"Processing {chr_name}, {anchor}")

            bed_file = Path(
                f"/home/nbrown62/data_mschatz1/nbrown62/panagram_data/tomato_sl4/introgression_analysis_v2/{anchor}_{chr_name}_{introgression_type}.bed"
            )
            paf_file = Path(
                f"/home/nbrown62/data_mschatz1/nbrown62/minimap2_data/tomato_sl4/{anchor}_edited_{reference}_edited.paf"
            )

            # load bed file - get sequences and merge centromeres
            bed_df, fasta_df = read_bed_file(bed_file, anchor, samples_file, index_dir)
            if bed_df is None:
                anchor_intro_df = get_intro_df_template(bin_size, chr_length)
            else:
                bed_df = merge_centromere_regions(bed_df, fasta_df, bin_size)

                # get bed file path and name without extension
                bed_file_merged_centromeres = bed_file.stem + "_with_centromeres.bed"
                bed_file_merged_centromeres = output_dir / bed_file_merged_centromeres

                # save bed file with merged centromeres
                bed_df[["Chromosome", "Start", "End", "Notes"]].to_csv(
                    bed_file_merged_centromeres, index=False, header=False, sep="\t"
                )

                # save bed file lifted over to SL4 coordinate space
                bed_file_liftover = bed_file.stem + f"_liftover_{reference}.bed"
                bed_file_liftover = output_dir / bed_file_liftover
                liftover_to_reference(str(bed_file_merged_centromeres), paf_file, bed_file_liftover)

                # save introgressions as bins
                anchor_intro_df = bed_to_bins(bed_file_liftover, bin_size, chr_length)

            anchor_intro_df.index = [anchor]
            anchor_intro_df = anchor_intro_df.rename_axis("Sample")

            if chr_intro_df is None:
                chr_intro_df = anchor_intro_df
            else:
                # add introgressions to df with other anchors
                chr_intro_df = pd.concat([chr_intro_df, anchor_intro_df.iloc[[0]]])

        # save off called introgressions
        chr_intro_file = output_dir / f"{chr_name}.{introgression_type}.txt"
        chr_intro_df.to_csv(chr_intro_file, sep="\t")

        metrics = score_introgressions(chr_intro_file, gt_intro_file, threshold)
        metrics.index = [chr_name]
        metrics_file = output_dir / f"metrics_{chr_name}_{introgression_type}.tsv"
        metrics.to_csv(metrics_file, sep="\t", index=False)

        if all_metrics_df is None:
            all_metrics_df = metrics
        else:
            all_metrics_df = pd.concat([all_metrics_df, metrics.iloc[[0]]])

    all_metrics_file = output_dir / f"all_metrics_{introgression_type}.tsv"
    all_metrics_df.to_csv(all_metrics_file, sep="\t")
    return


if __name__ == "__main__":
    score_all_introgressions()

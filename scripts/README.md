# Panagram's Auxiliary Scripts

This folder contains extra scripts used by some Panagram analyses. Most of these scripts are used by
Panagram directly and do not need to be called separately. Others can be used directly, such as the
association testing script described below.

## Association Testing

Panagram can be used to test for associations between user-provided variables and regions of the
bitmap. For each gene or exon in a GFF file, samples in the pangenome are clustered into two
haplotype groups based on the pattern of their bitmaps across that interval. Each interval is then
tested to determine whether haplotype group membership is associated with a phenotype. The model
accounts for relatedness among samples using a genome-wide relationship matrix. Intervals with
haplotype groups that are too small are excluded, and p-values are corrected for multiple testing
using false discovery rate correction.

Run the full Panagram association testing pipeline with `association_testing.py`.

This script:

1. Builds a genome-wide bitmap tree.
2. Creates gene or exon haplotype groups for all chromosomes/contigs.
3. Runs association testing per chromosome/contig.
4. Combines association results.
5. Applies FDR correction.
6. Makes volcano plots.

Chromosome-level haplotype querying and association testing are run in parallel. On a pangenome of
51 dog samples with 4 phenotypes to test, the script took a few hours to run with 40 threads. Note
that this script parses GFF files, which are not fully standardized; you may need to edit the script
to parse gene/exon names correctly. There are `NOTE` sections in the script to denote functions that
may need editing.

### Usage

```bash
python association_testing.py \
    -a ANCHORNAME \
    -p PANAGRAM_INDEX_DIR \
    -g GFF_FILE \
    -s SAMPLE_FILE \
    -m PHENOTYPE_FILE \
    -t TYPE \
    -o OUTPUT_DIR \
    -n N_THREADS
```

### Required arguments

| Short flag | Long flag          | Description                                                                                           |
| ---------- | ------------------ | ----------------------------------------------------------------------------------------------------- |
| `-a`       | `--anchor-name`    | Anchor genome name in Panagram's samples.tsv. The bitmap will be used from this anchor's perspective. |
| `-p`       | `--panagram-dir`   | Path to the Panagram index directory.                                                                 |
| `-g`       | `--gff`            | Path to the GFF annotation file for your anchor.                                                      |
| `-s`       | `--sample-file`    | Path to the Panagram samples.tsv file.                                                                |
| `-m`       | `--phenotype-file` | Path to the phenotype CSV file. See format below.                                                     |
| `-t`       | `--type`           | Feature type to test. Must be either `gene` or `exon`.                                                |
| `-o`       | `--output-dir`     | Output directory for all results.                                                                     |
| `-n`       | `--threads`        | Number of chromosomes/contigs to process in parallel.                                                 |

Example Phenotype file:

The phenotype file should be a CSV file with samples as rows and phenotypes as columns. The first
column should contain sample IDs. All phenotypes should be numeric. Categorical variables should be
encoded numerically, for example as 0, 1, 2, and so on.

```csv
sample_id,height,flowering_time
sample1,12.4,35
sample2,15.1,42
sample3,11.8,37
```

### Output structure

The output directory will contain:

```text
OUTPUT_DIR/
├── tree/
├── haplotypes/
├── associations/
└── plots/
```

#### Tree

```text
OUTPUT_DIR/tree/bitmap_tree.nwk
```

Genome-wide bitmap tree in Newick format.

#### Haplotype groups

```text
OUTPUT_DIR/haplotypes/{chromosome}_{type}_haplotypes.tsv
```

One haplotype file is produced per chromosome/contig.

Columns:

```text
gene_or_exon_id    sample    haplotype_group
```

Example:

```text
GeneA    sample1    0
GeneA    sample2    1
GeneA    sample3    0
```

For exon features, repeated exon entries from the same gene are numbered sequentially by appending
`_1`, `_2`, etc. to the gene name.

#### Association files

```text
OUTPUT_DIR/associations/{chromosome}_{type}_associations.tsv
```

One association result file is produced per chromosome/contig.

Columns:

```text
gene_name    phenotype_name    group1_size    group2_size    group1_mean    group2_mean    p_value
```

#### Final Outputs

Combined volcano plot:

```text
OUTPUT_DIR/plots/{type}_volcano.png
```

Sorted FDR-corrected association table:

```text
OUTPUT_DIR/plots/{type}_volcano_sorted.tsv
```

Individual phenotype volcano plots:

```text
OUTPUT_DIR/plots/{phenotype_name}_volcano.png
```

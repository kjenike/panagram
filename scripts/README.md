# Panagram's Auxiliary Scripts

This folder contains extra scripts used by some Panagram analyses. Most of these scripts are used by
Panagram directly and do not need to be called separately. Others can be used directly, such as the
association testing script described below.

## Association Testing

Panagram can be used to test for associations between user-provided variables and local regions of
the bitmap. For each tested feature or window, samples in the pangenome are clustered into two
haplotype groups based on the pattern of their bitmaps across that interval. Each interval is then
tested to determine whether haplotype group membership is associated with a phenotype. The model
accounts for relatedness among samples using a genome-wide relationship matrix. Intervals with
haplotype groups that are too small are excluded, and p-values are corrected for multiple testing
using false discovery rate correction.

Run the full Panagram association testing pipeline with `association_testing.py`.

This script:

1. Builds a genome-wide bitmap tree.
2. Creates haplotype groups for genes, exons, regulatory windows, or genome-wide windows.
3. Runs association testing in parallel across feature/window chunks.
4. Combines association results.
5. Applies FDR correction.
6. Makes volcano plots.

Feature/window-level haplotype querying and association testing are run in parallel. On a pangenome
of 51 dog samples with 4 phenotypes to test, the script took a few hours to run with 40 threads for
gene/exon-scale analyses. Whole-genome window analyses may take longer because they can involve tens
of millions of windows. Note that this script parses GFF files, which are not fully standardized;
you may need to edit the script to parse gene/exon names correctly. There are `NOTE` sections in the
script to denote functions that may need editing.

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
    -n N_THREADS \
    --chunk-size CHUNK_SIZE
```

### Required arguments

| Short flag | Long flag          | Description                                                                                           |
| ---------- | ------------------ | ----------------------------------------------------------------------------------------------------- |
| `-a`       | `--anchor-name`    | Anchor genome name in Panagram's samples.tsv. The bitmap will be used from this anchor's perspective. |
| `-p`       | `--panagram-dir`   | Path to the Panagram index directory.                                                                 |
| `-s`       | `--sample-file`    | Path to the Panagram samples.tsv file.                                                                |
| `-m`       | `--phenotype-file` | Path to the phenotype CSV file. See format below.                                                     |
| `-t`       | `--type`           | Feature type to test. Must be `gene`, `exon`, `reg`, or `all`.                                        |
| `-o`       | `--output-dir`     | Output directory for all results.                                                                     |
| `-n`       | `--threads`        | Number of feature/window chunks to process in parallel.                                               |

### Optional arguments

| Flag              | Description                                                                                                       |
| ----------------- | ----------------------------------------------------------------------------------------------------------------- |
| `-g`, `--gff`     | Path to the GFF annotation file for the anchor genome. Required for `gene`, `exon`, and `reg`; ignored for `all`. |
| `--chunk-size`    | Number of features/windows per parallel chunk. Default: `1000`. For whole-genome scans, `10000` is recommended.   |
| `--sort-temp-dir` | Temporary directory used by Unix `sort` during FDR correction. Useful for large runs on scratch storage.          |

### Feature types

| Type   | Description                                                                                                                                                         |
| ------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `gene` | Tests annotated gene intervals from the GFF file.                                                                                                                   |
| `exon` | Tests annotated exon intervals from the GFF file. Repeated exon entries from the same gene are numbered sequentially.                                               |
| `reg`  | Tests overlapping 200 bp windows in a region spanning 500 kb upstream of each gene to 500 kb downstream of each gene, excluding windows that overlap the gene body. |
| `all`  | Tests overlapping 200 bp windows across the entire anchor genome/contigs. This mode ignores the GFF file.                                                           |

For `reg` and `all`, windows are 200 bp long with a 100 bp step.

Window feature IDs are named as:

```text
chromosome_start
```

Example:

```text
chr1_RagTag_100000
```

### Example command for a gene analysis

```bash
python association_testing.py \
    -a Appa_hap1 \
    -p PANAGRAM_INDEX_DIR \
    -g anchor_annotations.gff \
    -s samples.tsv \
    -m phenotypes.csv \
    -t gene \
    -o association_output \
    -n 40
```

### Example phenotype file

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
OUTPUT_DIR/haplotypes/{chunk_id}_{type}_haplotypes.tsv
```

One haplotype file is produced per feature/window chunk. The file is stored in wide format, with one
row per feature/window and one column per sample. Values are binary haplotype group assignments.

Columns:

```text
feature_id    sample1    sample2    sample3    ...
```

Example:

```text
feature_id    Zeke_hap2    Zeke_hap1    Appa_hap1    Appa_hap2
LOC119876711  1            1            1            1
OR6C4B        0            0            0            0
```

#### Association files

```text
OUTPUT_DIR/associations/{chunk_id}_{type}_associations.tsv
```

One association result file is produced per feature/window chunk.

Columns:

```text
feature_id    phenotype_name    group1_size    group2_size    group1_mean    group2_mean    p_value
```

### Final outputs

Combined volcano plot:

```text
OUTPUT_DIR/plots/{type}_volcano.png
```

Sorted FDR-corrected association table:

```text
OUTPUT_DIR/associations/all_{type}_associations.fdr_sorted.tsv
```

Individual phenotype volcano plots:

```text
OUTPUT_DIR/plots/{phenotype_name}_volcano.png
```

The final FDR-corrected table contains:

```text
feature_id
phenotype_name
group1_size
group2_size
group1_mean
group2_mean
p_value
fdr
norm_mean_difference
plotted
```

Associations are conventionally considered significant at:

```text
fdr <= 0.05
```

## Annotating Window-level Associations

The `annotate_significant_windows.py` script annotates significant whole-genome window associations
with nearby genes. It is intended for use with the FDR-corrected output from the `all` mode of
`association_testing.py`.

### Usage

```bash
python annotate_significant_windows.py \
    -a ASSOCIATION_FILE \
    -g GFF_FILE \
    -o OUTPUT_FILE
```

### Required arguments

| Short flag | Long flag            | Description                                                                                                               |
| ---------- | -------------------- | ------------------------------------------------------------------------------------------------------------------------- |
| `-a`       | `--association-file` | FDR-sorted association TSV from the main pipeline. Usually `OUTPUT_DIR/associations/all_all_associations.fdr_sorted.tsv`. |
| `-g`       | `--gff`              | GFF annotation file for the anchor genome.                                                                                |
| `-o`       | `--output-file`      | Output annotated significant association TSV.                                                                             |

### Optional arguments

| Flag            | Description                                                         |
| --------------- | ------------------------------------------------------------------- |
| `--fdr-cutoff`  | FDR cutoff for significant associations. Default: `0.05`.           |
| `--window-size` | Window size used in `all` mode association testing. Default: `200`. |
| `--chunksize`   | Number of association rows to read at a time. Default: `500000`.    |

### Output

The output file contains one row per significant window-phenotype association.

Columns:

```text
feature_id
phenotype_name
chrom
start
end
p_value
fdr
intersecting_genes
upstream_gene
upstream_distance
downstream_gene
downstream_distance
```

If a window intersects one or more genes, the intersecting genes are reported in
`intersecting_genes`, and the upstream/downstream fields are left blank. If a window does not
intersect a gene, the nearest upstream and downstream genes are reported instead.

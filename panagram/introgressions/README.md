# Introgression Calling with Panagram

Introgressions are regions where one species has inherited a sequence from another related species.
Panagram's bitmap enables calculating the fraction of shared kmers between an anchor and all other
accessions in a pangenome. This allows one to call introgressions by looking for differences between
a reference anchor and another accession, or by looking for similarities between one accession and
another suspected introgression donor.

The introgression caller found here can be used to identify potential introgressed regions across
a pangenome and save them to a BED file for further analysis. The introgression caller works best on
a pangenome of accessions that are the same species with a closely related introgression donor and
reference anchor. Chromosome regions are divided into discrete bins and introgressions are called on
a per-bin basis.

There are additionally scripts to generate simulated pangenomes with various numbers of mutations
and introgressions.

## Usage

Once you have assembled a pangenome with Panagram, you will need an additional config.yaml file and
group.tsv file to control the introgression caller parameters. See below for the format for these
files. The introgression caller can be used as follows:

```
python introgression_runner.py <config.yaml> <--sweep>
```

Add the `--sweep` flag to try a range of kmer similarity thresholds from 0.1-0.95. Note that this
kicks off a number of threads equal to `18 * num. threads chosen in config file`.
[screen](https://www.gnu.org/software/screen/) must be installed on your system for `--sweep` to
work.

## Example

There is an example pangenome with introgressions called, scored, and visualized for a simulated
pangenome at [Coming Soon!]. The example also contains config and group files for reference.

## Group File

The group.tsv is a tab-separated file with 2 columns:
- *name*: the name given to each accession in the samples.tsv file for Panagram
- *group*: the group an accession is a part of. 'REF' denotes which accession is the reference. Other
groups can have any name; typically there are at least 2 other groups for introgression donors and
introgression recipients. The caller can be set to call all introgressions for one group based on
differences with 'REF' or call introgressions based on similarities between one group (typically
donors) and another (typically recipients).

## Config File Parameters

There are 3 parts to the introgression caller: calling, postprocessing, and scoring. Calling will
generate predicted introgression BED files for each accession you choose, for your chromosome(s) of
choice. Postprocessing helps smooth these BED files and remove potential false positives. Scoring
can only be used if comparing called introgressions to the SV-based approach
[here](https://github.com/malonge/CallIntrogressions) or to compare against introgressions generated
by the simulator. It will calculate precision, recall, and other metrics based on the number of
introgressions that match the results from another method.

Parameters are controlled by a YAML file. See `example_config.yaml` for reasonable defaults. The
parameters for each section are as follows:

### General Parameters

| Parameter   | Type   | Description                                                               |
|-------------|--------|---------------------------------------------------------------------------|
| output_dir  | string | path to output folder                                                     |
| index_dir   | string | path to folder where Panagram was run for pangenome                       |
| tsv         | string | path to group.tsv file                                                    |
| bin         | int    | determines size of bin when discretizing chromosome regions               |
| ref         | string | name of reference accession; should be the only accession in the REF group|
| threads     | int    | number of threads to use                                                  |

### Calling Parameters

| Parameter | Type                | Description                                                                 |
|----------|---------------------|-----------------------------------------------------------------------------|
| run      | boolean             | whether to run calling                                              |
| grp      | string              | groups of suspected introgression recipients                                |
| anc      | list[string]/null   | accessions to run introgression caller for (set grp to null to use this list)|
| chr      | list[string]/null   | chromosomes to check for introgressions (null = all)                         |
| cmp      | list[string]        | REF and/or groups of suspected introgression donors to compare against       |
| thr      | float               | bins below threshold are introgressions for REF; above for other groups      |
| stp      | int                 | kmer step size when sampling from bitmap                                     |
| gnm      | float/-1/null       | shift kmer sims to this mean; -1 auto-calc; null disables                     |
| trm      | int/null            | omit values outside this many stdevs from gnm calculations                      |
| sft      | string/null         | choose either 'mean' or 'median' smoothing filter; null disables               |
| ssz      | int/null            | number of bins in sft filter window                                           |
| urf      | boolean             | use reference coordinate system when cmp is REF; won't be used for other groups |
| rmf      | boolean             | remove fixed kmers shared by all accessions                                  |
| rmu      | list[string]/null   | remove kmers not in REF or ogrp for these noisy accessions; forces urf=false for these accessions |
| ogrp     | list[string]/null   | remove kmers not present in these groups (exclude REF)                       |
| edg      | boolean             | highlight edges and dampen center kmer similarities                           |
| isc      | boolean             | increase contrast by boosting kmer similarities                               |
| vis      | boolean             | visualize binned bitmap and detected introgressions                          |

### Postprocessing Parameters

| Parameter | Type         | Description                                                                                                                                                                      |
|----------|--------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| run      | boolean      | whether to run postprocessing                                                                                                                                                    |
| act      | list[string] | choose postprocessing steps: 'lift' (liftover to ref coordinate space), 'fgap' (fill gaps between introgressed bins), 'fcen' (fill larger gaps caused by centromere), 'rmbn' (remove small introgressions) |
| min      | int/null     | if rmbn, remove introgressions smaller than this number of bins                                                                                                                            |
| gap      | int/null     | if fgap, fill gaps between introgressions up to this number of bins                                                                                                                |
| map      | str/null     | if lift and not paf, minimap parameters for aligning to ref; must include -c flag; null to use '-x asm20 -c'                                                                 |
| paf      | str/null     | if lift already run, path to existing alignment folder to skip alignment                                                                                                 |

### Scoring Parameters

| Parameter | Type              | Description                                                                                                                     |
|----------|-------------------|---------------------------------------------------------------------------------------------------------------------------------|
| run      | boolean           | whether to run scoring                                                                                                          |
| gdt      | str               | path to folder containing Jaccard similarity text files or ground truth introgressions in the same format                      |
| act      | list[string]/null | apply postprocessing to ground truth; null to skip                                                                              |
| min      | int/null          | see postprocessing                                                                                                              |
| gap      | int/null          | see postprocessing                                                                                                              |
| thr      | float             | all values above this threshold are considered introgressed in the ground truth                                                 |
| cmp      | list[string]      | groups to check introgressions against; REF and merged results are compared against the set of introgressions in any group in this list |
| vis      | boolean           | visualize true/false positives and negatives per bin                                                                            |

## Introgression Simulator Usage

The simulator takes a reference FASTA and applies mutations to generate a wild relative. It then
takes the reference and replaces sections of it with the aligning section of the wild relative to
simulate a F1 offspring with introgressions. Subsequent generations of offspring are created by
applying mutations to the F1 offspring, up to the mutation rate chosen in the parameters.

```
python simulate_introgressions.py \
--ref <reference fasta to use as a base> \
--out-folder <output folder>

optional arguments:
  -h, --help            show this help message and exit
  --ref REF             Reference FASTA (input)
  --out-folder OUT_FOLDER
                        Output folder
  --offspring OFFSPRING
                        Number of offspring to simulate from created wild
                        relative (default: 1)
  --rounds ROUNDS       Number of rounds to apply mutations to offspring. Will
                        output offspring for each round (default: 5)
  --num-introgressions NUM_INTROGRESSIONS
                        Number of introgressions per generation (default: 48)
  --introgression-size-min INTROGRESSION_SIZE_MIN
                        Min introgression size (bp)
  --introgression-size-max INTROGRESSION_SIZE_MAX
                        Max introgression size (bp)
  --rel-sub-rate REL_SUB_RATE
                        Wild relative per-base substitution rate (default:
                        3.3e-3)
  --rel-ins-rate REL_INS_RATE
                        Wild relative per-base insertion rate (default:
                        3.3e-3)
  --rel-del-rate REL_DEL_RATE
                        Wild relative per-base deletion rate (default: 3.3e-3)
  --rel-ins-size-min REL_INS_SIZE_MIN
                        Wild relative insertion size min (bp)
  --rel-ins-size-max REL_INS_SIZE_MAX
                        Wild relative insertion size max (bp)
  --rel-del-size-min REL_DEL_SIZE_MIN
                        Wild relative deletion size min (bp)
  --rel-del-size-max REL_DEL_SIZE_MAX
                        Wild relative deletion size max (bp)
  --mut-sub-rate MUT_SUB_RATE
                        Offspring per-base substitution rate
  --mut-ins-rate MUT_INS_RATE
                        Offspring per-base insertion rate
  --mut-del-rate MUT_DEL_RATE
                        Offspring per-base deletion rate
  --mut-ins-size-min MUT_INS_SIZE_MIN
                        Offspring insertion size min (bp)
  --mut-ins-size-max MUT_INS_SIZE_MAX
                        Offspring insertion size max (bp)
  --mut-del-size-min MUT_DEL_SIZE_MIN
                        Offspring deletion size min (bp)
  --mut-del-size-max MUT_DEL_SIZE_MAX
                        Offspring deletion size max (bp)
  --seed SEED           Random seed for reproducibility (default: 42)
```


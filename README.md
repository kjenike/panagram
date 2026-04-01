<div align="center">
  <img src="./panagram/assets/panagram.png" alt="Panagram" width="300"/>
</div>

# Panagram: Interactive, alignment-free pan-genome browser

#### Katie Jenike, Nicole Brown, Sam Kovaka, Shujun Ou, Stephen Hwang, Srividya Ramakrishnan, Ben Langmead, Zach Lippman, Ian R Henderson, Michael C Schatz


Welcome to Panagram! Panagram is [an alignment-free pan-genome viewer](https://www.dropbox.com/s/g7snjgr8bs6c2uj/2023.01.17.Panagram.pdf).

# Installation

```
git clone --recursive https://github.com/kjenike/panagram.git
cd panagram
pip install .
```

The `--recursive` option is required to install the KMC dependency. If you forget to include it,
you can update the repository with the command `git submodule update --init`.

Installation may fail if pip is not up-to-date or if setuptools is not up-to-date. In order to
update pip and setuptools run:

```
pip install --upgrade pip
pip install --upgrade setuptools
```

## Dependencies

Requires python version >=3.11, pip, samtools, bgzip, and tabix. All other dependencies should be
automatically installed via pip.

Panagram relies on [KMC](https://github.com/refresh-bio/KMC) to build its kmer index.
This should be installed automatically, however it is possible that the KMC installation will fail
but panagram will successfully install. In this case `panagram view` can be run, but
`panagram index` will return an error. You may be able to debug the KMC installation by running
`make -C KMC py_kmc_api` and attempting to fix any errors, then re-run `pip install -v .` after the
 errors are fixed.

# Running
Panagram runs in two steps, the anchoring step (index command) and viewing (view command).

## Anchoring

To prepare Panagram for anchoring, run:
```
panagram index <samples.tsv> -k <k> --prepare
```
To choose which genomes will act as anchors, use `--anchor_genomes <one or more space-separated names>`.
A common choice is to use one or more reference genomes as anchors.

To run the indexing step,
start by preparing the panagram index. It is best to create an empty folder that will act as Panagram's
index folder. Within this folder, create a subfolder called FASTAS; this is where you can place any
FASTAS that you want to include in your pangenome. You can also create a folder called GFFs; if you
have any annotation files in GFF3 format, you can place them in here. Next, you will need to tell
Panagram where your FASTAS and GFFs are. For this, you will need a tsv file with a list of the samples.

The samples.tsv file should contain one sample per line. You need the name, fasta, and gff columns.
On each line include the sample name and minimally the fasta file location. Currently genome names
should only contain alphanumeric characters and underscores due to KMC requirements. If you have no
annotations for a file, leave the gff column blank. If you have multiple annotation files per sample,
you can concatenate them into one gff file. The id and anchor columns will be created by Panagram. See
below for an example.

```
name	fasta	gff	id	anchor
sample1	FASTAS/sample1.fasta	GFFS/sample1.gff	0	True
sample2	FASTAS/sample2.fasta	GFFS/sample2.gff	1	True
sample3	FASTAS/sample3.fasta	GFFS/sample3.gff	2	True
sample4	FASTAS/sample4.fasta	GFFS/sample4.gff	3	True
sample5	FASTAS/sample5.fasta	GFFS/sample5.gff	4	True
sample6	FASTAS/sample6.fasta	GFFS/sample6.gff	5	True
```

As another example, this samples.tsv would be a comparison with just two genomes.

```
name	fasta	gff	id	anchor
Col_0	FASTAS/wlod_Col-0.ragtag_scaffolds.fa	GFFS/wlod_Col-0.ragtag_scaffolds.gff	0	True
Tanz_1	FASTAS/wlod_Tanz-1.patch.scaffold.Chr.fa	GFFS/wlod_Tanz-1.patch.scaffold.Chr.gff	1	True
```

It is super important that any gff files are in the correct format.
We strongly suggest that if you run into any problems, you first check the format annotation format.
This can be done with command line tools like gff3validator or online
[here](https://genometools.org/cgi-bin/gff3validator.cgi).


Picking an acceptable k-mer length for the data set can be tricky. For samples that are very similar,
a larger k may be more approperiate. While samples that are more diverged may benefit from a smaller
k-mer length. The papers [here](https://www.cell.com/iscience/fulltext/S2589-0042(24)00275-X?uuid=uuid%3A8d061319-27f8-49ca-b7ee-0d33ec846225)
and [here](https://pubmed.ncbi.nlm.nih.gov/39890468/) give some detail on picking "good" k-mer length,
but if in doubt, k=21 usually works fine.

Once the preparation step is run, you can run Panagram's anchoring via snakemake and specify the
number of threads you want to use with this command:

```
snakemake --cores <num. threads> all
```

This step anchors KMC bitvectors to FASTA files to create a pan-kmer bitmap.

## View

Once anchoring is complete, navigate to the index folder and view your pangenome with
`panagram view .` This runs a local Dash server. The pan-genome browser can be viewed at
 http://127.0.0.1:8050/ by default.

Here is the full set of flags you can choose:

```
usage: panagram view [-h] <index_dir/> [genome] [chrom] [start] [end]
  index_dir           Panagram index directory
  genome              Initial anchor genome (optional)
  chrom               Initial chromosome (optional)
  start               Initial start coordinate (optional)
  end                 Initial end coordinate (optional)
  --ndebug            Run server in production mode (important for a public-
                      facing server)
  --port str          Server port (default: 8050)
  --host str          Server address (default: 127.0.0.1)
  --url_base str      A local URL prefix to use app-wide (passed to
                      Dash.dash(url_base_pathname=...)) (default: /)
  --bookmarks str     Bed file with bookmarked regions (default: None)

```

# Introgression Calling

<div align="center">
  <img src="./panagram/introgressions/assets/introgressions.svg" alt="Panagram" width="200"/>
</div>

Panagram's bitmap also enables calling introgressions between members of your pan-genome. For all
information on the introgression calling module, see the README [here](./panagram/introgressions).


# Bitdump

If you want to see the bitmap generated by Panagram, you can use the following:

```
usage: panagram bitdump [-h] [-v bool] index_dir coords step
  Query pan-kmer bitmap generated by "panagram index"/

  index_dir             Panagram index directory
  coords                Coordinates to query in chr:start-end format
  step                  Spacing between output kmers (optimized for multiples
                        of 100) (default: 1)
  -v bool, --verbose bool
                        Output the full bitmap (default: False)
```

# Example Run

First download the example_data.zip bacterial data from:
[http://data.schatz-lab.org/panagram/](http://data.schatz-lab.org/panagram/) or use this
[direct link](https://bx.bio.jhu.edu/data/panagram/example_data.zip).

Unzip the archive and you will find 5 bacterial genomes plus their annotations:
```
unzip example_data.zip
```

To run, first index the genomes:

```
cd example_data
panagram index samples.tsv -k 21 --prepare
snakemake --verbose --cores 30 all
```

Then you can panagram to visualize (from the example_data directory):
```
panagram view .
```

From there, you can view the results in your browser at [http://127.0.0.1:8050/](http://127.0.0.1:8050).

# Hosting and Proxies

Panagram uses [Dash](https://dash.plotly.com/introduction) to serve the plotly visualizations.
By default the dedicated webserver runs on localhost (127.0.0.1) on port 8050, but you can reverse
proxy to a different port and path using a web engine such as [nginx](https://www.nginx.com/).

For nginx, first reconfigure your nginx configuration file to add (note to be very careful
with the use of the slash ('/') character):

```
    location /panagram {
      proxy_pass http://127.0.0.1:8050;
    }
```

The retart nginx with:

```
systemctl stop nginx
systemctl start nginx
```

For a secure public-facing server, be sure to run with the option `panagram view --ndebug` to
disable debug mode. You may also wish to change the base URL path with the `--url_base` option,
for example to something like `--url_base /panagram/`. The port and host name can be specified
by the `--port` and `--host` options.

Finally you will need to run panagram using `panagram view <dir>`. You will probably want to run
this in a loop in case it needs to be restarted, such as:

```
until panagram view --ndebug .; do echo "restarting"; sleep 1; done
```

# Contributing and Support

Have a question or found a bug? [Open an issue](https://github.com/kjenike/panagram/issues)!


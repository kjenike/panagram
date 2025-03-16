# Installation

Requires OpenMP, htslib (shown below), and compiled KMC

```
#build KMC
make -C ../KMC

#download and build htslib
git clone https://github.com/samtools/htslib.git
cd htslib
autoreconf -i  # Build the configure script and install files it uses
./configure    # Optional but recommended, for choosing extra functionality
make
make install

#compile anchor binary
cd ..
make run_anchor
```

# Usage
Copy the `Snakefile` and `run_anchor` binary from this directory to build the indexing using the C++ implementation. Does not support annotations.
```
panagram index --prepare samples.tsv
cp ~/panagram/cpp/Snakefile ~/panagram/cpp/run_anchor . #copy Snakefile from panagram/cpp
vim Snakefile ... #edit to specify anchor binary path
snakemake -c <cores> .
```

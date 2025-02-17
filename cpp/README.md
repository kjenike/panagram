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
make anchor
```

# Usage
```
panagram index --prepare samples.tsv
cp ~/panagram/cpp/Snakefile . #copy Snakefile from panagram/cpp
snakemake -c <cores> .
```

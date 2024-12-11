import sys
import os
import os.path
from os import path
import subprocess
import numpy as np
import pandas as pd
import bgzip
import gzip
import csv
import glob
import pysam
from collections import defaultdict, Counter
from time import time
from Bio import bgzf, SeqIO
import yaml
import multiprocessing as mp
from types import SimpleNamespace
import shutil
import snakemake
import re
import logging
from panagram.index import Index


# let's try some toy examples of manipulating a bitmap

# import numpy as np
# arrs = dict()
# arr = np.array([3,4,5])
# step = 0

# # for step,arr in byte_arrs.items():
# arrs[step] = np.concatenate(arr, axis=1)
# arrs[step] = np.concatenate(np.array([4,5,6,7]), axis=1)


def better_dir(item):
    # don't print hidden functions
    methods = dir(item)
    return [method for method in methods if not method.startswith("_")]


# Reading in an index - you can do this after running snakemake
index_dir = "/home/nbrown62/data_mschatz1/nbrown62/panagram/example_data"
index = Index(index_dir)
print(index)

# Gives you something like Index(input='samples.tsv', mode=None, prefix='', k=21, cores=1, lowres_step=100, max_bin_kbp=200, min_bin_count=100, max_view_chrs=50, gff_gene_types=['gene'], gff_anno_types=None, gff_name='Name', anchor_genomes=['ecoli', 'ecoli_k12', 'klebsiella', 'salmonella', 'shigella'], prepare=False, kmc=KMC(memory=8, threads=1, use_existing=False), use_existing=1, threads=1, memory=1)
print(index.anchor_genomes)

# All available methods for Index class
print(better_dir(index))

# Get a genome class - they are listed by their names in a dictionary
print(index.genome_names[0])
ecoli_genome = index.genomes[index.genome_names[0]]

# All available methods for Genome class
print(better_dir(ecoli_genome))

# pandas df of chromosome names, sizes, and gene count
# might tell you which chrs are most interesting
chromosome_info = ecoli_genome.chrs
print(chromosome_info)

# bitmaps themselves are compressed objects until accessed
print(ecoli_genome.bitmaps)

# Get part of the bitmap from the genome
# Query a chromosome at a certain position to get a piece of the bitmap
# TODO: Not sure what step does yet? Sum of kmer occurrence across a window of size step?
ecoli_bitmap = ecoli_genome.query("NZ_CP015023.1", 0, 5506781, step=1)


# You can also query using the index itself
# print(index.query_bitmap(ecoli_genome, chrom, start=None, end=None, step=1))

# what if you wanted to look at intergressions
# set a threshold for what counts as an intergression
threshold = 0.75  # 75% co-occurrence across the pangenome
ecoli_bitmap["frac_co_occurrence"] = ecoli_bitmap.sum(axis=1)
ecoli_bitmap["frac_co_occurrence"] = ecoli_bitmap["frac_co_occurrence"] / len(
    index.genomes
)
print(ecoli_bitmap)

# the index on the left (+ starting point index) gives you the indices of intergressions
print(ecoli_bitmap[ecoli_bitmap["frac_co_occurrence"] >= threshold])

# in this example, there are 196262 indices (i.e., there is one at 34329 and one at 5465887)
# TODO: merge nearby intergression locations into one
# A run of length > run_threshold is a true intergression
# TODO: how long are intergressions typically?
# use a sliding window and take average co-occurance perhaps within a k-sized window

diff_threshold = 10000  # consider different intergression if more than n away?
# index_diff = df['index'].diff().fillna(0)

import sys
import os

# import os.path
# from os import path
from pathlib import Path
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
import plotly.express as px


def better_dir(item):
    # don't print hidden functions
    methods = dir(item)
    return [method for method in methods if not method.startswith("_")]


def visualize(pair, output_file, inverse=False):
    # take a look at what pair looks like after manipulation
    # pair[pair >= 1] = 10
    if inverse:
        fig = px.imshow(
            pair,
            color_continuous_scale=px.colors.sequential.Plasma[::-1],
            x=pair.columns,
            y=pair.index,
        )
    else:
        fig = px.imshow(pair, x=pair.columns, y=pair.index)
    fig.write_image(output_file)
    return


index_dir = "/home/nbrown62/data_mschatz1/nbrown62/panagram_data/tomato"
anchor = "SL5"  # "SL5"
# chr_name = "BGV006775_MAS2.0ch11"
output_dir = "/home/nbrown62/data_mschatz1/nbrown62/panagram_data/tomato/introgression_analysis_v1/"


output_dir = Path(output_dir)
output_dir.mkdir(parents=True, exist_ok=True)
index = Index(index_dir)
genome = index.genomes[anchor]
chrs = genome.sizes.keys()
# # print(index.genomes)
bitmap_step = 100
max_chr_bins = 350
k = 31

for chr_name in chrs:
    # get an entire chr's bitmap
    chr_size = genome.sizes[chr_name]
    chr_bitmap = genome.query(chr_name, 0, chr_size, step=bitmap_step)

    # get correlation matrix
    start_coord = 0
    end_coord = chr_size
    bin_size = ((end_coord - start_coord) // max_chr_bins) + 1
    num_kmers_in_bin = bin_size - k + 1

    pan, pair = index.bitmap_to_bins(chr_bitmap, bin_size)

    visualize(pair, output_dir / f"{anchor}_{chr_name}_original_heatmap.png", inverse=True)

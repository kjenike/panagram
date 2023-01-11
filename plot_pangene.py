#from bokeh.plotting import figure, output_file, save
import sys
import csv
import pysam
#import py_kmc_api as kmc
import bgzip
from Bio import bgzf
from plotly.subplots import make_subplots
import plotly.figure_factory as ff
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from upsetplot import generate_counts, plot
from lenspy import DynamicPlot
import dash
from mycolorpy import colorlist as mcp
import time
from Bio import Phylo
from scipy import signal
from dash import Dash, dcc, html, Input, Output, ctx, State, no_update
import math
from io import StringIO
from scipy.cluster.hierarchy import linkage
from scipy.cluster import hierarchy
import dash_bootstrap_components as dbc
from dash_bootstrap_templates import ThemeChangerAIO, template_from_url
from scipy.spatial.distance import pdist, squareform
import plotly.figure_factory as ff
#from cStringIO import StringIO
#plt.rcParams["figure.figsize"] = (40,40)
config_f = sys.argv[1] #"test_data/Maize_testing/maize_config.txt"

SG_window =53
poly_order = 3
SG_polynomial_order = 3
SG_check = [1]

#/Users/katiejenike/Desktop/Lab_notes/PanSol/PANGEN/DEBUG/QC_PLOTS/SQUI2_PANGEN_PLOTS
#python plot_chr_interactive.py result_SQUI2.chr1.txt SQUI2
#python plot_pangene.py test_data/result_SQUI2.chr1.txt SQUI2
#fasta_files = "test_data/Maize_testing/PRE_FILES/fastas.tsv"

#ref_genome = "Solqui2"
n_skips_start = 100
tree_file = ""#"test_data/all.alns.aug5.2022.fa.treefile"

fai_f = ""

#chrs_list = ["chr1","chr2","chr3","chr4","chr5",
#    "chr6","chr7","chr8","chr9","chr10",]
#chr_lens = [308452471, 243675191, 238017767, 250330460, 226353449, 
#181357234, 185808916, 182411202, 163004744, 152435371]


#Read in the config file:
with open(config_f, "r") as f:
    line1 = f.readline()
    while line1:
        num_samples = int(line1.strip().split(':')[1])
        line2 = f.readline()
        kmer_len = int(line2.strip().split(':')[1])

        line3 = f.readline()
        num_chrs = int(line3.strip().split(':')[1])

        line4 = f.readline()
        buff = int(line4.strip().split(':')[1])

        line5 = f.readline()
        x_start_init = int(line5.strip().split(':')[1])

        line6 = f.readline()
        x_stop_init = int(line6.strip().split(':')[1])

        line7 = f.readline()
        bins = int(line7.strip().split(':')[1])

        line8 = f.readline()
        gene_content_by_chr_f = line8.strip().split(':')[1]

        line9 = f.readline()
        gene_content_by_gene_f = line9.strip().split(':')[1]

        line10 = f.readline()
        fasta_files = line10.strip().split(':')[1]

        line11 = f.readline()
        genome_comp_file = line11.strip().split(':')[1]

        line12 = f.readline()
        rep_file = line12.strip().split(':')[1]

        line13 = f.readline()
        bins_file = line13.strip().split(':')[1]

        line14 = f.readline()
        rep_types_file = line14.strip().split(':')[1]

        line15 = f.readline()
        fai_f = line15.strip().split(':')[1]

        line16 = f.readline()
        genome_comp_pre = line16.strip().split(':')[1]

        line17 = f.readline()
        genome_comp_post = line17.strip().split(':')[1]

        line18 = f.readline()
        mash_filenames = line18.strip().split(':')[1]

        line19 = f.readline()
        mash_edges = line19.strip().split(':')[1]

        line20 = f.readline()
        genome_size_f = line20.strip().split(':')[1]

        line21 = f.readline()
        bit_file_prefix = line21.strip().split(':')[1]

        line22 = f.readline()
        anchor_name = line22.strip().split(':')[1]

        line1 = f.readline()

colors = mcp.gen_color(cmap="viridis_r",n=num_samples)
sns.set_palette('viridis', num_samples)
gene_file = gene_content_by_gene_f

#index_fasta = 
#genome_comp_file = "test_data/Maize_testing/PRE_FILES/totals.txt"
labels = []
with open(fasta_files, "r") as f:
    line = f.readline()
    while line:
        tmp = line.split('\t')[0]
        if tmp == "Saet_bc2058":
            tmp = "M82"
        elif tmp == "Saet_bc2059":
            tmp = "Sins2"
        labels.append(tmp)
        line = f.readline()
#anchor_name = labels[-1]


chrs_list = []
chr_lens = []
cntr = 0
#chr_nums = 12
with open(fai_f, "r") as f:
    line = f.readline()
    while line and cntr < num_chrs:
        #print(line)
        chrs_list.append(line.strip().split('\t')[0])
        chr_lens.append(int(line.strip().split('\t')[1]))
        cntr += 1
        line = f.readline()

BGZ_SUFFIX = ".pank"
ANN_SUFFIX = ".pana"

class KmerBitmap:
    def __init__(self, prefix, kmc=None, genome_tsv=None, anchors=None):
        self.prefix = prefix
        self.bgz_fname = prefix + BGZ_SUFFIX
        self.ann_fname = prefix + ANN_SUFFIX

        self.offsets = dict()
        self.ref_len = 0

        if kmc is None and genome_tsv is None:
            self._init_read()
        else:
            self._init_write(kmc, genome_tsv, anchors)

    def _init_read(self):
        self.write_mode = False

        with open(self.ann_fname) as self.ann:
            self.ngenomes, self.nbytes = map(int, self.ann.readline().split())
            #self._set_ngenomes(ngenomes)
            #self.ngenomes = len(self.fastas)

            self.genomes = list()

            nt_offs = 0
            for line in self.ann:
                nt_offs, genome, seq_name = line.split()[:3]
                self.offsets[seq_name] = int(nt_offs)
                if len(self.genomes) == 0 or genome != self.genomes[-1]:
                    self.genomes.append(genome)
            self.ref_len = nt_offs

        #t = time()
        #idx_in = open(self.prefix + ".panx", "rb")
        #nblocks = np.fromfile(idx_in, "uint64", 1)[0]
        #dtype = [("rstart", "uint64"), ("dstart", "uint64")]
        #blocks = np.zeros(int(nblocks)+1, dtype=dtype)
        #blocks[1:] = np.fromfile(idx_in, dtype, nblocks)
        #self.bgz_blocks = blocks.astype([("rstart", int), ("dstart", int)])


        #self.bgz = bgzf.BgzfReader(self.bgz_fname, "rb")
        self.bgz = open(self.bgz_fname, "rb")

    def query(self, name, start, end, step=1):
        byte_start = self.nbytes * (self.offsets[name] + start)
        length  = end - start

        #blk = np.searchsorted(self.bgz_blocks["dstart"], byte_start, side="right")-1
        #blk_offs = byte_start - self.bgz_blocks["dstart"][blk]
        #blk_start = self.bgz_blocks["rstart"][blk]


        #self.bgz.seek(bgzf.make_virtual_offset(blk_start, blk_offs))
        self.bgz.seek(byte_start)
        buf = self.bgz.read(length * self.nbytes)

        pac = np.frombuffer(buf, "uint8").reshape((length, self.nbytes))
        if step > 1:
            pac = pac[::step]
        ret = np.unpackbits(pac, bitorder="little", axis=1)[:,:self.ngenomes]

        #pac = np.frombuffer(buf, "uint32")
        #ret = np.zeros((len(pac), self.ngenomes), dtype="uint8")
        #for i in range(self.ngenomes):
        #    ret[:,i] = (pac >> i) & 1

        return ret
    
    @property
    def genome_names(self):
        return list(self.genomes)
    def _init_write(self, kmc_file, genome_tsv, anchors):
        self.write_mode = True
        #self.bgz = bgzip.BGZipWriter(open(self.bgz_fname, "wb"))
        self.bgz = open(self.bgz_fname, "wb")
        self.ann = open(self.ann_fname, "w")

        self.kmc_db = kmc.KMCFile()
        self.kmc_db.OpenForRA(kmc_file)

        self.genomes = list()
        self.fastas = dict()
        with open(genome_tsv) as genome_file:
            genomes_in = csv.reader(genome_file, delimiter="\t")
            for name, fasta in genomes_in:
                self.genomes.append(name)
                self.fastas[name] = fasta

        self.ngenomes = len(self.fastas)
        self.nbytes = int(np.ceil(self.ngenomes / 8))

        self.ann.write(f"{self.ngenomes}\t{self.nbytes}\n")

        for name, fasta in self.fastas.items():
            if anchors is None or name in anchors:
                self._load_fasta(name, fasta)
    def _load_fasta(self, name, fname):

        with pysam.FastaFile(fname) as fasta:
            vec = kmc.CountVec()
            t = time()
            for seq_name in fasta.references:
                seq = fasta.fetch(seq_name)#, 0, 10000000) 

                self.kmc_db.GetCountersForRead(seq, vec)

                arr = np.array(vec, dtype="uint32")
                arr = arr.view("uint8").reshape((len(arr),4))
                if self.nbytes < 4:
                    arr = arr[:,:self.nbytes]

                self.bgz.write(arr.tobytes())

                self.ann.write(f"{self.ref_len}\t{name}\t{seq_name}\n")
                self.ref_len += len(arr)

                print(f"Anchored {seq_name}", time()-t)
                t = time()

    def close(self):
        self.ann.close()
        self.bgz.close()


class KmerRef:
    def __init__(self, prefix):
        self.pac_fname = f"{prefix}.pac"
        ann_in = open(f"{prefix}.ann")

        self.ngenomes, self.nbits = map(int, ann_in.readline().split())
        self.nbytes = self.nbits // 8

        self.offsets = dict()
        for line in ann_in:
            nt_offs, name = line.split()[:2]
            self.offsets[name] = int(nt_offs)
        
        ann_in.close()
    def get_counts(self, name, start, end, step=1):
        byte_start = self.nbytes * (self.offsets[name] + start)
        length  = end - start
        pac = np.fromfile(self.pac_fname, offset=byte_start, count=length, dtype=f"int{self.nbits}")

        if step != 1:
            pac = pac[np.arange(0, len(pac), step)]

        ret = np.zeros((len(pac), self.ngenomes), dtype=bool)
        for i in range(self.ngenomes):
            ret[:,i] = (pac >> i) & 1

        return ret


def bitvec_to_mat(bitvecs, genome_count):
    #bitvec_to_mat(np.array([45]), 9).sum(axis=1)
    ret = np.zeros((len(bitvecs), genome_count), dtype=bool)
    for i in range(genome_count):
        ret[:,i] = (bitvecs >> i) & 1
    return ret


def get_newick(node, parent_dist, leaf_names, newick='') -> str:
    """
    Convert sciply.cluster.hierarchy.to_tree()-output to Newick format.

    :param node: output of sciply.cluster.hierarchy.to_tree()
    :param parent_dist: output of sciply.cluster.hierarchy.to_tree().dist
    :param leaf_names: list of leaf names
    :param newick: leave empty, this variable is used in recursion.
    :returns: tree in Newick format
    """
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parent_dist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parent_dist - node.dist, newick)
        else:
            newick = ");"
        newick = get_newick(node.get_left(), node.dist, leaf_names, newick=newick)
        newick = get_newick(node.get_right(), node.dist, leaf_names, newick=",%s" % (newick))
        newick = "(%s" % (newick)
        return newick

def find_shared(some_num):
    total = 0
    for i in some_num:
        total += name_idx_short[i]
    return total

def get_init_rep_types(rep_types_f=rep_types_file):
    rep_types = {}
    for i in range(0, num_chrs):
        rep_types["chr" + str(i+1)] = {}
    with open(rep_types_f, "r") as f:
        line = f.readline()
        while line:
            for i in range(0, num_chrs):
                rep_types["chr" + str(i+1)][line.strip()] = []
            line = f.readline()
    return rep_types
'''    
rep_types = {
        "18S_rRNA_gene":[],
        "25S_rRNA_gene":[],
        "5_8S_rRNA_gene":[],
        "5S_rRNA_gene":[],
        "CACTA_TIR_transposon":[],
        "Copia_LTR_retrotransposon":[],
        "Gypsy_LTR_retrotransposon":[],
        "hAT_TIR_transposon":[],
        "helitron":[],
        "long_terminal_repeat":[],
        "low_complexity":[],
        "LTR_retrotransposon":[],
        "Mutator_TIR_transposon":[],
        "PIF_Harbinger_TIR_transposon":[],
        "rDNA_intergenic_spacer_element":[],
        "repeat_region":[],
        "target_site_duplication":[],
        "Tc1_Mariner_TIR_transposon":[]
        }
'''


tree_fig = {}

def parse_counts_simple(cnts_tmp):
    names_simp = []
    for this in cnts_tmp:
        names_simp.append(name_idx_long[this])
    return names_simp

##### PLOTS

def get_index(j, bin_size, x_start):
    indx = int((j-x_start)/bin_size)
    return indx

def get_index_0_based(j, bin_size):
    indx = int(j/bin_size)
    return indx

def read_count_file(count_file, y_indx ):
    #this is reading the traditional counts. The counts that tell us how repetative stuff is 
    with open(count_file, "r") as f:
        line = f.readline()
        while line:
            zs_tmp = line.strip().split(":")[1].split(',')[:-1]
            line = f.readline()
    #full_counts = [int(i) for i in zs_tmp]
    #This will be an estimate of how reptative things are
    zs_tmp_2 = [len(i) for i in zs_tmp] #list(map(lambda x: len(x), zs_tmp))
    
    return zs_tmp_2 #, zs_tmp 
'''
def read_mini_count_files(x_start, x_stop):
    mini_counts = []
    #need to read in one or two files
    #find the files to read (this will be based on the x_start and x_stops)
    #Re-do the file division, so it isn't based on characters but is based on commas 
    #Then parse the counts to get the simple counts 
    tmpstop = 0
    while tmpstop < x_start:
        tmpstop += mini_file_size
    tmpstart = tmpstop - mini_file_size
    file_1 = "test_data/CNTR_FILES/chr1." + str(tmpstart) + "." + str(tmpstop) + ".txt"
    with open(file_1, "r") as f:
        line = f.readline()
        while line:
            mini_counts = line.strip().split(",")
            line = f.readline()
    if x_stop > tmpstop:
        #Need another file 
        tmpstart2 = tmpstop
        tmpstop = tmpstop + mini_file_size
        file_2 = "test_data/CNTR_FILES/chr1." + str(tmpstart2) + "." + str(tmpstop) + ".txt"
        with open(file_2, "r") as f:
            line = f.readline()
            while line:
                mini_counts += line.strip().split(",")
                line = f.readline()
    tmpstart3 = x_start - tmpstart 
    tmpstop3 = tmpstart3 + (x_stop-x_start)
    exact_mini_counts = [ int(x) for x in mini_counts[tmpstart3:tmpstop3] ] #mini_counts[tmpstart3:tmpstop3]
    simple_exact_mini_counts = bitvec_to_mat(np.array(exact_mini_counts), num_samples).sum(axis=1) #parse_counts_simple(exact_mini_counts)

    return simple_exact_mini_counts, exact_mini_counts
'''
def parse_counts(zs_tmp, zs, bin_size_x, y_indx, x_start, x_stop):
    bin_size_y = 1
    #bin_size_x = #len(zs_tmp)/bins
    cntr = 0
    x = 0
    for i in zs_tmp[x_start:x_stop] :#range(0, (x_stop-x_start)):#range(x_start, x_stop):
        #x = get_index(i, bin_size_x, x_start)
        #x = get_index_0_based(i, bin_size_x)
        #x = get_index_0_based(cntr, bin_size_x)
        if (cntr % bin_size_x) == 0 :
            x += 1
        zs[(i-1)][x] += 1
        #zs[(y-1)][x] += 1
        cntr += 1

    return zs

def read_annotations(ann_file, ann_types):
    with open(ann_file, "r") as f:
        line = f.readline()
        while line:
            tmp = line.split('\t')
            this_chr = tmp[0]
            this_type = tmp[3].strip()

            if this_type in rep_types["chr1"].keys() and this_chr in chrs_list: #and tmp[0]==chrs:
                
                #print(line.strip())
                ann_types[this_chr][this_type].append(int(tmp[1]))
                ann_types[this_chr][this_type].append(int(tmp[2]))
                ann_types[this_chr][this_type].append('None')
            line = f.readline()
    return ann_types

def perc_universal(start,stop, countme, n_skips, this_chr):
    univ = ((anchor.query(this_chr, start, stop, 1).sum(axis=1))==countme).sum()
    #float((names_simp[int(start/n_skips):int(stop/n_skips)]==countme).sum())
    #univ = float(names_simp[int(start/n_skips):int(stop/n_skips)].count(countme))
    return float(univ/(stop-start))*100
    #return float(univ/(((stop)-(start/n_skips))))*100

def read_gene_annotations(ann_file, gene_names, gene_locals, exon_locals, exon_names, cds_locals, n_skips, chrs):
    #rep_file = "Solqui2.repeats.chr1.100k.txt"
    gene_content = {}
    gene_anns = {}
    for i in range(1, num_chrs+1):
        gene_content["chr"+str(i)] = {"Names":[],"Universal":[], "Unique":[]}
        #gene_anns[chr_name[i-1]] = {}
    #gene_content = {"Names":[],"Universal":[], "Unique":[]}
    with open(ann_file, "r") as f:
        line = f.readline()
        while line:
            tmp = line.strip().split('\t')
            this_chr = tmp[0]
            if this_chr in gene_locals.keys():
                #this_name = tmp[8]#line.split('\t')[1].strip()
                this_type = tmp[3]
                this_start = int(tmp[1])
                this_stop = int(tmp[2])
                if this_type == "gene": #this_name.split(":")[0] == "ID=gene":
                    gene_locals[this_chr].append(this_start)
                    gene_locals[this_chr].append(this_stop)
                    gene_locals[this_chr].append('None')
                    this_name = tmp[6].split(';')[0].split("=")[1]
                    #tmp_name = this_name.split(';')[0].split("=")[1]
                    gene_names[this_chr].append(this_name)
                    gene_names[this_chr].append(this_name)
                    gene_names[this_chr].append('None')
                    gene_content[this_chr]["Names"].append(this_name)#perc_universal(int(tmp[2]), int(tmp[3]))
                    gene_content[this_chr]["Universal"].append(float(tmp[5]))
                    gene_content[this_chr]["Unique"].append(float(tmp[4]))
                    gene_anns[this_name] = tmp[6]
                elif this_type == "exon": #this_name.split(":")[0] == "ID=exon":
                    exon_locals[this_chr].append(this_start)
                    exon_locals[this_chr].append(this_stop)
                    exon_locals[this_chr].append('None')
                    tmp_name = tmp[4].split(';')[0].split("=")[1]#this_name.split(';')[0].split("=")[1]
                    exon_names[this_chr].append(tmp_name)
                    exon_names[this_chr].append(tmp_name)
                    exon_names[this_chr].append('None')
                #elif this_name.split(":")[0] == "ID=CDS":
                #    cds_locals.append(int(tmp[2]))
                #    cds_locals.append(int(tmp[3]))
                #    cds_locals.append('None')
            #print()
            line = f.readline()
    #print(gene_locals["chr1"][-100:])
    #print(len(gene_locals["chr1"]))
    #print("*********")
    return gene_locals, gene_names, exon_locals, exon_names, cds_locals, gene_content, gene_anns

def get_x_coordinates(tree):
    #Adapted from:
    #https://github.com/plotly/dash-phylogeny
    """Associates to  each clade an x-coord.
       returns dict {clade: x-coord}
    """
    xcoords = tree.depths()
    # tree.depth() maps tree clades to depths (by branch length).
    # returns a dict {clade: depth} where clade runs over all Clade instances of the tree, and depth
    # is the distance from root to clade
    #  If there are no branch lengths, assign unit branch lengths
    if not max(xcoords.values()):
        xcoords = tree.depths(unit_branch_lengths=True)
    return xcoords
def get_y_coordinates(tree, dist=1.3):
    #Adapted from:
    #https://github.com/plotly/dash-phylogeny
    """
       returns  dict {clade: y-coord}
       The y-coordinates are  (float) multiple of integers (i*dist below)
       dist depends on the number of tree leafs
    """
    maxheight = tree.count_terminals()  # Counts the number of tree leafs.
    # Rows are defined by the tips/leafs
    ycoords = dict((leaf, maxheight - i * dist) for i, leaf in enumerate(reversed(tree.get_terminals())))
    def calc_row(clade):
        for subclade in clade:
            if subclade not in ycoords:
                calc_row(subclade)
        ycoords[clade] = (ycoords[clade.clades[0]] +
                          ycoords[clade.clades[-1]]) / 2
    if tree.root.clades:
        calc_row(tree.root)
    return ycoords
def get_clade_lines(orientation='horizontal', y_curr=0, x_start=0, x_curr=0, y_bot=0, y_top=0,
                    line_color='rgb(25,25,25)', line_width=0.5):
    #Adapted from:
    #https://github.com/plotly/dash-phylogeny
    """define a shape of type 'line', for branch
    """
    branch_line = dict(type='line',
                       layer='below',
                       line=dict(color=line_color,
                                 width=line_width)
                       )
    if orientation == 'horizontal':
        branch_line.update(x0=x_start,
                           y0=y_curr,
                           x1=x_curr,
                           y1=y_curr)
    elif orientation == 'vertical':
        branch_line.update(x0=x_curr,
                           y0=y_bot,
                           x1=x_curr,
                           y1=y_top)
    else:
        raise ValueError("Line type can be 'horizontal' or 'vertical'")
    return branch_line
def biggest_num_in_clade(clade, kmer_num):
    if clade.clades:
        m = 0
        for c in clade.clades:
            tmp_m = biggest_num_in_clade(c, kmer_num)
            if tmp_m > m:
                m = tmp_m
        return m
    else:
        tmp_name = str(clade)
        #print(tmp_name)
        #if tmp_name.count("_") > 0 and tmp_name.count('Saet') == 0:
        #    tmp_name = tmp_name.split("_")[1]
        clade_num = kmer_num[tmp_name]
        return clade_num

def draw_clade(color_code, total_kmers, kmer_num, palette, clade, x_start, line_shapes, line_color="grey", line_width=5, x_coords=0, y_coords=0):
    #Adapted from:
    #https://github.com/plotly/dash-phylogeny
    """Recursively draw the tree branches, down from the given clade"""
    x_curr = x_coords[clade]
    y_curr = y_coords[clade]
    if str(clade) != "Clade" and str(clade) != "AT":
        tmp_name = str(clade)
        #if tmp_name.count("_") > 0:
        #    tmp_name = tmp_name.split("_")[1]
        #print(tmp_name)
        line_color = color_code[tmp_name]#palette[int(((kmer_num[tmp_name])/total_kmers)*100)-1]
    elif clade.clades:
        
        line_color = palette[int(biggest_num_in_clade(clade, kmer_num))+10]
        #Now we have to find the all of the children, and which one has the highest value   
    # Draw a horizontal line from start to here
    branch_line = get_clade_lines(orientation='horizontal', y_curr=y_curr, x_start=x_start, x_curr=x_curr,
                                  line_color=line_color, line_width=line_width)
    line_shapes.append(branch_line)
    if clade.clades:
        # Draw a vertical line connecting all children
        y_top = y_coords[clade.clades[0]]
        y_bot = y_coords[clade.clades[-1]]
        #Now we have to get the right color though. It should be the max value of any children 
        line_shapes.append(get_clade_lines(orientation='vertical', x_curr=x_curr, y_bot=y_bot, y_top=y_top,
                                           line_color=line_color, line_width=line_width))
        # Draw descendants
        for child in clade:
            draw_clade(color_code, total_kmers, kmer_num, palette, child, x_curr, line_shapes, x_coords=x_coords, y_coords=y_coords)

def parse_counts_complex(raw_counts):
    
    kmer_num = {}
    for l in labels:
        kmer_num[l] = 0

    for c in raw_counts:
        #Now we have to go through each digit 
        if c != "0":
            for t in name_idx_complex[c]:
                kmer_num[t] += 1

    return kmer_num

def parse_counts_complex_for_tree(raw_counts):
    cols = {}
    kmer_num = {}
    cntr = 0
    for l in labels:
        kmer_num[l] = 0
        cols[l] = cntr
        cntr += 1

    data = []
    for c in raw_counts:
        #Now we have to go through each digit 
        if c != "0":
            starting_idx = [0]*num_samples
            for t in name_idx_complex[c]:
                starting_idx[cols[t]] = 1
                kmer_num[t] += 1
            data.append(starting_idx)
    data = np.vstack(data)
    return data, kmer_num

def create_tree(tree_file, x_start_init, x_stop_init, raw_counts, n_skips):
    #Adapted from:
    #https://github.com/plotly/dash-phylogeny
    #adjusted_x_start = int(x_start_init/n_skips)
    #adjusted_x_stop = int(x_stop_init/n_skips)
    #tree_tmp1, kmer_num_tmp = parse_counts_complex_for_tree(raw_counts)
    
    tree_tmp1 = raw_counts #bitvec_to_mat(np.array(raw_counts), 9)
    #We need to sum up the number of times that kmers occur in each column (aka, each sample)
    kmer_num_tmp = tree_tmp1.sum(axis=0)
    
    #print(kmer_num_tmp[0])
    #print(treedata)
    matrix = linkage(
        tree_tmp1.transpose(),
        method='ward',
        metric='euclidean'
    )
    tree_tmp2 = hierarchy.to_tree(matrix, False)
    treedata = get_newick(tree_tmp2, tree_tmp2.dist, labels)
    
    palette = sns.color_palette("RdPu", 130).as_hex()
    total_kmers = max(kmer_num_tmp) #[-1] #kmer_num_tmp["Solqui2"]
    kmer_num = {}
    color_code = {}
    kmer_num_raw = {}
    #cntr = 0
    for k in range(0, len(kmer_num_tmp)):#kmer_num_tmp.keys():
        kmer_num[labels[k]] = float(((kmer_num_tmp[k])/total_kmers)*100)
        kmer_num_raw[labels[k]] = kmer_num_tmp[k]
        #print()
        #print(total_kmers)
        #print(kmer_num_tmp)
        #print(int(((kmer_num_tmp[k])/total_kmers)*100)+25)
        #print()
        color_code[labels[k]] = palette[int(((kmer_num_tmp[k])/total_kmers)*100)+10]
        #cntr += 1
    

    #treedata = "(A, (B, C), (D, E))"
    handle = StringIO(treedata)
    tree = Phylo.read(handle, "newick")
    #tree = Phylo.read(tree_file, "newick")
    x_coords = get_x_coordinates(tree)
    y_coords = get_y_coordinates(tree)
    line_shapes = []
    draw_clade(color_code, total_kmers, kmer_num, palette, tree.root, 0, line_shapes, line_color= "blue",#'rgb(25,25,25)',
                line_width=12, x_coords=x_coords,
                y_coords=y_coords)
    my_tree_clades = x_coords.keys()
    X = []
    Y = []
    text = []
    color = []
    sizes = []
    color_hist = {}
    for cl in my_tree_clades:
        X.append(x_coords[cl])
        Y.append(y_coords[cl])
        text.append(cl.name)
        tmp_name = str(cl.name)
        #if tmp_name.count("_") > 0:
        #    tmp_name = tmp_name.split("_")[1]
        if str(cl.name) == "None":
            this_color = "grey"
            color.append(this_color)
            sizes.append(5)
        else:
            this_color = color_code[tmp_name]
            color.append(this_color)
            sizes.append(10)
    #graph_title = "Pansol Phylgenetic Tree"#create_title(virus_name, nb_genome)
    intermediate_node_color = 'rgb(100,100,100)'
    axis = dict(showline=False,
                zeroline=False,
                showgrid=False,
                showticklabels=False,
                title=''  # y title
                )
    label_legend = []
    an_x = []
    an_y = []
    for t in range(0,len(text)):
        if text[t] != None:
            label_legend.append(text[t])
            an_x.append(X[t])
            an_y.append(Y[t])
    nodes = []
    node = dict(type='scatter',
                    x=X,
                    y=Y,
                    mode='markers',
                    marker=dict(color=color,
                                size=sizes),
                    text=text,  # vignet information of each node
                    #hoverinfo='',
                    showlegend=False
                    )
    nodes.append(node)
    def make_anns(x,y, text, kmer_num, i):
        tmp_txt = text#""
        #if text.count("_") > 0:
        #    tmp_txt = text.split("_")[1]
        #else:
        #    tmp_txt = text
        #print(kmer_num_tmp[tmp_txt])
        tmp_txt += " - " + str(kmer_num[tmp_txt])[:4] + "% " #(" + str(kmer_num_raw[tmp_txt]) + ")" 
        return dict(xref='x', yref='y', x=x, y=y, text="\t" + tmp_txt,
            showarrow=False,
            xanchor='left',
            yanchor='middle',)
    annotations = []
    for i in range(0, len(label_legend)):
        annotations.append(make_anns(an_x[i], an_y[i], label_legend[i], kmer_num, i))
    max_x = max(an_x)
    layout = dict(#title=graph_title,
                  paper_bgcolor='rgba(0,0,0,0)',
                  #dragmode="select",
                  font=dict(family='Balto', size=20),
                  # width=1000,
                  height=600,
                  autosize=True,
                  showlegend=True,
                  xaxis=dict(showline=True,
                             zeroline=False,
                             showgrid=True,  # To visualize the vertical lines
                             ticklen=4,
                             showticklabels=True,
                             title='Branch Length',
                             autorange=False,
                             #range=[0, 0.1]
                             ),
                  yaxis=axis,
                  hovermode='closest',
                  shapes=line_shapes,
                  plot_bgcolor='rgb(250,250,250)',
                  legend={'x': 0, 'y': 1},
                  annotations=annotations,
                  xaxis_range=[0,int(max_x*1.2)]
                  )
    fig = make_subplots(
        rows=1, cols=1,
        specs=[[{"type": "scatter", 'colspan':1}]], #, None, {"type": "bar"} ]], #, {"type": "pie"} , {"type": "pie"} ]],
        horizontal_spacing=0.01,
        #subplot_titles=("This region","CDS","Exons","Genes", "Whole chromosome")
    )
    fig.add_trace(node, row=1, col=1)
    fig.update_layout(layout)
    #Sort the hist bars 
    hist_y = []#list(kmer_num_tmp.values()).sort()
    hist_x = []
    #color_hist = []
    #use the color_codes 
    #for y in range(0, len(kmer_num_tmp)): #list(kmer_num_tmp.keys()):
    #    color_hist.append(color_code[labels[y]])
    #fig.add_trace(go.Bar(y=list(kmer_num_tmp.keys()), x=list(kmer_num_tmp.values()), showlegend=False,
    #    marker=dict(color=color_hist), orientation='h'), 
    #    row=1, col=3)
    fig.update_yaxes(visible=False, showticklabels=False)
    fig.update_layout(margin=dict(
            t=20,
            b=10,
            l=10,
            r=10))
    return fig

def get_local_info(x_all, gene_comp, bar_sum_regional, bar_sum_global):
    fig = make_subplots(
        rows=2, cols=2,
        specs=[[{"type": "bar", "colspan": 2}, None],
           [{"type": "bar"}, {"type": "bar"}]],
        subplot_titles=("Whole chromosome",  "This region", 
            "Genes", ), 
        vertical_spacing=0.1
    )
    x = []
    for i in range(1,num_samples+1):
        x.append(i)
    #x = [1,2,3,4,5,6,7,8,9]
    
    #colors = ["#fde725", "#addc30", "#5ec962", "#28ae80", "#21918c", "#2c728e", "#3b528b", "#472d7b", "#440154"]
    #This region
    #print(bar_sum_global)
    #print(bar_sum_regional)
    y=[(i/sum(bar_sum_regional[1:])*100) for i in bar_sum_regional[1:]]
    y_whole=[(i/sum(bar_sum_global)*100) for i in bar_sum_global]
    
    fig.add_trace(go.Bar(x=x, y=[a_i - b_i for a_i, b_i in zip(y, y_whole)], marker_color=colors, showlegend=False), row=2, col=1)
    fig.add_trace(go.Bar(x=x, y=y_whole, marker_color=colors, showlegend=False), row=1, col=1)
    #Genes
    #y=[(i/sum(gene_comp[1:])*100) for i in gene_comp[1:]]
    #print(gene_comp)
    #print(y_whole)
    #fig.add_trace(go.Bar(x=x, y=y_whole, marker_color="#7400b8", showlegend=False), row=1, col=4)
    fig.add_trace(go.Bar(x=x, y=[a_i - b_i for a_i, b_i in zip(gene_comp, y_whole)], marker_color=colors, showlegend=False), row=2, col=2)
    #fig.update_layout(xaxis_title_text="K-mers shared in X samples", yaxis_title_text='Frequency (log)')
    fig.update_xaxes(title_text="# of genomes", row=2, col=1)
    fig.update_yaxes(title_text="Difference from whole chromosome", row=2, col=1)
    fig.update_yaxes(title_text="Percent of k-mers", row=1, col=1)
    fig.update_layout(height=1000)
    return fig

def plot_interactive( n_skips, #window_filter, poly_order, shared_kmers, 
    layout, gene_comp, bins, names_simp, name, zs_tmp, 
    rep_types, plot_rep, plot_gene, x_start, x_stop, gene_locals, gene_names, exon_locals, exon_names):
    window_filter = SG_window
    poly_order = SG_polynomial_order
    shared_kmers = [1]
    tmp_lst = []
    fig = make_subplots(        
        rows=13, cols=1,
        shared_xaxes=True,
        #vertical_spacing=0.03,
        specs=[[{"type": "scatter"} ], #, None, None, None, None, None],
        [{"type": "scatter", 'rowspan':2} ],#, None, None, None, None, None],
        [None ], #,None, None, None, None, None],
        [{"type": "scatter", 'rowspan':2} ],#, None, None, None, None, None],
        [None ],#,None, None, None, None, None],
        [{"type": "bar", 'rowspan':8} ],#, None, None, None, None, None],
        [None ], #,None, None, None, None, None],
        [None ], #,None, None, None, None, None],
        [None ], #None, None, None, None, None],
        [None ], #,None, None, None, None, None],
        [None ], #,None, None, None, None, None],
        [None ], #,None, None, None, None, None],
        [None ], #,None, None, None, None, None],
        #[{"type": "heatmap"} ],#, None, None, None, None, None], 
        #[{"type": "pie", 'rowspan':2}, {"type": "pie", 'rowspan':2}, {"type": "pie", 'rowspan':2},{"type": "pie", 'rowspan':2}, {"type": "scatter", 'rowspan':2}, {"type": "scatter", 'rowspan':2}],
        #[None,None, None, None, None, None]
        ],
        subplot_titles=("Ref. Sequence Position","Gene Annotation", "Repeat Annotation",  "Conserved K-mers" )
    )

    t = time.time()
    #colors = ["grey","#fde725", "#addc30", "#5ec962", "#28ae80", "#21918c", "#2c728e", "#3b528b", "#472d7b", "#440154"]

    #We are adjusting the start and stop positions to account for the skipping. 
    #The adjusted value should be the index, whereas the x_start and x_stop are the real coordinates 
    adjusted_x_start = int(x_start/n_skips)
    adjusted_x_stop = int(x_stop/n_skips)

    #Get the bins
    bin_size = ((x_stop-x_start)/bins)
    adjusted_bin_size = (bin_size/n_skips)
    cntr = 0
    print("a1:", time.time()-t)
    t = time.time()
    
    #
    cats_tmp = [([0] * (bins+1)) for _ in range(num_samples+1)]
    cntr = 0
    x_tracker = 0
    #Here we are filling in the bins for the main figure.    
    bin_size_int = int(bin_size) + 1

    adjusted_bin_size_init = int(adjusted_bin_size) + 1

    cntr = x_start
    x = []
    while cntr < x_stop:
        x.append(cntr)
        cntr += bin_size_int
    cntr = 0

    for i in names_simp: #[adjusted_x_start:adjusted_x_stop] : #range(x_start, x_stop): #names_simp[x_start:x_stop]:#range(x_start, x_stop):#names_simp:
        #i tells us which y axis to use. 
        #cntr tells us what x axis to use 
        #print(x_tracker)  
        if (cntr % adjusted_bin_size_init) == 0 : #Is the x position within the same bin? Or do we need to move to the next bin? 
            x_tracker += 1 #X-tracker keeps track of which bin we are using (on the x-axis)
        
        #if (i) < len(cats_tmp) and x_tracker < len(cats_tmp[i]):
        if (i) < len(cats_tmp) and x_tracker < len(cats_tmp[i]):
            cats_tmp[i][x_tracker] += 1 
        cntr += 1 #n_skips
    #print(names_simp)

    print("b2:", time.time()-t)
    t = time.time()

    #Add a line plot that will cover the different repetative elements. 
    #plot_reps()
    if plot_rep == True:
        #print("About to add next trace")
        rep_colors = ["#f0f921", "#f8df25", "#fdc627", "#fdaf31", "#f99a3e", "#f3854b", "#e97257", "#de6164", "#d24f71", "#c43e7f", 
                "#b42e8d", "#a21d9a", "#8e0ca4", "#7801a8", "#6100a7", "#4903a0", "#2f0596", "#0d0887", "grey", "grey", "grey"]
        cntr = 0
        for i in rep_types.keys():
            rep_y = []
            rep_types_tmp = []

            j = 0
            while j < len(rep_types[i]) and (rep_types[i][j] < x_stop):
                if x_start < int(rep_types[i][j]): #< x_stop:
                    rep_types_tmp.append(rep_types[i][j])
                    rep_types_tmp.append(rep_types[i][j+1])
                    rep_types_tmp.append(rep_types[i][j+2])
                    rep_y.append(cntr)
                    rep_y.append(cntr)
                    rep_y.append(cntr)
                j += 3
            fig.add_trace(go.Scattergl(x=rep_types_tmp, y=rep_y, line=dict(color=rep_colors[cntr]), name=i, legendgroup="group2", 
                legendgrouptitle_text="Annotations"), row=4, col=1)
            cntr += 1
        fig.update_yaxes(visible=False, row=4, col=1)
        fig.update_xaxes(showticklabels=False, row=4, col=1)
        print("c:", time.time()-t)
        t = time.time()
    print("a2.2:", time.time()-t)
    t = time.time()
    if plot_gene == True:
        gene_locals_tmp = []
        gene_names_tmp = []
        intron_locals_tmp = []
        i = 0
        while i < len(gene_locals):
            if (x_start < int(gene_locals[i]) < x_stop) and (x_start < int(gene_locals[i+1]) < x_stop):
                gene_locals_tmp.append(gene_locals[i])
                gene_locals_tmp.append(gene_locals[i+1])
                gene_locals_tmp.append(gene_locals[i+2])

                gene_names_tmp.append(gene_names[i])
                gene_names_tmp.append(gene_names[i+1])
                gene_names_tmp.append(gene_names[i+2])

                intron_locals_tmp.append(gene_locals[i])
                intron_locals_tmp.append(gene_locals[i+1])
                intron_locals_tmp.append(gene_locals[i+2])
            i += 3

        gene_y = [2]*len(gene_locals_tmp)
        intron_y = [1]*len(intron_locals_tmp)        
        fig.add_trace(go.Scattergl(x=gene_locals_tmp, y=gene_y, line=dict(color="#3b528b", width=10), 
            text=gene_names_tmp, hovertemplate='<br>x:%{x}<br>m:%{text}', legendgroup="group2", 
            name="Gene"), row=2, col=1)
        fig.update_layout(clickmode='event+select')
        #Now add the exons: 
        exon_locals_tmp = []
        exon_names_tmp = []
        i = 0
        while i < len(exon_locals):
            if (x_start < int(exon_locals[i]) < x_stop) and (x_start < int(exon_locals[i+1]) < x_stop):
                exon_locals_tmp.append(exon_locals[i])
                exon_locals_tmp.append(exon_locals[i+1])
                exon_locals_tmp.append(exon_locals[i+2])

                exon_names_tmp.append(exon_names[i])
                exon_names_tmp.append(exon_names[i+1])
                exon_names_tmp.append(exon_names[i+2])
            if int(exon_locals[i+1]) > x_stop:
                i += 1000000000
            i += 3

        exon_y = [1]*len(exon_locals_tmp)
        fig.add_trace(go.Scattergl(x=exon_locals_tmp, y=exon_y, line=dict(color="#addc30", width=5), 
            text=exon_names_tmp, hovertemplate='<br>x:%{x}<br>m:%{text}', legendgroup="group2", 
            name="Exon"), row=2, col=1)
        
        #And now we add the dashed lines between the exons 
        fig.add_trace(go.Scattergl(x=intron_locals_tmp, y=intron_y, line=dict(color="#addc30", width=1, dash='dot'), 
            hovertemplate='<br>x:%{x}<br>m:%{text}', legendgroup="group2", 
            name="Intron"), row=2, col=1)

        fig.update_yaxes(visible=False, range=[-1,4], row=2, col=1)
        fig.update_xaxes(showticklabels=False, row=2, col=1)
    
    print("a2.1:", time.time()-t)
    t = time.time()
    #This is the conserved kmer plotting section
    bar_sum_regional = []
    bar_sum_names = []
    #bar_sum_global = []
    bar_sum_regional.append(sum(cats_tmp[0]))
    bar_sum_names.append(str(0))
    fig.add_trace(go.Bar(x=x, y=cats_tmp[0], name=str(0),
            legendgroup="group1", 
            legendgrouptitle_text="Conserved K-mers",
            marker=dict(color='grey'), 
            marker_line=dict(color='grey')
            ), 
            row=6, col=1 )
    for i in range(1, len(cats_tmp)):
        bar_sum_regional.append(sum(cats_tmp[i]))
        #bar_sum_global.append(names_simp.count(i))
        bar_sum_names.append(str(i))
        fig.add_trace(go.Bar(x=x, y=cats_tmp[i], name=str(i),
            legendgroup="group1", 
            legendgrouptitle_text="Conserved K-mers",
            marker=dict(color=colors[i-1]), 
            marker_line=dict(color=colors[i-1])
            ), 
            row=6, col=1 )
        #cntr += 1
    print("a2.3:", time.time()-t)
    t = time.time()
    fig.update_layout(barmode='stack', bargap=0.0)
    #fig.update_layout(clickmode='select')
    fig.update_xaxes(showticklabels=False, row=6, col=1)
    #Now we will add the smoothing line. There are three parameters that we adjust
    #window_filter = 53
    #poly_order = 3
    #shared_kmers = 1
    for sk in shared_kmers:
        y_tmp = cats_tmp[int(sk)][:-1]
        for i in range(0, int(sk)):
            y_tmp = [a + b for a, b in zip(y_tmp, cats_tmp[int(i)][:-1])] #cats_tmp[int(sk)][:-1]
        fig.add_trace(go.Scatter(x=x, y=signal.savgol_filter(y_tmp,window_filter,poly_order), 
            name="Savitzky-Golay - "+str(sk), marker=dict(color="grey"), mode='lines'), row=6, col=1)
    #fig.update_traces(visible=False, selector=dict(mode="markers"), row=6, col=1)
    print("a2:", time.time()-t)
    t = time.time()

    #Now we add the reference sequence:
    y_ref = [1, 1]
    x_ref = [x_start, x_stop]
    tickvals = []
    ticktxt = []
    cntr = x_start
    yvals = []
    interval = int((x_stop-x_start)/10)
    while cntr <= x_stop:
        tickvals.append(cntr)
        ticktxt.append(str(cntr))
        yvals.append(1)
        cntr += interval
    fig.add_trace(go.Scattergl(x=tickvals, y=yvals, text=ticktxt, 
        textposition='top center', showlegend=False, 
        mode='lines+markers+text', line=dict(color="grey"), 
        marker = dict(size=5, symbol='line-ns')), row=1, col=1)
    fig.update_yaxes(visible=False, range=[0.9,4], row=1, col=1)
    fig.update_xaxes(visible=False, title_text="Sequence position", row=1, col=1)
    fig.update_layout(template="simple_white" ) #,
    fig.update_xaxes(title_text="Sequence position", row=6, col=1)
    fig.update_yaxes(title_text="# of k-mers", row=6, col=1)
    fig.update_layout(height=1000, xaxis_range=[x_start,x_stop], font=dict(
        #family="Balto",
        size=16,
        ))
    print("a3:", time.time()-t)
    t = time.time()
    return fig, bar_sum_names, bar_sum_regional, colors, gene_names_tmp

def make_chr_whole(names_simp, whole_bins, gene_locals, x_start, x_stop, n_skips):
    #Let's use 1000 bins for the whole chromosome
    #whole_bins = 1000

    whole_bin_size = int((len(names_simp)*n_skips)/whole_bins)
    
    #X is the position in the chromosome 
    x = []
    for i in range(0, whole_bins):
        x.append(int(i*whole_bin_size))
    
    #Y will be the same for each 
    y = [1]*whole_bins
    #Z is what we care about. This will require binning. 
    z_1_tmp    = [0]*(whole_bins+1)
    z_9_tmp    = [0]*(whole_bins+1)
    z_cnts_tmp = [0.0]*(whole_bins+1)
    

    z_genes = [0]*(whole_bins+1)
    g = 0
    #for g in range(0,len(gene_locals)):
    while g < len(gene_locals):
        x_bin = int(int(gene_locals[g])/whole_bin_size)
        z_genes[x_bin] += 1
        g += 3

    cntr = 0
    for i in names_simp:
        tmp_x = int(cntr/whole_bin_size)
        if i == 1:
            #tmp_x = int(cntr/whole_bin_size)
            z_1_tmp[tmp_x] += 1
        elif i == num_samples:
            #tmp_x = int(cntr/whole_bin_size)
            z_9_tmp[tmp_x] += 1
            #And now we need to figure out the counts 
        #if cnts[cntr] != "0":
        #    z_cnts_tmp[tmp_x] += (i/int(cnts[cntr]))
        cntr += n_skips
    
    z_1    = [0]*(whole_bins+1)
    z_9    = [0]*(whole_bins+1)
    #z_cnts = [0.0]*(whole_bins+1)
    
    for i in range(0, (whole_bins+1)):
        z_1[i]    = (z_1_tmp[i]/whole_bin_size)*100
        z_9[i]    = (z_9_tmp[i]/whole_bin_size)*100
        #z_cnts[i] = z_cnts_tmp[i]/whole_bin_size #(z_cnts_tmp[i]/z_9_tmp[i])

    #x_all = x + x + x
    #y_all = [1]*whole_bins + [2]*whole_bins + [3]*whole_bins
    #z_all = z_genes + z_1 + z_9
    return x, y, z_1, z_9, z_genes #, z_cnts #x_all, y_all, z_all,  x, z_cnts

def read_chr_whole():
    with open(bins_file, "r") as f:
        line1 = f.readline()
        line2 = f.readline()
        line3 = f.readline()
        cntr = 0
        x = []
        z_9 = []
        z_1 = []
        y = []
        while line1:
            x_tmp = line1.strip().split(",")[:-1]
            z_9_tmp = line2.strip().split(",")[:-1]
            x_1_tmp = line3.strip().split(",")[:-1]
            x.append([])
            z_9.append([])
            z_1.append([])
            y.append([1]*len(x_tmp))
            for i in range(0, len(x_tmp)):
                x[cntr].append(int(x_tmp[i]))
                z_9[cntr].append(int(z_9_tmp[i]))
                z_1[cntr].append(int(x_1_tmp[i]))
            line1 = f.readline()
            line2 = f.readline()
            line3 = f.readline()
            cntr += 1
    
    return x, z_1, z_9, y

def make_gene_whole_chr(x, locs):
    z_genes = [0]*(len(x)+1)
    g = 0
    #print(locs[-100:])
    #print(len(locs))
    #y = []
    #for g in range(0,len(gene_locals)):
    while g < len(locs):
        #print(g)
        #print(locs[g])
        x_bin = int(int(locs[g])/200000)
        #print(x_bin)
        #print("***********")
        #if x_bin < len(z_genes):
        z_genes[x_bin] += 1
        g += 3
    #for i in z_genes:
    #    y.append("chr1_genes")
    return z_genes

def plot_chr_whole(x, z_1, z_9, z_genes, x_start, x_stop, y): 
    
    #print(z_9)
    chr_fig = make_subplots(rows=3, cols=1, 
        specs=[[{"type": "heatmap",}], [{"type": "heatmap",}], [{"type": "heatmap",}]],
        shared_xaxes=True,
        subplot_titles=("K-mer and gene density accross whole chromosome", "",
            ""),
        vertical_spacing = 0.0
        )

    chr_fig.add_trace(go.Heatmap(x=x, z=z_9, y=y, type = 'heatmap', colorscale='magma', showlegend=False, showscale=False), row=1, col=1)
    chr_fig.add_trace(go.Scatter(x=[x_start, x_start, None, x_stop, x_stop, None, x_start, x_stop, ], showlegend=False,
                   y=[0.5, 1.5, None, 0.5, 1.5, None, 1.45, 1.45 ],
                   mode='lines',
                   line_color='#1dd3b0', line_width=8), row=1, col=1)
    ###
    chr_fig.add_trace(go.Heatmap(x=x, z=z_1, y=y, type = 'heatmap', colorscale='magma', showscale=False), row=2, col=1)
    chr_fig.add_trace(go.Scatter(x=[x_start, x_start, None, x_stop, x_stop], showlegend=False,
                   y=[0.5, 1.5, None, 0.5, 1.5],
                   mode='lines',
                   line_color='#1dd3b0', line_width=8), row=2, col=1)
    ###
    chr_fig.add_trace(go.Heatmap(x=x, z=z_genes, y=y, type = 'heatmap', colorscale='magma', showscale=False ), row=3, col=1)
    chr_fig.add_trace(go.Scatter(x=[x_start, x_start, None, x_stop, x_stop, None, x_start, x_stop,], showlegend=False, 
                   y=[0.5, 1.5, None, 0.5, 1.5, None, 0.55, 0.55],
                   mode='lines',
                   line_color='#1dd3b0', line_width=8), row=3, col=1)
    ###
    #chr_fig.add_trace(go.Scattergl(x=x, y=z_cnts, showlegend=False ), row=4, col=1)
    #chr_fig.add_trace(go.Scatter(x=[x_start, x_start, None, x_stop, x_stop, None, x_start, x_stop, ],
    #               y=[0, 3, None, 0, 3, None, 0, 0 ], showlegend=False,
    #               mode='lines',
    #               line_color='#1dd3b0', line_width=4), row=4, col=1)
    #chr_fig.update_xaxes(visible=False, title_text="Sequence position",)
    #chr_fig.update_yaxes(visible=False, title_text="All", row=1, col=1)
    #chr_fig.update_yaxes(visible=False, title_text="Unique", row=2, col=1)
    #chr_fig.update_layout(yaxis_title="Title")
    
    chr_fig.update_yaxes( range=[0.5,1.5], showticklabels=False, row=1, col=1)
    chr_fig.update_yaxes( range=[0.5,1.5], showticklabels=False, row=2, col=1)
    chr_fig.update_yaxes( range=[0.5,1.5], showticklabels=False, row=3, col=1)
    
    #chr_fig.update_yaxes(visible=False, range=[0,3], row=4, col=1)

    chr_fig.update_xaxes(fixedrange=True, row=1, col=1)
    chr_fig.update_xaxes(fixedrange=True, row=2, col=1)
    chr_fig.update_xaxes(fixedrange=True, row=3, col=1)
    chr_fig.update_xaxes(title_text="Sequence position", row=3, col=1)
    chr_fig.update_yaxes(title_text="Univ.", row=1, col=1)
    chr_fig.update_yaxes(title_text="Uniq.", row=2, col=1)
    chr_fig.update_yaxes(title_text="Genes", row=3, col=1)
    #chr_fig.update_xaxes(fixedrange=True, row=4, col=1)
    chr_fig.update_layout(clickmode='event+select', dragmode="select", selectdirection='h')
    #chr_fig.update_traces(showscale=False)
    chr_fig.update_layout(height=350)
    chr_fig.update_layout(margin=dict(
            b=10,
            l=10,
            r=10),
        font=dict(
            #family="Balto",
            size=15,
        )
    )
    
    return chr_fig

def plot_whole_genome(x, z_1, z_9, y):
    spec = []
    sub_titles = []
    h = num_chrs*250
    for i in range(0, num_chrs):
        spec.append([{"type": "heatmap",}])
        spec.append([{"type": "heatmap",}])
        spec.append([{"type": "heatmap",}])
        sub_titles.append("Chromosome "+str(i+1))
        sub_titles.append("")
        sub_titles.append("")
    
    wg_fig = make_subplots(rows=(3*num_chrs), cols=1, 
        specs=spec, #[[{"type": "heatmap",}], [{"type": "heatmap",}], [{"type": "heatmap"}]],
        shared_xaxes=True,
        subplot_titles=sub_titles, #("K-mer and gene density accross whole chromosome", "", ""),
        vertical_spacing = 0.0
        )

    for i in range(1, len(x)+1):
        #cntr = 0
        #while cntr < num_chrs:
        wg_fig.add_trace(go.Heatmap(x=x[i-1], z=z_9[i-1], y=y[i-1], type = 'heatmap', colorscale='magma', showlegend=False,showscale=False), row=((i*3)-2), col=1)
        wg_fig.add_trace(go.Heatmap(x=x[i-1], z=z_1[i-1], y=y[i-1], type = 'heatmap', colorscale='magma', showscale=False), row=((i*3)-1), col=1)
        
        if i == 1:
            wg_fig.update_layout(xaxis={'side': 'top'}) 

        wg_fig.update_yaxes( range=[0.5,1.5], showticklabels=False, row=((i*3)-2), col=1)
        wg_fig.update_yaxes( range=[0.5,1.5], showticklabels=False, row=((i*3)-1), col=1)
        
        wg_fig.update_xaxes(fixedrange=True, row=((i*3)-2), col=1)
        wg_fig.update_xaxes(fixedrange=True, row=((i*3)-1), col=1)
        
        wg_fig.update_xaxes(title_text="Sequence position", row=((i*3)), col=1)
        wg_fig.update_yaxes(title_text="Universal", row=((i*3)-2), col=1)
        wg_fig.update_yaxes(title_text="Unique", row=((i*3)-1), col=1)
    
    wg_fig.update_layout(clickmode='event', plot_bgcolor='rgba(0,0,0,0)')
    wg_fig.update_layout(height=h)
    wg_fig.update_layout(margin=dict(
            b=10,
            l=10,
            r=10),
        font=dict(
            #family="Balto",
            size=14,
        ),
        paper_bgcolor='rgba(0,0,0,0)')
    return wg_fig

def plot_gene_content(gene_content, sort_by, colors, uniq_avg, univ_avg, local_gene_list):
    #print(gene_content)
    colors = ['#ffd60a', '#440154']
    df = pd.DataFrame(gene_content)
    x = []
    cntr = 0
    for i in range(0, len(df['Universal'])):
        x.append(cntr)
        cntr +=1
    cntr = 0
    df_sorted = df.sort_values(sort_by[-1])
    df_sorted['X'] = x
    #print(df_sorted)
    fig = go.Figure(data=[go.Scattergl(x=x, y=df_sorted[sort_by[cntr]], text=df_sorted['Names'], marker=dict(color=colors[cntr]),
            name="% "+sort_by[cntr], mode="markers")])
    cntr += 1
    while cntr < len(sort_by): #s in sort_by:  
        fig.add_trace(go.Scattergl(x=x, y=df_sorted[sort_by[cntr]], text=df_sorted['Names'],  marker=dict(color=colors[cntr]),
            name="% " + sort_by[cntr], mode="markers"))
        cntr += 1
        #fig.add_trace(go.Scattergl(x=x, y=df_sorted[cntr], text=df_sorted['Names'], marker=dict(color=colors[cntr]),
        #    name="% "+sort_by[cntr], mode="markers"))
    fig.update_layout(clickmode='event+select')    
    #fig.add_trace(go.Scattergl(x=x, y=df_sorted['Universal'], text=df_sorted['Names'],  marker=dict(color='#440154'),
    #    name="% Universal", mode="markers"))
    fig.update_layout(hovermode='x unified')

    df2 = df_sorted.loc[df_sorted['Names'].isin(local_gene_list)]
    #print(df2)
    fig.add_trace(go.Scattergl(x=df2['X'], y=df2['Universal'], marker=dict(color='#FF2192', size=10), 
        mode="markers", hoverinfo='skip', showlegend=False))
    fig.add_trace(go.Scattergl(x=df2['X'], y=df2['Unique'], marker=dict(color='#FF2192', size=10), 
        mode="markers", hoverinfo='skip', showlegend=False))

    fig.add_hline(y=uniq_avg, line_dash='dash', line_color='goldenrod')
    fig.add_hline(y=univ_avg, line_dash='dash', line_color='#440154')
    fig.update_layout(height=600)
    fig.update_layout(
        title={"text":"k-mer content in individual genes",
            "xanchor": "center",
            "x":0.5,
            "y":0.9,
            "yanchor": "top"
            },
        xaxis_title=sort_by[0] + " k-mer content rank",
        yaxis_title="% of k-mers",
        margin=dict(
            t=20,
            b=10,
            l=10,
            r=10),
        font=dict(
            #family="Balto",
            size=16,
        ),
        #plot_bgcolor='rgba(0,0,0,0)'
    )
    return fig

def update_output_div(input_value):
    return f'Output: {input_value}'

def read_pangenome_comp(ordered_labels):
    genome_comp_totals = {}
    genome_names = []
    cats = []
    for i in range(0,num_samples):
        cats.append([])
    #Each row here is a different k-mer comp 
    with open(genome_comp_file, "r") as f:
        line1 = f.readline()
        line2 = f.readline()
        while line1:
            if line1.strip() == "Saet_bc2059":
                thisname = "Sins2"
            elif line1.strip() == "Saet_bc2058":
                thisname = "M82"
            else:
                thisname = line1.strip()
            genome_names.append(thisname)
            genome_comp_totals[thisname] = line2.strip().split(':')[1].split(',')
            tmp = line2.strip().split(':')[1].split(',')
            tmp_total = 0
            for i in range(0,num_samples):
                tmp_total += int(tmp[i])
            for i in range(0,num_samples):
                cats[i].append(float(int(tmp[i])/tmp_total)*100)

            #genome_comp_totals[line1.strip()] = line2.strip().split(":")[1].split(',')
            line1 = f.readline()
            line2 = f.readline()
    fig = make_subplots(rows=1, cols=1)
    for g in range(0,len(ordered_labels)):    
        fig.add_trace(go.Bar(y=genome_names, x=cats[g], name=str(g+1), #orientation='h',
            #legendgroup="group1", 
            legendgrouptitle_text="# Samples",
            marker=dict(color=colors[g]), 
            marker_line=dict(color=colors[g]), orientation='h'
            #log_y=True
            ), row=1, col=1)
    #fig.update_xaxes(type="log")
    
    #print(genome_comp_totals)
    #fig = go.Figure([go.Bar(x=genome_comp_totals.keys(), y=genome_comp_totals.)])
    fig.update_layout(barmode='stack' )#,orientation='h')
    fig.update_layout(height=1500, font=dict(
        #family="Balto",
        size=26,
        ), plot_bgcolor='rgba(0,0,0,0)')
    return fig, genome_comp_totals

def read_genome_comp(anchor_name):
    chrs_comp = []
    #x = []
    cntr = 1
    with open(genome_comp_pre+anchor_name+genome_comp_post, "r") as f:
        line = f.readline()
        while line:
            #x.append(str(cntr))
            tmp = line.strip().split(':')[1].split(',')[:-1]
            tmp2 = []
            for i in tmp:
                tmp2.append(int(i))
            chrs_comp.append(tmp2)
            cntr += 1
            line = f.readline()
    fig = make_subplots(rows=len(chrs_comp), cols=1)
    x = []
    for i in range(1, num_samples+1):
        x.append(str(i))
    #y=[]
    #tots = sum(chrs_comp[0])
    #for j in chrs_comp[0]:
    #    y.append(int(j)/tots*100)
    #fig.add_trace(go.Bar(x=x, y=y, marker_color=colors, ), #marker_color=colors[1:len(x)]), 
    #    row=1,col=1)
    for i in range(0,len(chrs_comp)):
        y=[]
        tots = sum(chrs_comp[i])
        for j in chrs_comp[i]:
            y.append(int(j)/tots*100)
        fig.add_trace(go.Bar(x=x, y=y, marker_color=colors, showlegend=False,), #marker_color=colors[1:len(x)]), 
            row=i+1,col=1)
    fig.update_yaxes(type="log")
    fig.update_layout(paper_bgcolor='rgba(0,0,0,0)')
    return fig, chrs_comp

def make_all_genome_dend():
    file_names = {}
    with open(mash_filenames) as f:
        for line in f:
            line = line.rstrip()
            #print(line)
            if line == "Saet_bc2058.msh":
                thisname = "M82.msh"
            elif line == "Saet_bc2059.msh":
                thisname = "Sins2.msh"
            else:
                thisname = line
            file_names[thisname.split('.')[0]] = thisname.split('.')[0]
    sample_list = list(file_names.values())
    dim = len(sample_list)
    dist_mat = np.zeros((dim, dim), np.float64)

    with open(mash_edges) as f:
        for line in f:
            f_tmp, t_tmp, d, p, x = line.rstrip().split("\t")
            f_tmp2 = f_tmp.split('/')[6] #f_tmp.split('-')[1]
            t_tmp2 = t_tmp.split('/')[6]
            if f_tmp2 == "revio_test":
                f = "Saet_" + f_tmp.split('/')[-1].split('.')[0]
                if f == "Saet_bc2059":
                    f = "Sins2"
                elif f == "Saet_bc2058":
                    f = "M82"
            else:
                f = f_tmp2

            if t_tmp2 == "revio_test":
                t = "Saet_" + t_tmp.split('/')[-1].split('.')[0]
                if t == "Saet_bc2059":
                    t = "Sins2"
                elif t == "Saet_bc2058":
                    t = "M82"
            else:
                t = t_tmp2
            
            #print(f)
            if f in file_names:
                new_f = file_names[f]
                new_t = file_names[t]
                
                i = sample_list.index(new_f)
                j = sample_list.index(new_t)
                dist_mat[i][j] = d
                dist_mat[j][i] = d

    df = pd.DataFrame(dist_mat, columns=sample_list)
    labels = df.columns
    branch_colors = ['purple','purple','purple','purple','purple','purple']
    fig = ff.create_dendrogram(df, colorscale=branch_colors, labels=labels, orientation='bottom' ) #, orientation='bottom'
    fig.update_layout(paper_bgcolor='rgba(0,0,0,0)', font=dict(
            size=24,))
    for i in range(len(fig['data'])):
        fig['data'][i]['yaxis'] = 'y2'

    dendro_side = ff.create_dendrogram(df, orientation='right', colorscale=branch_colors, )
    dendro_side.update_layout(paper_bgcolor='rgba(0,0,0,0)')
    for i in range(len(dendro_side['data'])):
        dendro_side['data'][i]['xaxis'] = 'x2'
    # Add Side Dendrogram Data to Figure
    for data in dendro_side['data']:
        fig.add_trace(data)
    dendro_leaves = dendro_side['layout']['yaxis']['ticktext']
    dendro_leaves = list(map(int, dendro_leaves))
    data_dist = pdist(df)
    heat_data = squareform(data_dist)
    heat_data = heat_data[dendro_leaves,:]
    heat_data = heat_data[:,dendro_leaves]
    heatmap = [
        go.Heatmap(
            x = dendro_leaves,
            y = dendro_leaves,
            z = heat_data,
            #labels=labels
            colorscale = 'plasma_r'#'spectral'#'haline'#'portland_r'
        )
    ]
    heatmap[0]['x'] = fig['layout']['xaxis']['tickvals']
    heatmap[0]['y'] = dendro_side['layout']['yaxis']['tickvals']
    for data in heatmap:
        fig.add_trace(data)
    fig.update_layout({'width':1500, 'height':1500,
                             'showlegend':False, 'hovermode': 'closest',
                             'paper_bgcolor':'rgba(0,0,0,0)',
                             'plot_bgcolor':'rgba(0,0,0,0)'
                             })
    # Edit xaxis
    fig.update_layout(xaxis={'domain': [.2, 1],
        'mirror': True,
        'showgrid': False,
        'showline': False,
        'zeroline': False,
        'gridcolor':'rgba(0,0,0,0)',
        'ticks':""})
    # Edit xaxis2
    fig.update_layout(xaxis2={'domain': [0, .2],
        'mirror': True,
        'showgrid': False,
        'showline': False,
        'zeroline': False,
        'showticklabels': False,
        'gridcolor':'rgba(0,0,0,0)',
        'ticks':""})
    fig.update_layout(yaxis={'domain': [0, 0.825],#
        'mirror': True,
        'showgrid': False,
        'showline': False,
        'zeroline': False,
        'showticklabels': False,
        'gridcolor':'rgba(0,0,0,0)',
        'ticks': ""
        })
    # Edit yaxis2
    fig.update_layout(yaxis2={'domain':[.79, .975],
        'mirror': True,
        'showgrid': False,
        'showline': False,
        'zeroline': False,
        'showticklabels': False,
        'gridcolor':'rgba(0,0,0,0)',
        'ticks':""})
    return fig, labels

def read_genome_size_files():
    sizes = []
    seqs = []
    with open(genome_size_f, "r") as f:
        line_name = f.readline()
        line = f.readline()
        while line:
            print(line_name)
            if line_name.strip() == "Saet_bc2059":
                line_name = "Sins2"
            elif line_name.strip() == "Saet_bc2058":
                line_name = "M82"
            sizes.append([line_name.strip(), int(line.split('\t')[1])])
            seqs.append([line_name.strip(), int(line.split('\t')[0])])
            #seqs[line_name.strip()] = int(line.split('\t')[1])
            
            line_name = f.readline()
            line = f.readline()
    sizes_df = pd.DataFrame(sizes, columns=['Sample','Size'])
    seqs_df = pd.DataFrame(seqs, columns=['Sample','Seq'])
    
    sizes_df = sizes_df.sort_values(by=['Size'])
    seqs_df = seqs_df.sort_values(by=['Seq'])

    fig_sizes = make_subplots(rows=1, cols=1)
    fig_seqs = make_subplots(rows=1, cols=1)

    fig_sizes.add_trace(go.Scattergl(x=sizes_df['Sample'],y=sizes_df['Size']))
    fig_sizes.update_yaxes(title_text="Size of genome",)
    fig_sizes.add_vline(x=anchor_name, line_dash="dash", line_color="darkblue")

    fig_seqs.add_trace(go.Scattergl(x=seqs_df['Sample'],y=seqs_df['Seq']))
    fig_seqs.update_yaxes(title_text="# of scaffolds",)
    fig_seqs.add_vline(x=anchor_name, line_dash="dash", line_color="darkblue")
    fig_seqs.update_layout(font=dict(size=20,))
    fig_sizes.update_layout(font=dict(size=20,))
    return fig_sizes, fig_seqs

def plot_avgs(genome_comp_totals):
    #print(genome_comp_totals)
    fig = make_subplots(rows=1, cols=1)
    tmp = []
    for k in genome_comp_totals.keys():
        running_sum = 0
        total = 0
        for i in range(1, len(genome_comp_totals[k])):
            total += int(genome_comp_totals[k][i-1])
            running_sum += (int(genome_comp_totals[k][i-1])*i)
        #print("*******************")
        tmp.append([k,(running_sum/total)])
    #print(tmp)
    df = pd.DataFrame(tmp, columns=['Sample','Avg'])
    #print(df)
    df = df.sort_values(by=['Avg'])
    #print(df)
    fig.add_trace(go.Scattergl(x=df['Sample'],y=df['Avg']))
    fig.add_vline(x=anchor_name, line_dash="dash", line_color="darkblue")
    #fig.add_vline(x=anchor_name)
    fig.update_yaxes(title_text="Average k-mer",)
    fig.update_layout(font=dict(size=20,))
    #fig.add_trace(go.Scattergl(x=list(genome_comp_totals.keys()),y=tmp[]))
    return fig

def make_avg_kmer_fig(chrs_comp):
    fig = make_subplots(rows=1, cols=1)
    x = []
    y = []
    for c in range(0, num_chrs):#len(chrs_comp)):
        #print(chrs_comp[c])
        #print(c)
        x.append(chrs_list[c])
        total = 0
        running_sum = 0
        for i in range(1, len(chrs_comp[c])+1):
            total += int(chrs_comp[c][i-1])
            running_sum += (int(chrs_comp[c][i-1])*i)
        y.append(running_sum/total)
    fig.add_trace(go.Scattergl(x=x, y=y))
    fig.update_yaxes(title_text="Average k-mer",)
    fig.update_xaxes(title_text="Chromosome",)
    fig.update_layout(plot_bgcolor='rgba(0,0,0,0)',  font=dict(size=20))
    return fig

def make_genes_per_chr_fig(gene_names):
    fig = make_subplots(rows=1, cols=1)
    x = []
    y = []
    for c in range(0, len(gene_names)):
        x.append(chrs_list[c])
        y.append(int(len(gene_names[chrs_list[c]])/3)/chr_lens[c])
    fig.add_trace(go.Scattergl(x=x, y=y))
    fig.update_yaxes(title_text="Gene density",)
    fig.update_xaxes(title_text="Chromosome",)
    fig.update_layout(plot_bgcolor='rgba(0,0,0,0)',  font=dict(size=20))
    return fig

def make_gene_per_genome_fig(gene_content):
    fig = make_subplots(rows=1, cols=1)
    colors = ['#ffd60a', '#440154']
    
    df = pd.DataFrame(gene_content)
    #print(df)
    #df = df.transpose()
    x = []
    #print(df.loc['Universal'])
    d_new = {"Unique":[],"Universal":[],"Names":[]}
    for d in df.loc['Universal']:
        d_new['Universal'] += d
    for d in df.loc['Unique']:
        d_new['Unique'] += d
    for d in df.loc['Names']:
        d_new['Names'] += d
    df = pd.DataFrame(d_new)
    
    cntr = 0
    for i in range(0, len(df['Universal'])):
        x.append(cntr)
        cntr +=1
    #cntr = 0
    df_sorted = df.sort_values('Universal')
    df_sorted['X'] = x
    #print(df_sorted)
    fig.add_trace(go.Scattergl(x=x, y=df_sorted['Unique'], text=df_sorted['Names'], marker=dict(color=colors[0]),
            name="% "+"Unique", mode="markers"))
    fig.add_trace(go.Scattergl(x=x, y=df_sorted['Universal'], text=df_sorted['Names'], marker=dict(color=colors[1]),
            name="% "+"Universal", mode="markers"))
    fig.update_layout(plot_bgcolor='rgba(0,0,0,0)', font=dict(size=20) )
    fig.update_xaxes(title_text="Genes",)
    fig.update_yaxes(title_text="% of gene",)
    #cntr += 1
    #while cntr < len(sort_by): #s in sort_by:  
    #    fig.add_trace(go.Scattergl(x=x, y=df_sorted[sort_by[cntr]], text=df_sorted['Names'],  marker=dict(color=colors[cntr]),
    #        name="% " + sort_by[cntr], mode="markers"))
    #    cntr += 1
    fig.update_layout(hovermode='x unified')
    return fig

##### READ DATA
t = time.time()
cnts_tmp = []
x = []
all_chrs = {} #This is where we keep the raw, str, counts 
rep_types = get_init_rep_types()
#anchor = KmerRef(bit_file_prefix)

anchor = KmerBitmap(bit_file_prefix)

#name, start, end = args.coords
#bits = bitmap.query(anchor_name + "." + chrs, chr_start, chr_end, n_skips_start)

#anchor.query(anchor_name + ".chr1" , 75000000, 85000000, n_skips_start)

for i in range(1, num_chrs+1):
    chrs="chr" + str(i)
    chr_start = 0
    chr_end = chr_lens[i-1] #308452471
    print(anchor_name + "." + chrs)
    all_chrs[chrs] = anchor.query(chrs, chr_start, chr_end, n_skips_start)
    #all_chrs[chrs] = anchor.get_counts(chrs, chr_start, chr_end, n_skips_start)


chrs = "chr1"
print("a0:", time.time()-t)
t = time.time()
#chrs = "chr1"
#Parsing the counts will return the number of genomes present at each kmer position, based on the counts
#names_simp = parse_counts_simple(all_chrs[chrs])
#tmp = bitvec_to_mat(np.array(all_chrs[chrs]), 9)
names_simp = all_chrs[chrs].sum(axis=1)#tmp.sum(axis=1)
#names_simp = bits.sum(axis=1)
#print(len(names_simp))
bar_sum_global = {}

'''
with open(gene_content_by_chr_f, "r") as f:
    line = f.readline()
    while line:
        c = line.split(":")[0]
        tmp = line.strip().split(':')[1].split(',')[:-1]
        bar_sum_global[c] = [float(t) for t in tmp]
        line = f.readline()
'''

for k in range(1, num_chrs+1):
    c = "chr"+ str(k)
    bar_sum_global_tmp = []
    for i in range(1,(num_samples+1)):
        bar_sum_global_tmp.append((all_chrs[c].sum(axis=1)==i).sum())
    bar_sum_global[c] = bar_sum_global_tmp

#print(bar_sum_global)
#bar_sum_global = names_simp
tmp = 0
#We actually need to 
#print(bar_sum_global)
#for i in range(0, (num_samples+1)):
#    bar_sum_global.append(names_simp.count(i))
#print("Parsed first counts file")
print("a0.2:", time.time()-t)
t = time.time()

#This next step is reading in the traditional counts file, which will tell us how reptative stuff is 
#y=[1,10,100,1000,10000]
 #read_count_file(count_file, y)

#Read in the repeats file (this has the annotated repeats from Shujun) 
rep_types = read_annotations(rep_file, rep_types)
#print(rep_types.keys())
#print(rep_types["chr1"])
print("a0.4:", time.time()-t)
t = time.time()
#print("Read annotation file")
#grep -v "#" Solqui2_genes_1.0.0.gff | awk '{if($1=="chr1") print $0}' | awk '{if($3=="gene") print $1"\t"$9"\t"$4"\t"$5}' > Squi2.genes.chr1.txt
zs_tmp = []
gene_names = {}
gene_locals = {}
exon_locals = {}
exon_names = {}
cds_locals = {}
cds_names = {}
#exon_comp = {}
gene_comp = {}

for i in range(1, num_chrs+1):
    tmp_chr = "chr"+str(i)
    gene_names[tmp_chr] = []
    gene_locals[tmp_chr] = []
    exon_locals[tmp_chr] = []
    exon_names[tmp_chr] = []
    gene_comp[tmp_chr] = [0]*(num_samples+1)
    #exon_comp[tmp_chr] = [0]*(num_samples+1)

gene_locals, gene_names, exon_locals, exon_names, cds_locals, gene_content, gene_anns = read_gene_annotations(gene_file, 
    gene_names, gene_locals, exon_locals, exon_names, cds_locals, n_skips_start, "chr1")

#We need to get the genic kmer composition 

#gene_comp = [0]*(num_samples+1)
print("a0.45:", time.time()-t)
t = time.time()
#cntr = 0
#cntr2 = 0
#print(gene_locals)
def get_idx_content(gene_locals, gene_cntr, j):
    est_pos_start = int(gene_locals[gene_cntr]/n_skips_start)
    est_pos_stop = int(gene_locals[gene_cntr+1]/n_skips_start)
    tmp = (names_simp[est_pos_start:est_pos_stop]==j).sum()
    return tmp




#We are going through the gene locations (gene_locals) and characterizing their k-mer composition
'''
for i in range(1, num_chrs+1):
    tmp_chr = "chr"+str(i)
    cntr = 0
    cntr2 = 0
    while cntr < len(gene_locals[tmp_chr]):
        tmp1 = 0
        for j in range(0, 3):
            tmp2 = get_idx_content(gene_locals[tmp_chr], cntr, j) #names_simp[(gene_locals[cntr]/n_skips):(gene_locals[cntr+1]/n_skips)].count(j)
            tmp1 += tmp2
            gene_comp[tmp_chr][j] += tmp2 #names_simp[gene_locals[cntr]:gene_locals[cntr+1]].count(j)
        cntr2 += 1
        for j in range(3, (num_samples+1)):
            gene_comp[tmp_chr][j] += get_idx_content(gene_locals[tmp_chr], cntr, j) #names_simp[(gene_locals[cntr]*n_skips):(gene_locals[cntr+1]*n_skips)].count(j)
        cntr += 3
'''
with open(gene_content_by_chr_f, "r") as f:
    line = f.readline()
    while line:
        c = line.split(":")[0]
        tmp = line.strip().split(':')[1].split(',')[:-1]
        gene_comp[c] = [float(t) for t in tmp]
        line = f.readline()
    cntr = 0
print("a0.5:", time.time()-t)
t = time.time()
#This will have the figure componants that we need 


layout = go.Layout(
    margin=go.layout.Margin(
        l=10,  # left margin
        r=10,  # right margin
        b=10,  # bottom margin
        t=10  # top margin
    )
    )
##Set the smoothing variables. These will be changable later on. 
window_filter = 53
SG_window =53
poly_order = 3
SG_polynomial_order = 3
SG_check = [1]

shared_kmers = [1]
#print(len(names_simp))
x_start_init_adjust = int(x_start_init/n_skips_start)
x_stop_init_adjust = int(x_stop_init/n_skips_start)
#### NEED TO ADJUST THE START INIT AND STOP INIT SO THAT WE TAKE THE SKIPS INTO ACCOUNT
#Here is where we make the main plot 
#print("*********")
#print(len(gene_locals["chr1"]))
fig, bar_sum_names, bar_sum_regional, colors, gene_names_tmp = plot_interactive( n_skips_start, #window_filter, poly_order, shared_kmers,
    layout, #exon_comp["chr1"], 
    gene_comp["chr1"], bins, names_simp[x_start_init_adjust:x_stop_init_adjust], "chr1", zs_tmp, rep_types["chr1"], True, True, 
    x_start_init, x_stop_init, gene_locals["chr1"], gene_names["chr1"], exon_locals["chr1"], exon_names["chr1"],
    )
#print(gene_names_tmp)
gene_names_local = []
cntr = 0
while cntr < len(gene_names_tmp):
    gene_names_local.append(gene_names_tmp[cntr])
    cntr += 3
#print(gene_names_local)
#print(len(gene_locals["chr1"]))
sort_by = ["Unique","Universal"]
tmp = [(i/sum(bar_sum_global['chr1'])*100) for i in bar_sum_global['chr1']]
uniq_avg = tmp[1]
univ_avg = tmp[-1]
gene_content_plot = plot_gene_content(gene_content["chr1"], sort_by, colors, uniq_avg, univ_avg, gene_names_local)
#print("AFTER GENE CONTENT", len(gene_locals["chr1"]))
#sample_kmer_counts = parse_counts_complex(all_chrs[chrs][x_start_init:x_stop_init])
#x, y, z_1, z_9, z_genes = make_chr_whole(names_simp, 1000, gene_locals, x_start_init, x_stop_init)
x, z_1, z_9, y  = read_chr_whole()
#print("AFTER READ CHR", len(gene_locals["chr1"]))
z_genes = {}
for i in range(1, num_chrs+1):
    z_genes["chr"+str(i)] = make_gene_whole_chr(x[i-1], gene_locals["chr"+str(i)])
#print("AFTER MAKING GENE BINS", len(gene_locals["chr1"]))
chr_fig = plot_chr_whole(x[0], z_1[0], z_9[0], z_genes["chr1"], x_start_init, x_stop_init, y[0])
whole_genome_fig = plot_whole_genome(x, z_1, z_9, y)
#chr_fig = plot_chr_whole(names_simp, full_cnts, 1000, gene_locals, x_start_init, x_stop_init)
print("a4:", time.time()-t)
t = time.time()

#Now we need to make the pangenome plot
#This is going to be a bar chart 

whole_genome_hists, chrs_comp = read_genome_comp(anchor_name)
all_genomes_dend, ordered_labels = make_all_genome_dend()
pangenome_comp, genome_comp_totals = read_pangenome_comp(ordered_labels)
pangenome_avg = plot_avgs(genome_comp_totals)
pangenome_sizes, pangenome_num_seqs = read_genome_size_files()
avg_kmer_per_chr_fig = make_avg_kmer_fig(chrs_comp)
genes_per_chr_fig = make_genes_per_chr_fig(gene_names)
gene_content_per_genome_fig = make_gene_per_genome_fig(gene_content)
#pangenome_avg = 
#print(all_chrs[chrs])
#Set up the dash app 
chrs = "chr1"
#dbc_css = "https://cdn.jsdelivr.net/gh/AnnMarieW/dash-bootstrap-templates/dbc.min.css"
app = dash.Dash(external_stylesheets=["https://www.w3schools.com/w3css/4/w3.css"]) 
#app = Dash(__name__, external_stylesheets=[dbc.themes.PULSE, dbc_css])
config = {"toImageButtonOptions" : {"format" : "svg", "width" : None, "height" : None},}
app.layout = html.Div([
    #html.Div(id = 'parent', children = [
    html.H1(id = 'H1', children = 'Panagram', style = {'textAlign':'center', "font-size": 48}), 
    dcc.Tabs(id="tabs-graph", style={"font-size": 36}, children=[
        dcc.Tab(label='Pangenome', style={"background-color": "lightgrey", "font-size": 36,}, 
            children=[
                html.Div(id="PanInfo1", children=[
                    html.I("Panagram of " + str(num_samples) + " genomes"),
                    html.Br(),
                    html.Label('Select a genome: '),
                    dcc.Dropdown(labels, anchor_name, style=dict(width='40%', height='110%', verticalAlign="middle")),
                    html.Br(),
                    ],style={'padding':'1%', #'flex':1, #'padding-top' : '1%', 'padding-left' : '1%', 'padding-bottom' : '1%', 'padding-right' : '1%', 
                    'font-size':'24px',
                    'textAlign': 'left', "border":"2px grey solid"}),
                html.Div(className="w3-half", children=[
                    dcc.Graph(id="all_genomes",style={"font-size": 20, "height" : 1500},
                        figure = all_genomes_dend, config=config
                        ),
                    ]),
                html.Div(className="w3-half", children=[
                    dcc.Graph(id="all_genomes_kmer_comp",
                        figure = pangenome_comp, config=config, style={"font-size": 30, "height" : 1500}
                        )
                    ]),
                html.Div(className="w3-third", children=[
                    dcc.Graph(id="all_genome_sizes",
                        figure = pangenome_sizes, config=config, style={"font-size": 30, "height" : 500}
                        )
                    ]),
                html.Div(className="w3-third", children=[
                    dcc.Graph(id="average_kmer_content",
                        figure = pangenome_avg, config=config, style={"font-size": 30, "height" : 500}
                        )
                    ]),
                html.Div(className="w3-third", children=[
                    dcc.Graph(id="pangenome_num_seqs",
                        figure = pangenome_num_seqs, config=config, style={"font-size": 30, "height" : 500}
                        )
                    ]),
            ]),
        dcc.Tab(label='Anchor genome', value='homepage', style={"background-color": "lightgrey", "font-size": 36,}, 
            children=[
                html.Div(children=[
                    html.I("Anchor genome " + anchor_name),
                    html.Br(),
                    html.Label('Select a chromosome: '),
                    dcc.Dropdown(chrs_list, chrs, style=dict(width='40%', height='110%', verticalAlign="middle", ), id="Anchor_tab_dropdown"),
                    #html.Div(id="Anchor_tab_dropdown"),
                    html.Br(),
                    ], style={'padding-top' : '1%', 'padding-left' : '1%', 'padding-bottom' : '1%', 
                    'padding-right' : '1%', 'font-size':'24px',
                    'textAlign': 'left', "border":"2px grey solid", }),
                    #html.Div(className="w3-threequarter",
                    #children=[
                html.Div(className="w3-third",children=[
                    dcc.Graph(id="gene_content",
                        figure = gene_content_per_genome_fig, config=config), #gene_content_fig)
                    ],), #style={"height" : 700}),
                html.Div(className="w3-third",children=[
                    dcc.Graph(id="genes_per_chr",
                        figure = genes_per_chr_fig, config=config)
                    ], ),#style={"height" : 700}),
                html.Div(className="w3-third",children=[
                    dcc.Graph(id="avg_kmer_chr",
                        figure = avg_kmer_per_chr_fig, config=config)
                    ], ),#style={"height" : 700}),

                html.Div(className="w3-threequarter",children=[
                            dcc.Graph(id="all_chromosomes",
                                figure = whole_genome_fig, config=config, style={"height" : num_chrs*250})
                    ]),
                html.Div(className="w3-quarter",children=[
                            dcc.Graph(id="all_chromosomes_hists", style={"height" : num_chrs*250},
                                figure = whole_genome_hists, config=config
                                )
                            ]),
            ]),
            #]),
            #),
        dcc.Tab(label='Chromosome', value='chr_view', style={"background-color": "lightgrey", "font-size": 36,},
            children=[
                #html.Br(),
                html.Div(children=[ #id="Chrs_Info", children=[
                    #children=[
                    html.Div(children=[
                    html.I(anchor_name + "."), # + chrs + ":" ), #+ str(x_start_init) + "-" + str(x_stop_init)), #"Chromosome " + chrs,
                    #html.I()
                    html.I(chrs + ":", id="chr_name"),
                    dcc.Input(id="Chrs_Info", placeholder=str(x_start_init) + "-" + str(x_stop_init), ),
                    html.Br(),
                    html.Br(),
                    html.I("Genes in this region: " + str(len(gene_names_local)), id='regional_genes'),
                    html.Br(),
                    ],style={"display": "inline-block"}),
                    
                    html.Div(children=[
                    
                    html.I("K-mer length: " + str(kmer_len),style={"display": "inline-block",}),
                    html.Br(),
                    html.I("Number of bins: " + str(bins), style={'display': 'inline-block'}),
                    html.Br(),
                    html.I("Step size: " + str(n_skips_start), style={"display": "inline-block",}, id='step_size_out'),
                    html.Br(),
                    ],style={"display": "inline-block", 'padding-left' : '10%'}),
                    

                    ],
                    
                #],
                #value="Chromosome ", 
                #html.Div(id='chr_info', value=chrs, type='text')
                style={'padding-top' : '1%', 'padding-left' : '1%', 'padding-bottom' : '1%', 'padding-right' : '1%', 
                'font-size':'24px',
                #"display": "inline-block",
                            #'textAlign': 'left', 
                            "border":"2px grey solid"}
                            ),
            #]
            html.Div(children=[
                html.Div(className="w3-container", children=[
                        #left figure
                        dcc.Graph(id="chromosome",
                            figure = chr_fig, config=config, style={"font-size": 20, "height" : 350})
                ])
            ]),
            html.Div(children=[
                html.Div(className="w3-container", children=[
                    html.Div(className="w3-threequarter", children=[
                        #left figure - calling this the "Main" figure
                        dcc.Graph(id="primary",figure = fig, config=config, 
                            style={"height": 1000,  "font-size": 20})
                    ]), #style={'padding-top' : '1%', 'padding-left' : '1%', 'padding-bottom' : '1%', 'padding-right' : '1%', 
                    html.Div(className="w3-quarter", children=[
                        dcc.Graph(id="Secondary", 
                            #Now we have the phylogenetic tree
                            figure = get_local_info(names_simp[x_start_init:x_stop_init], #exon_comp[chrs], 
                                gene_comp[chrs], 
                                bar_sum_regional, bar_sum_global[chrs]), 
                            config=config,
                            style={"height": 1000, "font-size": 20}),
                    ])
                ])
            ]),
            html.Div(children=[
                html.Div(className="w3-container", children=[
                    html.Div(className="w3-third", children=[
                            dcc.Graph(id="Genes", 
                                figure=gene_content_plot, 
                                config=config,
                                style={"font-size": 20}),
                        ]),
                    html.Div(className="w3-twothird", children=[
                        dcc.Graph(id="Third", 
                            #This is the histogram section
                            
                            figure = create_tree(tree_file, x_start_init,x_stop_init, all_chrs[chrs][int(x_start_init/n_skips_start):int(x_stop_init/n_skips_start)], n_skips_start),
                            config=config,
                            style={"font-size": 20}),
                    ]),
                    
                    
                    ])
                ])
        ]),  
    #]),
    ]),

    html.Div("chr1", style={"display" : "none"}, id="selected-chrom-state"),
    #html.Div("chr1", style={"display" : "none"}, id="prev_selected-state"),
    
    html.Div(x_start_init,id='num-start',style={"display" : "none"} ),
    html.Div(x_stop_init,id='num-stop',style={"display" : "none"} ),
    #html.Div( [], style={"display" : "none"}, id="selected-chrom-data") #style={"display" : "none"}, 
])

def get_buffer(tmp_start, tmp_stop, n_skips):
    min_dist = n_skips*bins
    act_dist = tmp_stop-tmp_start
    if act_dist >= min_dist:
        return 0 #The distance between start and stop are good enough 
    else: #
        return 1 #int((min_dist - act_dist)/2)+1 
#app.config.requests_pathname_prefix = app.config.routes_pathname_prefix.split('/')[-1]
last_reset_click = 0

@app.callback(
    Output('chromosome','figure'),
    Output('primary','figure'),
    Output('Secondary', 'figure'),
    Output('Third', 'figure'),
    Output('Genes', 'figure'),
    Output('chromosome', 'relayoutData'),
    Output('selected-chrom-state', 'children'),
    Output('Chrs_Info', 'value'),
    Output('chr_name', 'children'),
    Output('regional_genes', 'children'),
    #Output('step_size_out', 'children'),

    Input('all_chromosomes','relayoutData'),
    Input('chromosome','figure'),
    Input('primary','figure'),
    Input('Secondary', 'figure'),
    Input('Third', 'figure'),
    Input('Genes', 'figure'),
    #Input('my_check','value'),
    Input('num-start','children'),
    Input('num-stop','children'),
    #Input('gene_jump', 'value'),
    #Input('bins', 'value'),
    #Input('textarea-state-button', 'n_clicks'),
    #Input('save-button', 'n_clicks'),
    #Input('SG_window','value'),
    #Input('SG_polynomial_order','value'),
    #Input('SG_check','value'),
    Input('primary', 'clickData'),
    Input('primary', 'relayoutData'),
    Input('chromosome', 'selectedData'),
    Input('Genes', 'clickData'),
    Input('Chrs_Info', 'value'),
    Input('Anchor_tab_dropdown', 'value'),

    #Input('Genes', 'clickData'),

    State('selected-chrom-state', 'children'),
    #Sent 11
    #State('prev_selected-state', 'children'),
    )
def update_figure(wgr, chr_fig, fig1, fig2, fig3, fig4, #clicker, 
    x_start, x_stop, #gene_jump, 
    #bins, #n_clicks, 
    #save_clicks, 
    #SG_window, SG_polynomial_order, SG_check, 
    clickData, relayoutData, chr_relayoutData, gene_jump_bottom, user_chr_coords, anchor_tab_dropdown, chrs, ):#clicker_update, gene_jump, bins):
    #Given 12 
    triggered_id = ctx.triggered_id
    print(triggered_id)
    #print(wgr)
    
    #print()
    n_skips = 100
    click_me_genes = True
    click_me_rep = True
    if triggered_id == 'Anchor_tab_dropdown':
        print(anchor_tab_dropdown)
        chrs = anchor_tab_dropdown
        chr_num = int(chrs.split('r')[1])
        return update_all_figs(chr_num, chr_relayoutData, click_me_rep, click_me_genes, chrs)
    elif triggered_id == 'Chrs_Info':
        #print(user_chr_coords)
        if user_chr_coords.count(' ') > 0:
            new_x_start = int(user_chr_coords.strip().split('-')[0]) 
            new_x_stop  = int(user_chr_coords.strip().split('-')[1])
            return user_chr_coords_triggered(new_x_start, new_x_stop, click_me_genes, click_me_rep, chr_relayoutData, 
                chrs, fig4)
            #print("woohoo!!!!!!!!")
        else:
            return no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update
    elif triggered_id == 'all_chromosomes':
        #Get the chromosome number 
        if wgr != None and len(wgr)>1:
            chr_num_tmp = list(wgr.keys())[0].split('.')[0].split('axis')[1] #if wgr != None and len(wgr)>1:
            if len(chr_num_tmp) == 0:
                chr_num = 1
            else:
                chr_num = int(int(chr_num_tmp)/3)+1
            #chr_fig = #This is the top pane, showing the whole chromosome. 
            #fig1 = #The main figure
            #fig2 = #The "local info" figure. Showing the k-mer composition in different regions 
            #fig3 = #Tree figure 
            #fig4 = #Sorted gene figure 
            return update_all_figs(chr_num, chr_relayoutData, click_me_rep, click_me_genes, #SG_window, SG_polynomial_order,
                #SG_check, 
                chrs) #chr_fig, fig1, fig2, fig3, fig4, chr_relayoutData
    elif triggered_id == 'Genes':
        return sorted_gene_fig(chr_fig, fig1, fig2, fig3, fig4, gene_jump_bottom, #SG_window, SG_polynomial_order,
            #SG_check, 
            n_skips, click_me_genes, click_me_rep, chr_relayoutData, chrs)
    elif triggered_id == 'chromosome':
        return chromosome_gene_triggered(chr_fig, fig1, fig2, fig3, fig4, gene_jump_bottom, #SG_window, SG_polynomial_order,
            #SG_check, 
            n_skips, click_me_genes, click_me_rep, chr_relayoutData, chrs)
    elif triggered_id == 'primary':
        return primary_fig_triggered(chr_fig, fig1, fig2, fig3, fig4, gene_jump_bottom, #SG_window, SG_polynomial_order,
            #SG_check, 
            n_skips, click_me_genes, click_me_rep, chr_relayoutData, relayoutData, clickData, x_start, x_stop, chrs)
    local_gene_list = []
    cntr = 0
    while cntr < len(gene_names_tmp):
        local_gene_list.append(gene_names_tmp[cntr])
        cntr += 3
    return chr_fig, fig1, fig2, fig3, fig4, chr_relayoutData, chrs, update_output_div(chrs,x_start,x_stop,anchor_name), update_out_chr(chrs) , update_gene_locals(local_gene_list) #, clickData

def user_chr_coords_triggered(new_x_start, new_x_stop, click_me_genes, click_me_rep, chr_relayoutData, chrs, fig4):

    exact_mini_counts = anchor.query(chrs, new_x_start, new_x_stop, n_skips_start)
    simple_cnts_for_plots = exact_mini_counts.sum(axis=1)
    chr_num = int(chrs.split('r')[1])

    fig1, bar_sum_names, bar_sum_regional, colors, gene_names_tmp = plot_interactive( n_skips_start, 
        layout, 
        gene_comp[chrs], bins, simple_cnts_for_plots, 
        chrs, zs_tmp, rep_types[chrs], click_me_rep, click_me_genes,  
        int(new_x_start), int(new_x_stop), gene_locals[chrs], gene_names[chrs], exon_locals[chrs], exon_names[chrs],
        )
    fig2 = get_local_info(
        names_simp[int(new_x_start/n_skips_start):int(new_x_stop/n_skips_start)], 
        gene_comp[chrs], bar_sum_regional, bar_sum_global[chrs])
    fig3 = create_tree(tree_file, new_x_start, new_x_stop, 
        exact_mini_counts, n_skips_start)
    chr_fig = plot_chr_whole(x[chr_num-1], z_1[chr_num-1], z_9[chr_num-1], z_genes[chrs], new_x_start, new_x_stop,  y[chr_num-1])
    local_gene_list = []
    cntr = 0
    while cntr < len(gene_names_tmp):
        local_gene_list.append(gene_names_tmp[cntr])
        cntr += 3
    fig4 = plot_gene_content(gene_content[chrs], sort_by, colors, uniq_avg, univ_avg, local_gene_list)
    return chr_fig, fig1, fig2, fig3, fig4, chr_relayoutData, chrs, update_output_div(chrs,new_x_start,new_x_stop,anchor_name), update_out_chr(chrs), update_gene_locals(local_gene_list)

def sorted_gene_fig(chr_fig, fig1, fig2, fig3, fig4, gene_jump_bottom, #SG_window, SG_polynomial_order,
            #SG_check, 
            n_skips, click_me_genes, click_me_rep, chr_relayoutData, chrs, ):
    print("****Sorted gene fig******")
    print(gene_jump_bottom)
    print(chrs)
    print("**********")
    names_simp = all_chrs[chrs].sum(axis=1)
    chr_num = int(chrs.split('r')[1])
    if gene_jump_bottom != None and gene_jump_bottom['points'][0]['text'] in gene_names[chrs]:
        print("first elif")
        #chr_num = int(chrs.split('r')[1])
        tmp_idx = gene_names[chrs].index(gene_jump_bottom['points'][0]['text'])
        if get_buffer(int(gene_locals[chrs][tmp_idx]), int(gene_locals[chrs][tmp_idx+1]), n_skips_start) == 1:
            #we need to move to single nucleotide resoltion 
            x_start = int(gene_locals[chrs][tmp_idx])-10
            x_stop = int(gene_locals[chrs][tmp_idx+1])+10
            exact_mini_counts = anchor.query(chrs, x_start, x_stop, 1)
            simple_cnts_for_plots = exact_mini_counts.sum(axis=1)
            #simple_cnts_for_plots, exact_mini_counts = read_mini_count_files(x_start, x_stop)
            n_skips = 1
        else: 
            x_start = int(gene_locals[chrs][tmp_idx]) - buff #int(1000*n_skips/2)
            x_stop = int(gene_locals[chrs][tmp_idx+1]) + buff #int(1000*n_skips/2)
            n_skips = n_skips_start
            simple_cnts_for_plots = names_simp[int(x_start/n_skips):int(x_stop/n_skips)]
            exact_mini_counts = all_chrs[chrs][int(x_start/n_skips):int(x_stop/n_skips)]
        fig1, bar_sum_names, bar_sum_regional, colors, gene_names_tmp = plot_interactive( n_skips, #int(SG_window), 
            #int(SG_polynomial_order), SG_check,
            layout, #exon_comp[chrs], 
            gene_comp[chrs], bins, simple_cnts_for_plots, chrs, zs_tmp, rep_types[chrs], click_me_rep, click_me_genes,  
            int(x_start), int(x_stop), gene_locals[chrs], gene_names[chrs], exon_locals[chrs], exon_names[chrs],
            )
        fig2 = get_local_info(
            simple_cnts_for_plots, #exon_comp[chrs], 
            gene_comp[chrs], bar_sum_regional, bar_sum_global[chrs])
        fig3 = create_tree(tree_file, x_start, x_stop, exact_mini_counts, n_skips)
        gene_jump = ""
        if x_start_init != x_start:
            chr_fig = plot_chr_whole(x[chr_num-1], z_1[chr_num-1], z_9[chr_num-1], z_genes[chrs], x_start, x_stop,  y[chr_num-1])
    local_gene_list = []
    cntr = 0
    while cntr < len(gene_names_tmp):
        local_gene_list.append(gene_names_tmp[cntr])
        cntr += 3
    fig4 = plot_gene_content(gene_content[chrs], sort_by, colors, uniq_avg, univ_avg, local_gene_list)
    return chr_fig, fig1, fig2, fig3, fig4, chr_relayoutData, chrs, update_output_div(chrs,x_start,x_stop,anchor_name), update_out_chr(chrs), update_gene_locals(local_gene_list)
def chromosome_gene_triggered(chr_fig, fig1, fig2, fig3, fig4, gene_jump_bottom, #SG_window, SG_polynomial_order,
    #SG_check, 
    n_skips, click_me_genes, click_me_rep, chr_relayoutData, chrs, ):
    print(chrs)
    chr_num = int(chrs.split('r')[1])
    names_simp = all_chrs[chrs].sum(axis=1)
    if "x2" in chr_relayoutData['range'].keys():
            x_start = int(chr_relayoutData['range']['x2'][0])
            x_stop = int(chr_relayoutData['range']['x2'][1])
    elif "x" in chr_relayoutData['range'].keys():
            x_start = int(chr_relayoutData['range']['x'][0])
            x_stop = int(chr_relayoutData['range']['x'][1])
    else:
            x_start = int(chr_relayoutData['range']['x3'][0])
            x_stop = int(chr_relayoutData['range']['x3'][1])
    if get_buffer(x_start, x_stop, n_skips_start) == 1:
        #we need to move to single nucleotide resoltion
        exact_mini_counts = anchor.query(chrs, x_start, x_stop, 1)
        simple_cnts_for_plots = exact_mini_counts.sum(axis=1)
        n_skips = 1
    else:
        exact_mini_counts = all_chrs[chrs][int(x_start/n_skips_start):int(x_stop/n_skips_start)]
        simple_cnts_for_plots = names_simp[int(x_start/n_skips_start):int(x_stop/n_skips_start)]
        n_skips = n_skips_start
    fig1, bar_sum_names, bar_sum_regional, colors, gene_names_tmp = plot_interactive( n_skips, #int(SG_window), 
            #int(SG_polynomial_order), SG_check,
            layout, #exon_comp[chrs], 
            gene_comp[chrs], bins, simple_cnts_for_plots, chrs, zs_tmp, rep_types[chrs], click_me_rep, click_me_genes,  
            int(x_start), int(x_stop), gene_locals[chrs], gene_names[chrs], exon_locals[chrs], exon_names[chrs],
            )
    fig2 = get_local_info(
            simple_cnts_for_plots, #exon_comp[chrs], 
            gene_comp[chrs], bar_sum_regional, bar_sum_global[chrs])
    #And now we update the histograms 
    fig3 = create_tree(tree_file, x_start, x_stop, exact_mini_counts, n_skips)
    chr_relayoutData = None#{'autosize': True}
    if x_start_init != x_start:
        chr_fig = plot_chr_whole(x[chr_num-1], z_1[chr_num-1], z_9[chr_num-1], z_genes[chrs], x_start, x_stop,  y[chr_num-1])
    local_gene_list = []
    cntr = 0
    while cntr < len(gene_names_tmp):
        local_gene_list.append(gene_names_tmp[cntr])
        cntr += 3
    fig4 = plot_gene_content(gene_content[chrs], sort_by, colors, uniq_avg, univ_avg, local_gene_list)
    return chr_fig, fig1, fig2, fig3, fig4, chr_relayoutData, chrs, update_output_div(chrs,x_start,x_stop,anchor_name), update_out_chr(chrs), update_gene_locals(local_gene_list)
def primary_fig_triggered(chr_fig, fig1, fig2, fig3, fig4, gene_jump_bottom, #SG_window, SG_polynomial_order,
    #SG_check, 
    n_skips, click_me_genes, click_me_rep, chr_relayoutData, relayoutData, clickData, x_start, x_stop, chrs, ):
    print(chrs)
    local_gene_list = gene_names_local
    names_simp = all_chrs[chrs].sum(axis=1)
    chr_num = int(chrs.split('r')[1])
    if relayoutData != None and 'xaxis4.range[0]' in relayoutData.keys():
        print("fourth elif")
        x_start = int(relayoutData['xaxis4.range[0]'])
        x_stop = int(relayoutData['xaxis4.range[1]'])
        if get_buffer(x_start, x_stop, n_skips_start) == 1:
            exact_mini_counts = anchor.query(chrs, x_start, x_stop, 1)
            simple_cnts_for_plots = exact_mini_counts.sum(axis=1)
            n_skips = 1
        else:
            exact_mini_counts = all_chrs[chrs][int(x_start/n_skips_start):int((x_stop)/n_skips_start)]
            simple_cnts_for_plots = names_simp[int(x_start/n_skips_start):int((x_stop)/n_skips_start)]
            n_skips = n_skips_start 
        #fig1, fig2, fig3 = update_all_figs(x_start, x_stop, click_me_rep, click_me_genes)
        fig1, bar_sum_names, bar_sum_regional, colors, gene_names_tmp = plot_interactive( n_skips_start, #int(SG_window), 
            #int(SG_polynomial_order), SG_check,
            layout, #exon_comp[chrs], 
            gene_comp[chrs], bins, 
            simple_cnts_for_plots,
            #names_simp[int(x_start/n_skips_start):int(x_stop/n_skips_start)], 
            chrs, zs_tmp, rep_types[chrs], click_me_rep, click_me_genes,  
            int(x_start), int(x_stop), gene_locals[chrs], gene_names[chrs], exon_locals[chrs], exon_names[chrs],
            )
        fig2 = get_local_info(
            #full_cnts[x_start:x_stop],
            simple_cnts_for_plots, 
            #names_simp[int(x_start/n_skips_start):int(x_stop/n_skips_start)], #exon_comp[chrs], 
            gene_comp[chrs], bar_sum_regional, bar_sum_global[chrs])
        #And now we update the histograms 
        fig3 = create_tree(tree_file, x_start, x_stop, 
            exact_mini_counts,
            #all_chrs[chrs][int(x_start/n_skips_start):int(x_stop/n_skips_start)], 
            n_skips)
        chr_fig = plot_chr_whole(x[chr_num-1], z_1[chr_num-1], z_9[chr_num-1], z_genes[chrs], x_start, x_stop,  y[chr_num-1])
        local_gene_list = []
        cntr = 0
        while cntr < len(gene_names_tmp):
            local_gene_list.append(gene_names_tmp[cntr])
            cntr += 3
        fig4 = plot_gene_content(gene_content[chrs], sort_by, colors, uniq_avg, univ_avg, local_gene_list)
    elif clickData != None:#len(print(clickData['points'])) > 0:
        print("third elif")
        #print(clickData)
        tmp_idx = gene_names[chrs].index(clickData['points'][0]['text'])
        #gene_buffer = get_buffer(int(gene_locals[tmp_idx]), int(gene_locals[tmp_idx+1]), n_skips_start)
        x_start = int(gene_locals[chrs][tmp_idx])-10
        x_stop = int(gene_locals[chrs][tmp_idx+1])+10
        if get_buffer(x_start, x_stop, n_skips_start) == 1:
            #we need to move to single nucleotide resoltion
            exact_mini_counts = anchor.query(chrs, x_start, x_stop, 1)
            simple_cnts_for_plots = exact_mini_counts.sum(axis=1)
            n_skips = 1
        else:
            exact_mini_counts = all_chrs[chrs][int(x_start/n_skips_start):int((x_stop)/n_skips_start)]
            simple_cnts_for_plots = names_simp[int(x_start/n_skips_start):int((x_stop)/n_skips_start)]
            n_skips = n_skips_start 
        #fig1, fig2, fig3 = update_all_figs(x_start, x_stop, click_me_rep, click_me_genes)
        fig1, bar_sum_names, bar_sum_regional, colors, gene_names_tmp = plot_interactive( n_skips, #int(SG_window), 
            #int(SG_polynomial_order), SG_check,
            layout, #exon_comp[chrs], 
            gene_comp[chrs], bins, simple_cnts_for_plots, chrs, zs_tmp, rep_types[chrs], 
            click_me_rep, click_me_genes,  
            int(x_start), int(x_stop), gene_locals[chrs], gene_names[chrs], exon_locals[chrs], exon_names[chrs],
            )
        fig2 = get_local_info(
            #full_cnts[x_start:x_stop], 
            simple_cnts_for_plots, 
            #exon_comp[chrs], 
            gene_comp[chrs], bar_sum_regional, bar_sum_global[chrs])
        #And now we update the histograms 
        fig3 = create_tree(tree_file, x_start, x_stop, 
            exact_mini_counts, n_skips)
        chr_fig = plot_chr_whole(x[chr_num-1], z_1[chr_num-1], z_9[chr_num-1], z_genes[chrs], x_start, x_stop,  y[chr_num-1])
        local_gene_list = []
        cntr = 0
        while cntr < len(gene_names_tmp):
            local_gene_list.append(gene_names_tmp[cntr])
            cntr += 3 

        fig4 = plot_gene_content(gene_content[chrs], sort_by, colors, uniq_avg, univ_avg, local_gene_list)   
    return chr_fig, fig1, fig2, fig3, fig4, chr_relayoutData, chrs, update_output_div(chrs,x_start,x_stop,anchor_name), update_out_chr(chrs),update_gene_locals(local_gene_list)
def render_tab_content(active_tab, chrs):
    return chr_fig, fig1, fig2, fig3, fig4, chr_relayoutData, chrs, update_output_div(chrs,x_start_init,x_stop_init,anchor_name), update_out_chr(chrs), update_gene_locals(local_gene_list)
def update_all_figs(chr_num, chr_relayoutData, click_me_rep, click_me_genes, #SG_window, SG_polynomial_order,
    #SG_check, 
    chrs, ):
    #chr_fig = This is the top pane, showing the whole chromosome. 
    #fig1    = The main figure
    #fig2    = The "local info" figure. Showing the k-mer composition in different regions 
    #fig3    = Tree figure 
    #fig4    = Sorted gene figure
    chrs = "chr" + str(chr_num)
    names_simp = all_chrs[chrs].sum(axis=1) #all_chrs["chr"+str(chr_num)].sum(axis=1)#tmp.sum(axis=1)
    #print(names_simp[0:100])
    #chrs = "chr" + str(chr_num)
    print(len(names_simp))
    print(x_start_init_adjust)
    print(x_stop_init_adjust)
    print(names_simp[x_start_init_adjust:x_stop_init_adjust])
    fig1, bar_sum_names, bar_sum_regional, colors, gene_names_tmp = plot_interactive( n_skips_start, #int(SG_window), 
            #int(SG_polynomial_order), SG_check,
            layout, #exon_comp[chrs], 
            gene_comp[chrs], bins, names_simp[x_start_init_adjust:x_stop_init_adjust], 
            chrs, zs_tmp, rep_types[chrs], click_me_rep, click_me_genes,  
            int(x_start_init), int(x_stop_init), gene_locals[chrs], gene_names[chrs], exon_locals[chrs], 
            exon_names[chrs])
    local_gene_list = []
    cntr = 0
    while cntr < len(gene_names_tmp):
        local_gene_list.append(gene_names_tmp[cntr])
        cntr += 3
    sort_by = ["Unique","Universal"]
    tmp = [(i/sum(bar_sum_global['chr1'])*100) for i in bar_sum_global['chr1']]
    uniq_avg = tmp[1]
    univ_avg = tmp[-1]
    fig4 = plot_gene_content(gene_content[chrs], sort_by, colors, uniq_avg, univ_avg, local_gene_list)
    
    #x, z_1, z_9, y  = read_chr_whole()
    #z_genes = make_gene_whole_chr(x[chr_num-1], gene_locals)
    chr_fig = plot_chr_whole(x[chr_num-1], z_1[chr_num-1], z_9[chr_num-1], z_genes[chrs], x_start_init, 
        x_stop_init, y[chr_num-1])

    fig2 = get_local_info(
            names_simp[x_start_init_adjust:x_stop_init_adjust], #exon_comp[chrs], 
            gene_comp[chrs], 
            bar_sum_regional, bar_sum_global[chrs])
    
    fig3 = create_tree(tree_file, x_start_init_adjust, x_stop_init_adjust, 
        all_chrs[chrs][x_start_init_adjust:x_stop_init_adjust], n_skips_start)
    
    #chr_fig = plot_chr_whole(x[chr_num-1], z_1[chr_num-1], z_9[chr_num-1], z_genes[chr_num-1], x_start_init, x_stop_init, y[chr_num-1])
    return chr_fig, fig1, fig2, fig3, fig4, chr_relayoutData, chrs, update_output_div(chrs,x_start_init,x_stop_init, anchor_name), update_out_chr(chrs), update_gene_locals(local_gene_list)
def update_output_div(chrs, start, stop, anchor_name):
    return f'{start}-{stop}'
def update_out_chr(chrs):
    return f'{chrs}:'
def update_gene_locals(local_gene_list):
    printme = ""
    if len(local_gene_list)==1:
        printme += "Genes: "
        printme += local_gene_list[0] + ": "
        printme += gene_anns[local_gene_list[0]]
    elif len(local_gene_list)<=10:
        printme += "Genes: "
        for i in local_gene_list:
            printme += i + ", "
    else:
        printme = "Genes in this region: " + str(len(local_gene_list))
    return f'{printme}'

#@app.callback(Output('textarea-state-button', 'n_clicks'),
#    Output('num-start','value'),
#    Output('num-stop','value'),
#    [Input('reset_clicks','n_clicks')])
#def update(reset):
#    return 0, x_start_init, x_stop_init

app.run_server(host='127.0.0.1', debug=True)





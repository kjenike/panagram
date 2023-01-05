#from bokeh.plotting import figure, output_file, save
import sys
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
import time
from Bio import Phylo
from scipy import signal
from dash import Dash, dcc, html, Input, Output, ctx, State
import math
from io import StringIO
from scipy.cluster.hierarchy import linkage
from scipy.cluster import hierarchy
#from cStringIO import StringIO
#plt.rcParams["figure.figsize"] = (40,40)
num_samples = 9
num_chrs = 10
buff = 1000
mini_file_size = 1000000
bit_file_prefix = "test_data/b73.anchor.ids"
x_start_init = 75000000#0
x_stop_init = 85000000 #len(zs_tmp)
bins=400
sns.set_palette('viridis', num_samples)
#plt.rcParams.update({'font.size': 22})

#/Users/katiejenike/Desktop/Lab_notes/PanSol/PANGEN/DEBUG/QC_PLOTS/SQUI2_PANGEN_PLOTS
#python plot_chr_interactive.py result_SQUI2.chr1.txt SQUI2
#python plot_pangene.py test_data/result_SQUI2.chr1.txt SQUI2
labels = ["B73","P39","M162W","Il14H","B97",
    "HP301","CML69","CML228","M37W"]
ref_genome = "Solqui2"
n_skips_start = 100
in_file=sys.argv[1]
sample = sys.argv[2]
count_file = "test_data/SQUI2.chr1.txt"
rep_file = "test_data/Zm-B73-REFERENCE-NAM-5.0.TE.gff3"#"test_data/Solqui2_v2.0.fasta.mod.EDTA.TEanno.gff3"#"test_data/Solqui2.repeats.chr1.txt"
gene_file = "test_data/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3"#"test_data/Solqui2.liftoff_gene_models.gff3"#"test_data/Squi2.chr1.exons.genes.cds.txt"
tree_file = "test_data/all.alns.aug5.2022.fa.treefile"
bins_file = "test_data/b73.chr1.txt"#"test_data/squi_chr1_bins.txt"
rep_types_file = "test_data/maize_te.txt"

chr_lens = [308452471, 243675191, 238017767, 250330460, 226353449, 
181357234, 185808916, 182411202, 163004744, 152435371]

sample_labels_1s   =  {
    "0": [""],
    "1": ["Solqui2"],
    "3": ["Solabu2"],
    "5": ["Solcan1"],
    "6": ["Solqui2","Solcan1"],
    "4": ["Solqui2","Solabu2"],
    "8": ["Solcan1","Solabu2"],
    "9": ["Solqui2","Solcan1","Solabu2"]
}
sample_labels_10s  =  {
    "0": [],
    "1": ["Solaet3"],
    "3": ["Solyc"],
    "5": ["Solmac3"],
    "6": ["Solaet3", "Solmac3"],
    "4": ["Solaet3", "Solyc"],
    "8": ["Solmac3","Solyc"],
    "9": ["Solaet3","Solyc","Solmac3"]
}
sample_labels_100s =  {
    "0": [],
    "1": ["Solmur2hap1"],
    "3": ["Solmur2hap2"],
    "5": ["Solpri1"],
    "6": ["Solpri1","Solmur2hap1"],
    "4": ["Solmur2hap1","Solmur2hap2"],
    "8": ["Solmur2hap2","Solpri1"],
    "9": ["Solpri1","Solmur2hap1","Solmur2hap2"]
}

name_idx_short = {"0": 0, "1": 1, "3": 1, "5": 1, "4": 2, "6": 2, "8":2, "9":3}

name_idx_long = {}
name_idx_complex = {}
for i in ["0","1","3","4","5","6","8","9"]:
    name_idx_long[i] = name_idx_short[i]
    name_idx_complex[i] = sample_labels_1s[i]
    for j in ["0","1","3","4","5","6","8","9"]:
        tmp_str = i+j
        tmp_total = name_idx_short[i] + name_idx_short[j]
        name_idx_long[tmp_str] = tmp_total
        tmp_lst = sample_labels_10s[i] + sample_labels_1s[j]
        name_idx_complex[tmp_str] = tmp_lst
        for k in ["0","1","3","4","5","6","8","9"]:
            tmp_str = i+j+k
            tmp_total = name_idx_short[i] + name_idx_short[j] + name_idx_short[k]
            name_idx_long[tmp_str] = tmp_total
            tmp_lst = sample_labels_100s[i] + sample_labels_10s[j] + sample_labels_1s[k]
            name_idx_complex[tmp_str] = tmp_lst

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

        ret = np.zeros((len(pac), self.ngenomes), dtype=bool)
        for i in range(self.ngenomes):
            ret[:,i] = (pac >> i) & 1

        if step != 1:
            ret = ret[np.arange(0, len(pac), step)]

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
    with open(rep_types_f, "r") as f:
        line = f.readline()
        while line:
            rep_types[line.strip()] = []
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
sample_counts = {
        1: 1,
        3: 1, 
        5: 1,
        
        4: 2,
        6: 2,
        8: 2,

        9: 3
        }

sample_colors = {
        1: "#febd2a",
        3: "#f48849",
        5: "#db5c68",

        4: "#b83289",
        6: "#8b0aa5",
        8: "#5302a3",

        9: "#0d0887"
        }

sample_labels = {
        1: "Solqui2",
        3: "Solabu2",
        5: "Solcan1",

        10: "Solaet3",
        30: "Solyc",
        50: "Solmac3",

        100: "Solmur2hap1",
        300: "Solmur2hap2",
        500: "Solpri1"
        }
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
    simple_exact_mini_counts = bitvec_to_mat(np.array(exact_mini_counts), 9).sum(axis=1) #parse_counts_simple(exact_mini_counts)

    return simple_exact_mini_counts, exact_mini_counts

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

def read_annotations(ann_file, ann_types, chrs):
    with open(ann_file, "r") as f:
        line = f.readline()
        while line:
            tmp = line.split('\t')
            this_type = tmp[2].strip()

            if this_type in rep_types and tmp[0]==chrs:
                #print(line.strip())
                ann_types[this_type].append(int(tmp[3]))
                ann_types[this_type].append(int(tmp[4]))
                ann_types[this_type].append('None')
            line = f.readline()
    return ann_types

def perc_universal(start,stop, countme, n_skips):
    univ = ((anchor.get_counts(chrs, start, stop, n_skips).sum(axis=1))==countme).sum()
    #float((names_simp[int(start/n_skips):int(stop/n_skips)]==countme).sum())
    #univ = float(names_simp[int(start/n_skips):int(stop/n_skips)].count(countme))
    return float(univ/(int(stop/n_skips)-int(start/n_skips)))*100
    #return float(univ/(((stop)-(start/n_skips))))*100

def read_gene_annotations(ann_file, gene_names, gene_locals, exon_locals, exon_names, cds_locals, n_skips, chrs):
    #rep_file = "Solqui2.repeats.chr1.100k.txt"
    gene_content = {}
    for i in range(1, num_chrs+1):
        gene_content["chr"+str(i)] = {"Names":[],"Universal":[], "Unique":[]}
    #gene_content = {"Names":[],"Universal":[], "Unique":[]}
    with open(ann_file, "r") as f:
        line = f.readline()
        while line:
            #print(line.strip())
            tmp = line.strip().split('\t')
            this_chr = tmp[0]
            #print(this_chr)
            if this_chr in gene_locals.keys():
                this_name = tmp[8]#line.split('\t')[1].strip()
                this_type = tmp[2]
                #if tmp[0] == chrs:
                #print(this_type)
                if this_type == "gene": #this_name.split(":")[0] == "ID=gene":
                    gene_locals[this_chr].append(int(tmp[3]))
                    gene_locals[this_chr].append(int(tmp[4]))
                    gene_locals[this_chr].append('None')
                    tmp_name = this_name.split(';')[0].split("=")[1]
                    gene_names[this_chr].append(tmp_name)
                    gene_names[this_chr].append(tmp_name)
                    gene_names[this_chr].append('None')
                    gene_content[this_chr]["Names"].append(tmp_name)#perc_universal(int(tmp[2]), int(tmp[3]))
                    gene_content[this_chr]["Universal"].append(perc_universal(int(tmp[3]), int(tmp[4]), num_samples, 1))
                    gene_content[this_chr]["Unique"].append(perc_universal(int(tmp[3]), int(tmp[4]), 1, 1))
                elif this_type == "exon": #this_name.split(":")[0] == "ID=exon":
                    exon_locals[this_chr].append(int(tmp[3]))
                    exon_locals[this_chr].append(int(tmp[4]))
                    exon_locals[this_chr].append('None')
                    tmp_name = this_name.split(';')[0].split("=")[1]
                    exon_names[this_chr].append(tmp_name)
                    exon_names[this_chr].append(tmp_name)
                    exon_names[this_chr].append('None')
                #elif this_name.split(":")[0] == "ID=CDS":
                #    cds_locals.append(int(tmp[2]))
                #    cds_locals.append(int(tmp[3]))
                #    cds_locals.append('None')
            #print()
            line = f.readline()
    print(gene_locals["chr1"][-100:])
    print(len(gene_locals["chr1"]))
    print("*********")
    return gene_locals, gene_names, exon_locals, exon_names, cds_locals, gene_content

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
        if tmp_name.count("_") > 0:
            tmp_name = tmp_name.split("_")[1]
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
        if tmp_name.count("_") > 0:
            tmp_name = tmp_name.split("_")[1]
        #print(tmp_name)
        line_color = color_code[tmp_name]#palette[int(((kmer_num[tmp_name])/total_kmers)*100)-1]
    elif clade.clades:
        
        line_color = palette[int(biggest_num_in_clade(clade, kmer_num))+25]
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
            starting_idx = [0]*9
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
    total_kmers = kmer_num_tmp[0] #kmer_num_tmp["Solqui2"]
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
        color_code[labels[k]] = palette[int(((kmer_num_tmp[k])/total_kmers)*100)+25]
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
        if tmp_name.count("_") > 0:
            tmp_name = tmp_name.split("_")[1]
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
        tmp_txt = ""
        if text.count("_") > 0:
            tmp_txt = text.split("_")[1]
        else:
            tmp_txt = text
        #print(kmer_num_tmp[tmp_txt])
        tmp_txt += " - " + str(kmer_num[tmp_txt])[:4] + "% (" + str(kmer_num_raw[tmp_txt]) + ")" 
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
                  font=dict(family='Balto', size=16),
                  # width=1000,
                  height=500,
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

def get_local_info(x_all, exon_comp, gene_comp, bar_sum_regional, bar_sum_global):
    fig = make_subplots(
        rows=2, cols=2,
        #shared_yaxes=True,
        specs=[[{"type": "bar"}, {"type": "bar"}], 
            [{"type": "bar"}, {"type": "bar"} ]],
        subplot_titles=("Whole chromosome", "This region", "Exons", 
            "Genes", ), 
        vertical_spacing=0.1
    )
    x = []
    for i in range(1,num_samples+1):
        x.append(i)
    #x = [1,2,3,4,5,6,7,8,9]
    colors = ["#fde725", "#addc30", "#5ec962", "#28ae80", "#21918c", "#2c728e", "#3b528b", "#472d7b", "#440154"]
    #This region
    #print(bar_sum_global)
    y=[(i/sum(bar_sum_regional[1:])*100) for i in bar_sum_regional[1:]]
    y_whole=[(i/sum(bar_sum_global)*100) for i in bar_sum_global]
    #print("vvvvvv")
    #print(y)
    #print(y_whole)
    #print("^^^^^^^^")
    #fig.add_trace(go.Bar(x=x, y=y_whole, marker_color="#7400b8", showlegend=False), row=1, col=1)
    fig.add_trace(go.Bar(x=x, y=[a_i - b_i for a_i, b_i in zip(y, y_whole)], marker_color=colors, showlegend=False), row=1, col=2)
    #fig.update_layout(xaxis_title_text="K-mers shared in X samples", yaxis_title_text='Frequency (log)', barmode='overlay')

    #CDS
    #y=[(i/sum(cds_comp[1:])*100) for i in cds_comp[1:]]
    #fig.add_trace(go.Bar(x=x, y=y_whole, marker_color="#7400b8", showlegend=False), row=1, col=2)
    fig.add_trace(go.Bar(x=x, y=y_whole, marker_color=colors, showlegend=False), row=1, col=1)
    #fig.update_layout(xaxis_title_text="K-mers shared in X samples", yaxis_title_text='Frequency (log)')

    #Exons
    y=[(i/sum(exon_comp[1:])*100) for i in exon_comp[1:]]
    #fig.add_trace(go.Bar(x=x, y=y_whole, marker_color="#7400b8", showlegend=False), row=1, col=3)
    fig.add_trace(go.Bar(x=x, y=[a_i - b_i for a_i, b_i in zip(y, y_whole)], marker_color=colors, showlegend=False), row=2, col=1)
    #fig.update_layout(xaxis_title_text="K-mers shared in X samples", yaxis_title_text='Frequency (log)')
    #fig.update_layout(xaxis_title_text="K-mers shared in X samples", yaxis_title_text='Frequency (log)', barmode='overlay')

    #Genes
    y=[(i/sum(gene_comp[1:])*100) for i in gene_comp[1:]]
    #fig.add_trace(go.Bar(x=x, y=y_whole, marker_color="#7400b8", showlegend=False), row=1, col=4)
    fig.add_trace(go.Bar(x=x, y=[a_i - b_i for a_i, b_i in zip(y, y_whole)], marker_color=colors, showlegend=False), row=2, col=2)
    #fig.update_layout(xaxis_title_text="K-mers shared in X samples", yaxis_title_text='Frequency (log)')
    fig.update_xaxes(title_text="# of genomes", row=2, col=1)
    fig.update_yaxes(title_text="Difference from whole chromosome", row=2, col=1)
    fig.update_yaxes(title_text="Perent of k-mers", row=1, col=1)
    fig.update_layout(height=1000)
    return fig

def plot_interactive( n_skips, window_filter, poly_order, shared_kmers, layout, exon_comp, gene_comp, bins, names_simp, name, zs_tmp, 
    rep_types, plot_rep, plot_gene, x_start, x_stop, gene_locals, gene_names):

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
    colors = ["grey","#fde725", "#addc30", "#5ec962", "#28ae80", "#21918c", "#2c728e", "#3b528b", "#472d7b", "#440154"]

    #We are adjusting the start and stop positions to account for the skipping. 
    #The adjusted value should be the index, whereas the x_start and x_stop are the real coordinates 
    adjusted_x_start = int(x_start/n_skips)
    adjusted_x_stop = int(x_stop/n_skips)

    #Get the bins
    bin_size = ((x_stop-x_start)/bins)
    adjusted_bin_size = (bin_size/n_skips)
    #print()
    #print(bin_size)
    #print()
    #Setting up our x axis scale 
    #cntr = x_start
    #x = []
    #while cntr < x_stop:
    #    x.append(cntr)
    #    cntr += bin_size
    cntr = 0
    print("a1:", time.time()-t)
    t = time.time()
    
    #
    cats_tmp = [([0] * (bins+1)) for _ in range(10)]
    cntr = 0
    x_tracker = 0
    
    #Here we are filling in the bins for the main figure.    
    bin_size_int = int(bin_size) + 1

    adjusted_bin_size_init = int(adjusted_bin_size) + 1
    #print('**********')
    #print(adjusted_bin_size_init)
    #print(bin_size_int)

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
        if x_tracker < len(cats_tmp[i]):
            cats_tmp[i][x_tracker] += 1 
        cntr += 1 #n_skips

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
        while i < len(exon_locals[chrs]):
            if (x_start < int(exon_locals[chrs][i]) < x_stop) and (x_start < int(exon_locals[chrs][i+1]) < x_stop):
                exon_locals_tmp.append(exon_locals[chrs][i])
                exon_locals_tmp.append(exon_locals[chrs][i+1])
                exon_locals_tmp.append(exon_locals[chrs][i+2])

                exon_names_tmp.append(exon_names[chrs][i])
                exon_names_tmp.append(exon_names[chrs][i+1])
                exon_names_tmp.append(exon_names[chrs][i+2])
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
    for i in range(0, len(cats_tmp)):
        bar_sum_regional.append(sum(cats_tmp[i]))
        #bar_sum_global.append(names_simp.count(i))
        bar_sum_names.append("Present in "+str(i) + " samples")
        fig.add_trace(go.Bar(x=x, y=cats_tmp[i], name="Present in "+str(i) + " samples",
            legendgroup="group1", 
            legendgrouptitle_text="Conserved K-mers",
            marker=dict(color=colors[i]), 
            marker_line=dict(color=colors[i])
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
    fig.update_layout(height=1000, xaxis_range=[x_start,x_stop])
    print("a3:", time.time()-t)
    t = time.time()
    return fig, bar_sum_names, bar_sum_regional, colors

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
    print(locs[-100:])
    print(len(locs))
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

    chr_fig.add_trace(go.Heatmap(x=x, z=z_9, y=y, type = 'heatmap', colorscale='plasma', showlegend=False, showscale=False), row=1, col=1)
    chr_fig.add_trace(go.Scatter(x=[x_start, x_start, None, x_stop, x_stop, None, x_start, x_stop, ], showlegend=False,
                   y=[0.5, 1.5, None, 0.5, 1.5, None, 1.5, 1.5 ],
                   mode='lines',
                   line_color='#1dd3b0', line_width=4), row=1, col=1)
    ###
    chr_fig.add_trace(go.Heatmap(x=x, z=z_1, y=y, type = 'heatmap', colorscale='plasma', showscale=False), row=2, col=1)
    chr_fig.add_trace(go.Scatter(x=[x_start, x_start, None, x_stop, x_stop], showlegend=False,
                   y=[0.5, 1.5, None, 0.5, 1.5],
                   mode='lines',
                   line_color='#1dd3b0', line_width=4), row=2, col=1)
    ###
    chr_fig.add_trace(go.Heatmap(x=x, z=z_genes, y=y, type = 'heatmap', colorscale='plasma', showscale=False ), row=3, col=1)
    chr_fig.add_trace(go.Scatter(x=[x_start, x_start, None, x_stop, x_stop], showlegend=False, 
                   y=[0.5, 1.5, None, 0.5, 1.5],
                   mode='lines',
                   line_color='#1dd3b0', line_width=4), row=3, col=1)
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
    chr_fig.update_yaxes(title_text="Universal", row=1, col=1)
    chr_fig.update_yaxes(title_text="Unique", row=2, col=1)
    chr_fig.update_yaxes(title_text="Genes", row=3, col=1)
    #chr_fig.update_xaxes(fixedrange=True, row=4, col=1)
    chr_fig.update_layout(clickmode='event+select', dragmode="select", selectdirection='h')
    #chr_fig.update_traces(showscale=False)
    chr_fig.update_layout(height=350)
    chr_fig.update_layout(margin=dict(
            b=10,
            l=10,
            r=10))
    
    return chr_fig

def plot_whole_genome(x, z_1, z_9, y):
    spec = []
    sub_titles = []
    h = num_chrs*350
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
        wg_fig.add_trace(go.Heatmap(x=x[i-1], z=z_9[i-1], y=y[i-1], type = 'heatmap', colorscale='plasma', showlegend=False, showscale=False), row=((i*3)-2), col=1)
        wg_fig.add_trace(go.Heatmap(x=x[i-1], z=z_1[i-1], y=y[i-1], type = 'heatmap', colorscale='plasma', showscale=False), row=((i*3)-1), col=1)
        
        wg_fig.update_yaxes( range=[0.5,1.5], showticklabels=False, row=((i*3)-2), col=1)
        wg_fig.update_yaxes( range=[0.5,1.5], showticklabels=False, row=((i*3)-1), col=1)
        
        wg_fig.update_xaxes(fixedrange=True, row=((i*3)-2), col=1)
        wg_fig.update_xaxes(fixedrange=True, row=((i*3)-1), col=1)
        
        wg_fig.update_xaxes(title_text="Sequence position", row=((i*3)), col=1)
        wg_fig.update_yaxes(title_text="Universal", row=((i*3)-2), col=1)
        wg_fig.update_yaxes(title_text="Unique", row=((i*3)-1), col=1)
    
    wg_fig.update_layout(clickmode='event')
    wg_fig.update_layout(height=h)
    wg_fig.update_layout(margin=dict(
            b=10,
            l=10,
            r=10))
    return wg_fig

def plot_gene_content(gene_content, sort_by, colors):
    #print(gene_content)
    colors = ['#fde725', '#440154']
    df = pd.DataFrame(gene_content)
    x = []
    cntr = 0
    for i in range(0, len(df['Universal'])):
        x.append(cntr)
        cntr +=1
    cntr = 0
    df_sorted = df.sort_values(sort_by[-1])
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
    fig.update_layout(height=500)
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
            r=10)
    )
    return fig

def update_output_div(input_value):
    return f'Output: {input_value}'


##### READ DATA
t = time.time()
cnts_tmp = []
x = []
all_chrs = {} #This is where we keep the raw, str, counts 
rep_types = get_init_rep_types()
anchor = KmerRef(bit_file_prefix)

for i in range(1, num_chrs+1):
    chrs="chr" + str(i)
    chr_start = 0
    chr_end = chr_lens[i-1] #308452471
    all_chrs[chrs] = anchor.get_counts(chrs, chr_start, chr_end, n_skips_start)

print("a0:", time.time()-t)
t = time.time()
#chrs = "chr1"
#Parsing the counts will return the number of genomes present at each kmer position, based on the counts
#names_simp = parse_counts_simple(all_chrs[chrs])
#tmp = bitvec_to_mat(np.array(all_chrs[chrs]), 9)
names_simp = all_chrs[chrs].sum(axis=1)#tmp.sum(axis=1)
#print(len(names_simp))
bar_sum_global = {}
for k in range(1, num_chrs+1):
    c = "chr"+ str(k)
    bar_sum_global_tmp = []
    for i in range(1,(num_samples+1)):
        bar_sum_global_tmp.append((names_simp==i).sum())
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
rep_types = read_annotations(rep_file, rep_types, chrs)
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
exon_comp = {}
gene_comp = {}

for i in range(1, num_chrs+1):
    tmp_chr = "chr"+str(i)
    gene_names[tmp_chr] = []
    gene_locals[tmp_chr] = []
    exon_locals[tmp_chr] = []
    exon_names[tmp_chr] = []
    gene_comp[tmp_chr] = [0]*(num_samples+1)
    exon_comp[tmp_chr] = [0]*(num_samples+1)

gene_locals, gene_names, exon_locals, exon_names, cds_locals, gene_content = read_gene_annotations(gene_file, 
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

    
    cntr = 0
    while cntr < len(exon_locals[tmp_chr]):
        for j in range(0, (num_samples+1)):
            exon_comp[tmp_chr][j] += get_idx_content(exon_locals[tmp_chr], cntr, j) #names_simp[(exon_locals[cntr]*n_skips):(exon_locals[cntr+1]*n_skips)].count(j)
        cntr += 3
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

##Make the gene content plot
# Plot the % of each gene that is unique
# Plot the % of each gene that is universal 
# The yellow and purple plot that appears in the lower right hand order
#gene_content_plot = plot_gene_content(gene_content, sort_by, colors)

##Set the smoothing variables. These will be changable later on. 
window_filter = 53
poly_order = 3
shared_kmers = [1]
#print(len(names_simp))
x_start_init_adjust = int(x_start_init/n_skips_start)
x_stop_init_adjust = int(x_stop_init/n_skips_start)
#### NEED TO ADJUST THE START INIT AND STOP INIT SO THAT WE TAKE THE SKIPS INTO ACCOUNT
#Here is where we make the main plot 
#print("*********")
#print(len(gene_locals["chr1"]))
fig, bar_sum_names, bar_sum_regional, colors = plot_interactive( n_skips_start, window_filter, poly_order, shared_kmers,
    layout, exon_comp["chr1"], 
    gene_comp["chr1"], bins, names_simp[x_start_init_adjust:x_stop_init_adjust], "chr1", zs_tmp, rep_types, True, True, 
    x_start_init, x_stop_init, gene_locals["chr1"], gene_names["chr1"])
#print(len(gene_locals["chr1"]))
sort_by = ["Unique","Universal"]
gene_content_plot = plot_gene_content(gene_content["chr1"], sort_by, colors)
#print("AFTER GENE CONTENT", len(gene_locals["chr1"]))
#sample_kmer_counts = parse_counts_complex(all_chrs[chrs][x_start_init:x_stop_init])
#x, y, z_1, z_9, z_genes = make_chr_whole(names_simp, 1000, gene_locals, x_start_init, x_stop_init)
x, z_1, z_9, y  = read_chr_whole()
#print("AFTER READ CHR", len(gene_locals["chr1"]))
z_genes = make_gene_whole_chr(x[0], gene_locals["chr1"])
#print("AFTER MAKING GENE BINS", len(gene_locals["chr1"]))
chr_fig = plot_chr_whole(x[0], z_1[0], z_9[0], z_genes, x_start_init, x_stop_init, y[0])
whole_genome_fig = plot_whole_genome(x, z_1, z_9, y)
#chr_fig = plot_chr_whole(names_simp, full_cnts, 1000, gene_locals, x_start_init, x_stop_init)
print("a4:", time.time()-t)
t = time.time()
#print(all_chrs[chrs])
#Set up the dash app 
chrs = "chr1"
app = dash.Dash(external_stylesheets=["https://www.w3schools.com/w3css/4/w3.css"]) 

app.layout = html.Div([
    #html.Div(id = 'parent', children = [
    html.H1(id = 'H1', children = 'Pan-Genus k-mers', style = {'textAlign':'center'}), 
    dcc.Tabs(id="tabs-graph", value='homepage', style={"font-size": 26}, children=[
        dcc.Tab(label='Whole genome', value='homepage', style={"background-color": "lightgrey", "font-size": 26}, children=[
            html.Div(children=[
                dcc.Graph(id="all_chromosomes",
                    figure = whole_genome_fig)
            ])
        ]),
        dcc.Tab(label='Chromosome scale', value='chr_view', style={"background-color": "lightgrey", "font-size": 26},  children=[
            html.Div(children=[
                html.Div(className="w3-container", children=[
                        #left figure
                        dcc.Graph(id="chromosome",
                            figure = chr_fig)
                ])
            ]),
            html.Div(children=[
                html.Div(className="w3-container", children=[
                    html.Div(className="w3-threequarter", children=[
                        #left figure - I'm calling this the "Main" figure
                        dcc.Graph(id="primary",figure = fig)
                    ]), #style={'padding-top' : '1%', 'padding-left' : '1%', 'padding-bottom' : '1%', 'padding-right' : '1%', 
                    html.Div(className="w3-quarter", children=[
                        dcc.Graph(id="Secondary", 
                            #Now we have the phylogenetic tree
                            figure = get_local_info(names_simp[x_start_init:x_stop_init], exon_comp[chrs], gene_comp[chrs], 
                                bar_sum_regional, bar_sum_global[chrs])),
                    ])
                ])
            ]),
            html.Div(children=[
                html.Div(className="w3-container", children=[
                    html.Div(className="w3-quarter", children=[
                            dcc.Graph(id="Genes", 
                                figure=gene_content_plot),
                        ]),
                    html.Div(className="w3-half", children=[
                        dcc.Graph(id="Third", 
                            #This is the histogram section
                            
                            figure = create_tree(tree_file, x_start_init,x_stop_init, all_chrs[chrs][int(x_start_init/n_skips_start):int(x_stop_init/n_skips_start)], n_skips_start)),
                    ]),
                    html.Div(className="w3-quarter", children=[
                        #right figure (This is where we will have user input)
                    #    #right figure
                            html.I("Change x limits or jump to a gene of interest", style={ 
                                'font-size':'18px', 'fontWeight': 'bold'}),
                            html.Br(),
                            html.Br(),
                            html.I("X start: "),
                            dcc.Input(
                                id='num-start',
                                type='number',
                                value=x_start_init,
                                placeholder="X start",
                                debounce=True
                            ),
                            html.I("\t\tGene name: "),
                            dcc.Input(
                                id='gene_jump',
                                type='text',
                                placeholder="Type gene name here",
                                debounce=True
                            ),
                            html.Br(),
                            html.I("X stop: "),
                            dcc.Input(
                                id='num-stop',
                                type='number',
                                placeholder="X stop",
                                value=x_stop_init,
                                debounce=True
                            ),
                            html.I("\t\tBins: "),
                            dcc.Input(
                                id='bins',
                                type='number',
                                placeholder="# of bins",
                                value=bins
                            ),
                            html.Br(),
                            html.Br(),

                            dcc.Checklist(["Genes", "Repeats"],[], id="my_check"),
                            html.Br(),
                            
                            html.I("Smoothing: "),
                            html.Br(),
                            dcc.Checklist(['1', '2', '3', '4', '5', '6', '7', '8', '9'],
                                ['1'],
                                inline=True,
                                id="SG_check",
                            ),
                            #html.Br(),
                            html.I("SG window: "),
                            dcc.Input(
                                id='SG_window',
                                type='number',
                                value=window_filter
                            ),
                            html.Br(),
                            html.I("\t\tSG polynomial order: "),
                            dcc.Input(
                                id='SG_polynomial_order',
                                type='number',
                                value=poly_order
                            ),

                            html.Br(),

                            html.Pre(id='tester'),

                            html.Button('Submit', id='textarea-state-button', n_clicks=0),
                            html.Button('Save', id='save-button', n_clicks=0),
                            html.Button('Reset', id='reset_clicks', n_clicks=0),
                        ], style={'padding-top' : '1%', 'padding-left' : '1%', 'padding-bottom' : '1%', 'padding-right' : '1%', 'font-size':'15px',
                            'textAlign': 'left', "border":"2px grey solid"})
                    ])
                ])
        ]),
        
    #]),
    ]),
    html.Div("chr1", id="selected-chrom-state"),
    html.Div("chr1", id="prev_selected-state"),
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
    Output('prev_selected-state', 'children'),

    Input('all_chromosomes','relayoutData'),
    Input('chromosome','figure'),
    Input('primary','figure'),
    Input('Secondary', 'figure'),
    Input('Third', 'figure'),
    Input('Genes', 'figure'),
    Input('my_check','value'),
    Input('num-start','value'),
    Input('num-stop','value'),
    Input('gene_jump', 'value'),
    Input('bins', 'value'),
    Input('textarea-state-button', 'n_clicks'),
    Input('save-button', 'n_clicks'),
    Input('SG_window','value'),
    Input('SG_polynomial_order','value'),
    Input('SG_check','value'),
    Input('primary', 'clickData'),
    Input('primary', 'relayoutData'),
    Input('chromosome', 'selectedData'),
    Input('Genes', 'clickData'),

    State('selected-chrom-state', 'children'),
    State('prev_selected-state', 'children'),
    )
def update_figure(wgr, chr_fig, fig1, fig2, fig3, fig4, clicker, x_start, x_stop, gene_jump, bins, n_clicks, save_clicks, SG_window, 
    SG_polynomial_order, SG_check, clickData, relayoutData, chr_relayoutData, gene_jump_bottom, chrs, pchrs):#clicker_update, gene_jump, bins):
    
    triggered_id = ctx.triggered_id
    print(triggered_id)
    print(wgr)
    
    print()
    n_skips = 100
    click_me_genes = True
    click_me_rep = True
    if triggered_id == 'all_chromosomes':
        #Get the chromosome number 
        if wgr != None and len(wgr)>1:
            chr_num_tmp = list(wgr.keys())[0].split('.')[0].split('axis')[1] #if wgr != None and len(wgr)>1:
            if len(chr_num_tmp) == 0:
                chr_num = 1
            else:
                chr_num = int(int(chr_num_tmp)/2)
            print("SET IT", chr_num)
            #print(list(wgr.keys())[0].split('.')[0].split('axis')[1])
            #chr_fig = #This is the top pane, showing the whole chromosome. 
            #fig1 = #The main figure
            #fig2 = #The "local info" figure. Showing the k-mer composition in different regions 
            #fig3 = #Tree figure 
            #fig4 = #Sorted gene figure 
            return update_all_figs(chr_num, chr_relayoutData, click_me_rep, click_me_genes, SG_window, SG_polynomial_order,
                SG_check, chrs, pchrs) #chr_fig, fig1, fig2, fig3, fig4, chr_relayoutData
    elif triggered_id == 'Genes':
        return sorted_gene_fig(chr_fig, fig1, fig2, fig3, fig4, gene_jump_bottom, SG_window, SG_polynomial_order,
            SG_check, n_skips, click_me_genes, click_me_rep, chr_relayoutData, chrs, pchrs)
    elif triggered_id == 'chromosome':
        return chromosome_gene_triggered(chr_fig, fig1, fig2, fig3, fig4, gene_jump_bottom, SG_window, SG_polynomial_order,
            SG_check, n_skips, click_me_genes, click_me_rep, chr_relayoutData, chrs, pchrs)
    elif triggered_id == 'primary':
        return primary_fig_triggered(chr_fig, fig1, fig2, fig3, fig4, gene_jump_bottom, SG_window, SG_polynomial_order,
            SG_check, n_skips, click_me_genes, click_me_rep, chr_relayoutData, relayoutData, clickData, x_start, x_stop, chrs, pchrs)
    return chr_fig, fig1, fig2, fig3, fig4, chr_relayoutData, chrs, pchrs #, clickData

def sorted_gene_fig(chr_fig, fig1, fig2, fig3, fig4, gene_jump_bottom, SG_window, SG_polynomial_order,
            SG_check, n_skips, click_me_genes, click_me_rep, chr_relayoutData, chrs, pchrs):
    print("****Sorted gene fig******")
    print(gene_jump_bottom)
    print(chrs)
    print("**********")
    #chr_start = 0
    #chr_num = int(chrs.split('r')[1])
    #chr_end = chr_lens[chr_num-1]
    #all_chrs_data = anchor.get_counts(chrs, chr_start, chr_end, n_skips_start)
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
            exact_mini_counts = anchor.get_counts(chrs, x_start, x_stop, 1)
            simple_cnts_for_plots = exact_mini_counts.sum(axis=1)
            #simple_cnts_for_plots, exact_mini_counts = read_mini_count_files(x_start, x_stop)
            n_skips = 1
        else: 
            x_start = int(gene_locals[chrs][tmp_idx]) - buff #int(1000*n_skips/2)
            x_stop = int(gene_locals[chrs][tmp_idx+1]) + buff #int(1000*n_skips/2)
            n_skips = n_skips_start
            simple_cnts_for_plots = names_simp[int(x_start/n_skips):int(x_stop/n_skips)]
            exact_mini_counts = all_chrs[chrs][int(x_start/n_skips):int(x_stop/n_skips)]
        fig1, bar_sum_names, bar_sum_regional, colors = plot_interactive( n_skips, int(SG_window), 
            int(SG_polynomial_order), SG_check,
            layout, exon_comp[chrs], gene_comp[chrs], bins, simple_cnts_for_plots, chrs, zs_tmp, rep_types, click_me_rep, click_me_genes,  
            int(x_start), int(x_stop), gene_locals[chrs], gene_names[chrs])
        fig2 = get_local_info(
            simple_cnts_for_plots, exon_comp[chrs], gene_comp[chrs], bar_sum_regional, bar_sum_global[chrs])
        fig3 = create_tree(tree_file, x_start, x_stop, exact_mini_counts, n_skips)
        gene_jump = ""
        if x_start_init != x_start:
            chr_fig = plot_chr_whole(x[chr_num-1], z_1[chr_num-1], z_9[chr_num-1], z_genes, x_start, x_stop,  y[chr_num-1])
    return chr_fig, fig1, fig2, fig3, fig4, chr_relayoutData, chrs, pchrs
def chromosome_gene_triggered(chr_fig, fig1, fig2, fig3, fig4, gene_jump_bottom, SG_window, SG_polynomial_order,
    SG_check, n_skips, click_me_genes, click_me_rep, chr_relayoutData, chrs, pchrs):
    print(chrs)
    #chr_start = 0
    #chr_num = int(chrs.split('r')[1])
    #chr_end = chr_lens[chr_num-1]
    #all_chrs_data = anchor.get_counts(chrs, chr_start, chr_end, n_skips_start)
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
        exact_mini_counts = anchor.get_counts(chrs, x_start, x_stop, 1)
        simple_cnts_for_plots = exact_mini_counts.sum(axis=1)
        n_skips = 1
    else:
        exact_mini_counts = all_chrs[chrs][int(x_start/n_skips_start):int(x_stop/n_skips_start)]
        simple_cnts_for_plots = names_simp[int(x_start/n_skips_start):int(x_stop/n_skips_start)]
        n_skips = n_skips_start
    fig1, bar_sum_names, bar_sum_regional, colors = plot_interactive( n_skips, int(SG_window), 
            int(SG_polynomial_order), SG_check,
            layout, exon_comp[chrs], gene_comp[chrs], bins, simple_cnts_for_plots, chrs, zs_tmp, rep_types, click_me_rep, click_me_genes,  
            int(x_start), int(x_stop), gene_locals[chrs], gene_names[chrs])
    fig2 = get_local_info(
            simple_cnts_for_plots, exon_comp[chrs], gene_comp[chrs], bar_sum_regional, bar_sum_global[chrs])
    #And now we update the histograms 
    fig3 = create_tree(tree_file, x_start, x_stop, exact_mini_counts, n_skips)
    chr_relayoutData = None#{'autosize': True}
    if x_start_init != x_start:
        chr_fig = plot_chr_whole(x[chr_num-1], z_1[chr_num-1], z_9[chr_num-1], z_genes, x_start, x_stop,  y[chr_num-1])
    return chr_fig, fig1, fig2, fig3, fig4, chr_relayoutData, chrs, pchrs
def primary_fig_triggered(chr_fig, fig1, fig2, fig3, fig4, gene_jump_bottom, SG_window, SG_polynomial_order,
    SG_check, n_skips, click_me_genes, click_me_rep, chr_relayoutData, relayoutData, clickData, x_start, x_stop, chrs, pchrs):
    print(chrs)
    #chr_start = 0
    #chr_num = int(chrs.split('r')[1])
    #chr_end = chr_lens[chr_num-1]
    #all_chrs_data = anchor.get_counts(chrs, chr_start, chr_end, n_skips_start)
    names_simp = all_chrs[chrs].sum(axis=1)
    if relayoutData != None and 'xaxis4.range[0]' in relayoutData.keys():
        print("fourth elif")
        x_start = int(relayoutData['xaxis4.range[0]'])
        x_stop = int(relayoutData['xaxis4.range[1]'])
        #fig1, fig2, fig3 = update_all_figs(x_start, x_stop, click_me_rep, click_me_genes)
        fig1, bar_sum_names, bar_sum_regional, colors = plot_interactive( n_skips_start, int(SG_window), 
            int(SG_polynomial_order), SG_check,
            layout, exon_comp[chrs], gene_comp[chrs], bins, names_simp[int(x_start/n_skips_start):int(x_stop/n_skips_start)], chrs, zs_tmp, rep_types, click_me_rep, click_me_genes,  
            int(x_start), int(x_stop), gene_locals[chrs], gene_names[chrs])
        fig2 = get_local_info(
            #full_cnts[x_start:x_stop], 
            names_simp[int(x_start/n_skips_start):int(x_stop/n_skips_start)], exon_comp[chrs], gene_comp[chrs], bar_sum_regional, bar_sum_global[chrs])
        #And now we update the histograms 
        fig3 = create_tree(tree_file, x_start, x_stop, 
            all_chrs[chrs][int(x_start/n_skips_start):int(x_stop/n_skips_start)], n_skips_start)
        chr_fig = plot_chr_whole(x[chr_num-1], z_1[chr_num-1], z_9[chr_num-1], z_genes, x_start, x_stop,  y[chr_num-1])
    elif clickData != None:#len(print(clickData['points'])) > 0:
        print("third elif")
        tmp_idx = gene_names.index(clickData['points'][0]['text'])
        #gene_buffer = get_buffer(int(gene_locals[tmp_idx]), int(gene_locals[tmp_idx+1]), n_skips_start)
        x_start = int(gene_locals[tmp_idx])-10
        x_stop = int(gene_locals[tmp_idx+1])+10
        if get_buffer(x_start, x_stop, n_skips_start) == 1:
            #we need to move to single nucleotide resoltion
            exact_mini_counts = anchor.get_counts(chrs, x_start, x_stop, 1)
            simple_cnts_for_plots = exact_mini_counts.sum(axis=1)
            n_skips = 1
        else:
            exact_mini_counts = all_chrs[chrs][int(x_start/n_skips_start):int((x_stop)/n_skips_start)]
            simple_cnts_for_plots = names_simp[int(x_start/n_skips_start):int((x_stop)/n_skips_start)]
            n_skips = n_skips_start 
        #fig1, fig2, fig3 = update_all_figs(x_start, x_stop, click_me_rep, click_me_genes)
        fig1, bar_sum_names, bar_sum_regional, colors = plot_interactive( n_skips, int(SG_window), 
            int(SG_polynomial_order), SG_check,
            layout, exon_comp[chrs], gene_comp[chrs], bins, simple_cnts_for_plots, chrs, zs_tmp, rep_types, 
            click_me_rep, click_me_genes,  
            int(x_start), int(x_stop), gene_locals[chrs], gene_names[chrs])
        fig2 = get_local_info(
            #full_cnts[x_start:x_stop], 
            simple_cnts_for_plots, 
            exon_comp[chrs], gene_comp[chrs], bar_sum_regional, bar_sum_global[chrs])
        #And now we update the histograms 
        fig3 = create_tree(tree_file, x_start, x_stop, 
            exact_mini_counts, n_skips)
        chr_fig = plot_chr_whole(x[chr_num-1], z_1[chr_num-1], z_9[chr_num-1], z_genes, x_start, x_stop,  y[chr_num-1])
    return chr_fig, fig1, fig2, fig3, fig4, chr_relayoutData, chrs, pchrs
def render_tab_content(active_tab, chrs):
    return chr_fig, fig1, fig2, fig3, fig4, chr_relayoutData, chrs, pchrs
def update_all_figs(chr_num, chr_relayoutData, click_me_rep, click_me_genes, SG_window, SG_polynomial_order,
    SG_check, chrs, pchrs):
    #chr_fig = This is the top pane, showing the whole chromosome. 
    #fig1    = The main figure
    #fig2    = The "local info" figure. Showing the k-mer composition in different regions 
    #fig3    = Tree figure 
    #fig4    = Sorted gene figure
    print(chrs)
    chrs = "chr" + str(chr_num)
    print(chrs)
    #chr_start = 0
    #chr_end = chr_lens[chr_num-1]
    #all_chrs[chrs] = anchor.get_counts(chrs, chr_start, chr_end, n_skips_start)
    names_simp = all_chrs[chrs].sum(axis=1) #all_chrs["chr"+str(chr_num)].sum(axis=1)#tmp.sum(axis=1)
    print(names_simp[0:100])
    
    '''bar_sum_global = []
    for i in range(1,(num_samples+1)):
        bar_sum_global.append((names_simp==i).sum())
    gene_names = []
    gene_locals = []
    exon_locals = []
    exon_names = []
    cds_locals = []
    cds_names = []
    #rep_types = []
    zs_tmp = []
    '''
    #gene_file = "test_data/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3" #"test_data/Squi2.chr"+str(chr_num)+".exons.genes.cds.txt"
    #rep_file = "test_data/Zm-B73-REFERENCE-NAM-5.0.TE.gff3" #"test_data/Solqui2.repeats.chr"+str(chr_num)+".txt"
    #gene_locals, gene_names, exon_locals, exon_names, cds_locals, gene_content = read_gene_annotations(gene_file, 
    #    gene_names, gene_locals, exon_locals, exon_names, cds_locals, n_skips_start, chrs)
    #rep_types = get_init_rep_types()
    #print(rep_file)
    #rep_types = read_annotations(rep_file, rep_types, chrs)

    #Get the different compositions 
    '''
    cntr = 0
    cntr2 = 0
    while cntr < len(gene_locals):
        tmp1 = 0
        for j in range(0, 3):
            tmp2 = get_idx_content(gene_locals, cntr, j) #names_simp[(gene_locals[cntr]/n_skips):(gene_locals[cntr+1]/n_skips)].count(j)
            tmp1 += tmp2
            gene_comp[j] += tmp2 #names_simp[gene_locals[cntr]:gene_locals[cntr+1]].count(j)
        cntr2 += 1
        for j in range(3, (num_samples+1)):
            gene_comp[j] += get_idx_content(gene_locals, cntr, j) #names_simp[(gene_locals[cntr]*n_skips):(gene_locals[cntr+1]*n_skips)].count(j)
        cntr += 3

    exon_comp = [0]*(num_samples+1)
    cntr = 0
    while cntr < len(exon_locals):
        for j in range(0, (num_samples+1)):
            exon_comp[j] += get_idx_content(exon_locals, cntr, j) #names_simp[(exon_locals[cntr]*n_skips):(exon_locals[cntr+1]*n_skips)].count(j)
        cntr += 3

    cds_comp = [0]*(num_samples+1)
    '''
    #Make the names_sim data 
    #Find zs_tmp

    #chrs = "chr" + str(chr_num)
    fig1, bar_sum_names, bar_sum_regional, colors = plot_interactive( n_skips_start, int(SG_window), 
            int(SG_polynomial_order), SG_check,
            layout, exon_comp[chrs], gene_comp[chrs], bins, names_simp[x_start_init_adjust:x_stop_init_adjust], 
            chrs, zs_tmp, rep_types, click_me_rep, click_me_genes,  
            int(x_start_init), int(x_stop_init), gene_locals[chrs], gene_names[chrs])
    
    sort_by = ["Unique","Universal"]
    fig4 = plot_gene_content(gene_content[chrs], sort_by, colors)
    
    #x, z_1, z_9, y  = read_chr_whole()
    #z_genes = make_gene_whole_chr(x[chr_num-1], gene_locals)
    chr_fig = plot_chr_whole(x[chr_num-1], z_1[chr_num-1], z_9[chr_num-1], z_genes, x_start_init, 
        x_stop_init, y[chr_num-1])

    fig2 = get_local_info(
            names_simp[x_start_init_adjust:x_stop_init_adjust], exon_comp[chrs], gene_comp[chrs], 
            bar_sum_regional, bar_sum_global[chrs])
    
    fig3 = create_tree(tree_file, x_start_init_adjust, x_stop_init_adjust, 
        all_chrs[chrs][x_start_init_adjust:x_stop_init_adjust], n_skips_start)
    
    #chr_fig = plot_chr_whole(x[chr_num-1], z_1[chr_num-1], z_9[chr_num-1], z_genes[chr_num-1], x_start_init, x_stop_init, y[chr_num-1])
    return chr_fig, fig1, fig2, fig3, fig4, chr_relayoutData, chrs, pchrs


@app.callback(Output('textarea-state-button', 'n_clicks'),
    Output('num-start','value'),
    Output('num-stop','value'),
    [Input('reset_clicks','n_clicks')])
def update(reset):
    return 0, x_start_init, x_stop_init

app.run_server(host='127.0.0.1', debug=True)





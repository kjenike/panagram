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
#from Bio import Phylo
#import dash_cytoscape
#import itertools
#import dash_core_components as dcc
#import dash_html_components as html
from dash import Dash, dcc, html, Input, Output, ctx
import math
#plt.rcParams["figure.figsize"] = (40,40)
sns.set_palette('viridis', 9)
#plt.rcParams.update({'font.size': 22})

#bins = 500

#/Users/katiejenike/Desktop/Lab_notes/PanSol/PANGEN/DEBUG/QC_PLOTS/SQUI2_PANGEN_PLOTS
#python plot_chr_interactive.py result_SQUI2.chr1.txt SQUI2

in_file=sys.argv[1]
sample = sys.argv[2]
count_file = "SQUI2.chr1.txt"
rep_file = "Solqui2.repeats.chr1.txt"
gene_file = "Squi2.exons.genes.cds.txt"
tree_file = "all.alns.aug5.2022.fa.treefile"

##### Make the naming index
#names = {}
#allowed_nums = [0, 1, 3, 4, 5, 6, 8, 9]
#for i in allowed_nums:
#    for j in allowed_nums:
#        for k in allowed_nums:
#            this_name = str(int((str(i) + str(j) + str(k))))
#            names["pos_"+this_name] = []

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

def find_shared(some_num):
    total = 0
    for i in some_num:
        total += name_idx_short[i]
    return total

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

gene_types = {
    
}

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

#print(len(cnts_tmp))
def parse_counts_1():
    cntr = 0
    cnts = []
    colors = []
    y = []
    last = 0
    for this in cnts_tmp:
        #if cntr % 1000 == 1:
        if this > 1:    
            if len(x) > 1 :
                #last = #cnts[-1]
                if this != last:
                    x.append(cntr)
                    #y.append(1)
                    cnts.append(sample_counts[this])
                    colors.append(sample_colors[this])
                #last = this
            else:
                x.append(cntr)
                #y.append(1)
                cnts.append(sample_counts[this])
                colors.append(sample_colors[this])
        cntr += 1
        last = this

def parse_counts_2():
    
    cntr = 0
    for this in cnts_tmp:
        if this == 1:
            ones.append(cntr)
            ones_n.append("One")
        elif this == 4:
            fours.append(cntr)
            fours_n.append("Four")
        elif this == 6:
            sixs.append(cntr)
            sixs_n.append("Six")
        elif this == 9:
            nines.append(cntr)
            nines_n.append("Nine")
        cntr += 1

def parse_counts_3():
    cntr = 0
    for this in cnts_tmp:
        name = "pos_"+this
        #name_n = "name_" + this
        names[name].append(cntr)
        #name_n.append(this)
        cntr += 1

def parse_counts_simple(cnts_tmp):
    names_simp = []
    for this in cnts_tmp:
        names_simp.append(name_idx_long[this])
    return names_simp

print("About to plot")

##### PLOTS

def get_index(j, bin_size, x_start):
    indx = int((j-x_start)/bin_size)
    return indx

def get_index_0_based(j, bin_size):
    indx = int(j/bin_size)
    return indx

def read_count_file_old(count_file):
    xs = []
    ys = []
    ys_tmp = []
    with open(count_file, "r") as f:
        line = f.readline()
        while line:
            ys_tmp = line.strip().split(":")[1].split(',')[:-1]
            line = f.readline()
    avg_y = 0
    for i in range(0, len(ys_tmp)):
        xs.append(i)
        ys.append(int(ys_tmp[i]))

    return xs, ys

def read_count_file(count_file, y_indx ):
    with open(count_file, "r") as f:
        line = f.readline()
        while line:
            zs_tmp = line.strip().split(":")[1].split(',')[:-1]
            line = f.readline()
    #full_counts = [int(i) for i in zs_tmp]
    zs_tmp_2 = [len(i) for i in zs_tmp] #list(map(lambda x: len(x), zs_tmp))
    return zs_tmp_2, zs_tmp 

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

tree_fig = {}

def get_x_coordinates(tree):
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
    """Recursively draw the tree branches, down from the given clade"""
    x_curr = x_coords[clade]
    y_curr = y_coords[clade]
    if str(clade) != "Clade" and str(clade) != "AT":
        tmp_name = str(clade)
        if tmp_name.count("_") > 0:
            tmp_name = tmp_name.split("_")[1]
        print(tmp_name)
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
    
    kmer_num = { "Solyc":0, "Solabu2":0, "Solaet3":0, "Solmac3":0, 
        "Solpri1":0,"Solcan1":0,"Solqui2":0,"Solmur2hap1":0,"Solmur2hap2":0,
    }

    for c in raw_counts:
        #Now we have to go through each digit 
        if c != "0":
            for t in name_idx_complex[c]:
                kmer_num[t] += 1

    return kmer_num

def create_tree(tree_file, x_start_init, x_stop_init, raw_counts):
    kmer_num_tmp = parse_counts_complex(raw_counts[x_start_init:x_stop_init])
    palette = sns.color_palette("RdPu_r", 160).as_hex()
    total_kmers = kmer_num_tmp["Solqui2"]
    kmer_num = {}
    color_code = {}
    for k in kmer_num_tmp.keys():
        kmer_num[k] = float(((kmer_num_tmp[k])/total_kmers)*100)
        color_code[k] = palette[int(((kmer_num_tmp[k])/total_kmers)*100)+10]
    #print(kmer_num)
    #print(color_code)
    tree = Phylo.read(tree_file, "newick")
    x_coords = get_x_coordinates(tree)
    y_coords = get_y_coordinates(tree)
    line_shapes = []
    draw_clade(color_code, total_kmers, kmer_num, palette, tree.root, 0, line_shapes, line_color= "blue",#'rgb(25,25,25)',
                line_width=1, x_coords=x_coords,
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
    def make_anns(x,y, text, kmer_num):
        tmp_txt = ""
        if text.count("_") > 0:
            tmp_txt = text.split("_")[1]
        else:
            tmp_txt = text
        tmp_txt += " - " + str(kmer_num[tmp_txt])[:4] + "% (" + str(kmer_num_tmp[tmp_txt]) + ")" 
        return dict(xref='x', yref='y', x=x, y=y, text="\t" + tmp_txt,
            showarrow=False,
            xanchor='left',
            yanchor='middle',)
    annotations = []
    for i in range(0, len(label_legend)):
        annotations.append(make_anns(an_x[i], an_y[i], label_legend[i], kmer_num))
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
                             range=[0, 0.1]
                             ),
                  yaxis=axis,
                  hovermode='closest',
                  shapes=line_shapes,
                  plot_bgcolor='rgb(250,250,250)',
                  legend={'x': 0, 'y': 1},
                  annotations=annotations,
                  #xaxis_range=[0,0.1]
                  )
    fig = make_subplots(
        rows=1, cols=3,
        specs=[[{"type": "scatter", 'colspan':2}, None, {"type": "bar"} ]], #, {"type": "pie"} , {"type": "pie"} ]],
        horizontal_spacing=0.01,
        #subplot_titles=("This region","CDS","Exons","Genes", "Whole chromosome")
    )
    fig.add_trace(node, row=1, col=1)
    fig.update_layout(layout)
    #Sort the hist bars 
    hist_y = []#list(kmer_num_tmp.values()).sort()
    hist_x = []
    color_hist = []
    #use the color_codes 
    for y in list(kmer_num_tmp.keys()):
        color_hist.append(color_code[y])
    fig.add_trace(go.Bar(y=list(kmer_num_tmp.keys()), x=list(kmer_num_tmp.values()), showlegend=False,
        marker=dict(color=color_hist), orientation='h'), 
        row=1, col=3)
    fig.update_yaxes(visible=False, showticklabels=False)
    fig.update_layout(margin=dict(
            t=20,
            b=10,
            l=10,
            r=10))
    return fig

def get_global_info(colors, exon_comp, gene_comp, bar_sum_names, bar_sum_regional, bar_sum_global):

    fig = make_subplots(
        rows=1, cols=5,
        specs=[[{"type": "pie"}, {"type": "pie"}, {"type": "pie"} , {"type": "pie"} ]],
        subplot_titles=("This region","CDS","Exons","Genes", "Whole chromosome")
    )
    fig.add_trace(go.Pie(labels=bar_sum_names, values=bar_sum_regional, showlegend=False, marker_colors=colors, 
        name="Regional"), row=1, col=1)
    fig.add_trace(go.Pie(labels=bar_sum_names, values=cds_comp, showlegend=False, marker_colors=colors, 
        name="CDS"), row=1, col=2)
    fig.add_trace(go.Pie(labels=bar_sum_names, values=exon_comp, showlegend=False, marker_colors=colors, 
        name="Exons"), row=1, col=3)
    fig.add_trace(go.Pie(labels=bar_sum_names, values=gene_comp, showlegend=False, marker_colors=colors, 
        name="Genes"), row=1, col=4)
    fig.add_trace(go.Pie(labels=bar_sum_names, values=bar_sum_global, showlegend=False, marker_colors=colors, 
        name="Chromosome"), row=1, col=5)
    
    fig.update_layout(height=1000)
    return fig

def get_local_info(full_counts_ref, x_all, exon_comp, gene_comp, bar_sum_regional, bar_sum_global, cds_comp):
    #fig = make_subplots(
    #    rows=1, cols=4,
    #    specs=[[{"type": "pie"}, {"type": "pie"}, {"type": "pie"}, {"type": "pie"} ]]
    #)#go.Figure()
    x_ref = [int(i) for i in full_counts_ref]
    fig = make_subplots(
        rows=2, cols=2,
        #shared_yaxes=True,
        specs=[[{"type": "bar"}, {"type": "bar"}], 
            [{"type": "bar"}, {"type": "bar"} ]],
        subplot_titles=("Whole chromosome", "This region", "Exons", 
            "Genes", ), 
        vertical_spacing=0.1
    )
    x = [1,2,3,4,5,6,7,8,9]
    colors = ["#fde725", "#addc30", "#5ec962", "#28ae80", "#21918c", "#2c728e", "#3b528b", "#472d7b", "#440154"]
    #Repeats
    #fig.add_trace(go.Histogram(x=x_ref, marker_color="#21918c", opacity=0.7, showlegend=False), row=1, col=1)
    #fig.add_trace(go.Bar(x=x_ref_whole, marker_color="#330C73", opacity=0.7), row=1, col=1)
    #fig.update_yaxes(type='log')
    #fig.update_layout(xaxis_range=[0,25000], xaxis_title_text="K-mer count in reference", yaxis_title_text='Frequency (log)')

    #This region
    y=[(i/sum(bar_sum_regional[1:])*100) for i in bar_sum_regional[1:]]
    y_whole=[(i/sum(bar_sum_global[1:])*100) for i in bar_sum_global[1:]]
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

    #Whole chromosome
    #y=[(i/sum(bar_sum_global[1:])*100) for i in bar_sum_global[1:]]
    #fig.add_trace(go.Bar(x=x, y=y_whole, marker_color="#7400b8",  showlegend=False), row=1, col=5)
    #fig.update_layout(xaxis_title_text="K-mers shared in X samples", yaxis_title_text='Frequency (log)')
    #fig.add_trace(go.Histogram(x=x_all_whole, marker_color="#330C73", opacity=0.7), row=2, col=1)

    fig.update_layout(height=1000)
    return fig

def plot_interactive( window_filter, poly_order, shared_kmers, layout, exon_comp, gene_comp, bins, names_simp, name, zs_tmp, 
    rep_types, plot_rep, plot_gene, x_start, x_stop, gene_locals, gene_names):

    tmp_lst = []
    fig = make_subplots(        
        rows=14, cols=1,
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
        [{"type": "heatmap"} ],#, None, None, None, None, None], 
        #[{"type": "pie", 'rowspan':2}, {"type": "pie", 'rowspan':2}, {"type": "pie", 'rowspan':2},{"type": "pie", 'rowspan':2}, {"type": "scatter", 'rowspan':2}, {"type": "scatter", 'rowspan':2}],
        #[None,None, None, None, None, None]
        ],
        subplot_titles=("Ref. Sequence Position","Gene Annotation", "Repeat Annotation",  "Conserved K-mers" , "K-mer repeats")
    )

    t = time.time()
    colors = ["grey","#fde725", "#addc30", "#5ec962", "#28ae80", "#21918c", "#2c728e", "#3b528b", "#472d7b", "#440154"]

    
    #Get the bins
    bin_size = (x_stop-x_start)/bins
    cntr = x_start
    x = []
    while cntr < x_stop:
        x.append(cntr)
        cntr += bin_size
    cntr = 0
    print("a1:", time.time()-t)
    t = time.time()
    
    cats_tmp = [([0] * (bins+1)) for _ in range(10)]
    cntr = 0
    x_tracker = 0
    print("a1.4:", time.time()-t)
    t = time.time()
    
    #cntr_2 = x_start
    bin_size_int = int(bin_size)+1
    for i in names_simp[x_start:x_stop] : #range(x_start, x_stop): #names_simp[x_start:x_stop]:#range(x_start, x_stop):#names_simp:
        #i tells us which y axis to use. 
        #cntr tells us what x axis to use 

        if (cntr % bin_size_int) == 0 :
            x_tracker += 1

        if i < len(cats_tmp) and x_tracker < len(cats_tmp[i]):
            cats_tmp[i][x_tracker] += 1 

        cntr += 1
    print("a1.5:", time.time()-t)
    t = time.time()
    #cntr = 0
    y=[1,10,100,1000,10000,100000]
    zs = [[0]*(bins+1) for _ in range(len(y))]
    zs_2 = parse_counts(zs_tmp, zs, bin_size_int, y,  x_start, x_stop)
    y=["1","10","100","1000","10000", "100000"]
    print("b1:", time.time()-t)
    t = time.time()
    #y=[1,2,10,20,100,200,1000,2000,10000] #,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
    fig.add_trace(go.Heatmap(x=x, y=y, z=zs_2, type = 'heatmap', colorscale='plasma'), row=14, col=1)
    fig.update_xaxes(showticklabels=False, row=14, col=1)
    fig.update_traces(colorbar_len=0.25, colorbar_y=0.1, selector=dict(type='heatmap'))
    #fig.update_traces(colorbar_yanchor="top", selector=dict(type='heatmap'))
    
    #colorbar_yanchor="bottom", 
    #fig.update_traces(colorbar_orientation='h', selector=dict(type='heatmap'))
    #fig.update_xaxes(type='category', row=1, col=1)
    #fig.update_yaxes(type="log", row=1, col=1)
    print("b2:", time.time()-t)
    t = time.time()

    #Add a line plot that will cover the different repetative elements. 
    #plot_reps()
    if plot_rep == True:
        #print("About to add next trace")
        rep_colors = ["#f0f921", "#f8df25", "#fdc627", "#fdaf31", "#f99a3e", "#f3854b", "#e97257", "#de6164", "#d24f71", "#c43e7f", 
                "#b42e8d", "#a21d9a", "#8e0ca4", "#7801a8", "#6100a7", "#4903a0", "#2f0596", "#0d0887"]
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
        #print("e:", time.time()-t)
        #t = time.time()
        #fig.update_layout(xaxis_range=[x_start,x_stop])
    
    #This is the conserved kmer plotting section
    bar_sum_regional = []
    bar_sum_names = []
    bar_sum_global = []
    for i in range(0, len(cats_tmp)):
        bar_sum_regional.append(sum(cats_tmp[i]))
        bar_sum_global.append(names_simp.count(i))
        bar_sum_names.append("Present in "+str(i) + " samples")
        fig.add_trace(go.Bar(x=x, y=cats_tmp[i], name="Present in "+str(i) + " samples",
            legendgroup="group1", 
            legendgrouptitle_text="Conserved K-mers",
            marker=dict(color=colors[i]), 
            marker_line=dict(color=colors[i])), row=6, col=1 )
        #cntr += 1

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

    fig.update_layout(height=1000, xaxis_range=[x_start,x_stop])
    print("a3:", time.time()-t)
    t = time.time()
    return fig, bar_sum_names, bar_sum_regional, bar_sum_global, colors

def make_chr_whole(names_simp, cnts, whole_bins, gene_locals, x_start, x_stop):
    #Let's use 1000 bins for the whole chromosome
    #whole_bins = 1000

    whole_bin_size = int(len(names_simp)/whole_bins)
    
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
        elif i == 9:
            #tmp_x = int(cntr/whole_bin_size)
            z_9_tmp[tmp_x] += 1
            #And now we need to figure out the counts 
        if cnts[cntr] != "0":
            z_cnts_tmp[tmp_x] += (i/int(cnts[cntr]))
        cntr += 1
    
    z_1    = [0]*(whole_bins+1)
    z_9    = [0]*(whole_bins+1)
    z_cnts = [0.0]*(whole_bins+1)
    
    for i in range(0, (whole_bins+1)):
        z_1[i]    = (z_1_tmp[i]/whole_bin_size)*100
        z_9[i]    = (z_9_tmp[i]/whole_bin_size)*100
        z_cnts[i] = z_cnts_tmp[i]/whole_bin_size #(z_cnts_tmp[i]/z_9_tmp[i])

    #x_all = x + x + x
    #y_all = [1]*whole_bins + [2]*whole_bins + [3]*whole_bins
    #z_all = z_genes + z_1 + z_9
    return x, y, z_1, z_9, z_genes, z_cnts #x_all, y_all, z_all,  x, z_cnts
    
def plot_chr_whole(x, y, z_1, z_9, z_genes, z_cnts, x_start, x_stop): 
    
    chr_fig = make_subplots(rows=4, cols=1, 
        specs=[[{"type": "heatmap",}], [{"type": "heatmap",}], [{"type": "heatmap",}], [{"type": "scatter",}]],
        shared_xaxes=True,
        #subplot_titles=("", 
        #    "(inter k-mers frequency) / (intra k-mers frequency)"),
        vertical_spacing = 0.0
        )

    chr_fig.add_trace(go.Heatmap(x=x, y=y, z=z_9, type = 'heatmap', colorscale='plasma', showlegend=False, showscale=False), row=1, col=1)
    chr_fig.add_trace(go.Scatter(x=[x_start, x_start, None, x_stop, x_stop, None, x_start, x_stop, ], showlegend=False,
                   y=[0.5, 1.5, None, 0.5, 1.5, None, 1.5, 1.5 ],
                   mode='lines',
                   line_color='#1dd3b0', line_width=4), row=1, col=1)
    ###
    chr_fig.add_trace(go.Heatmap(x=x, y=y, z=z_1, type = 'heatmap', colorscale='plasma', showscale=False), row=2, col=1)
    chr_fig.add_trace(go.Scatter(x=[x_start, x_start, None, x_stop, x_stop], showlegend=False,
                   y=[0.5, 1.5, None, 0.5, 1.5],
                   mode='lines',
                   line_color='#1dd3b0', line_width=4), row=2, col=1)
    ###
    chr_fig.add_trace(go.Heatmap(x=x, y=y, z=z_genes, type = 'heatmap', colorscale='plasma', showscale=False ), row=3, col=1)
    chr_fig.add_trace(go.Scatter(x=[x_start, x_start, None, x_stop, x_stop], showlegend=False, 
                   y=[0.5, 1.5, None, 0.5, 1.5],
                   mode='lines',
                   line_color='#1dd3b0', line_width=4), row=3, col=1)
    ###
    chr_fig.add_trace(go.Scattergl(x=x, y=z_cnts, showlegend=False ), row=4, col=1)
    chr_fig.add_trace(go.Scatter(x=[x_start, x_start, None, x_stop, x_stop, None, x_start, x_stop, ],
                   y=[0, 3, None, 0, 3, None, 0, 0 ], showlegend=False,
                   mode='lines',
                   line_color='#1dd3b0', line_width=4), row=4, col=1)
    #chr_fig.update_xaxes(visible=False, title_text="Sequence position",)
    #chr_fig.update_yaxes(visible=False, title_text="All", row=1, col=1)
    #chr_fig.update_yaxes(visible=False, title_text="Unique", row=2, col=1)
    #chr_fig.update_layout(yaxis_title="Title")
    chr_fig.update_yaxes( range=[0.5,1.5], showticklabels=False, row=1, col=1)
    chr_fig.update_yaxes( range=[0.5,1.5], showticklabels=False, row=2, col=1)
    chr_fig.update_yaxes( range=[0.5,1.5], showticklabels=False, row=3, col=1)
    chr_fig.update_yaxes(visible=False, range=[0,3], row=4, col=1)

    chr_fig.update_xaxes(fixedrange=True, row=1, col=1)
    chr_fig.update_xaxes(fixedrange=True, row=2, col=1)
    chr_fig.update_xaxes(fixedrange=True, row=3, col=1)
    chr_fig.update_xaxes(fixedrange=True, row=4, col=1)
    chr_fig.update_layout(clickmode='event+select', dragmode="select", selectdirection='h')
    #chr_fig.update_traces(showscale=False)
    chr_fig.update_layout(height=350)
    chr_fig.update_layout(margin=dict(
            b=10,
            l=10,
            r=10))
    
    return chr_fig

def plot_gene_content(gene_content):
    df = pd.DataFrame(gene_content)
    x = []
    cntr = 0
    for i in range(0, len(df['Universal'])):
        x.append(cntr)
        cntr +=1
    df_sorted = df.sort_values('Universal')
    fig = go.Figure(data=[go.Scattergl(x=x, y=df_sorted['Unique'], text=df_sorted['Names'], marker=dict(color='#FDD000'),
        name="% Unique", mode="markers")])
    fig.add_trace(go.Scattergl(x=x, y=df_sorted['Universal'], text=df_sorted['Names'],  marker=dict(color='#440154'),
        name="% Universal", mode="markers"))
    fig.update_layout(hovermode='x unified')
    fig.update_layout(height=500)
    fig.update_layout(
        title={"text":"k-mer content in individual genes",
            "xanchor": "center",
            "x":0.5,
            "y":0.9,
            "yanchor": "top"
            },
        xaxis_title="Universal k-mer content rank",
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
chrs="chr1"
x = []
all_chrs = {}
with open(in_file, "r") as f:
    line = f.readline()
    while line:
        all_chrs[chrs] = line.strip().split(":")[1].split(',')[:-1]
        line = f.readline()

print("a0:", time.time()-t)
t = time.time()

#for chrs in all_chrs.keys():
print("a0.1:", time.time()-t)
t = time.time()
names_simp = parse_counts_simple(all_chrs[chrs])
print("Parsed first counts file")
print("a0.2:", time.time()-t)
t = time.time()

y=[1,10,100,1000,10000]

zs_tmp, full_cnts = read_count_file(count_file, y)
print("a0.3:", time.time()-t)
t = time.time()
print("Read in count file")
#Read in the repeats file (this has the annotated repeats from Shujun) 
def read_annotations(ann_file, ann_types):
    with open(ann_file, "r") as f:
        line = f.readline()
        while line:
            tmp = line.split('\t')
            this_type = tmp[1].strip()
            ann_types[this_type].append(int(tmp[2]))
            ann_types[this_type].append(int(tmp[3]))
            ann_types[this_type].append('None')
            line = f.readline()
    return ann_types

def perc_universal(start,stop, countme):
    univ = float(names_simp[start:stop].count(countme))
    #print(names_simp[start:stop])
    #if univ != 0:
    #    print(univ)
    return float(univ/(stop-start))*100

def read_gene_annotations(ann_file, gene_names, gene_locals, exon_locals, exon_names, cds_locals):
    #rep_file = "Solqui2.repeats.chr1.100k.txt"
    gene_content = {"Names":[],"Universal":[], "Unique":[]}
    with open(ann_file, "r") as f:
        line = f.readline()
        while line:
            tmp = line.strip().split('\t')
            this_name = tmp[1]#line.split('\t')[1].strip()
            if this_name.split(":")[0] == "ID=gene":
                gene_locals.append(int(tmp[2]))
                gene_locals.append(int(tmp[3]))
                gene_locals.append('None')
                tmp_name = this_name.split(';')[0].split(":")[1]
                gene_names.append(tmp_name)
                gene_names.append(tmp_name)
                gene_names.append('None')
                gene_content["Names"].append(tmp_name)#perc_universal(int(tmp[2]), int(tmp[3]))
                gene_content["Universal"].append(perc_universal(int(tmp[2]), int(tmp[3]), 9))
                gene_content["Unique"].append(perc_universal(int(tmp[2]), int(tmp[3]), 1))
            elif this_name.split(":")[0] == "ID=exon":
                exon_locals.append(int(tmp[2]))
                exon_locals.append(int(tmp[3]))
                exon_locals.append('None')
                tmp_name = this_name.split(';')[0].split(":")[1]
                exon_names.append(tmp_name)
                exon_names.append(tmp_name)
                exon_names.append('None')
            elif this_name.split(":")[0] == "ID=CDS":
                cds_locals.append(int(tmp[2]))
                cds_locals.append(int(tmp[3]))
                cds_locals.append('None')

            line = f.readline()
    return gene_locals, gene_names, exon_locals, exon_names, cds_locals, gene_content


rep_types = read_annotations(rep_file, rep_types)
print("a0.4:", time.time()-t)
t = time.time()
print("Read annotation file")
#grep -v "#" Solqui2_genes_1.0.0.gff | awk '{if($1=="chr1") print $0}' | awk '{if($3=="gene") print $1"\t"$9"\t"$4"\t"$5}' > Squi2.genes.chr1.txt

gene_names = []
gene_locals = []
exon_locals = []
exon_names = []
cds_locals = []
cds_names = []

gene_locals, gene_names, exon_locals, exon_names, cds_locals, gene_content = read_gene_annotations(gene_file, 
    gene_names, gene_locals, exon_locals, exon_names, cds_locals)

#Get the uniqueness within the genes 
#sorted_gene_uniqness_x, gene_uniqness = [], []

#We need to get the genic kmer composition 
gene_comp = [0]*10
cntr = 0
cntr2 = 0
while cntr < len(gene_locals):
    #print(names_simp[gene_locals[cntr]:gene_locals[cntr+1]])
    tmp1 = 0
    for j in range(0, 3):
        tmp2 = names_simp[gene_locals[cntr]:gene_locals[cntr+1]].count(j)
        tmp1 += tmp2
        gene_comp[j] += tmp2 #names_simp[gene_locals[cntr]:gene_locals[cntr+1]].count(j)
    #print(tmp1/(gene_locals[cntr+1]-gene_locals[cntr]))
    #gene_uniqness.append(tmp1/(gene_locals[cntr+1]-gene_locals[cntr]))
    #sorted_gene_uniqness_x.append(cntr2)
    cntr2 += 1
    for j in range(3, 10):
        gene_comp[j] += names_simp[gene_locals[cntr]:gene_locals[cntr+1]].count(j)
    cntr += 3

#sorted_gene_uniqness = gene_uniqness.sort()

#print(gene_uniqness)

exon_comp = [0]*10
cntr = 0
while cntr < len(exon_locals):
    #print(names_simp[gene_locals[cntr]:gene_locals[cntr+1]])
    for j in range(0, 10):
        exon_comp[j] += names_simp[exon_locals[cntr]:exon_locals[cntr+1]].count(j)
    cntr += 3

cds_comp = [0]*10
cntr = 0
'''while cntr < len(cds_locals):
    #print(names_simp[gene_locals[cntr]:gene_locals[cntr+1]])
    for j in range(0, 10):
        cds_comp[j] += names_simp[cds_locals[cntr]:cds_locals[cntr+1]].count(j)
    cntr += 3
'''
print("a0.5:", time.time()-t)
t = time.time()
#This will have the figure componants that we need 
x_start_init = 75000000#0
x_stop_init = 85000000 #len(zs_tmp)
bins=400

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
gene_content_plot = plot_gene_content(gene_content)

##Add the smoothing
window_filter = 53
poly_order = 3
shared_kmers = [1]
fig, bar_sum_names, bar_sum_regional, bar_sum_global, colors = plot_interactive( window_filter, poly_order, shared_kmers,
    layout, exon_comp, 
    gene_comp, bins, names_simp, chrs, zs_tmp, rep_types, True, True, 
    x_start_init, x_stop_init, gene_locals, gene_names)

#sample_kmer_counts = parse_counts_complex(all_chrs[chrs][x_start_init:x_stop_init])
x, y, z_1, z_9, z_genes, z_cnts = make_chr_whole(names_simp, full_cnts, 1000, gene_locals, x_start_init, x_stop_init)
chr_fig = plot_chr_whole(x, y, z_1, z_9, z_genes, z_cnts, x_start_init, x_stop_init)
#chr_fig = plot_chr_whole(names_simp, full_cnts, 1000, gene_locals, x_start_init, x_stop_init)
print("a4:", time.time()-t)
t = time.time()
#Set up the dash app 
app = dash.Dash(external_stylesheets=["https://www.w3schools.com/w3css/4/w3.css"]) 

app.layout = html.Div([
    #html.Div(id = 'parent', children = [
    html.H1(id = 'H1', children = 'Pan-Genus k-mers', style = {'textAlign':'center'}), 
    
    html.Div(children=[
        html.Div(className="w3-container", children=[
            #html.Div(className="w3-threequarter", children=[
                #left figure
                dcc.Graph(id="chromosome",
                    figure = chr_fig)
            #])
        ])
    ]),
    
    html.Div(children=[
        html.Div(className="w3-container", children=[
            html.Div(className="w3-threequarter", children=[
                #left figure - I'm calling this the "Main" figure
                dcc.Graph(id="primary",figure = fig)
            ]), #style={'padding-top' : '1%', 'padding-left' : '1%', 'padding-bottom' : '1%', 'padding-right' : '1%', 
               #     "border":"6px #1dd3b0 solid"}),
            html.Div(className="w3-quarter", children=[
                dcc.Graph(id="Secondary", 
                    #Now we have the phylogenetic tree
                    figure = get_local_info(full_cnts[x_start_init:x_stop_init], names_simp[x_start_init:x_stop_init], exon_comp, gene_comp, bar_sum_regional, bar_sum_global, cds_comp)),
                    #create_tree(tree_file, sample_kmer_counts)),
                    #figure=get_global_info(colors, exon_comp, gene_comp, bar_sum_names, bar_sum_regional, bar_sum_global)),
                    #figure=get_local_info(full_cnts[x_start_init:x_stop_init], names_simp[x_start_init:x_stop_init]))
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
                #html.H1('Percent of k-mers in each category'),
                dcc.Graph(id="Third", 
                    #This is the histogram section
                    figure= create_tree(tree_file, x_start_init,x_stop_init, all_chrs[chrs])),
                    #get_local_info(full_cnts[x_start_init:x_stop_init], names_simp[x_start_init:x_stop_init], exon_comp, gene_comp, bar_sum_regional, bar_sum_global, cds_comp)),
                    #figure=get_global_info(colors, exon_comp, gene_comp, bar_sum_names, bar_sum_regional, bar_sum_global)),
            ]#, style={'padding-top' : '1%', 'padding-left' : '1%', 'padding-bottom' : '1%', 'padding-right' : '1%', 'font-size':'15px',
            #        'textAlign': 'center', "border":"2px grey solid"}
            ),
            html.Div(className="w3-quarter", children=[
                #right figure (This is where we will have user input)
                #html.Div(className="w3-third", children=[
            #    #right figure
            #    dcc.Graph(id="Secondary", figure=get_local_info(colors, exon_comp, gene_comp, bar_sum_names, bar_sum_regional, bar_sum_global))
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
                        value=500
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
                #dcc.Graph(id="Secondary", figure=get_local_info(colors, exon_comp, gene_comp, bar_sum_names, bar_sum_regional, bar_sum_global))
            #])
        ])
    ]),
])

def update_all_figs(x_start, x_stop, click_me_rep, click_me_genes):
    fig1, bar_sum_names, bar_sum_regional, bar_sum_global, colors = plot_interactive( int(SG_window), 
            int(SG_polynomial_order), SG_check,
            layout, exon_comp, gene_comp, bins, names_simp, chrs, zs_tmp, rep_types, click_me_rep, click_me_genes,  
            int(x_start), int(x_stop), gene_locals, gene_names)
    fig2 = get_local_info(full_cnts[x_start:x_stop], 
            names_simp[x_start:x_stop], exon_comp, gene_comp, bar_sum_regional, bar_sum_global, cds_comp)
    #And now we update the histograms 
    fig3 = create_tree(tree_file, x_start, x_stop, all_chrs[chrs])
    return fig1, fig2, fig3

#app.config.requests_pathname_prefix = app.config.routes_pathname_prefix.split('/')[-1]
last_reset_click = 0
@app.callback(
    Output('chromosome','figure'),
    Output('primary','figure'),
    Output('Secondary', 'figure'),
    Output('Third', 'figure'),
    Output('chromosome', 'relayoutData'),

    Input('chromosome','figure'),
    Input('primary','figure'),
    Input('Secondary', 'figure'),
    Input('Third', 'figure'),
    
    #Output()
    #Output('primary', 'clickData'),
    Input('my_check','value'),
    Input('num-start','value'),
    Input('num-stop','value'),
    #Input('bins','value'),
    #Input('update_check', 'value'),
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
    #Input('reset_clicks', 'n_clicks')
    )
def update_figure(chr_fig, fig1, fig2, fig3, clicker, x_start, x_stop, gene_jump, bins, n_clicks, save_clicks, SG_window, 
    SG_polynomial_order, SG_check, clickData, relayoutData, chr_relayoutData):#clicker_update, gene_jump, bins):

    click_me_genes = True
    click_me_rep = True
    if gene_jump in gene_names:
        tmp_idx = gene_names.index(gene_jump)
        x_start = int(gene_locals[tmp_idx]) - 1000
        x_stop = int(gene_locals[tmp_idx+1]) + 1000
        #fig1, fig2, fig3 = update_all_figs(x_start, x_stop, click_me_rep, click_me_genes)
        fig1, bar_sum_names, bar_sum_regional, bar_sum_global, colors = plot_interactive( int(SG_window), 
            int(SG_polynomial_order), SG_check,
            layout, exon_comp, gene_comp, bins, names_simp, chrs, zs_tmp, rep_types, click_me_rep, click_me_genes,  
            int(x_start), int(x_stop), gene_locals, gene_names)
        fig2 = get_local_info(full_cnts[x_start:x_stop], 
            names_simp[x_start:x_stop], exon_comp, gene_comp, bar_sum_regional, bar_sum_global, cds_comp)
        #And now we update the histograms 
        fig3 = create_tree(tree_file, x_start, x_stop, all_chrs[chrs])
        gene_jump = ""
    elif x_start != x_start_init and x_stop != x_stop_init:
        x_start = x_start
        x_stop = x_stop
        #fig1, fig2, fig3 = update_all_figs(x_start, x_stop, click_me_rep, click_me_genes)
        fig1, bar_sum_names, bar_sum_regional, bar_sum_global, colors = plot_interactive( int(SG_window), 
            int(SG_polynomial_order), SG_check,
            layout, exon_comp, gene_comp, bins, names_simp, chrs, zs_tmp, rep_types, click_me_rep, click_me_genes,  
            int(x_start), int(x_stop), gene_locals, gene_names)
        fig2 = get_local_info(full_cnts[x_start:x_stop], 
            names_simp[x_start:x_stop], exon_comp, gene_comp, bar_sum_regional, bar_sum_global, cds_comp)
        #And now we update the histograms 
        fig3 = create_tree(tree_file, x_start, x_stop, all_chrs[chrs])
    elif clickData != None:#len(print(clickData['points'])) > 0:
        print()
        print("Click data!!!")
        print(clickData)
        print()
        tmp_idx = gene_names.index(clickData['points'][0]['text'])
        x_start = int(gene_locals[tmp_idx]) - 1000
        x_stop = int(gene_locals[tmp_idx+1]) + 1000
        #fig1, fig2, fig3 = update_all_figs(x_start, x_stop, click_me_rep, click_me_genes)
        fig1, bar_sum_names, bar_sum_regional, bar_sum_global, colors = plot_interactive( int(SG_window), 
            int(SG_polynomial_order), SG_check,
            layout, exon_comp, gene_comp, bins, names_simp, chrs, zs_tmp, rep_types, click_me_rep, click_me_genes,  
            int(x_start), int(x_stop), gene_locals, gene_names)
        fig2 = get_local_info(full_cnts[x_start:x_stop], 
            names_simp[x_start:x_stop], exon_comp, gene_comp, bar_sum_regional, bar_sum_global, cds_comp)
        #And now we update the histograms 
        fig3 = create_tree(tree_file, x_start, x_stop, all_chrs[chrs])
    elif relayoutData != None and 'xaxis4.range[0]' in relayoutData.keys():
        print("*********")
        #print(relayoutData)
        print(relayoutData['xaxis4.range[0]'])
        print(relayoutData['xaxis4.range[1]'])
        print("*********")
        x_start = int(relayoutData['xaxis4.range[0]'])
        x_stop = int(relayoutData['xaxis4.range[1]'])
        #fig1, fig2, fig3 = update_all_figs(x_start, x_stop, click_me_rep, click_me_genes)
        fig1, bar_sum_names, bar_sum_regional, bar_sum_global, colors = plot_interactive( int(SG_window), 
            int(SG_polynomial_order), SG_check,
            layout, exon_comp, gene_comp, bins, names_simp, chrs, zs_tmp, rep_types, click_me_rep, click_me_genes,  
            int(x_start), int(x_stop), gene_locals, gene_names)
        fig2 = get_local_info(full_cnts[x_start:x_stop], 
            names_simp[x_start:x_stop], exon_comp, gene_comp, bar_sum_regional, bar_sum_global, cds_comp)
        #And now we update the histograms 
        fig3 = create_tree(tree_file, x_start, x_stop, all_chrs[chrs])
    #print(chr_relayoutData)
    elif chr_relayoutData != None and "range" in chr_relayoutData.keys():
        #print(chr_relayoutData)
        if "x2" in chr_relayoutData['range'].keys():
            x_start = int(chr_relayoutData['range']['x2'][0])
            x_stop = int(chr_relayoutData['range']['x2'][1])
        elif "x" in chr_relayoutData['range'].keys():
            x_start = int(chr_relayoutData['range']['x'][0])
            x_stop = int(chr_relayoutData['range']['x'][1])
        else:
            x_start = int(chr_relayoutData['range']['x3'][0])
            x_stop = int(chr_relayoutData['range']['x3'][1])
        #fig1, fig2, fig3 = update_all_figs(x_start, x_stop, click_me_rep, click_me_genes)
        fig1, bar_sum_names, bar_sum_regional, bar_sum_global, colors = plot_interactive( int(SG_window), 
            int(SG_polynomial_order), SG_check,
            layout, exon_comp, gene_comp, bins, names_simp, chrs, zs_tmp, rep_types, click_me_rep, click_me_genes,  
            int(x_start), int(x_stop), gene_locals, gene_names)
        fig2 = get_local_info(full_cnts[x_start:x_stop], 
            names_simp[x_start:x_stop], exon_comp, gene_comp, bar_sum_regional, bar_sum_global, cds_comp)
        #And now we update the histograms 
        fig3 = create_tree(tree_file, x_start, x_stop, all_chrs[chrs])
        chr_relayoutData = None#{'autosize': True}
    print(chr_relayoutData)
    if x_start_init != x_start:
        chr_fig = plot_chr_whole(x, y, z_1, z_9, z_genes, z_cnts, x_start, x_stop)
    if save_clicks > 0:
        fig.write_image("fig1.pdf")
    return chr_fig, fig1, fig2, fig3, chr_relayoutData #, clickData

@app.callback(Output('textarea-state-button', 'n_clicks'),
    Output('num-start','value'),
    Output('num-stop','value'),
    [Input('reset_clicks','n_clicks')])
def update(reset):
    print("Reset:")
    print(reset)
    return 0, x_start_init, x_stop_init

app.run_server(host='127.0.0.1', debug=True)





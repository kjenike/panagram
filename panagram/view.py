import sys
import csv
import pysam
import os.path
from os import path
from mycolorpy import colorlist as mcp
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
import dash
import time
from Bio import Phylo
from scipy import signal
from dash import Dash, dcc, html, Input, Output, ctx, State, no_update
import math
from io import StringIO
from scipy.cluster.hierarchy import linkage
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist, squareform
import plotly.figure_factory as ff
from .index import Index 

def view(params):
    SG_window =53
    poly_order = 3
    SG_polynomial_order = 3
    SG_check = [1]
    #window_size = 200000

    #n_skips_start = 100
    buff = 1000

    index = Index(params.index_dir) #Directory that contains the anchor direcotry

    anchor_name, chrs = index.chrs.index[0]
    if params.genome is not None:
        anchor_name = params.genome
    if params.chrom is not None:
        chrs = params.chrom
    x_start_init = 0 if params.start is None else params.start
    x_stop_init  = 1000000 if params.end is None else params.end
    bins = params.max_chr_bins
    opt_bed_file = params.bookmarks
    window_size = index.chr_bin_kbp*1000
    n_skips_start = index.lowres_step

    num_samples = len(index.genomes)
    kmer_len = index.k
    rep_list = index.gff_anno_types
    labels = list(index.genomes)

    sns.set_palette('viridis', num_samples)
    colors = mcp.gen_color(cmap="viridis_r",n=num_samples)

    #opt_bed_file = "None"#sys.argv[3]
    #Read in the config file:
    #with open(config_f, "r") as f:
    #    line1 = f.readline()
    #    while line1:
    #        num_samples = int(line1.strip().split(':')[1])
    #        line2 = f.readline()
    #        kmer_len = int(line2.strip().split(':')[1])

            #line4 = f.readline()
            #buff = int(line4.strip().split(':')[1])

    #        line5 = f.readline()
    #        x_start_init = int(line5.strip().split(':')[1])

    #        line6 = f.readline()
    #        x_stop_init = int(line6.strip().split(':')[1])

    #        line7 = f.readline()
    #        bins = int(line7.strip().split(':')[1])

    #        line8 = f.readline()
    #        rep_types_file = line8.strip().split(':')[1]

            #line9 = f.readline()
            #mash_filenames = line9.strip().split(':')[1]

    #        line10 = f.readline()
    #        mash_edges = line10.strip().split(':')[1]

    #        line11 = f.readline()
    #        bit_file_prefix = line11.strip().split(':')[1]

    #        line12 = f.readline()
    #        anchor_name = line12.strip().split(':')[1]

    #        line13 = f.readline()
    #        opt_bed_file = line13.strip()

    #        line1 = f.readline()


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

    def get_init_rep_types(rep_types_f):
        #rep_types = {}
        rep_list = []
        #for i in range(0, num_chrs):
        #    rep_types["chr" + str(i+1)] = {}
        with open(rep_types_f, "r") as f:
            line = f.readline()
            while line:
                rep_list.append(line.strip())
                #for i in range(0, num_chrs):
                #    rep_types["chr" + str(i+1)][line.strip()] = []
                line = f.readline()
        return  rep_list

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

    def perc_universal(start,stop, countme, n_skips, this_chr):
        univ = ((index.query_bitmap(anchor_name, this_chr, start, stop, 1).sum(axis=1))==countme).sum()
        return float(univ/(stop-start))*100

    def read_gene_annotations(anchor_name):
        gene_content = {}
        gene_anns = {}
        for chrom in index.chrs.loc[anchor_name].index:
            gene_content[chrom] = {"Names":[],"Universal":[], "Unique":[]}
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

                line = f.readline()
        #genes = idx.query_genes("shigella", "NC_004337.2", 0, 10000)
        #bounds = genes.loc[:, ["start", "end"]]
        #bounds["break"] = None #pd.NA
        #print(bounds.to_numpy().flatten())
        return gene_locals, gene_names, exon_locals, exon_names, gene_content, gene_anns

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
            #line_color = "grey"
            line_color = color_code[tmp_name]#palette[int(((kmer_num[tmp_name])/total_kmers)*100)-1]
        elif clade.clades:
            
            line_color = palette[int(biggest_num_in_clade(clade, kmer_num))+10]
            #line_color = "grey"
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

    def create_tree(x_start_init, x_stop_init, raw_counts, n_skips):
        #Adapted from:
        #https://github.com/plotly/dash-phylogeny

        tree_tmp1 = raw_counts 
        #We need to sum up the number of times that kmers occur in each column (aka, each sample)
        kmer_num_tmp = tree_tmp1.sum(axis=0)
        
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
                #this_color = "grey"
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

    def get_local_info(x_all, bar_sum_regional, bar_sum_global_tmp, anchor_name, chrs):
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
        y=[(i/sum(bar_sum_regional[1:])*100) for i in bar_sum_regional[1:]]
        y_whole=[(i/sum(bar_sum_global_tmp)*100) for i in bar_sum_global_tmp]
        #.sum(axis=1)
        fig.add_trace(go.Bar(x=x, y=[a_i - b_i for a_i, b_i in zip(y, y_whole)], marker_color=colors, showlegend=False), row=2, col=1)
        fig.add_trace(go.Bar(x=x, y=y_whole, marker_color=colors, showlegend=False), row=1, col=1)
        #Genes
        #y=[(i/sum(gene_comp[1:])*100) for i in gene_comp[1:]]
        #fig.add_trace(go.Bar(x=x, y=y_whole, marker_color="#7400b8", showlegend=False), row=1, col=4)
        totals = 0
        gene_comp_tmp = []
        gene_comp = []
        for i in range(1, num_samples+1):
            totals += index.chrs.loc[(anchor_name, chrs), "gene_occ_"+str(i)]
            gene_comp_tmp.append(index.chrs.loc[(anchor_name, chrs), "gene_occ_"+str(i)])
        
        for t in gene_comp_tmp:
            gene_comp.append(t/totals*100)
        fig.add_trace(go.Bar(x=x, y=[a_i - b_i for a_i, b_i in zip(gene_comp, y_whole)], marker_color=colors, showlegend=False), row=2, col=2)
        #fig.update_layout(xaxis_title_text="K-mers shared in X samples", yaxis_title_text='Frequency (log)')
        fig.update_xaxes(title_text="# of genomes", row=2, col=1)
        fig.update_yaxes(title_text="Difference from whole chromosome", row=2, col=1)
        fig.update_yaxes(title_text="Percent of k-mers", row=1, col=1)
        fig.update_layout(height=1000)
        return fig

    def plot_interactive( n_skips, #window_filter, poly_order, shared_kmers, 
        layout, bins, names_simp, name, zs_tmp, 
        plot_rep, plot_gene, x_start, x_stop, gene_locals, gene_names, anchor_name, chrs): #exon_locals, exon_names):
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
            ],
            subplot_titles=("Ref. Sequence Position","Gene Annotation", "Annotation",  "Conserved K-mers" )
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
        #print("a1:", time.time()-t)
        #t = time.time()
        
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
            if (cntr % adjusted_bin_size_init) == 0 : #Is the x position within the same bin? Or do we need to move to the next bin? 
                x_tracker += 1 #X-tracker keeps track of which bin we are using (on the x-axis)
            
            #if (i) < len(cats_tmp) and x_tracker < len(cats_tmp[i]):
            if (i) < len(cats_tmp) and x_tracker < len(cats_tmp[i]):
                cats_tmp[i][x_tracker] += 1 
            cntr += 1 #n_skips

        #print("b2:", time.time()-t)
        #t = time.time()
        plot_rep = True
        #Add a line plot that will cover the different repetative elements. 
        if plot_rep == True:
            #print("About to add repeats trace")
            rep_colors = ["#f0f921", "#f8df25", "#fdc627", "#fdaf31", "#f99a3e", "#f3854b", "#e97257", "#de6164", "#d24f71", "#c43e7f", 
                    "#b42e8d", "#a21d9a", "#8e0ca4", "#7801a8", "#6100a7", "#4903a0", "#2f0596", "#0d0887", "grey", "grey", "grey"]
            cntr = 0
            #df = index.query_anno(anchor_name, chrs, x_start, x_stop)
            #print(rep_list)
            for i in rep_list: #rep_types.keys():
                df = index.query_anno(anchor_name, chrs, x_start, x_stop)
                #print(df)
                if i == "exon":
                    #print("EXON!")
                    exon_y   = []
                    exon_tmp = []
                    df = df[df["type"]==i]
                    bounds = df.loc[:, ["start", "end"]]
                    bounds["break"] = None #pd.NA
                    exon_tmp = bounds.to_numpy().flatten()
                    e = 0
                    while e < len(exon_tmp):
                        exon_y.append(0.5)#[cntr]*len(rep_types_tmp)
                        exon_y.append(0.5)
                        exon_y.append(None)
                        e += 3
                    fig.add_trace(go.Scattergl(x=exon_tmp, y=exon_y, 
                        line=dict(color="#a0da39", width=5), name=i), row=2, col=1)
                else:
                    rep_y = []
                    rep_types_tmp = []
                    #df = index.query_anno(anchor_name, chrs, x_start, x_stop)
                    df = df[df["type"]==i]
                    bounds = df.loc[:, ["start", "end"]]
                    bounds["break"] = None #pd.NA
                    anno_locals_tmp = bounds.to_numpy().flatten()
                    r = 0
                    while r < len(anno_locals_tmp):
                        rep_y.append(cntr)#[cntr]*len(rep_types_tmp)
                        rep_y.append(cntr)
                        rep_y.append(None)
                        r += 3
                    fig.add_trace(go.Scattergl(x=anno_locals_tmp, y=rep_y, 
                        line=dict(color=rep_colors[cntr]), name=i, legendgroup="group2", 
                        legendgrouptitle_text="Annotations"), row=4, col=1)
                    cntr += 1
            fig.update_yaxes(visible=False, row=4, col=1)
            fig.update_xaxes(showticklabels=False, row=4, col=1)
        #print("a2.2:", time.time()-t)
        #t = time.time()
        if plot_gene == True:
            gene_locals_tmp = []
            gene_names_tmp = []
            intron_locals_tmp = []
            i = 0
            while i < len(gene_names):
                gene_names_tmp.append(gene_names[i])
                gene_names_tmp.append(gene_names[i])
                gene_names_tmp.append(None)
                i += 1
            gene_y = [2]*len(gene_names_tmp)

            fig.add_trace(go.Scattergl(x=gene_locals, y=gene_y, 
                line=dict(color="#3b528b", width=10), 
                text=gene_names_tmp, hovertemplate='<br>x:%{x}<br>m:%{text}', legendgroup="group2", 
                name="Gene"), row=2, col=1)
            fig.update_layout(clickmode='event+select')
            #Now add the exons: 
            i = 0

            fig.update_yaxes(visible=False, range=[-1,4], row=2, col=1)
            fig.update_xaxes(showticklabels=False, row=2, col=1)
        
        #print("a2.1:", time.time()-t)
        #t = time.time()
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
        #print("a2.3:", time.time()-t)
        #t = time.time()
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
        #print("a2:", time.time()-t)
        #t = time.time()

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
        fig.update_yaxes(title_text="# of k-mers", range=[0,adjusted_bin_size_init]  , row=6, col=1)
        fig.update_layout(height=1000, xaxis_range=[x_start,x_stop], font=dict(
            #family="Balto",

            size=16,
            ))
        #print("a3:", time.time()-t)
        #t = time.time()
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
        return x, y, z_1, z_9, z_genes #, z_cnts #x_all, y_all, z_all,  x, z_cnts

    def read_chr_whole(this_bin_file):
        with open(this_bin_file, "r") as f:
            line1 = f.readline()
            line2 = f.readline()
            line3 = f.readline()
            cntr = 0
            x = {}
            z_9 = {}
            z_1 = {}
            y = {}
            while line1:
                x_tmp = line1.strip().split(":")[1].split(",")[:-1]
                z_9_tmp = line2.strip().split(",")[:-1]
                this_chr_name = line1.strip().split(":")[0]
                x_1_tmp = line3.strip().split(",")[:-1]
                x[this_chr_name] = []
                z_9[this_chr_name] = [] #.append([])
                z_1[this_chr_name] = [] #.append([])
                y[this_chr_name] = ([1]*len(x_tmp))
                for i in range(0, len(x_tmp)):
                    x[this_chr_name].append(int(x_tmp[i]))
                    z_9[this_chr_name].append(int(z_9_tmp[i]))
                    z_1[this_chr_name].append(int(x_1_tmp[i]))
                line1 = f.readline()
                line2 = f.readline()
                line3 = f.readline()
                cntr += 1
        
        return x, z_1, z_9, y

    def make_gene_whole_chr(x, locs):
        z_genes = [0]*(len(x)+2)
        g = 0
        while g < len(locs):
            #print(locs[g])
            x_bin = int(int(locs[g])/200000)
            if x_bin < len(z_genes):
                z_genes[x_bin] += 1
            g += 1
        return z_genes

    def plot_chr_whole( x_start, x_stop, anchor_name, this_chr): 
        
        #print(index.chr_bins.loc[(anchor_name, "total_occ_1")][this_chr])

        #x = index.chr_bins.loc[(anchor, this_chr , "total_occ_1")]
        z_1 = index.chr_bins.loc[(anchor_name, "total_occ_1")][this_chr]
        z_9 = index.chr_bins.loc[(anchor_name, "total_occ_"+str(num_samples))][this_chr]
        y, x = [], []
        cntr = 0
        for xtmp in z_1:
            x.append(cntr)
            y.append(1)
            cntr += window_size
        #print(z_9)
        z_genes = [0]*len(x)
        if index.gene_tabix[anchor_name] is not None:
            try:
                genes = index.query_genes(anchor_name, this_chr, 0, index.chrs.loc[anchor_name, this_chr]["size"])
            except ValueError:
                print("Skipping " + l + "\t" + c)
                genes = None

            if genes is not None:
                bounds = genes["start"].to_numpy() #genes.loc[:, ["start"]].to_numpy()
                #print(bounds)
                z_genes = make_gene_whole_chr(x, bounds) #[g.split(';')[1].split('=')[1] for g in genes['attr']]) 
                
                #print(z_genes)
                #bounds["break"] = None
        #print(genes)
        #print(this_chr)
        #print(anchor_name)
        chr_fig = make_subplots(rows=3, cols=1, 
            specs=[[{"type": "heatmap",}], [{"type": "heatmap",}], [{"type": "heatmap",}]],
            shared_xaxes=True,
            subplot_titles=("K-mer and gene density accross whole chromosome", "",
                ""),
            vertical_spacing = 0.0
            )

        chr_fig.add_trace(go.Heatmap(x=x, z=z_9, y=y, type = 'heatmap', colorscale='magma_r', showlegend=False, showscale=False), row=1, col=1)
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
        chr_fig.add_trace(go.Heatmap(x=x, z=z_genes, y=y, type = 'heatmap', colorscale='magma_r', showscale=False ), row=3, col=1)
        chr_fig.add_trace(go.Scatter(x=[x_start, x_start, None, x_stop, x_stop, None, x_start, x_stop,], showlegend=False, 
                       y=[0.5, 1.5, None, 0.5, 1.5, None, 0.55, 0.55],
                       mode='lines',
                       line_color='#1dd3b0', line_width=8), row=3, col=1)
        
        chr_fig.update_yaxes( range=[0.5,1.5], showticklabels=False, row=1, col=1)
        chr_fig.update_yaxes( range=[0.5,1.5], showticklabels=False, row=2, col=1)
        chr_fig.update_yaxes( range=[0.5,1.5], showticklabels=False, row=3, col=1)
        
        #chr_fig.update_yaxes(visible=False, range=[0,3], row=4, col=1)

        chr_fig.update_xaxes(fixedrange=True, range=[0,index.chrs.loc[anchor_name, this_chr]["size"]], row=1, col=1)
        chr_fig.update_xaxes(fixedrange=True, range=[0,index.chrs.loc[anchor_name, this_chr]["size"]], row=2, col=1)
        chr_fig.update_xaxes(fixedrange=True, row=3, col=1)
        chr_fig.update_xaxes(title_text="Sequence position", range=[0,index.chrs.loc[anchor_name, this_chr]["size"]], row=3, col=1)
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

    def plot_whole_genome(anchor_name):
        spec = []
        sub_titles = []
        h = num_chrs[anchor_name]*250
        for chrom in index.chrs.loc[anchor_name].index: #i in range(0, num_chrs[anchor_name]):
            spec.append([{"type": "heatmap",}])
            spec.append([{"type": "heatmap",}])
            spec.append([{"type": "heatmap",}])
            sub_titles.append(chrom)
            sub_titles.append("")
            sub_titles.append("")
        
        wg_fig = make_subplots(rows=(3*num_chrs[anchor_name]), cols=1, 
            specs=spec, #[[{"type": "heatmap",}], [{"type": "heatmap",}], [{"type": "heatmap"}]],
            shared_xaxes=True,
            subplot_titles=sub_titles, #("K-mer and gene density accross whole chromosome", "", ""),
            vertical_spacing = 0.0
            )
        cntr = 1
        for chrom in index.chrs.loc[anchor_name].index:
            x = list(index.chr_bins.loc[(anchor_name, chrom)].index)
            #x.append(1)
            #print(chrom)
            #print("x:")
            #print(x)
            #print(len(x))
            #print("y:")
            #print([1]*len(x))
            #print("z:")
            #z = list(index.chr_bins.loc[(anchor_name, "total_occ_"+str(num_samples))][chrom])
            #z.append(0)
            #print(z)
            if len(x)!=1:
                wg_fig.add_trace(go.Heatmap(x=x, z=index.chr_bins.loc[(anchor_name, "total_occ_"+str(num_samples))][chrom], 
                    y=[1]*(len(x)), type = 'heatmap', colorscale='magma_r', showlegend=False,showscale=False), row=((cntr*3)-2), col=1)
                wg_fig.add_trace(go.Heatmap(x=x, z=index.chr_bins.loc[(anchor_name, "total_occ_1")][chrom], 
                    y=[1]*(len(x)), type = 'heatmap', colorscale='magma', showscale=False), row=((cntr*3)-1), col=1)
            
            if cntr == 1:
                wg_fig.update_layout(xaxis={'side': 'top'}) 

            #wg_fig.update_yaxes( range=[0.5,1.5], showticklabels=False, row=((cntr*3)-2), col=1)
            #wg_fig.update_yaxes( range=[0.5,1.5], showticklabels=False, row=((cntr*3)-1), col=1)
            
            #wg_fig.update_xaxes(fixedrange=True, row=((cntr*3)-2), col=1)
            #wg_fig.update_xaxes(fixedrange=True, row=((cntr*3)-1), col=1)
            
            #wg_fig.update_xaxes(title_text="Sequence position", row=((cntr*3)), col=1)
            #wg_fig.update_yaxes(title_text="Universal", row=((cntr*3)-2), col=1)
            #wg_fig.update_yaxes(title_text="Unique", row=((cntr*3)-1), col=1)
            cntr += 1
        
        wg_fig.update_layout(clickmode='event', plot_bgcolor='rgba(0,0,0,0)')
        wg_fig.update_layout(height=h)
        #wg_fig.update_layout(margin=dict(
        #        b=10,
        #        l=10,
        #        r=10),
        #    font=dict(
        #        #family="Balto",
        #        size=14,
        #    ),
        #    paper_bgcolor='rgba(0,0,0,0)')
        return wg_fig

    def plot_gene_content(gene_names,universals,uniques,sizes, sort_by, colors, uniq_avg, univ_avg, x_start, x_stop, anchor_name, chrs):
        #print(gene_content)
        colors = ['#ffd60a', '#440154']
        d = {
            "Names":gene_names,
            "Universal": universals/sizes,
            "Unique":uniques/sizes
        }
        uniq_avg = uniq_avg/100
        univ_avg = univ_avg/100
        df = pd.DataFrame(d)
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
        local_genes = index.query_genes(anchor_name, chrs, x_start, x_stop)
        local_gene_list = [g.split(';')[0].split("=")[1] for g in local_genes['attr']]
        df2 = df_sorted.loc[df_sorted['Names'].isin(local_gene_list)]

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

    #Tab 1 bar chart of each genome 
    def read_pangenome_comp():
        genome_comp_totals = {}
        genome_names = []
        cats = []
        for i in range(0,num_samples):
            cats.append([])

        cols = ["total_occ_"+str(i) for i in range(1,num_samples+1)]
        counts = index.chrs[cols].groupby(level="genome").sum()
        fracs = counts.divide(counts.sum(axis=1), axis=0) 

        for l in labels: 
            genome_comp_totals[l] = []
            for c in cols:
                genome_comp_totals[l].append(fracs.loc[l,c])

        for i in range(0, num_samples):
            for l in labels:
                cats[i].append(genome_comp_totals[l][i])
        #print(cats)
        fig = make_subplots(rows=1, cols=1)
        for g in range(0,len(labels)):
            fig.add_trace(go.Bar(y=labels, x=cats[g], name=str(g+1), #orientation='h',
                #legendgroup="group1", 
                legendgrouptitle_text="# Samples",
                marker=dict(color=colors[g]), 
                marker_line=dict(color=colors[g]), orientation='h',
                #log_x=True
                ), row=1, col=1)
        fig.update_layout(barmode='stack' )#,orientation='h')
        #fig.update_xaxes(type="log")
        fig.update_layout(height=1500, font=dict(
            size=26,
            ), plot_bgcolor='rgba(0,0,0,0)')
        return fig, genome_comp_totals

    def read_genome_comp(anchor_name):
        chrs_comp = []
        for c in index.chrs.loc[anchor_name].index:
            #print(c)
            tmp = []
            for n in range(1, num_samples+1):
                tmp.append(index.chrs.loc[anchor_name, c]["total_occ_"+str(n)])
            chrs_comp.append(tmp)
        fig = make_subplots(rows=len(chrs_comp), cols=1)
        x = []
        for i in range(1, num_samples+1):
            x.append(str(i))
        
        for i in range(0,len(chrs_comp)):
            #for chrom in 
            y=[]
            tots = sum(chrs_comp[i])
            if tots > (window_size/n_skips_start):
                #print("WHAT?!!!")
                for j in chrs_comp[i]:
                    y.append(int(j)/tots*100)
                fig.add_trace(go.Bar(x=x, y=y, marker_color=colors, showlegend=False,), #marker_color=colors[1:len(x)]), 
                    row=i+1,col=1)
        fig.update_yaxes(type="log")
        fig.update_layout(paper_bgcolor='rgba(0,0,0,0)')
        return fig, chrs_comp

    #Tab 1 dendogram 
    def make_all_genome_dend(labels_og):
        #file_names = {}
        #with open(mash_filenames) as f:
        #    for line in f:
        #        line = line.rstrip()
                
        #        file_names[line.split('\t')[0]] = line.split('\t')[0]
        sample_list = labels_og #list(file_names.values())
        dim = len(sample_list)
        dist_mat = np.zeros((dim, dim), np.float64)

        with open(index.genome_dist_fname) as f:
            for line in f:
                f, t, d, p, x = line.rstrip().split("\t")
                #print(f)
                #print(t)
                i = sample_list.index(f)
                j = sample_list.index(t)
                dist_mat[i][j] = d
                dist_mat[j][i] = d

                #print(f)
                #if f in file_names:
                #    new_f = file_names[f]
                #    new_t = file_names[t]
                #    
                #    i = sample_list.index(new_f)
                #    j = sample_list.index(new_t)
                #    dist_mat[i][j] = d
                #    dist_mat[j][i] = d

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
        fig.update_layout({"autosize":True,#'width':1500, 'height':1500,
                                 'showlegend':False, 'hovermode': 'closest',
                                 'paper_bgcolor':'rgba(0,0,0,0)',
                                 'plot_bgcolor':'rgba(0,0,0,0)'
                                 })
        fig.update_yaxes(
            scaleanchor = "x",
            scaleratio = 1,
        ) 
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
        fig.update_layout(yaxis={'domain': [0, 0.8],#
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
        seqs = []
        sizes = []
        for l in labels:
            seqs.append([l, len(index.chrs.loc[l].index)])
            sizes.append([l, index.chrs.loc[l]["size"].sum()])

        sizes_df = pd.DataFrame(sizes, columns=['Sample','Size'])
        seqs_df = pd.DataFrame(seqs, columns=['Sample','Seq'])
        
        sizes_df = sizes_df.sort_values(by=['Size'])
        seqs_df = seqs_df.sort_values(by=['Seq'])
        return sizes_df, seqs_df

    def make_genome_size_plots(anchor_name, sizes_df, seqs_df):
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

    def find_avgs(genome_comp_totals):
        #print(genome_comp_totals)
        
        tmp = []
        for k in genome_comp_totals.keys():
            running_sum = 0
            total = 0
            for i in range(1, len(genome_comp_totals[k])):
                total += float(genome_comp_totals[k][i-1])
                running_sum += (float(genome_comp_totals[k][i-1])*i)
            #print("*******************")
            tmp.append([k,(running_sum/total)])
        #print(tmp)
        df = pd.DataFrame(tmp, columns=['Sample','Avg'])
        #print(df)
        df = df.sort_values(by=['Avg'])
        #print(df)
        return df


    def plot_avgs(anchor_name, df):
        fig = make_subplots(rows=1, cols=1)
        fig.add_trace(go.Scattergl(x=df['Sample'],y=df['Avg']))
        fig.add_vline(x=anchor_name, line_dash="dash", line_color="darkblue")
        #fig.add_vline(x=anchor_name)
        fig.update_yaxes(title_text="Average k-mer",)
        fig.update_layout(font=dict(size=20,))
        #fig.add_trace(go.Scattergl(x=list(genome_comp_totals.keys()),y=tmp[]))
        return fig

    def make_avg_kmer_fig(chrs_comp, this_anchor):
        fig = make_subplots(rows=1, cols=1)
        x = []
        y = []
        for c in range(0, num_chrs[this_anchor]):#len(chrs_comp)):
            #print(c)
            x.append(chrs_list[this_anchor][c])
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

    def make_genes_per_chr_fig(anchor_name):
        fig = make_subplots(rows=1, cols=1)
        x = []
        y = []
        for c in range(0, len(chrs_list[anchor_name])):
            this_chrom = chrs_list[anchor_name][c]
            x.append(this_chrom)
            try:
                y.append(len(index.query_genes(anchor_name, this_chrom, 0, index.chrs.loc[anchor_name, this_chrom]["size"]))/index.chrs.loc[anchor_name, this_chrom]["size"])
            except ValueError:
                print(anchor_name + "\t" + this_chrom)
        fig.add_trace(go.Scattergl(x=x, y=y))
        fig.update_yaxes(title_text="Gene density",)
        fig.update_xaxes(title_text="Chromosome",)
        fig.update_layout(plot_bgcolor='rgba(0,0,0,0)',  font=dict(size=20))
        return fig

    def make_gene_per_genome_fig(anchor_name):
        
        fig = make_subplots(rows=1, cols=1)
        colors = ['#ffd60a', '#440154']
        gene_names = []
        universals = []
        uniques = []
        sizes = []
        genes = list()
        #print("About to iterate over chromosomes")
        for chrom in index.chrs.loc[anchor_name].index:
            #print(chrom)
            try:
                g = index.query_genes(anchor_name, chrom, 0, index.chrs.loc[anchor_name, chrom]["size"])
            except ValueError: 
                #print("Theres a value error?!")
                continue
            genes.append(g)
        #print("About to concat")
        genes = pd.concat(genes)
        #print(genes)
        genes["name"] = genes["attr"].str.extract("ID=([^;]+)")
        genes["size"] = genes["end"] - genes["start"]
        #print("Finished the name and size stuff")

            #gene_names += [g.split(';')[0].split("=")[1] for g in genes['attr']]
            #universals += genes["universal"]
            #uniques += genes["unique"]
            #print(genes["end"]-genes["start"])
            #sizes += genes["end"]-genes["start"]
        #d = {
        #        "Names":gene_names,
        #        "Universal":universals/sizes,
        #        "Unique":uniques/sizes,
        #        #"Sizes": sizes
        #}
        
        #print(genes)
        #df = pd.DataFrame(d)
        x = [i for i in range(0, len(genes['universal']))]
        #print("Made x axis")
        genes['universal'] = genes['universal']/genes["size"]
        #print("Updated universal")
        genes['unique'] = genes['unique']/genes["size"]
        #print("Updated unique")
        df_sorted = genes.sort_values('universal')
        #print("sorted the dataframe")
        df_sorted['X'] = x
        #print(df_sorted)
        fig.add_trace(go.Scattergl(x=x, y=df_sorted['unique'], text=df_sorted["chr"] + ":" + df_sorted['name'], marker=dict(color=colors[0]),
                name="Proportion "+"Unique", mode="markers"))
        fig.add_trace(go.Scattergl(x=x, y=df_sorted['universal'], text=df_sorted["chr"] + ":" +df_sorted['name'], marker=dict(color=colors[1]),
                name="Proportion "+"Universal", mode="markers"))
        fig.update_layout(plot_bgcolor='rgba(0,0,0,0)', font=dict(size=20) )
        fig.update_xaxes(title_text="Genes",)
        fig.update_yaxes(title_text="Proportion of gene",)
        fig.update_layout(hovermode='x unified')
        return fig

    ##### READ DATA
    #print("Reading data!")


    #print(labels)
    chrs_list = {}
    chr_lens = {}
    num_chrs = {}
    cntr = 0
    pre_bed_info = {}
    #chr_nums = 12
    for l in labels:
        num_chrs[l] = len(index.chrs.loc[l])
        chrs_list[l] = index.chrs.loc[l].index
        #chr_lens[l] = index.chrs.loc[l]["size"]
        pre_bed_info[l] = []
        #print(index.chrs.loc[l]["size"])
    t = time.time()
    cnts_tmp = []
    x = []
    all_chrs = {} #This is where we keep the raw, str, counts 
    #rep_list = get_init_rep_types()
    
    if opt_bed_file != None: 
        #The optional bed file was provided 
        with open(opt_bed_file, "r") as f:
            line = f.readline()
            while line:
                tmp = line.strip().split("\t")
                pre_bed_info[tmp[0]].append(tmp[1] + ":" + tmp[2] + "-" + tmp[3])
                line = f.readline()


    for i in range(1, num_chrs[anchor_name]+1):
        chrs=index.chrs.loc[anchor_name].index[i-1] #"chr" + str(i)
        chr_start = 0

        chr_end = index.chrs.loc[anchor_name, chrs]["size"]#chr_lens[anchor_name][i-1] #308452471
        #print(chr_end)
        #print(anchor_name + "." + chrs)
        all_chrs[chrs] = index.query_bitmap(anchor_name, chrs, chr_start, chr_end, n_skips_start)
        #all_chrs[chrs] = anchor.get_counts(chrs, chr_start, chr_end, n_skips_start)
    chrs = index.chrs.loc[anchor_name].index[0] #"chr1"
    #print("a0:", time.time()-t)
    #t = time.time()
    #chrs = "chr1"
    #Parsing the counts will return the number of genomes present at each kmer position, based on the counts
    names_simp = all_chrs[chrs].sum(axis=1)#tmp.sum(axis=1)
    bar_sum_global = {}

    for l in labels:
        bar_sum_global[l] = {}
        cntr = 0

        for c in index.chrs.loc[l].index: 
            counts = index.chr_bins.loc[(l,c)].sum(axis=0)
            bar_sum_global[l][c] = counts.to_numpy()#bar_sum_global_tmp
            cntr +=1

    tmp = 0
    
    #print("a0.2:", time.time()-t)
    #t = time.time()
    
    #Read in the repeats file (this has the annotated repeats from Shujun)     
    #print("a0.4:", time.time()-t)
    #t = time.time()
    #grep -v "#" Solqui2_genes_1.0.0.gff | awk '{if($1=="chr1") print $0}' | awk '{if($3=="gene") print $1"\t"$9"\t"$4"\t"$5}' > Squi2.genes.chr1.txt
    zs_tmp = []
    gene_names = {}
    gene_locals = {}
    exon_locals = {}
    exon_names = {}
    all_rep_types = {}
    #gene_comp = {}
    gene_content = {}
    #gene_anns = {}

    #print("a0.45:", time.time()-t)
    #t = time.time()
    def get_idx_content(gene_locals, gene_cntr, j):
        est_pos_start = int(gene_locals[gene_cntr]/n_skips_start)
        est_pos_stop = int(gene_locals[gene_cntr+1]/n_skips_start)
        tmp = (names_simp[est_pos_start:est_pos_stop]==j).sum()
        return tmp

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
    x_start_init_adjust = int(x_start_init/n_skips_start)
    x_stop_init_adjust = int(x_stop_init/n_skips_start)
    #### NEED TO ADJUST THE START INIT AND STOP INIT SO THAT WE TAKE THE SKIPS INTO ACCOUNT
    #Here is where we make the main plot 
    genes = index.query_genes(anchor_name, chrs, 0, index.chrs.loc[anchor_name, chrs]["size"])
    #[g.split(';')[1].split('=')[1] for g in genes['attr']]
    bounds = genes.loc[:, ["start", "end"]]
    bounds["break"] = None #pd.NA
    #gene_locals[l][chrom] = bounds.to_numpy().flatten()
    gene_names = [g.split(';')[0].split("=")[1] for g in genes['attr']]
    main_fig, bar_sum_names, bar_sum_regional, colors, gene_names_tmp = plot_interactive( n_skips_start, #window_filter, poly_order, shared_kmers,
        layout, #exon_comp["chr1"], 
        bins, names_simp[x_start_init_adjust:x_stop_init_adjust], chrs, zs_tmp, 
        True, True, 
        x_start_init, x_stop_init, bounds.to_numpy().flatten(), #gene_locals[anchor_name][chrs], 
        gene_names, anchor_name, chrs#gene_names[anchor_name][chrs], 
        #exon_locals[anchor_name][chrs], exon_names[anchor_name][chrs],
        )
    local_genes = len(index.query_genes(anchor_name, chrs, x_start_init, x_stop_init))
    #local_bounds = local_genes.loc[:, ["start", "end"]]
    #local_bounds["break"] = None
    sort_by = ["Unique","Universal"]
    tmp = [(i/sum(bar_sum_global[anchor_name][chrs])*100) for i in bar_sum_global[anchor_name][chrs]]
    uniq_avg = tmp[1]
    univ_avg = tmp[-1]
    #genes = index.query_genes(anchor_name, chrs, 0, index.chrs.loc[anchor_name, chrs]["size"])
    #gene_names = [g.split(';')[1].split("=")[1] for g in genes['attr']]
    universals = genes["universal"]
    uniques = genes["unique"]
    #print(genes["end"]-genes["start"])
    sizes = genes["end"]-genes["start"] #genes["size"]
    gene_content_plot = plot_gene_content(gene_names,universals,uniques,sizes,sort_by, 
        colors, uniq_avg, univ_avg, x_start_init, x_stop_init, anchor_name, chrs)
    #print("AFTER GENE CONTENT", len(gene_locals["chr1"]))
    #sample_kmer_counts = parse_counts_complex(all_chrs[chrs][x_start_init:x_stop_init])
    #x, y, z_1, z_9, z_genes = make_chr_whole(names_simp, 1000, gene_locals, x_start_init, x_stop_init)
    x, z_1, z_9, y = {}, {}, {}, {}
    chr_fig = plot_chr_whole(x_start_init, x_stop_init, anchor_name, chrs)
    whole_genome_fig = plot_whole_genome(anchor_name)
    #chr_fig = plot_chr_whole(names_simp, full_cnts, 1000, gene_locals, x_start_init, x_stop_init)
    #print("a4:", time.time()-t)
    #t = time.time()

    #Now we need to make the pangenome plot
    #This is going to be a bar chart 

    whole_genome_hists, chrs_comp = read_genome_comp(anchor_name)
    all_genomes_dend, ordered_labels = make_all_genome_dend(labels)
    pangenome_comp, genome_comp_totals = read_pangenome_comp()
    tmp_df = find_avgs(genome_comp_totals)
    pangenome_avg = plot_avgs(anchor_name, tmp_df)
    pg_sizes, pg_num_seqs = read_genome_size_files()
    pangenome_sizes, pangenome_num_seqs = make_genome_size_plots(anchor_name, pg_sizes, pg_num_seqs)
    avg_kmer_per_chr_fig = make_avg_kmer_fig(chrs_comp, anchor_name)
    genes_per_chr_fig = make_genes_per_chr_fig(anchor_name)
    gene_content_per_genome_fig = make_gene_per_genome_fig(anchor_name)
    #pangenome_avg = 
    #print(all_chrs[chrs])
    #Set up the dash app 
    chrs = index.chrs.loc[anchor_name].index[0] #"chr1"
    #dbc_css = "https://cdn.jsdelivr.net/gh/AnnMarieW/dash-bootstrap-templates/dbc.min.css"
    app = dash.Dash(
        external_stylesheets=["https://www.w3schools.com/w3css/4/w3.css"],
        url_base_pathname=params.url_base
    ) 
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
                        dcc.Dropdown(labels, anchor_name, style=dict(width='40%', height='110%', verticalAlign="middle"), id="pangenome_tab_dropdown"),
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
                        html.I("Anchor genome " + anchor_name, id="anchor_labels"),
                        html.Br(),
                        html.Label('Select a chromosome: '),
                        dcc.Dropdown(chrs_list[anchor_name], chrs, style=dict(width='40%', height='110%', verticalAlign="middle", ), id="Anchor_tab_dropdown"),
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
                                    figure = whole_genome_fig, config=config, style={"height" : num_chrs[anchor_name]*250})
                        ]),
                    html.Div(className="w3-quarter",children=[
                                dcc.Graph(id="all_chromosomes_hists", style={"height" : num_chrs[anchor_name]*250},
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
                        #html.I(anchor_name + "."), # + chrs + ":" ), #+ str(x_start_init) + "-" + str(x_stop_init)), #"Chromosome " + chrs,
                        #html.I()
                        html.I(anchor_name + "." + chrs + ":", id="chr_name"),
                        dcc.Input(id="Chrs_Info", placeholder=str(x_start_init) + "-" + str(x_stop_init), ),
                        html.Br(),
                        #html.Br(),
                        html.I("Genes in this region: " + str(local_genes), id='regional_genes'),
                        html.Br(),
                        #html.Br(),
                        html.Label('Pre-selected regions: '),
                        dcc.Dropdown(pre_bed_info[anchor_name], chrs + ":" + str(x_start_init) + "-" + str(x_stop_init) , style=dict(width='100%', height='110%', verticalAlign="middle", ), id="pre_regions_dropdown"),
                        #html.Div(id="Anchor_tab_dropdown"),
                        #html.Br(),
                        ],style={"display": "inline-block"}),
                        
                        html.Div(children=[
                        
                        html.I("K-mer length: " + str(kmer_len),style={"display": "inline-block",}),
                        html.Br(),
                        html.I("Number of bins: " + str(bins), style={'display': 'inline-block'}),
                        html.Br(),
                        html.I("Step size: " + str(n_skips_start), style={"display": "inline-block",}, id='step_size_out'),
                        html.Br(),
                        #html.Br(),
                        #html.Br(),
                        #html.Br(),
                        #html.Br(),
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
                            dcc.Graph(id="primary",figure = main_fig, config=config, 
                                style={"height": 1000,  "font-size": 20})
                        ]), #style={'padding-top' : '1%', 'padding-left' : '1%', 'padding-bottom' : '1%', 'padding-right' : '1%', 
                        html.Div(className="w3-quarter", children=[
                            dcc.Graph(id="Secondary", 
                                #Now we have the phylogenetic tree
                                figure = get_local_info(names_simp[x_start_init:x_stop_init], #exon_comp[chrs], 
                                    #gene_comp[anchor_name][chrs], 
                                    bar_sum_regional, bar_sum_global[anchor_name][chrs], anchor_name, chrs), 
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
                                
                                figure = create_tree( x_start_init,x_stop_init, all_chrs[chrs][int(x_start_init/n_skips_start):int(x_stop_init/n_skips_start)], n_skips_start),
                                config=config,
                                style={"font-size": 20}),
                        ]),
                        
                        
                        ])
                    ])
            ]),  
        #]),
        ]),

        html.Div(chrs, style={"display" : "none"}, id="selected-chrom-state"),
        html.Div(anchor_name, style={"display" : "none"}, id="selected-anchor-state"),
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
        Output('selected-anchor-state', 'children'),
        Output('all_chromosomes', 'figure'),
        Output('anchor_labels', 'children'),
        Output('Anchor_tab_dropdown', 'options'),
        Output('pangenome_num_seqs', 'figure'),
        Output('average_kmer_content', 'figure'),
        Output('all_genome_sizes', 'figure'),
        Output('gene_content', 'figure'),
        Output('genes_per_chr', 'figure'),
        Output('avg_kmer_chr', 'figure'),
        #Output('step_size_out', 'children'),

        Input('all_chromosomes','relayoutData'),
        Input('chromosome','figure'),
        Input('primary','figure'),
        Input('Secondary', 'figure'),
        Input('Third', 'figure'),
        Input('Genes', 'figure'),
        Input('num-start','children'),
        Input('num-stop','children'),
        Input('primary', 'clickData'),
        Input('primary', 'relayoutData'),
        Input('chromosome', 'selectedData'),
        Input('Genes', 'clickData'),
        Input('Chrs_Info', 'value'),
        Input('Anchor_tab_dropdown', 'value'),
        Input('pangenome_tab_dropdown', 'value'),
        Input('pre_regions_dropdown', 'value'),
        #Input('Anchor_name', 'children'),
        #Input('Anchor_name', 'children'),

        State('selected-chrom-state', 'children'),
        State('selected-anchor-state', 'children'),
        #State('', 'children'),
        #Sent 11
        #State('prev_selected-state', 'children'),
        )
    def update_figure(wgr, chr_fig, fig1, fig2, fig3, fig4, #clicker, 
        x_start, x_stop, #gene_jump, 
        clickData, relayoutData, chr_relayoutData, gene_jump_bottom, user_chr_coords, 
        anchor_tab_dropdown, select_anchor_dropdown, pre_selected_region, chrs, anchor_name):#clicker_update, gene_jump, bins):
        #Given 12 
        triggered_id = ctx.triggered_id
        print(triggered_id)
        n_skips = 100
        click_me_genes = True
        click_me_rep = True
        if triggered_id == 'Anchor_tab_dropdown':
            #print(anchor_tab_dropdown)
            chrs = anchor_tab_dropdown
            chr_num = chrs_list[anchor_name].get_loc(chrs)
            #chr_num = chrs_list[anchor_name].index(chrs) #int(chrs.split('r')[1])
            return update_all_figs(chr_num, chr_relayoutData, click_me_rep, click_me_genes, chrs, anchor_name, 0, x_start_init, x_stop_init)
        elif triggered_id == "pre_regions_dropdown":
            #print(pre_selected_region)
            if pre_selected_region != None:
                chrs = pre_selected_region.split(":")[0]

                x_start = int(pre_selected_region.split(":")[1].split("-")[0])
                x_stop = int(pre_selected_region.split(":")[1].split("-")[1])
                chr_num = chrs_list[anchor_name].get_loc(chrs)+1
                return update_all_figs(chr_num, chr_relayoutData, click_me_rep, click_me_genes, chrs, anchor_name, 0, x_start, x_stop)
        elif triggered_id == 'pangenome_tab_dropdown':
            anchor_name = select_anchor_dropdown
            #print(chrs_list)
            chrs = chrs_list[anchor_name][0]

            chr_num = 1 #chrs_list[anchor_name].index(chrs) #int(chrs.split('r')[1])
            
            return update_all_figs(chr_num, chr_relayoutData, click_me_rep, click_me_genes, chrs, anchor_name, 1, x_start_init, x_stop_init)
        elif triggered_id == 'Chrs_Info':
            #print(user_chr_coords)
            if user_chr_coords.count(' ') > 0:
                new_x_start = int(user_chr_coords.strip().split('-')[0]) 
                new_x_stop  = int(user_chr_coords.strip().split('-')[1])
                return user_chr_coords_triggered(new_x_start, new_x_stop, click_me_genes, click_me_rep, chr_relayoutData, 
                    chrs, fig4, n_skips_start, anchor_name), no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update 
            else:
                return no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update 
        elif triggered_id == 'all_chromosomes':
            #Get the chromosome number 
            if wgr != None and len(wgr)>1:
                chr_num_tmp = list(wgr.keys())[0].split('.')[0].split('axis')[1] #if wgr != None and len(wgr)>1:
                if len(chr_num_tmp) == 0:
                    chr_num = 1
                else:
                    chr_num = int(int(chr_num_tmp)/3)+1
                return update_all_figs(chr_num, chr_relayoutData, click_me_rep, click_me_genes, chrs, anchor_name, 0, x_start_init, x_stop_init) #chr_fig, fig1, fig2, fig3, fig4, chr_relayoutData
        elif triggered_id == 'Genes':
            #print("ANCHOR")
            #print(anchor_name)
            return sorted_gene_fig(chr_fig, fig1, fig2, fig3, fig4, gene_jump_bottom, #SG_window, SG_polynomial_order,
                #SG_check, 
                n_skips, click_me_genes, click_me_rep, chr_relayoutData, chrs, anchor_name)

        elif triggered_id == 'chromosome':
            return chromosome_gene_triggered(chr_fig, fig1, fig2, fig3, fig4, gene_jump_bottom, n_skips, click_me_genes, click_me_rep, chr_relayoutData, chrs, anchor_name)
        elif triggered_id == 'primary':
            #print("Is the issue here?")
            return primary_fig_triggered(chr_fig, fig1, fig2, fig3, fig4, gene_jump_bottom, n_skips, click_me_genes, click_me_rep, chr_relayoutData, relayoutData, clickData, x_start, x_stop, chrs, anchor_name)
        local_gene_list = index.query_genes(anchor_name, chrs, int(x_start), int(x_stop))
        local_gene_list["name"] = local_gene_list["attr"].str.extract("ID=([^;]+)")
        #local_gene_list = []
        #cntr = 0
        #while cntr < len(gene_names_tmp):
        #    local_gene_list.append(gene_names_tmp[cntr])
        #    cntr += 3
        return chr_fig, fig1, fig2, fig3, fig4, chr_relayoutData, chrs, update_output_div(chrs,x_start,x_stop,anchor_name), update_out_chr(chrs, anchor_name) , update_gene_locals(local_gene_list, chrs, x_start, x_stop, anchor_name), anchor_name, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update   #, clickData

    def user_chr_coords_triggered(new_x_start, new_x_stop, click_me_genes, click_me_rep, chr_relayoutData, chrs, 
        fig4, n_skips_start, anchor_name):

        if get_buffer(new_x_start, new_x_stop, n_skips_start) == 1:
            if (new_x_stop - new_x_start) <= bins:
                new_x_start = new_x_start - 175
                new_x_stop = new_x_stop + 175
            exact_mini_counts = index.query_bitmap(anchor_name, chrs, new_x_start, new_x_stop, 1)
            simple_cnts_for_plots = exact_mini_counts.sum(axis=1)
            n_skips = 1
        else:
            exact_mini_counts = index.query_bitmap(anchor_name, chrs, new_x_start, new_x_stop, n_skips)
            simple_cnts_for_plots = exact_mini_counts.sum(axis=1)
        chr_num = chrs_list[anchor_name].get_loc(chrs)
        #chr_num = chrs_list[anchor_name].index(chrs)#int(chrs.split('r')[1])
        genes = index.query_genes(l, chrom, 0, index.chrs.loc[l, chrom]["size"])
        #[g.split(';')[1].split('=')[1] for g in genes['attr']]
        bounds = genes.loc[:, ["start", "end"]]
        bounds["break"] = None
        #gene_locals[l][chrom] = bounds.to_numpy().flatten()
        gene_names = [g.split(';')[0].split("=")[1] for g in genes['attr']]
        fig1, bar_sum_names, bar_sum_regional, colors, gene_names_tmp = plot_interactive( n_skips, 
            layout, 
            bins, simple_cnts_for_plots, 
            chrs, zs_tmp, click_me_rep, click_me_genes,  
            int(new_x_start), int(new_x_stop), 
            bounds.to_numpy().flatten(),#gene_locals[anchor_name][chrs], 
            gene_names, #[g.split(';')[0].split('=')[1] for g in genes['attr']],#gene_names[anchor_name][chrs], 
            anchor_name, chrs
            #exon_locals[anchor_name][chrs], exon_names[anchor_name][chrs],
            )
        fig2 = get_local_info(
            names_simp[int(new_x_start/n_skips):int(new_x_stop/n_skips)], 
            #gene_comp[anchor_name][chrs], 
            bar_sum_regional, bar_sum_global[anchor_name][chrs], anchor_name, chrs)
        fig3 = create_tree(new_x_start, new_x_stop, 
            exact_mini_counts, n_skips)
        chr_fig = plot_chr_whole(
            new_x_start, new_x_stop, anchor_name, chrs)
        #local_gene_list = []
        #cntr = 0
        #while cntr < len(gene_names_tmp):
        #    local_gene_list.append(gene_names_tmp[cntr])
        #    cntr += 3
        local_gene_list = index.query_genes(anchor_name, chrs, int(x_start), int(x_stop))
        local_gene_list["name"] = local_gene_list["attr"].str.extract("ID=([^;]+)")
        #genes = index.query_genes(anchor_name, chrs, 0, index.chrs.loc[anchor_name, chrs]["size"])
        
        universals = genes["universal"]
        uniques = genes["unique"]
        #print(genes["end"]-genes["start"])
        sizes = genes["end"]-genes["start"]
        fig4 = plot_gene_content(gene_names,universals,uniques,sizes, sort_by, colors, uniq_avg, 
            univ_avg, int(new_x_start), int(new_x_stop), anchor_name, chrs)
        return chr_fig, fig1, fig2, fig3, fig4, chr_relayoutData, chrs, update_output_div(chrs,new_x_start,new_x_stop,anchor_name), update_out_chr(chrs, anchor_name), update_gene_locals(local_gene_list, chrs, x_start, x_stop, anchor_name), anchor_name, update_anchor_name(anchor_name), no_update, no_update, no_update, no_update, no_update, no_update, no_update 

    def sorted_gene_fig(chr_fig, fig1, fig2, fig3, fig4, gene_jump_bottom, #SG_window, SG_polynomial_order,
                n_skips, click_me_genes, click_me_rep, chr_relayoutData, chrs, anchor_name ):
        chr_num = chrs_list[anchor_name].get_loc(chrs)
        #chr_num = chrs_list[anchor_name].index(chrs)
        names_simp = index.query_bitmap(anchor_name, chrs, 0, index.chrs.loc[anchor_name, chrs]["size"], n_skips_start).sum(axis=1) #all_chrs[chrs].sum(axis=1)
        #int(chrs.split('r')[1])
        if gene_jump_bottom != None: #and gene_jump_bottom['points'][0]['text'] in gene_names[anchor_name][chrs]:
            #print("first elif")
            
            genes = index.query_genes(anchor_name, chrs, 0, index.chrs.loc[anchor_name, chrs]["size"])

            #[g.split(';')[1].split('=')[1] for g in genes['attr']]
            bounds = genes.loc[:, ["start", "end"]]
            bounds["break"] = None
            gene_locals = bounds.to_numpy().flatten()
            #print(gene_jump_bottom)
            this_gene_name = gene_jump_bottom['points'][0]['text']
            #print(this_gene_name)
            genes["name"] = genes["attr"].str.extract("ID=([^;]+)")
            this_gene = genes[genes["name"] == this_gene_name]
            this_start = int(this_gene['start'])-1
            this_stop = int(this_gene['end'])+1
            if get_buffer(this_start,this_stop, n_skips_start) == 1:
                #we need to move to single nucleotide resoltion               
                if (this_stop - this_start) <= bins:
                    x_start = this_start - 175
                    x_stop = this_stop + 175
                else:
                    x_start = this_start-1
                    x_stop = this_stop+1
                exact_mini_counts = index.query_bitmap(anchor_name, chrs, x_start, x_stop, 1)
                simple_cnts_for_plots = exact_mini_counts.sum(axis=1)
                #simple_cnts_for_plots, exact_mini_counts = read_mini_count_files(x_start, x_stop)
                n_skips = 1
            else: 
                x_start = this_start - buff #int(1000*n_skips/2)
                x_stop = this_stop + buff #int(1000*n_skips/2)
                n_skips = n_skips_start
                simple_cnts_for_plots = names_simp[int(x_start/n_skips):int(x_stop/n_skips)]
                exact_mini_counts = all_chrs[chrs][int(x_start/n_skips):int(x_stop/n_skips)]
            
            fig1, bar_sum_names, bar_sum_regional, colors, gene_names_tmp = plot_interactive( n_skips, #int(SG_window), 
                #int(SG_polynomial_order), SG_check,
                layout, #exon_comp[chrs], 
                bins, simple_cnts_for_plots, chrs, zs_tmp, 
                click_me_rep, click_me_genes,  
                int(x_start), int(x_stop), gene_locals, 
                [g.split(';')[0].split('=')[1] for g in genes['attr']], #gene_names[anchor_name][chrs], 
                anchor_name, chrs
                #exon_locals[anchor_name][chrs], exon_names[anchor_name][chrs],
                )
            fig2 = get_local_info(
                simple_cnts_for_plots, #exon_comp[chrs], 
                #gene_comp[anchor_name][chrs], 
                bar_sum_regional, bar_sum_global[anchor_name][chrs], anchor_name, chrs)
            fig3 = create_tree(x_start, x_stop, exact_mini_counts, n_skips)
            gene_jump = ""
            if x_start_init != x_start:
                chr_fig = plot_chr_whole(
                    x_start, x_stop, anchor_name, chrs)
        #local_gene_list = []
        #cntr = 0
        #while cntr < len(gene_names_tmp):
        #    local_gene_list.append(gene_names_tmp[cntr])
        #    cntr += 3
        #print("Sorted gene list. Almost finished ")
        local_gene_list = index.query_genes(anchor_name, chrs, int(x_start), int(x_stop))
        local_gene_list["name"] = local_gene_list["attr"].str.extract("ID=([^;]+)")

        genes = index.query_genes(anchor_name, chrs, 0, index.chrs.loc[anchor_name, chrs]["size"])
        gene_names = [g.split(';')[0].split("=")[1] for g in genes['attr']]
        universals = genes["universal"]
        uniques = genes["unique"]
        #print(genes["end"]-genes["start"])
        sizes = genes["end"]-genes["start"]
        fig4 = plot_gene_content(gene_names,universals,uniques,sizes,
            sort_by, colors, uniq_avg, univ_avg, int(x_start), int(x_stop),anchor_name, chrs)
        return chr_fig, fig1, fig2, fig3, fig4, chr_relayoutData, chrs, update_output_div(chrs,x_start,x_stop,anchor_name), update_out_chr(chrs, anchor_name), update_gene_locals(local_gene_list, chrs, x_start, x_stop, anchor_name), anchor_name, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update 
    
    def chromosome_gene_triggered(chr_fig, fig1, fig2, fig3, fig4, gene_jump_bottom, #SG_window, SG_polynomial_order,
        n_skips, click_me_genes, click_me_rep, chr_relayoutData, chrs, anchor_name):
        #print("HERE!")
        chr_num = chrs_list[anchor_name].get_loc(chrs)
        names_simp = index.query_bitmap(anchor_name, chrs, 0, index.chrs.loc[anchor_name, chrs]["size"], n_skips_start).sum(axis=1) #all_chrs[chrs].sum(axis=1)
        
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
            if (x_stop - x_start) <= bins:
                x_start = x_start - 175
                x_stop = x_stop + 175
            exact_mini_counts = index.query_bitmap(anchor_name, chrs, x_start, x_stop, 1)
            simple_cnts_for_plots = exact_mini_counts.sum(axis=1)
            n_skips = 1
        else:
            exact_mini_counts = index.query_bitmap(anchor_name, chrs, int(x_start/n_skips_start), int(x_stop/n_skips_start), n_skips_start)#all_chrs[chrs][int(x_start/n_skips_start):int(x_stop/n_skips_start)]
            simple_cnts_for_plots = names_simp[int(x_start/n_skips_start):int(x_stop/n_skips_start)]
            n_skips = n_skips_start
        
        genes = index.query_genes(anchor_name, chrs, x_start, x_stop)
        #[g.split(';')[1].split('=')[1] for g in genes['attr']]
        bounds = genes.loc[:, ["start", "end"]]
        bounds["break"] = None
        #print([g.split(';')[1].split('=')[1] for g in genes['attr']])
        fig1, bar_sum_names, bar_sum_regional, colors, gene_names_tmp = plot_interactive( n_skips, #int(SG_window), 
                layout, #exon_comp[chrs], 
                bins, simple_cnts_for_plots, chrs, zs_tmp,  
                click_me_rep, click_me_genes,  
                int(x_start), int(x_stop), 
                bounds.to_numpy().flatten(),#gene_locals[anchor_name][chrs], 
                [g.split(';')[0].split('=')[1] for g in genes['attr']],#gene_names[anchor_name][chrs], 
                anchor_name, chrs
                #exon_locals[anchor_name][chrs], exon_names[anchor_name][chrs],
                )
        fig2 = get_local_info(
                simple_cnts_for_plots, #exon_comp[chrs], 
                #gene_comp[anchor_name][chrs], 
                bar_sum_regional, bar_sum_global[anchor_name][chrs], anchor_name, chrs)
        #And now we update the histograms 
        fig3 = create_tree( x_start, x_stop, exact_mini_counts, n_skips)
        chr_relayoutData = None#{'autosize': True}
        if x_start_init != x_start:
            chr_fig = plot_chr_whole(
                x_start, x_stop, anchor_name, chrs)
        #local_gene_list = []
        #cntr = 0
        #while cntr < len(gene_names_tmp):
        #    local_gene_list.append(gene_names_tmp[cntr])
        #    cntr += 3
        local_gene_list = genes #index.query_genes(anchor_name, chrs, int(x_start), int(x_stop))
        local_gene_list["name"] = local_gene_list["attr"].str.extract("ID=([^;]+)")
        genes = index.query_genes(anchor_name, chrs, 0, index.chrs.loc[anchor_name, chrs]["size"])
        gene_names = [g.split(';')[0].split("=")[1] for g in genes['attr']]
        universals = genes["universal"]
        uniques = genes["unique"]
        #print(genes["end"]-genes["start"])
        sizes = genes["end"]-genes["start"]
        fig4 = plot_gene_content(gene_names,universals,uniques,sizes, 
            sort_by, colors, uniq_avg, univ_avg, x_start, x_stop, anchor_name, chrs)
        return chr_fig, fig1, fig2, fig3, fig4, chr_relayoutData, chrs, update_output_div(chrs,x_start,x_stop,anchor_name), update_out_chr(chrs, anchor_name), update_gene_locals(local_gene_list, chrs, x_start, x_stop, anchor_name), anchor_name, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update 

    def primary_fig_triggered(chr_fig, fig1, fig2, fig3, fig4, gene_jump_bottom, #SG_window, SG_polynomial_order,
        n_skips, click_me_genes, click_me_rep, chr_relayoutData, relayoutData, clickData, 
        x_start, x_stop, chrs, anchor_name ):
        chr_num = chrs_list[anchor_name].get_loc(chrs)
        chr_size = index.chrs.loc[(anchor_name, chrs), "size"]
        names_simp = index.query_bitmap(anchor_name, chrs, 0, chr_size, n_skips_start).sum(axis=1)#all_chrs[chrs].sum(axis=1)
        genes = index.query_genes(anchor_name, chrs, 0, index.chrs.loc[anchor_name, chrs]["size"])
        #[g.split(';')[1].split('=')[1] for g in genes['attr']]
        bounds = genes.loc[:, ["start", "end"]]
        bounds["break"] = None
        gene_locals = bounds.to_numpy().flatten()
        #int(chrs.split('r')[1])
        if relayoutData != None and 'xaxis4.range[0]' in relayoutData.keys():
            #print("fourth elif")
            x_start = int(relayoutData['xaxis4.range[0]'])
            x_stop = int(relayoutData['xaxis4.range[1]'])
            if get_buffer(x_start, x_stop, n_skips_start) == 1:
                if (x_stop - x_start) <= bins:
                    x_start = x_start - 175
                    x_stop = x_stop + 175
                exact_mini_counts = index.query_bitmap(anchor_name, chrs, x_start, x_stop, 1)
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
                bins, 
                simple_cnts_for_plots,
                #names_simp[int(x_start/n_skips_start):int(x_stop/n_skips_start)], 
                chrs, zs_tmp, 
                click_me_rep, click_me_genes,  
                int(x_start), int(x_stop), 
                gene_locals,#gene_locals[anchor_name][chrs], 
                [g.split(';')[0].split('=')[1] for g in genes['attr']],#gene_names[anchor_name][chrs], 
                anchor_name, chrs
                #exon_locals[anchor_name][chrs], exon_names[anchor_name][chrs],
                )

            fig2 = get_local_info(
                #full_cnts[x_start:x_stop],
                simple_cnts_for_plots, 
                bar_sum_regional, bar_sum_global[anchor_name][chrs], anchor_name, chrs)
            #And now we update the histograms 
            fig3 = create_tree(x_start, x_stop, 
                exact_mini_counts,
                #all_chrs[chrs][int(x_start/n_skips_start):int(x_stop/n_skips_start)], 
                n_skips)
            chr_fig = plot_chr_whole(
                x_start, x_stop, anchor_name, chrs)
            #local_gene_list = []
            #cntr = 0
            #while cntr < len(gene_names_tmp):
            #    local_gene_list.append(gene_names_tmp[cntr])
            #    cntr += 3
            local_gene_list = index.query_genes(anchor_name, chrs, int(x_start), int(x_stop))
            local_gene_list["name"] = local_gene_list["attr"].str.extract("ID=([^;]+)")
            genes = index.query_genes(anchor_name, chrs, 0, index.chrs.loc[anchor_name, chrs]["size"])
            gene_names = [g.split(';')[0].split("=")[1] for g in genes['attr']]
            universals = genes["universal"]
            uniques = genes["unique"]
            #print(genes["end"]-genes["start"])
            sizes = genes["end"]-genes["start"]
            fig4 = plot_gene_content(gene_names,universals,uniques,sizes, 
                sort_by, colors, uniq_avg, univ_avg, int(x_start), int(x_stop), anchor_name, chrs)
        elif clickData != None:#len(print(clickData['points'])) > 0:
            #print("third elif")
            #print(clickData)
            #tmp_idx = gene_names[anchor_name][chrs].index(clickData['points'][0]['text'])
            #gene_buffer = get_buffer(int(gene_locals[tmp_idx]), int(gene_locals[tmp_idx+1]), n_skips_start)
            this_gene_name = clickData['points'][0]['text']
            
            genes["name"] = genes["attr"].str.extract("ID=([^;]+)")
            this_gene = genes[genes["name"] == this_gene_name]

            x_start = int(this_gene['start'])
            x_stop = int(this_gene['end'])
            if get_buffer(x_start, x_stop, n_skips_start) == 1:
                #we need to move to single nucleotide resoltion
                if (x_stop - x_start) <= bins:
                    x_start = x_start - 175
                    x_stop = x_stop + 175
                exact_mini_counts = index.query_bitmap(anchor_name, chrs, x_start, x_stop, 1)
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
                bins, simple_cnts_for_plots, chrs, zs_tmp, 
                click_me_rep, click_me_genes,  
                int(x_start), int(x_stop), gene_locals, 
                [g.split(';')[0].split('=')[1] for g in genes['attr']], #gene_names[anchor_name][chrs], 
                anchor_name, chrs
                )
            fig2 = get_local_info(
                #full_cnts[x_start:x_stop], 
                simple_cnts_for_plots, 
                #exon_comp[chrs], 
                #gene_comp[anchor_name][chrs], 
                bar_sum_regional, bar_sum_global[anchor_name][chrs], anchor_name, chrs)
            #And now we update the histograms 
            fig3 = create_tree(x_start, x_stop, 
                exact_mini_counts, n_skips)
            chr_fig = plot_chr_whole(x_start, x_stop, anchor_name, chrs)
            local_gene_list = index.query_genes(anchor_name, chrs, int(x_start), int(x_stop))
            local_gene_list["name"] = local_gene_list["attr"].str.extract("ID=([^;]+)")
            #local_gene_list = []
            #cntr = 0
            #while cntr < len(gene_names_tmp):
            #    local_gene_list.append(gene_names_tmp[cntr])
            #    cntr += 3 
            genes = index.query_genes(anchor_name, chrs, 0, index.chrs.loc[anchor_name, chrs]["size"])
            gene_names = [g.split(';')[0].split("=")[1] for g in genes['attr']]
            universals = genes["universal"]
            uniques = genes["unique"]
            sizes = genes["end"]-genes["start"]
            #local_gene_list = index.query_genes(anchor_name, chrs, int(x_start), int(x_stop))
            #local_gene_list = [g.split(';')[1].split("=")[1] for g in local_gene_list['attr']]
            fig4 = plot_gene_content(gene_names,universals,uniques,sizes, 
                sort_by, colors, uniq_avg, univ_avg, int(x_start), int(x_stop), anchor_name, chrs)   
        else:
            local_gene_list = index.query_genes(anchor_name, chrs, x_start_init,x_stop_init)
            local_gene_list["name"] = local_gene_list["attr"].str.extract("ID=([^;]+)")
            #local_gene_list = [g.split(';')[0].split("=")[1] for g in local_gene_list['attr']]
        return chr_fig, fig1, fig2, fig3, fig4, chr_relayoutData, chrs, update_output_div(chrs,x_start,x_stop,anchor_name), update_out_chr(chrs, anchor_name), update_gene_locals(local_gene_list, chrs, x_start, x_stop, anchor_name), anchor_name, no_update, no_update, no_update,no_update, no_update, no_update, no_update, no_update, no_update 
    def render_tab_content(active_tab, chrs):
        return chr_fig, fig1, fig2, fig3, fig4, chr_relayoutData, chrs, update_output_div(chrs,x_start_init,x_stop_init,anchor_name), update_out_chr(chrs, anchor_name), update_gene_locals(local_gene_list, chrs, x_start, x_stop, anchor_name), anchor_name, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update 
    def update_all_figs(chr_num, chr_relayoutData, click_me_rep, click_me_genes, #SG_window, SG_polynomial_order,
        chrs, anchor_name, redo_wg, x_start, x_stop):
        chrs = index.chrs.loc[anchor_name].index[chr_num-1]#"chr" + str(chr_num)
        names_simp = index.query_bitmap(anchor_name, chrs, x_start, x_stop, n_skips_start).sum(axis=1)
        genes = index.query_genes(anchor_name, chrs, 0, index.chrs.loc[anchor_name, chrs]["size"])
        #[g.split(';')[1].split('=')[1] for g in genes['attr']]
        bounds = genes.loc[:, ["start", "end"]]
        bounds["break"] = None
        #gene_locals = bounds.to_numpy().flatten()
        gene_names = [g.split(';')[0].split("=")[1] for g in genes['attr']]
        fig1, bar_sum_names, bar_sum_regional, colors, gene_names_tmp = plot_interactive( n_skips_start, #int(SG_window), 
                #int(SG_polynomial_order), SG_check,
                layout, #exon_comp[chrs], 
                bins, names_simp, #names_simp[x_start_init_adjust:x_stop_init_adjust], 
                chrs, zs_tmp, click_me_rep, click_me_genes,  
                int(x_start), int(x_stop), 
                bounds.to_numpy().flatten(),#gene_locals[anchor_name][chrs], 
                gene_names,
                #[g.split(';')[0].split('=')[1] for g in genes['attr']], #gene_names[anchor_name][chrs], 
                anchor_name, chrs
                #exon_locals[anchor_name][chrs], 
                #exon_names[anchor_name][chrs]
                )
        #local_gene_list = []
        local_gene_list = index.query_genes(anchor_name, chrs, int(x_start), int(x_stop))
        local_gene_list["name"] = local_gene_list["attr"].str.extract("ID=([^;]+)")
        #cntr = 0
        #while cntr < len(gene_names_tmp):
        #    local_gene_list.append(gene_names_tmp[cntr])
        #    cntr += 3
        sort_by = ["Unique","Universal"]
        tmp = [(i/sum(bar_sum_global[anchor_name][chrs])*100) for i in bar_sum_global[anchor_name][chrs]]
        uniq_avg = tmp[1]
        univ_avg = tmp[-1]
        #genes = index.query_genes(anchor_name, chrs, 0, index.chrs.loc[anchor_name, chrs]["size"])
        #gene_names = [g.split(';')[0].split("=")[1] for g in genes['attr']]
        universals = genes["universal"]
        uniques = genes["unique"]
        #print(genes["end"]-genes["start"])
        sizes = genes["end"]-genes["start"]
        fig4 = plot_gene_content(gene_names,universals,uniques,sizes, 
            sort_by, colors, uniq_avg, univ_avg, int(x_start), int(x_stop), anchor_name, chrs)
        
        #x, z_1, z_9, y  = read_chr_whole()
        #z_genes = make_gene_whole_chr(x[chr_num-1], gene_locals)
        #print(z_genes[anchor_name][chrs])
        chr_fig = plot_chr_whole(x_start, 
            x_stop, anchor_name, chrs)

        fig2 = get_local_info(names_simp,
                #names_simp[x_start_init_adjust:x_stop_init_adjust], #exon_comp[chrs], 
                #gene_comp[anchor_name][chrs], 
                bar_sum_regional, bar_sum_global[anchor_name][chrs], anchor_name, chrs)
        
        fig3 = create_tree(x_start, x_stop, 
            index.query_bitmap(anchor_name, chrs, x_start, x_stop, n_skips_start),
            #all_chrs[chrs][x_start_init_adjust:x_stop_init_adjust], 
            n_skips_start)
        #chr_fig = plot_chr_whole(x[chr_num-1], z_1[chr_num-1], z_9[chr_num-1], z_genes[chr_num-1], x_start_init, x_stop_init, y[chr_num-1])
        if redo_wg == 1:
            big_plot = plot_whole_genome(anchor_name)
            pg_sizes_fig, pg_scaffolds_fig = make_genome_size_plots(anchor_name, pg_sizes, pg_num_seqs)
            pg_avgs = plot_avgs(anchor_name, tmp_df)
            tab2_gene_density_fig = make_genes_per_chr_fig(anchor_name)
            #print("Finished tab 2 gene density fig")
            chrs_comp = []
            for c in index.chrs.loc[anchor_name].index:
                tmp = []
                for n in range(1, num_samples+1):
                    tmp.append(index.chrs.loc[anchor_name, c]["total_occ_"+str(n)])
                chrs_comp.append(tmp)
            #print("Made it to redowg")
            tab2_avg_fig = make_avg_kmer_fig(chrs_comp, anchor_name)
            #print("Finished tab 2 avg fig")
            tab2_sorted_genes = make_gene_per_genome_fig(anchor_name)
            #print("Finished tab 2 sorted genes fig")
        else:
            big_plot = no_update
            pg_sizes_fig = no_update
            pg_avgs = no_update
            pg_scaffolds_fig = no_update

            tab2_gene_density_fig = no_update
            tab2_avg_fig = no_update
            tab2_sorted_genes = no_update
        #print("Made it here!")
        return chr_fig, fig1, fig2, fig3, fig4, chr_relayoutData, chrs, update_output_div(chrs,x_start_init,x_stop_init, anchor_name), update_out_chr(chrs,anchor_name), update_gene_locals(local_gene_list, chrs, x_start, x_stop, anchor_name), anchor_name, big_plot, update_anchor_name(anchor_name), update_chromosome_list(anchor_name), pg_scaffolds_fig, pg_avgs, pg_sizes_fig, tab2_sorted_genes, tab2_gene_density_fig, tab2_avg_fig
    def update_output_div(chrs, start, stop, anchor_name):
        return f'{start}-{stop}'
    def update_out_chr(chrs, anchor_name):
        return f'{anchor_name}.{chrs}:'
    def update_gene_locals(local_gene_list, chrs, x_start, x_stop, anchor_name):
        printme = ""
        if len(local_gene_list)==1:
            #print("In updating gene locals section")
            #print(index.query_anno(anchor_name, chrs, x_start, x_stop))
            #print(local_gene_list['name'][0])
            #print(local_gene_list['attr'][0])
            #print("In updating gene locals section")
            printme += "Genes: "
            printme += local_gene_list['name'][0] + ": "
            printme += local_gene_list['attr'][0] #index.query_anno(anchor_name, chrs, x_start, x_stop)['attr']
            #gene_anns[anchor_name][local_gene_list[0]]
        elif len(local_gene_list)<=10:
            printme += "Genes: "
            for i in local_gene_list['name']:
                printme += i + ", "
        else:
            printme = "Genes in this region: " + str(len(local_gene_list))
        return f'{printme}'
    def update_anchor_name(anchor_name):
        return f'{anchor_name}'
    def update_chromosome_list(anchor_name):
        return_me = [{'label': i, 'value': i} for i in index.chrs.loc[anchor_name].index]
        return return_me 

    app.run_server(host=params.host, port=params.port, debug=not params.ndebug)





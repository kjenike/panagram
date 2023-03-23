import sys
import csv
import pysam
import os.path
from os import path
from mycolorpy import colorlist as mcp
from plotly.subplots import make_subplots
import plotly.figure_factory as ff
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
from scipy.cluster import hierarchy
from io import StringIO
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist, squareform
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
    #anchor_name = "kakapo"
    #chrs = "S1"
    if params.genome is not None:
        anchor_name = params.genome
    if params.chrom is not None:
        chrs = params.chrom

    x_start_init = 0 if params.start is None else params.start
    x_stop_init  =  index.chrs.loc[(anchor_name,chrs),"size"] if params.end is None else params.end
    bins = 200#params.max_chr_bins
    opt_bed_file = params.bookmarks
    window_size = index.chr_bin_kbp*1000
    n_skips_start = index.lowres_step

    num_samples = len(index.genomes)
    kmer_len = index.k
    rep_list = index.gff_anno_types
    labels = list(index.genomes)

    sns.set_palette('viridis', num_samples)
    colors = mcp.gen_color(cmap="viridis_r",n=num_samples)


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

    def get_x_coordinates(tree):
        #Adapted from:
        #https://github.com/plotly/dash-phylogeny
        """Associates to  each clade an x-coord.
           returns dict {clade: x-coord}
        """
        xcoords = tree.depths()
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

    def get_clade_lines(orientation='horizontal', y_curr=0, start_coord=0, x_curr=0, y_bot=0, y_top=0,
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
            branch_line.update(x0=start_coord,
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
            clade_num = kmer_num[tmp_name]
            return clade_num

    def draw_clade(color_code, total_kmers, kmer_num, palette, clade, start_coord, line_shapes, line_color="grey", line_width=5, x_coords=0, y_coords=0):
        #Adapted from:
        #https://github.com/plotly/dash-phylogeny
        """Recursively draw the tree branches, down from the given clade"""
        x_curr = x_coords[clade]
        y_curr = y_coords[clade]
        if str(clade) != "Clade" and str(clade) != "AT":
            tmp_name = str(clade)
            line_color = color_code[tmp_name]#palette[int(((kmer_num[tmp_name])/total_kmers)*100)-1]
        elif clade.clades:
            
            line_color = palette[int(biggest_num_in_clade(clade, kmer_num))+10]
            #line_color = "grey"
            #Now we have to find the all of the children, and which one has the highest value   
        # Draw a horizontal line from start to here
        branch_line = get_clade_lines(orientation='horizontal', y_curr=y_curr, start_coord=start_coord, x_curr=x_curr,
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
                      #font=dict(family='Balto', size=22),
                      # width=1000,
                      #height=1500,
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
            vertical_spacing=0.1,
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
        
        #fig.add_trace(go.Bar(x=x, y=[a_i - b_i for a_i, b_i in zip(y, y_whole)], marker_color=colors, showlegend=False), row=2, col=1)
        fig.add_trace(go.Bar(x=x, y=y, marker_color=colors, showlegend=False), row=2, col=1)
        fig.add_trace(go.Bar(x=x, y=y_whole, marker_color=colors, showlegend=False), row=1, col=1)
        #Genes
        #y=[(i/sum(gene_comp[1:])*100) for i in gene_comp[1:]]
        #fig.add_trace(go.Bar(x=x, y=y_whole, marker_color="#7400b8", showlegend=False), row=1, col=4)
        totals = 0
        gene_comp_tmp = []
        gene_comp = []
        for i in range(1, num_samples+1):
            totals += index.chr_occ.loc[(anchor_name, chrs), ("gene",i)]
            gene_comp_tmp.append(index.chr_occ.loc[(anchor_name, chrs), ("gene",i)])
        
        for t in gene_comp_tmp:
            gene_comp.append(t/totals*100)
        
        fig.add_trace(go.Bar(x=x, y=[a_i - b_i for a_i, b_i in zip(gene_comp, y_whole)], marker_color=colors, showlegend=False), row=2, col=2)
        #fig.update_layout(xaxis_title_text="K-mers shared in X samples", yaxis_title_text='Frequency (log)')
        fig.update_xaxes(title_text="# of genomes", row=2, col=1)
        fig.update_yaxes(title_text="Difference from whole chromosome", row=2, col=1)
        fig.update_yaxes(title_text="Percent of k-mers", row=1, col=1)
        fig.update_yaxes(type="log", row=1, col=1)
        fig.update_yaxes(type="log", row=2, col=1)
        fig.update_yaxes(type="log", row=2, col=2)
        fig.update_layout(height=1000)
        
        return fig

    def plot_interactive(n_skips, layout, bins, names_simp, name, zs_tmp, plot_rep, plot_gene, start_coord, end_coord, gene_locals, gene_names, anchor_name, chrs):
        window_filter = SG_window
        poly_order = SG_polynomial_order
        shared_kmers = [1]
        tmp_lst = []
        fig = make_subplots(        
            rows=13, cols=1,
            shared_xaxes=True,
            specs=[
                [{"type": "scatter"}], 
                [{"type": "scatter", 'rowspan':2}], [None ], 
                [{"type": "scatter", 'rowspan':2}], [None ],
                [{"type": "bar", 'rowspan':8}] ] + [[None]]*7,
            subplot_titles=("Ref. Sequence Position","Gene Annotation", "Annotation",  "Conserved K-mers" )
        )

        t = time.time()

        #We are adjusting the start and stop positions to account for the skipping. 
        #The adjusted value should be the index, whereas the start_coord and end_coord are the real coordinates 
        adjusted_x_start = int(start_coord/n_skips)
        adjusted_x_stop = int(end_coord/n_skips)

        #Get the bins
        bin_size = ((end_coord-start_coord)/bins)
        adjusted_bin_size = (bin_size/n_skips)
        cntr = 0
        
        cats_tmp = [([0] * (bins+1)) for _ in range(num_samples+1)]
        cntr = 0
        x_tracker = 0
        #Here we are filling in the bins for the main figure.    
        bin_size_int = int(bin_size) + 1

        adjusted_bin_size_init = int(adjusted_bin_size) + 1

        cntr = start_coord
        x = []
        while cntr < end_coord:
            x.append(cntr)
            cntr += bin_size_int
        cntr = 0
        for i in names_simp: #[adjusted_x_start:adjusted_x_stop] : #range(start_coord, end_coord): #names_simp[start_coord:end_coord]:#range(start_coord, end_coord):#names_simp:
            #i tells us which y axis to use. 
            #cntr tells us what x axis to use 
            if (cntr % adjusted_bin_size_init) == 0 : #Is the x position within the same bin? Or do we need to move to the next bin? 
                x_tracker += 1 #X-tracker keeps track of which bin we are using (on the x-axis)
            
            if (i) < len(cats_tmp) and x_tracker < len(cats_tmp[i]):
                cats_tmp[i][x_tracker] += 1 
            cntr += 1 #n_skips

        plot_rep = True
        #Add a line plot that will cover the different repetative elements. 
        if plot_rep == True:
            rep_colors = mcp.gen_color(cmap="plasma",n=len(rep_list ))
            cntr = 0

            for i in rep_list: #rep_types.keys():
                df = index.query_anno(anchor_name, chrs, start_coord, end_coord)
                if i == "exon":
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
                    #df = index.query_anno(anchor_name, chrs, start_coord, end_coord)
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
        
        #This is the conserved kmer plotting section
        bar_sum_regional = []
        bar_sum_names = []
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

        fig.update_layout(barmode='stack', bargap=0.0)
        fig.update_xaxes(showticklabels=False, row=6, col=1)

        #Now we will add the smoothing line. There are three parameters that we adjust
        for sk in shared_kmers:
            y_tmp = cats_tmp[int(sk)][:-1]
            for i in range(0, int(sk)):
                y_tmp = [a + b for a, b in zip(y_tmp, cats_tmp[int(i)][:-1])] #cats_tmp[int(sk)][:-1]
            fig.add_trace(go.Scatter(x=x, y=signal.savgol_filter(y_tmp,window_filter,poly_order), 
                name="Savitzky-Golay - "+str(sk), marker=dict(color="grey"), mode='lines'), row=6, col=1)

        #Now we add the reference sequence:
        y_ref = [1, 1]
        x_ref = [start_coord, end_coord]
        tickvals = []
        ticktxt = []
        cntr = start_coord
        yvals = []
        interval = int((end_coord-start_coord)/10)
        while cntr <= end_coord:
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
        fig.update_layout(height=1000, xaxis_range=[start_coord,end_coord], font=dict(size=16))
        return fig, bar_sum_names, bar_sum_regional, colors, gene_names_tmp

    def make_gene_whole_chr(x, locs):
        z_genes = [0]*(len(x)+2)
        g = 0
        while g < len(locs):
            x_bin = int(int(locs[g])/window_size)
            if x_bin < len(z_genes):
                z_genes[x_bin] += 1
            g += 1
        return z_genes

    def plot_chr_whole( start_coord, end_coord, anchor_name, this_chr): 
        z_1 = index.chr_bins.loc[anchor_name, ("total",1)][this_chr]
        z_9 = index.chr_bins.loc[anchor_name, ("total", num_samples)][this_chr]
        y, x = [], []
        cntr = 0
        for xtmp in z_1:
            x.append(cntr)
            y.append(1)
            cntr += window_size

        z_genes = [0]*len(x)
        if index.gene_tabix[anchor_name] is not None:
            try:
                genes = index.query_genes(anchor_name, this_chr, 0, index.chrs.loc[anchor_name, this_chr]["size"])
            except ValueError:
                genes = None

            if genes is not None:
                bounds = genes["start"].to_numpy() #genes.loc[:, ["start"]].to_numpy()
                z_genes = make_gene_whole_chr(x, bounds) #[g.split(';')[1].split('=')[1] for g in genes['attr']]) 
                
        chr_fig = make_subplots(rows=3, cols=1, 
            specs=[[{"type": "heatmap",}], [{"type": "heatmap",}], [{"type": "heatmap",}]],
            shared_xaxes=True,
            subplot_titles=("K-mer and gene density accross whole chromosome", "",
                ""),
            vertical_spacing = 0.0
            )

        chr_fig.add_trace(go.Heatmap(x=x, z=z_9, y=y, type = 'heatmap', colorscale='magma_r', showlegend=False, showscale=False), row=1, col=1)
        chr_fig.add_trace(go.Scatter(x=[start_coord, start_coord, None, end_coord, end_coord, None, start_coord, end_coord, ], showlegend=False,
                       y=[0.5, 1.5, None, 0.5, 1.5, None, 1.45, 1.45 ],
                       mode='lines',
                       line_color='#1dd3b0', line_width=8), row=1, col=1)

        chr_fig.add_trace(go.Heatmap(x=x, z=z_1, y=y, type = 'heatmap', colorscale='magma', showscale=False), row=2, col=1)
        chr_fig.add_trace(go.Scatter(x=[start_coord, start_coord, None, end_coord, end_coord], showlegend=False,
                       y=[0.5, 1.5, None, 0.5, 1.5],
                       mode='lines',
                       line_color='#1dd3b0', line_width=8), row=2, col=1)

        chr_fig.add_trace(go.Heatmap(x=x, z=z_genes, y=y, type = 'heatmap', colorscale='magma_r', showscale=False ), row=3, col=1)
        chr_fig.add_trace(go.Scatter(x=[start_coord, start_coord, None, end_coord, end_coord, None, start_coord, end_coord,], showlegend=False, 
                       y=[0.5, 1.5, None, 0.5, 1.5, None, 0.55, 0.55],
                       mode='lines',
                       line_color='#1dd3b0', line_width=8), row=3, col=1)
        
        chr_fig.update_yaxes( range=[0.5,1.5], showticklabels=False, row=1, col=1)
        chr_fig.update_yaxes( range=[0.5,1.5], showticklabels=False, row=2, col=1)
        chr_fig.update_yaxes( range=[0.5,1.5], showticklabels=False, row=3, col=1)

        chr_fig.update_xaxes(fixedrange=True, range=[0,index.chrs.loc[anchor_name, this_chr]["size"]], row=1, col=1)
        chr_fig.update_xaxes(fixedrange=True, range=[0,index.chrs.loc[anchor_name, this_chr]["size"]], row=2, col=1)
        chr_fig.update_xaxes(fixedrange=True, row=3, col=1)
        chr_fig.update_xaxes(title_text="Sequence position", range=[0,index.chrs.loc[anchor_name, this_chr]["size"]], row=3, col=1)
        chr_fig.update_yaxes(title_text="Univ.", row=1, col=1)
        chr_fig.update_yaxes(title_text="Uniq.", row=2, col=1)
        chr_fig.update_yaxes(title_text="Genes", row=3, col=1)
        chr_fig.update_layout(clickmode='event+select', dragmode="select", selectdirection='h')
        chr_fig.update_layout(height=350)
        chr_fig.update_layout(margin=dict(b=10, l=10, r=10), font=dict(size=15))
        
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
            if len(x)!=1:
                wg_fig.add_trace(go.Heatmap(x=x, z=index.chr_bins.loc[anchor_name, ("total",num_samples)][chrom], 
                    y=[1]*(len(x)), type = 'heatmap', colorscale='magma_r', showlegend=False,showscale=False), row=((cntr*3)-2), col=1)
                wg_fig.add_trace(go.Heatmap(x=x, z=index.chr_bins.loc[anchor_name, ("total",1)][chrom], 
                    y=[1]*(len(x)), type = 'heatmap', colorscale='magma', showscale=False), row=((cntr*3)-1), col=1)
            
            if cntr == 1:
                wg_fig.update_layout(xaxis={'side': 'top'}) 
            cntr += 1
        
        wg_fig.update_layout(clickmode='event', plot_bgcolor='rgba(0,0,0,0)')
        wg_fig.update_layout(height=h)
        return wg_fig

    def plot_gene_content(gene_names,universals,uniques,sizes, sort_by, colors, uniq_avg, univ_avg, start_coord, end_coord, anchor_name, chrs):
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
        fig = go.Figure(data=[go.Scattergl(x=x, y=df_sorted[sort_by[cntr]], text=df_sorted['Names'], marker=dict(color=colors[cntr]),
                name="% "+sort_by[cntr], mode="markers")])
        cntr += 1
        while cntr < len(sort_by): #s in sort_by:  
            fig.add_trace(go.Scattergl(x=x, y=df_sorted[sort_by[cntr]], text=df_sorted['Names'],  marker=dict(color=colors[cntr]),
                name="% " + sort_by[cntr], mode="markers"))
            cntr += 1
        fig.update_layout(clickmode='event+select')    
        fig.update_layout(hovermode='x unified')
        print("query", anchor_name, chrs, start_coord, end_coord)
        local_genes = index.query_genes(anchor_name, chrs, start_coord, end_coord)
        local_gene_list = [g.split(';')[0].split("=")[1] for g in local_genes['attr']]
        df2 = df_sorted.loc[df_sorted['Names'].isin(local_gene_list)]

        fig.add_trace(go.Scattergl(x=df2['X'], y=df2['Universal'], marker=dict(color='#FF2192', size=10), 
            mode="markers", hoverinfo='skip', showlegend=False))
        fig.add_trace(go.Scattergl(x=df2['X'], y=df2['Unique'], marker=dict(color='#FF2192', size=10), 
            mode="markers", hoverinfo='skip', showlegend=False))

        fig.add_hline(y=uniq_avg, line_dash='dash', line_color='goldenrod')
        fig.add_hline(y=univ_avg, line_dash='dash', line_color='#440154')
        #fig.update_layout(height=500)
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
                size=16,
            ),
        )
        return fig

    #Tab 1 bar chart of each genome 
    def read_pangenome_comp():

        fig = make_subplots(rows=1, cols=1)
        for i,g in enumerate(labels):
            fig.add_trace(go.Bar(y=labels, x=index.genome_occ_freq["total", i+1], name=g, #orientation='h',
                legendgrouptitle_text="# Samples",
                marker=dict(color=colors[i]), 
                marker_line=dict(color=colors[i]), orientation='h',
                ), row=1, col=1)
        fig.update_layout(barmode='stack' )#,orientation='h')
        fig.update_layout(font=dict( #height=1500, 
            size=26,
            ), plot_bgcolor='rgba(0,0,0,0)')
        fig.update_layout(
            updatemenus = [
            dict(type="buttons",direction="left", #showactive=True, x=0, y=0, 
                pad={"r": 10, "t": 10},
                showactive=True,
                x=0,#0.11,
                xanchor="left",
                y=1.05,
                yanchor="top",
                buttons=list([
                dict(args=[{'xaxis.type': 'linear'}],label="Linear", method="relayout"
                    ),
                dict(args=[{'xaxis.type': 'log'}],label="Log",method="relayout"
                    )
                ])
            ),])
        return fig#, genome_comp_totals

    def read_genome_comp(anchor_name):
        totals = index.chr_occ_freq.loc[anchor_name,"total"]

        fig = make_subplots(rows=len(totals), cols=1)

        for i,(c,freqs) in enumerate(totals.iterrows()):
            perc = freqs*100
            fig.add_trace(go.Bar(x=perc.index, y=perc, marker_color=colors, showlegend=False,), row=i+1,col=1)
        fig.update_yaxes(type="log")
        fig.update_layout(paper_bgcolor='rgba(0,0,0,0)')
        return fig

    #Tab 1 dendogram 
    def make_all_genome_dend():
        dist_mat = np.zeros((index.ngenomes, index.ngenomes), np.float64)

        with open(index.genome_dist_fname) as f:
            for line in f:
                f, t, d, p, x = line.rstrip().split("\t")
                i = index.genomes.get_loc(f)
                j = index.genomes.get_loc(t)
                dist_mat[i][j] = d
                dist_mat[j][i] = d

        df = pd.DataFrame(dist_mat, columns=index.genomes)
        branch_colors = ['purple','purple','purple','purple','purple','purple']
        fig = ff.create_dendrogram(df, colorscale=branch_colors, labels=df.columns, orientation='bottom' ) #, orientation='bottom'
        fig.update_layout(paper_bgcolor='rgba(0,0,0,0)', font=dict(
                size=10,))
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
        return fig

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

    def make_genome_count_plot(anchor_name):
        counts = index.genome_sizes["chr_count"].sort_values()
        fig = make_subplots(rows=1, cols=1)
        fig.add_trace(go.Scattergl(x=counts.index,y=counts))
        fig.update_yaxes(title_text="# of scaffolds",)
        fig.add_vline(x=anchor_name, line_dash="dash", line_color="darkblue")
        fig.update_layout(font=dict(size=20,))
        return fig

    def make_genome_size_plot(anchor_name):
        lens = index.genome_sizes["length"].sort_values()
        fig = make_subplots(rows=1, cols=1)
        fig.add_trace(go.Scattergl(x=lens.index,y=lens))
        fig.update_yaxes(title_text="Size of genome",)
        fig.add_vline(x=anchor_name, line_dash="dash", line_color="darkblue")
        fig.update_layout(font=dict(size=20,))

        return fig

    def plot_avgs(anchor_name):
        fig = make_subplots(rows=1, cols=1)
        fig.add_trace(go.Scattergl(x=index.genome_occ_avg.index,y=index.genome_occ_avg))
        fig.add_vline(x=anchor_name, line_dash="dash", line_color="darkblue")
        fig.update_yaxes(title_text="Average k-mer",)
        fig.update_layout(font=dict(size=20,))
        return fig

    def make_avg_kmer_fig(this_anchor):
        fig = make_subplots(rows=1, cols=1)

        avgs = index.chr_occ_avg[this_anchor]

        fig.add_trace(go.Scattergl(x=avgs.index, y=avgs))
        fig.update_yaxes(title_text="Average k-mer",)
        fig.update_xaxes(title_text="Chromosome",)
        fig.update_layout(plot_bgcolor='rgba(0,0,0,0)',  font=dict(size=20))

        fig.update_layout(
            updatemenus = [
            dict(type="buttons",direction="left", #showactive=True, x=0, y=0, 
                pad={"r": 10, "t": 10},
                showactive=True,
                x=0,#0.11,
                xanchor="left",
                y=1.2,
                yanchor="top",
                buttons=list([
                dict(args=[{'yaxis.type': 'linear'}],label="Linear", method="relayout"
                    ),
                dict(args=[{'yaxis.type': 'log'}],label="Log",method="relayout"
                    )
                ])
            ),])

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
        fig.update_layout(
            updatemenus = [
            dict(type="buttons",direction="left", #showactive=True, x=0, y=0, 
                pad={"r": 10, "t": 10},
                showactive=True,
                x=0,#0.11,
                xanchor="left",
                y=1.2,
                yanchor="top",
                buttons=list([
                dict(args=[{'yaxis.type': 'linear'}],label="Linear", method="relayout"
                    ),
                dict(args=[{'yaxis.type': 'log'}],label="Log",method="relayout"
                    )
                ])
            ),])
        return fig

    def make_gene_per_genome_fig(anchor_name):
        
        fig = make_subplots(rows=1, cols=1)
        colors = ['#ffd60a', '#440154']
        gene_names = []
        universals = []
        uniques = []
        sizes = []
        genes = list()
        for chrom in index.chrs.loc[anchor_name].index:
            try:
                g = index.query_genes(anchor_name, chrom, 0, index.chrs.loc[anchor_name, chrom]["size"])
            except ValueError: 
                continue
            genes.append(g)
        genes = pd.concat(genes)
        genes["name"] = genes["attr"].str.extract("Name=([^;]+)")
        genes["size"] = genes["end"] - genes["start"]
        
        x = [i for i in range(0, len(genes['universal']))]
        genes['universal'] = genes['universal']/genes["size"]
        genes['unique'] = genes['unique']/genes["size"]
        df_sorted = genes.sort_values('universal')
        df_sorted['X'] = x
        fig.add_trace(go.Scattergl(x=x, y=df_sorted['unique'], text=df_sorted["chr"] + ":" + df_sorted['name'], marker=dict(color=colors[0]),
                name="Proportion "+"Unique", mode="markers"))
        fig.add_trace(go.Scattergl(x=x, y=df_sorted['universal'], text=df_sorted["chr"] + ":" +df_sorted['name'], marker=dict(color=colors[1]),
                name="Proportion "+"Universal", mode="markers"))
        fig.update_layout(plot_bgcolor='rgba(0,0,0,0)', font=dict(size=20) )
        fig.update_xaxes(title_text="Genes",)
        fig.update_yaxes(title_text="Proportion of gene",)
        fig.update_layout(hovermode='x unified')

        return fig
    
    def plot_anno_conserve(anchor_name):
        file_name="gene_var.txt"
        all_avgs = []
        all_stds = []
        all_names = []
        all_lens = []

        d_anno = {"Standard_dev":all_stds, #np.log10(all_stds),#all_stds,
            "Average":all_avgs,
            "Name":all_names,
            "Length":all_lens,
            #"Sweep":is_sweep
        }
        df_anno = pd.DataFrame(d_anno)
        fig = px.scatter(df_anno,x="Average", y="Standard_dev", marginal_y='histogram', 
            marginal_x='histogram',
            hover_name="Name", 
            #color="Sweep",
            opacity=0.5,
            color=np.log10(df_anno["Length"])
        )
        return fig

    ##### READ DATA
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
    t = time.time()
    cnts_tmp = []
    x = []
    all_chrs = {} #This is where we keep the raw, str, counts 
    
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
        all_chrs[chrs] = index.query_bitmap(anchor_name, chrs, chr_start, chr_end, n_skips_start)
    chrs = index.chrs.loc[anchor_name].index[0] #"chr1"

    #Parsing the counts will return the number of genomes present at each kmer position, based on the counts
    names_simp = all_chrs[chrs].sum(axis=1)#tmp.sum(axis=1)
    bar_sum_global = {}

    for l in labels:
        bar_sum_global[l] = {}
        cntr = 0
        if l in index.anchor_genomes:
            for c in index.chrs.loc[l].index:
                counts = index.chr_bins.loc[(l,c)].sum(axis=0)
                bar_sum_global[l][c] = counts.to_numpy()#bar_sum_global_tmp
                cntr +=1

    tmp = 0
    
    #Read in the repeats file (this has the annotated repeats from Shujun)     
    zs_tmp = []
    gene_names = {}
    gene_locals = {}
    exon_locals = {}
    exon_names = {}
    all_rep_types = {}
    gene_content = {}

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
    #if len(index.chr_bins.loc[(anchor_name,chrs)]) > 1:
    genes = index.query_genes(anchor_name, chrs, 0, index.chrs.loc[anchor_name, chrs]["size"])
    #[g.split(';')[1].split('=')[1] for g in genes['attr']]
    bounds = genes.loc[:, ["start", "end"]]
    bounds["break"] = None #pd.NA
    #gene_locals[l][chrom] = bounds.to_numpy().flatten()
    gene_names = [g.split(';')[0].split("=")[1] for g in genes['attr']]

    main_fig, bar_sum_names, bar_sum_regional, colors, gene_names_tmp = plot_interactive(n_skips_start,
        layout, bins, names_simp[x_start_init_adjust:x_stop_init_adjust], chrs, zs_tmp, True, True, 
        x_start_init, x_stop_init, bounds.to_numpy().flatten(), gene_names, anchor_name, chrs
    )

    local_genes = len(index.query_genes(anchor_name, chrs, x_start_init, x_stop_init))

    sort_by = ["Unique","Universal"]
    tmp = [(i/sum(bar_sum_global[anchor_name][chrs])*100) for i in bar_sum_global[anchor_name][chrs]]
    uniq_avg = tmp[1]
    univ_avg = tmp[-1]

    universals = genes["universal"]
    uniques = genes["unique"]

    sizes = genes["end"]-genes["start"] #genes["size"]
    gene_content_plot = plot_gene_content(gene_names,universals,uniques,sizes,sort_by, 
        colors, uniq_avg, univ_avg, x_start_init, x_stop_init, anchor_name, chrs)

    x, z_1, z_9, y = {}, {}, {}, {}

    #Now we need to make the pangenome plot
    #This is going to be a bar chart 

    chr_fig                     = plot_chr_whole(x_start_init, x_stop_init, anchor_name, chrs)
    whole_genome_fig            = plot_whole_genome(anchor_name)
    whole_genome_hists_fig      = read_genome_comp(anchor_name)
    all_genomes_dend_fig        = make_all_genome_dend()
    pangenome_comp_fig          = read_pangenome_comp()
    pangenome_avg_fig           = plot_avgs(anchor_name)
    pg_size_fig                 = make_genome_size_plot(anchor_name)
    pg_count_fig                = make_genome_count_plot(anchor_name)
    avg_kmer_per_chr_fig        = make_avg_kmer_fig(anchor_name)
    genes_per_chr_fig           = make_genes_per_chr_fig(anchor_name)
    gene_content_per_genome_fig = make_gene_per_genome_fig(anchor_name)

    chrs = index.chrs.loc[anchor_name].index[0] #"chr1"

    config = {"toImageButtonOptions" : {"format" : "svg", "width" : None, "height" : None} , "scrollZoom": False}
    tab_style = {"background-color": "lightgrey", "font-size": 36}

    PANGENOME_TAB = [
            html.Div(id="PanInfo1", children=[
                html.I("Panagram of " + str(num_samples) + " genomes"),
                html.Br(),
                html.Label('Select a genome: '),
                dcc.Dropdown(labels, anchor_name, style=dict(width='40%', height='110%', verticalAlign="middle"), id="genome_select_dropdown"),
                html.Br(),
            ], style={'padding':'1%', 'font-size':'24px', 'textAlign': 'left', "border":"2px grey solid"}),
            html.Div(className="w3-half", children=[
                dcc.Graph(id="all_genomes",style={"font-size": 20, "height" : 1500},
                    config=config #, figure = all_genomes_dend_fig
                ),
            ]),
            html.Div(className="w3-half", children=[
                dcc.Graph(id="all_genomes_kmer_comp",
                    config=config, style={"font-size": 30, "height" : 1500}#, figure = pangenome_comp_fig
                )
            ]),
            html.Div(className="w3-third", children=[
                dcc.Graph(id="all_genome_sizes",
                    config=config, style={"font-size": 30, "height" : 500}#,figure = pg_size_fig, 
                )
            ]),
            html.Div(className="w3-third", children=[
                dcc.Graph(id="average_kmer_content",
                    config=config, style={"font-size": 30, "height" : 500}#figure = pangenome_avg_fig,
                )
            ]),
            html.Div(className="w3-third", children=[
                dcc.Graph(id="pangenome_num_seqs",
                    config=config, style={"font-size": 30, "height" : 500}#,figure = pg_count_fig
                )
            ])]

    #def render_anchor_tab():
    #    return 
    ANCHOR_TAB = [
            html.Div(children=[
                    html.I("Anchor genome " + anchor_name, id="anchor_labels"),
                    html.Br(),
                    html.Label('Select a chromosome: '),
                    dcc.Dropdown(chrs_list[anchor_name], chrs, style=dict(width='40%', height='110%', verticalAlign="middle", ), id="chr_select_dropdown"),
                    html.Br(),
                ], style={'padding-top' : '1%', 'padding-left' : '1%', 'padding-bottom' : '1%', 'padding-right' : '1%', 'font-size':'24px', 'textAlign': 'left', "border":"2px grey solid", }),
            html.Div(className="w3-third",children=[
                dcc.Graph(id="gene_content",
                #figure = gene_content_per_genome_fig, 
                config=config), #gene_content_fig)
            ],), 
            html.Div(className="w3-third",children=[
                dcc.Graph(id="genes_per_chr",
                #figure = genes_per_chr_fig, 
                config=config)
            ],),
            html.Div(className="w3-third",children=[
                dcc.Graph(id="avg_kmer_chr", config=config) #, figure = avg_kmer_per_chr_fig
            ],),
            html.Div(className="w3-threequarter",children=[
                dcc.Graph(id="all_chromosomes",
                    #figure = whole_genome_fig, 
                    config=config, style={"height" : num_chrs[anchor_name]*250})
            ]),
            html.Div(className="w3-quarter",children=[
                dcc.Graph(id="all_chromosomes_hists", style={"height" : num_chrs[anchor_name]*250},
                    figure = whole_genome_hists_fig, 
                    config=config
                )
            ]),
        ]

    #def render_chromosome_tab():
    #    return [
    CHROMOSOME_TAB = [
            html.Div(children=[ # info panel
                html.Div(children=[ #left side
                    html.I(anchor_name + "." + chrs + ":", id="chr_name"),
                    dcc.Input(id="Chrs_Info", placeholder=str(x_start_init) + "-" + str(x_stop_init), debounce=True),
                    html.Br(),
                    html.I("Genes in this region: " + str(local_genes), id='regional_genes'),
                    html.Br(),
                    html.Label('Pre-selected regions: '),
                    dcc.Dropdown(pre_bed_info[anchor_name], chrs + ":" + str(x_start_init) + "-" + str(x_stop_init) , style=dict(width='100%', height='110%', verticalAlign="middle", ), id="pre_regions_dropdown"),
                ],style={"display": "inline-block"}),
                
                html.Div(children=[#right side
                    html.I("K-mer length: " + str(kmer_len),style={"display": "inline-block",}),
                    html.Br(),
                    html.I("Number of bins: " + str(bins), style={'display': 'inline-block'}),
                    html.Br(),
                    html.I("Step size: " + str(n_skips_start), style={"display": "inline-block",}, id='step_size_out'),
                    html.Br(),
                ], style={"display": "inline-block", 'padding-left' : '10%'}),

            ], style = {'padding-top' : '1%', 'padding-left' : '1%', 'padding-bottom' : '1%', 'padding-right' : '1%', 'font-size':'24px', "border":"2px grey solid"}),

            html.Div(children=[ #summary 
                html.Div(className="w3-container", children=[
                        #left figure
                        dcc.Graph(id="chromosome",
                            config=config, style={"font-size": 20, "height" : 350}) #figure = chr_fig, 
                ])
            ]),

            html.Div(children=[ 
                html.Div(className="w3-container", children=[
                    html.Div(className="w3-threequarter", children=[
                        #left figure - calling this the "Main" figure
                        dcc.Graph(id="primary", config=config,  #,figure = main_fig
                            style={"height": 1000,  "font-size": 20})
                    ]), #style={'padding-top' : '1%', 'padding-left' : '1%', 'padding-bottom' : '1%', 'padding-right' : '1%', 
                    html.Div(className="w3-quarter", children=[
                        dcc.Graph(id="Secondary", 
                            #Now we have the phylogenetic tree
                            #figure = get_local_info(names_simp[x_start_init:x_stop_init], #exon_comp[chrs], 
                            #gene_comp[anchor_name][chrs], 
                            #bar_sum_regional, bar_sum_global[anchor_name][chrs], anchor_name, chrs), 
                            config=config,
                            style={"height": 1000, "font-size": 20}),
                    ])
                ])
            ]),
            html.Div(children=[
                html.Div(className="w3-container", children=[
                    html.Div(className="w3-third", children=[
                        dcc.Graph(id="Genes", 
                            #figure=gene_content_plot, 
                            config=config,
                            style={"font-size": 20, "height":750}),
                    ]),
                    html.Div(className="w3-twothird", children=[
                        dcc.Graph(id="Third", 
                            #This is the histogram section
                            #figure = create_tree( x_start_init,x_stop_init, all_chrs[chrs][int(x_start_init/n_skips_start):int(x_stop_init/n_skips_start)], n_skips_start),
                            config=config,
                            style={"font-size": 20, "height":750}),
                    ]),
                ])
            ])
        ]

    #def render_annotation_tab():
    #    return [
    ANNOTATION_TAB = html.Div(children=[
            dcc.Graph(id="annotation_conservation",
                figure=plot_anno_conserve(anchor_name),
                config=config,
                style={"height": 1000, "font-size": 20}
            )
        ])

    #TAB_RENDERERS = {
    #    "pangenome"  :   render_pangenome_tab,
    #    "anchor"      :  render_anchor_tab,
    #    "chromosome" :   render_chromosome_tab,
    #    "annotation"  :  render_annotation_tab
    #}

    app = dash.Dash(
        external_stylesheets=["https://www.w3schools.com/w3css/4/w3.css"],
        url_base_pathname=params.url_base,
        suppress_callback_exceptions=True
    ) 
    #app = Dash(__name__, external_stylesheets=[dbc.themes.PULSE, dbc_css])
    app.layout = html.Div([
        #html.Div(id = 'parent', children = [
        html.H1(id = 'H1', children = 'Panagram', style = {'textAlign':'center', "font-size": 64}), 
        dcc.Tabs(id="tabs", value="pangenome", style={"font-size": 36}, children=[
            dcc.Tab(value="pangenome", label='Pangenome', style=  tab_style, children=PANGENOME_TAB),
            dcc.Tab(value="anchor", label='Anchor genome', style= tab_style, children=ANCHOR_TAB),
            dcc.Tab(value="chromosome", label='Chromosome', style=tab_style, children=CHROMOSOME_TAB), 
            dcc.Tab(value="annotation", label='Annotation', style=tab_style, children=ANNOTATION_TAB), 
        ]),

        html.Div(id="tab-content"),

        html.Div(chrs, style={"display" : "none"}, id="selected-chrom-state"),
        html.Div(anchor_name, style={"display" : "none"}, id="selected-anchor-state"),
        
        html.Div(x_start_init,id='start-coord-state',style={"display" : "none"} ),
        html.Div(x_stop_init,id='end-coord-state',style={"display" : "none"} ),
    ])

    def get_buffer(tmp_start, tmp_stop, n_skips):
        min_dist = n_skips*bins*8
        act_dist = tmp_stop-tmp_start
        if act_dist >= min_dist:
            return 0 #The distance between start and stop are good enough 
        else: #
            return 1 #int((min_dist - act_dist)/2)+1 
    last_reset_click = 0

    def set_coords(anchor, chrom=None, start=0, end=None):
        anchor_chrs = index.chrs.loc[anchor].index
        if chrom not in anchor_chrs:
            chrom = anchor_chrs[0]
        if end is None:
            end = index.chrs.loc[(anchor, chrom), "size"]
        return anchor, chrom, start, end

    @app.callback(
        Output('selected-anchor-state', 'children'), #pangenome
        Output('selected-chrom-state', 'children'),  #anchor
        Output('start-coord-state','children'),   #start_coord div (constant?)
        Output('end-coord-state','children'),   #x_end div (constant?)
        Output('chromosome', 'relayoutData'),        #chromosome
        Output('Chrs_Info', 'value'),                #chromosome
        Output('chr_name', 'children'),              #chromosome

        Input('genome_select_dropdown', 'value'),     #pangenome, anchor, chromosome
        Input('chr_select_dropdown', 'value'),   #anchor, chromosome
        Input('all_chromosomes','relayoutData'), #anchor -> chromosome
        Input('primary', 'clickData'),           #chromosome
        Input('primary', 'relayoutData'),        #chromosome
        Input('chromosome', 'selectedData'),     #chromosome
        Input('Genes', 'clickData'),             #chromosome
        Input('Chrs_Info', 'value'),             #chromosome
        Input('pre_regions_dropdown', 'value'),  #chromosome

        #Input('all_chromosomes','relayoutData'), #anchor -> chromosome

        State('selected-anchor-state', 'children'),
        State('selected-chrom-state', 'children'),
        State('start-coord-state', 'children'),
        State('end-coord-state', 'children'),
    )
    def nav_callback(anchor_dropdown, chr_dropdown, anctab_chrs_relayout, chrtab_primary_click, chrtab_primary_relayout, chrtab_chr_select, chrtab_gene_click, user_chr_coords, pre_selected_region, anchor, chrom, start_coord, end_coord):
        triggered_id = ctx.triggered_id

        print("TRIG", triggered_id)
        n_skips = 100
        click_me_genes = True
        click_me_rep = True

        if triggered_id == "genome_select_dropdown":
            anchor, chrom, start_coord, end_coord = set_coords(anchor_dropdown)

        elif triggered_id == "chr_select_dropdown":
            anchor, chrom, start_coord, end_coord = set_coords(anchor, chr_dropdown)

        elif triggered_id == "pre_regions_dropdown":
            if pre_selected_region != None:
                chrom = pre_selected_region.split(":")[0]

                start_coord = int(pre_selected_region.split(":")[1].split("-")[0])
                end_coord = int(pre_selected_region.split(":")[1].split("-")[1])
                n_skips = 1

        #Select start/end coordinates, triggers CHROMOSOME
        elif triggered_id == 'Chrs_Info':
            start_coord = int(user_chr_coords.strip().split('-')[0]) 
            end_coord  = int(user_chr_coords.strip().split('-')[1])

        #Select chromosome (region) in anchors tab, triggers CHROMOSOME
        elif triggered_id == 'all_chromosomes':
            #Get the chromosome number 
            if anctab_chrs_relayout != None and len(anctab_chrs_relayout)>1:
                chr_num_tmp = list(anctab_chrs_relayout.keys())[0].split('.')[0].split('axis')[1] #if anctab_chrs_relayout != None and len(anctab_chrs_relayout)>1:
                if len(chr_num_tmp) == 0:
                    chr_num = 1
                else:
                    chr_num = int(int(chr_num_tmp)/3)+1
                
                chrom = index.chrs.index[chr_num-1]
                anchor, chrom, start_coord, end_coord = set_coords(anchor, chrom)
                print("set", anchor, chrom, start_coord, end_coord)

        #Chromosome gene plot, triggers CHROMOSOME
        elif triggered_id == 'Genes':
            if chrtab_gene_click != None:
                this_gene_name = chrtab_gene_click['points'][0]['text']

                genes = index.query_genes(anchor, chrom, int(start_coord), int(end_coord))
                genes["name"] = genes["attr"].str.extract("Name=([^;]+)")
                
                genes["name"] = genes["attr"].str.extract("Name=([^;]+)")
                this_gene = genes[genes["name"] == this_gene_name]
                print(this_gene)

                start_coord = int(this_gene['start'])
                end_coord = int(this_gene['end'])

        #Chromosome top plot, triggers CHROMOSOME
        elif triggered_id == 'chromosome':
            if not "range" in chrtab_chr_select:
                return (no_update,)*7
            #chromosome_gene_triggered(gene_jump_bottom, n_skips, click_me_genes, click_me_rep, chrtab_chr_select, chrom, anchor_name)
            print("CRHOM", chrtab_chr_select)
            if "x2" in chrtab_chr_select['range'].keys():
                start_coord = int(chrtab_chr_select['range']['x2'][0])
                end_coord = int(chrtab_chr_select['range']['x2'][1])
            elif "x" in chrtab_chr_select['range'].keys():
                start_coord = int(chrtab_chr_select['range']['x'][0])
                end_coord = int(chrtab_chr_select['range']['x'][1])
            else:
                start_coord = int(chrtab_chr_select['range']['x3'][0])
                end_coord = int(chrtab_chr_select['range']['x3'][1])
        

        #Chromosome main plot, triggers CHROMOSOME
        elif triggered_id == 'primary':
            if chrtab_primary_relayout != None and 'xaxis4.range[0]' in chrtab_primary_relayout.keys():
                start_coord = int(chrtab_primary_relayout['xaxis4.range[0]'])
                end_coord = int(chrtab_primary_relayout['xaxis4.range[1]'])
            else:
                print("OTHER CLICK", chrtab_primary_click)
            #elif chrtab_primary_click != None:

        if (end_coord - start_coord) <= bins:
            start_coord = start_coord - 175
            end_coord = end_coord + 175

        print("coordinates are",  anchor, chrom, start_coord, end_coord)

        return (
            anchor, chrom, start_coord, end_coord,
            chrtab_chr_select,
            update_output_div(chrom, start_coord, end_coord, anchor), 
            update_out_chr(chrom, anchor), 
        )


    @app.callback(
        Output('all_genomes', 'figure'),
        Output('all_genomes_kmer_comp', 'figure'),
        Output('pangenome_num_seqs', 'figure'),      #pangenome
        Output('average_kmer_content', 'figure'),    #pangenome
        Output('all_genome_sizes', 'figure'),        #pangenome

        Input('tabs', 'value'),
        Input('selected-anchor-state', 'children')
    )
    def pangenome_callback(tab, anchor_name):
        if tab != "pangenome":
            return ({},)*5 #+ (no_update,)
        return all_genomes_dend_fig, pangenome_comp_fig, make_genome_count_plot(anchor_name), plot_avgs(anchor_name), make_genome_size_plot(anchor_name)#, anchor_name

    @app.callback(
        Output('all_chromosomes', 'figure'),         #anchor
        Output('gene_content', 'figure'),            #anchor
        Output('genes_per_chr', 'figure'),           #anchor
        Output('avg_kmer_chr', 'figure'),            #anchor
        Output('chr_select_dropdown', 'options'),    #anchor
        Output('anchor_labels', 'children'),         #anchor
        #Output('selected-chrom-state', 'children'),  #anchor

        Input('tabs', 'value'),
        Input('selected-anchor-state', 'children')
    )
    def anchor_callback(tab, anchor_name):
        if tab != "anchor":
            return ({},)*4 + (no_update,)*2
        triggered_id = ctx.triggered_id

        return (
            plot_whole_genome(anchor_name),
            make_gene_per_genome_fig(anchor_name),
            make_genes_per_chr_fig(anchor_name),
            make_avg_kmer_fig(anchor_name),
            update_chromosome_list(anchor_name), 
            anchor_name, 
        )
        
    @app.callback(
        Output('chromosome','figure'),               #chromosome
        Output('primary','figure'),                  #chromosome
        Output('Secondary', 'figure'),               #chromosome
        Output('Third', 'figure'),                   #chromosome
        Output('Genes', 'figure'),                   #chromosome
        Output('regional_genes', 'children'),        #chromosome

        Input('tabs', 'value'),                  #all
        Input('start-coord-state','children'),   #start_coord div (constant?)
        Input('end-coord-state','children'),   #x_end div (constant?)

        State('selected-chrom-state', 'children'),
        State('selected-anchor-state', 'children'),
    )
    def chromosome_callback(tab, start_coord, end_coord, chrs, anchor_name):

        if tab != "chromosome":
            return ({},)*5 + (no_update,)#+ (no_update,)*4

        triggered_id = ctx.triggered_id

        n_skips = 100
        click_me_genes = True
        click_me_rep = True
        chr_num = chrs_list[anchor_name].get_loc(chrs)+1

        #Bookmarks, select from chromosome, triggers CHROMOSOME

        return update_all_figs(chr_num, click_me_rep, click_me_genes, chrs, anchor_name, 0, start_coord, end_coord,n_skips) 

    def user_chr_coords_triggered(new_x_start, new_x_stop, click_me_genes, click_me_rep, chr_relayoutData, chrs, n_skips_start, anchor_name):

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
            n_skips = n_skips_start

        chr_num = chrs_list[anchor_name].get_loc(chrs)
        genes = index.query_genes(anchor_name, chrs, 0, index.chrs.loc[anchor_name, chrs]["size"])

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
        local_gene_list = index.query_genes(anchor_name, chrs, int(new_x_start), int(new_x_stop))
        local_gene_list["name"] = local_gene_list["attr"].str.extract("Name=([^;]+)")
        #genes = index.query_genes(anchor_name, chrs, 0, index.chrs.loc[anchor_name, chrs]["size"])
        
        universals = genes["universal"]
        uniques = genes["unique"]
        sizes = genes["end"]-genes["start"]
        fig4 = plot_gene_content(gene_names,universals,uniques,sizes, sort_by, colors, uniq_avg, 
                                 univ_avg, int(new_x_start), int(new_x_stop), anchor_name, chrs)
        return (
            chr_fig, fig1, fig2, fig3, fig4, chr_relayoutData, 
        )

    def sorted_gene_fig(gene_jump_bottom, n_skips, click_me_genes, click_me_rep, chr_relayoutData, chrs, anchor_name):
        chr_num = chrs_list[anchor_name].get_loc(chrs)
        #chr_num = chrs_list[anchor_name].index(chrs)
        names_simp = index.query_bitmap(anchor_name, chrs, 0, index.chrs.loc[anchor_name, chrs]["size"], n_skips_start).sum(axis=1) #all_chrs[chrs].sum(axis=1)
        #int(chrs.split('r')[1])
        if gene_jump_bottom != None: #and gene_jump_bottom['points'][0]['text'] in gene_names[anchor_name][chrs]:
            
            genes = index.query_genes(anchor_name, chrs, 0, index.chrs.loc[anchor_name, chrs]["size"])

            #[g.split(';')[1].split('=')[1] for g in genes['attr']]
            bounds = genes.loc[:, ["start", "end"]]
            bounds["break"] = None
            gene_locals = bounds.to_numpy().flatten()
            genes["name"] = genes["attr"].str.extract("Name=([^;]+)")
            
            this_gene = genes[genes["name"] == this_gene_name]
            this_start = int(this_gene['start'])-1
            this_stop = int(this_gene['end'])+1
            if get_buffer(this_start,this_stop, n_skips_start) == 1:
                #we need to move to single nucleotide resoltion               
                if (this_stop - this_start) <= bins:
                    start_coord = this_start - 175
                    end_coord = this_stop + 175
                else:
                    start_coord = this_start-1
                    end_coord = this_stop+1
                exact_mini_counts = index.query_bitmap(anchor_name, chrs, start_coord, end_coord, 1)
                simple_cnts_for_plots = exact_mini_counts.sum(axis=1)
                #simple_cnts_for_plots, exact_mini_counts = read_mini_count_files(start_coord, end_coord)
                n_skips = 1
            else: 
                start_coord = this_start - buff #int(1000*n_skips/2)
                end_coord = this_stop + buff #int(1000*n_skips/2)
                n_skips = n_skips_start
                simple_cnts_for_plots = names_simp[int(start_coord/n_skips):int(end_coord/n_skips)]
                exact_mini_counts = all_chrs[chrs][int(start_coord/n_skips):int(end_coord/n_skips)]
            
            fig1, bar_sum_names, bar_sum_regional, colors, gene_names_tmp = plot_interactive( n_skips, #int(SG_window), 
                #int(SG_polynomial_order), SG_check,
                layout, #exon_comp[chrs], 
                bins, simple_cnts_for_plots, chrs, zs_tmp, 
                click_me_rep, click_me_genes,  
                int(start_coord), int(end_coord), gene_locals, 
                [g.split(';')[0].split('=')[1] for g in genes['attr']], #gene_names[anchor_name][chrs], 
                anchor_name, chrs
                #exon_locals[anchor_name][chrs], exon_names[anchor_name][chrs],
                )
            fig2 = get_local_info(
                simple_cnts_for_plots, #exon_comp[chrs], 
                #gene_comp[anchor_name][chrs], 
                bar_sum_regional, bar_sum_global[anchor_name][chrs], anchor_name, chrs)
            fig3 = create_tree(start_coord, end_coord, exact_mini_counts, n_skips)
            gene_jump = ""
            if x_start_init != start_coord:
                chr_fig = plot_chr_whole(
                    start_coord, end_coord, anchor_name, chrs)
        #local_gene_list = []
        #cntr = 0
        #while cntr < len(gene_names_tmp):
        #    local_gene_list.append(gene_names_tmp[cntr])
        #    cntr += 3
        local_gene_list = index.query_genes(anchor_name, chrs, int(start_coord), int(end_coord))
        local_gene_list["name"] = local_gene_list["attr"].str.extract("Name=([^;]+)")

        genes = index.query_genes(anchor_name, chrs, 0, index.chrs.loc[anchor_name, chrs]["size"])
        gene_names = [g.split(';')[0].split("=")[1] for g in genes['attr']]
        universals = genes["universal"]
        uniques = genes["unique"]
        sizes = genes["end"]-genes["start"]
        fig4 = plot_gene_content(gene_names,universals,uniques,sizes,
            sort_by, colors, uniq_avg, univ_avg, int(start_coord), int(end_coord),anchor_name, chrs)
        return (
            chr_fig, fig1, fig2, fig3, fig4
        )
    
    def chromosome_gene_triggered(gene_jump_bottom,n_skips, click_me_genes, click_me_rep, chr_relayoutData, chrs, anchor_name):
        chr_num = chrs_list[anchor_name].get_loc(chrs)
        names_simp = index.query_bitmap(anchor_name, chrs, 0, index.chrs.loc[anchor_name, chrs]["size"], n_skips_start).sum(axis=1) #all_chrs[chrs].sum(axis=1)
        
        if "x2" in chr_relayoutData['range'].keys():
            start_coord = int(chr_relayoutData['range']['x2'][0])
            end_coord = int(chr_relayoutData['range']['x2'][1])
        elif "x" in chr_relayoutData['range'].keys():
            start_coord = int(chr_relayoutData['range']['x'][0])
            end_coord = int(chr_relayoutData['range']['x'][1])
        else:
            start_coord = int(chr_relayoutData['range']['x3'][0])
            end_coord = int(chr_relayoutData['range']['x3'][1])
        if get_buffer(start_coord, end_coord, n_skips_start) == 1:
            #we need to move to single nucleotide resoltion
            if (end_coord - start_coord) <= bins:
                start_coord = start_coord - 175
                end_coord = end_coord + 175
            exact_mini_counts = index.query_bitmap(anchor_name, chrs, start_coord, end_coord, 1)
            simple_cnts_for_plots = exact_mini_counts.sum(axis=1)
            n_skips = 1
        else:
            exact_mini_counts = index.query_bitmap(anchor_name, chrs, int(start_coord/n_skips_start), int(end_coord/n_skips_start), n_skips_start)#all_chrs[chrs][int(start_coord/n_skips_start):int(end_coord/n_skips_start)]
            simple_cnts_for_plots = names_simp[int(start_coord/n_skips_start):int(end_coord/n_skips_start)]
            n_skips = n_skips_start
        
        genes = index.query_genes(anchor_name, chrs, start_coord, end_coord)
        #[g.split(';')[1].split('=')[1] for g in genes['attr']]
        bounds = genes.loc[:, ["start", "end"]]
        bounds["break"] = None
        fig1, bar_sum_names, bar_sum_regional, colors, gene_names_tmp = plot_interactive( n_skips, #int(SG_window), 
                layout, #exon_comp[chrs], 
                bins, simple_cnts_for_plots, chrs, zs_tmp,  
                click_me_rep, click_me_genes,  
                int(start_coord), int(end_coord), 
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
        fig3 = create_tree( start_coord, end_coord, exact_mini_counts, n_skips)
        chr_relayoutData = None#{'autosize': True}
        if x_start_init != start_coord:
            chr_fig = plot_chr_whole(
                start_coord, end_coord, anchor_name, chrs)
        #local_gene_list = []
        #cntr = 0
        #while cntr < len(gene_names_tmp):
        #    local_gene_list.append(gene_names_tmp[cntr])
        #    cntr += 3
        local_gene_list = genes #index.query_genes(anchor_name, chrs, int(start_coord), int(end_coord))
        local_gene_list["name"] = local_gene_list["attr"].str.extract("Name=([^;]+)")
        genes = index.query_genes(anchor_name, chrs, 0, index.chrs.loc[anchor_name, chrs]["size"])
        gene_names = [g.split(';')[0].split("=")[1] for g in genes['attr']]
        universals = genes["universal"]
        uniques = genes["unique"]
        sizes = genes["end"]-genes["start"]
        fig4 = plot_gene_content(gene_names,universals,uniques,sizes, 
            sort_by, colors, uniq_avg, univ_avg, start_coord, end_coord, anchor_name, chrs)
        return (
            chr_fig, fig1, fig2, fig3, fig4
        )

    def primary_fig_triggered(gene_jump_bottom, n_skips, click_me_genes, click_me_rep, chr_relayoutData, relayoutData, clickData, start_coord, end_coord, chrs, anchor_name):
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
            start_coord = int(relayoutData['xaxis4.range[0]'])
            end_coord = int(relayoutData['xaxis4.range[1]'])

        elif clickData != None:
            this_gene_name = clickData['points'][0]['text']
            
            genes["name"] = genes["attr"].str.extract("Name=([^;]+)")
            this_gene = genes[genes["name"] == this_gene_name]

            start_coord = int(this_gene['start'])
            end_coord = int(this_gene['end'])

        local_gene_list = index.query_genes(anchor_name, chrs, int(start_coord), int(end_coord))
        local_gene_list["name"] = local_gene_list["attr"].str.extract("Name=([^;]+)")

        if get_buffer(start_coord, end_coord, n_skips_start) == 1:
            if (end_coord - start_coord) <= bins:
                start_coord = start_coord - 175
                end_coord = end_coord + 175
            exact_mini_counts = index.query_bitmap(anchor_name, chrs, start_coord, end_coord, 1)
            simple_cnts_for_plots = exact_mini_counts.sum(axis=1)
            n_skips = 1
        else:
            exact_mini_counts = all_chrs[chrs][int(start_coord/n_skips_start):int((end_coord)/n_skips_start)]
            simple_cnts_for_plots = names_simp[int(start_coord/n_skips_start):int((end_coord)/n_skips_start)]
            n_skips = n_skips_start 

        genes = index.query_genes(anchor_name, chrs, 0, index.chrs.loc[anchor_name, chrs]["size"])
        gene_names = [g.split(';')[0].split("=")[1] for g in genes['attr']]
        universals = genes["universal"]
        uniques = genes["unique"]
        sizes = genes["end"]-genes["start"]

        chr_fig = plot_chr_whole(start_coord, end_coord, anchor_name, chrs)
        fig1, bar_sum_names, bar_sum_regional, colors, gene_names_tmp = plot_interactive(
            n_skips, layout, bins, simple_cnts_for_plots, chrs, zs_tmp, click_me_rep, click_me_genes,  
            int(start_coord), int(end_coord), gene_locals, [g.split(';')[0].split('=')[1] for g in genes['attr']],
            anchor_name, chrs
        )
        fig2 = get_local_info(simple_cnts_for_plots, bar_sum_regional, bar_sum_global[anchor_name][chrs], anchor_name, chrs)
        fig3 = create_tree(start_coord, end_coord, exact_mini_counts, n_skips)
        fig4 = plot_gene_content(gene_names,universals,uniques,sizes, 
            sort_by, colors, uniq_avg, univ_avg, int(start_coord), int(end_coord), anchor_name, chrs)   
        return (
            chr_fig, fig1, fig2, fig3, fig4, 
        )

    def update_all_figs(chr_num, click_me_rep, click_me_genes, chrom, anchor_name, redo_wg, start_coord, end_coord, n_skips):
        chrom = index.chrs.loc[anchor_name].index[chr_num-1]#"chr" + str(chr_num)

        genes = index.query_genes(anchor_name, chrom, 0, index.chrs.loc[anchor_name, chrom]["size"])
        bounds = genes.loc[:, ["start", "end"]]
        bounds["break"] = None
        gene_names = [g.split(';')[0].split("=")[1] for g in genes['attr']]

        if get_buffer(start_coord, end_coord, n_skips_start) == 1:
            exact_mini_counts = index.query_bitmap(anchor_name, chrs, start_coord, end_coord, 1)
            simple_cnts_for_plots = exact_mini_counts.sum(axis=1)
            n_skips = 1

        names_simp = index.query_bitmap(anchor_name, chrom, start_coord, end_coord, n_skips).sum(axis=1)

        fig1, bar_sum_names, bar_sum_regional, colors, gene_names_tmp = plot_interactive(
            n_skips, layout, bins, names_simp, 
            chrom, zs_tmp, click_me_rep, click_me_genes,  
            int(start_coord), int(end_coord), 
            bounds.to_numpy().flatten(), 
            gene_names, anchor_name, chrom
        )

        local_gene_list = index.query_genes(anchor_name, chrom, int(start_coord), int(end_coord))
        local_gene_list["name"] = local_gene_list["attr"].str.extract("Name=([^;]+)")

        sort_by = ["Unique","Universal"]
        tmp = [(i/sum(bar_sum_global[anchor_name][chrom])*100) for i in bar_sum_global[anchor_name][chrom]]
        uniq_avg = tmp[1]
        univ_avg = tmp[-1]

        universals = genes["universal"]
        uniques = genes["unique"]
        sizes = genes["end"]-genes["start"]

        fig4 = plot_gene_content(gene_names,universals,uniques,sizes,sort_by, colors, uniq_avg, univ_avg, int(start_coord), int(end_coord), anchor_name, chrom)
        
        chr_fig = plot_chr_whole(start_coord, end_coord, anchor_name, chrom)

        fig2 = get_local_info(names_simp, bar_sum_regional, bar_sum_global[anchor_name][chrom], anchor_name, chrom)
        
        fig3 = create_tree(
            start_coord, end_coord, 
            index.query_bitmap(anchor_name, chrom, start_coord, end_coord, n_skips),
            n_skips)

        if redo_wg == 1:
            big_plot = plot_whole_genome(anchor_name)
            pg_sizes_fig = make_genome_size_plot(anchor_name)
            pg_scaffolds_fig = make_genome_count_plot(anchor_name)
            pg_avgs = plot_avgs(anchor_name)
            tab2_gene_density_fig = make_genes_per_chr_fig(anchor_name)
            tab2_avg_fig = make_avg_kmer_fig(anchor_name)
            tab2_sorted_genes = make_gene_per_genome_fig(anchor_name)
        else:
            big_plot = no_update
            pg_sizes_fig = no_update
            pg_avgs = no_update
            pg_scaffolds_fig = no_update

            tab2_gene_density_fig = no_update
            tab2_avg_fig = no_update
            tab2_sorted_genes = no_update

        return (
            chr_fig, fig1, fig2, fig3, fig4,  
            update_gene_locals(local_gene_list, chrom, start_coord, end_coord, anchor_name)
        )

    def update_output_div(chrs, start, stop, anchor_name):
        print("SET", start, stop)
        return f'{start}-{stop}'

    def update_out_chr(chrs, anchor_name):
        return f'{anchor_name}.{chrs}:'

    def update_gene_locals(local_gene_list, chrs, start_coord, end_coord, anchor_name):
        printme = ""
        if len(local_gene_list)==1:
            printme += "Genes: "
            printme += local_gene_list['name'][0] + ": "
            printme += local_gene_list['attr'][0] #index.query_anno(anchor_name, chrs, start_coord, end_coord)['attr']
        elif len(local_gene_list)<=25:
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

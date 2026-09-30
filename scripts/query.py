import sys
import csv
import pysam
import os.path
from os import path
from mycolorpy import colorlist as mcp
from matplotlib.colors import hex2color
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
from panagram.index import Index
from panagram import figs

max_bins = 10
anchor_name = sys.argv[1]
chrom = sys.argv[2]
start_coord = int(sys.argv[3])
end_coord = int(sys.argv[4])#cluster_f = sys.argv[4]

working_dir = sys.argv[5]#"/home/kjenike1/data-mschatz1/kjenike/PANAGRAM/SOL/INTROGRESS/"
step = 1
#gff_file = working_dir+anchor_name+"_scaffold.gff3"

index  = Index(working_dir)
n_skips = 100
bitmap = index.query_bitmap(anchor_name, chrom, start_coord, end_coord, n_skips)
bin_size = ((end_coord - start_coord) // 200) + 1
pancounts, paircounts = index.bitmap_to_bins(bitmap, bin_size)
genes = []
colors = mcp.gen_color(cmap="viridis_r",n=index.ngenomes)
colors_heatmap = mcp.gen_color(cmap="viridis_r",n=index.ngenomes)


def read_gff(chrom,start_coord,end_coord):
    genes = []
    exons = []
    to_highlight = []
    if highlight == 1:
        highlightme = read_gene_highlights()
    else:
        highlightme = []

    with open(gff_file,"r") as f:
        line = f.readline()
        while line:
            tmp = line.strip().split("\t")
            if tmp[0] == chrom and tmp[2]=="gene" and int(tmp[3]) <= end_coord and int(tmp[4]) >= start_coord:
                if int(tmp[3]) < start_coord:
                    start_here = start_coord
                else:
                    start_here = int(tmp[3])
                if int(tmp[4]) > end_coord:
                    end_here = end_coord
                else:
                    end_here = int(tmp[4])
                genes.append(start_here)
                genes.append(end_here)
                genes.append(None)

                this_id = tmp[8].split("\"")[1]#.split(";")[0]
                if this_id in highlightme:
                    print(this_id)
                    to_highlight.append(start_here)
                    to_highlight.append(end_here)
                    to_highlight.append(None)
            elif tmp[0] == chrom and tmp[2]=="exon" and int(tmp[3]) <= end_coord and int(tmp[4]) >= start_coord:
                if int(tmp[3]) < start_coord:
                    start_here = start_coord
                else:
                    start_here = int(tmp[3])
                if int(tmp[4]) > end_coord:
                    end_here = end_coord
                else:
                    end_here = int(tmp[4])
                exons.append(start_here)
                exons.append(end_here)
                exons.append(None)
            line = f.readline()
    #print(genes)
    #print(exons)
    return genes, exons, to_highlight

def read_clusters(cluster_f):
    cluster_xs = []
    cluster_ys = []
    cluster_zs = []
    with open(cluster_f,"r") as f:
        line = f.readline()
        while line:
            tmp = line.strip().split("\t")
            if tmp[0] == chrom:
                cluster_xs.append(int(tmp[1]))
                cluster_zs.append(int(tmp[3]))
                cluster_ys.append(1)
            line = f.readline()
    return cluster_xs,cluster_ys,cluster_zs
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

def get_clade_lines(orientation='horizontal', y_curr=0, start_coord=0, x_curr=0, y_bot=0, y_top=0,
                        line_color='rgb(25,25,25)', line_width=0.5):
        #Adapted from:
        #https://github.com/plotly/dash-phylogeny
        """define a shape of type 'line', for branch
        """
        branch_line = dict(type='line',
                           layer='above',
                           xref="x",yref="y7",
                           line=dict(color=line_color,
                                     width=line_width)
                           )
        if orientation == 'horizontal':
            branch_line.update(x0=start_coord,
                               y0=y_curr,
                               x1=x_curr,
                               y1=y_curr,
                            )
        elif orientation == 'vertical':
            branch_line.update(x0=x_curr,
                               y0=y_bot,
                               x1=x_curr,
                               y1=y_top,
                )
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
            line_color = color_code[tmp_name]#palette[int(((kmer_num[tmp_name])/total_kmers100)-1]
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
def create_tree(bitmap):
        #Adapted from:
        #https://github.com/plotly/dash-phylogeny

        #We need to sum up the number of times that kmers occur in each column (aka, each sample)
        kmer_num_tmp = bitmap.sum(axis=0)

        matrix = linkage(
            bitmap.transpose(),
            method='ward',
            metric='euclidean'
        )
        tree_tmp2 = hierarchy.to_tree(matrix, False)
        treedata = get_newick(tree_tmp2, tree_tmp2.dist, index.genome_names)

        palette = sns.color_palette("RdPu", 130).as_hex()
        total_kmers = max(kmer_num_tmp) #[-1] #kmer_num_tmp["Solqui2"]
        kmer_num = {}
        color_code = {}
        kmer_num_raw = {}
        #cntr = 0
        for k in range(0, len(kmer_num_tmp)):#kmer_num_tmp.keys():
            kmer_num[index.genome_names[k]] = float(((kmer_num_tmp.iloc[k])/total_kmers)*100)
            kmer_num_raw[index.genome_names[k]] = kmer_num_tmp.iloc[k]
            color_code[index.genome_names[k]] = palette[int(((kmer_num_tmp.iloc[k])/total_kmers)*100)+10]
            #cntr += 1
        #treedata = "(A, (B, C), (D, E))"
        handle = StringIO(treedata)
        tree = Phylo.read(handle, "newick")
        #tree = Phylo.read(tree_file, "newick")
        x_coords = get_x_coordinates(tree)
        y_coords = get_y_coordinates(tree)
        #print(y_coords)
        line_shapes = []
        draw_clade(color_code, total_kmers, kmer_num, palette, tree.root, 0, line_shapes, line_color= "blue",#'rgb(25,25,25)',
                    line_width=6, x_coords=x_coords,
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
                sizes.append(0)
            else:
                #this_color = "grey"
                this_color = color_code[tmp_name]
                color.append(this_color)
                sizes.append(0)
        #print(y_coords)
        #print(max(y_coords))
        #print(min(y_coords))
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
                    x=X, y=Y,
                    xaxis="x",yaxis="y7",
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
            return dict(xref='x', yref='y7', x=x, y=y, text="\t" + tmp_txt,
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
                      xaxis=dict(showline=False,#True,
                                 zeroline=False,
                                 showgrid=False,#True,  # To visualize the vertical lines
                                 ticklen=4,
                                 showticklabels=False,#True,
                                 title='Branch Length',
                                 autorange=False,
                                 #range=[0, 0.1]
                                 ),
                      yaxis=axis,
                      #yaxis_title=None,
                      hovermode='closest',
                      shapes=line_shapes,
                      plot_bgcolor='rgb(250,250,250)',
                      legend={'x': 0, 'y': 1},
                      annotations=annotations,
                      xaxis_range=[0,int(max_x*1.2)]
                  )
        return line_shapes,annotations,node,label_legend,kmer_num

def read_sample_order():
    order = []
    with open("sample_order.txt","r") as f:
        line = f.readline()
        while line:
            order.append(line.strip())
            line = f.readline()
    return order

def update_biodigs(df_tmp,var,tree_order):
    
    ordered = {"Accession_rename":[],"var_order":[]}
    for t in range(0,len(tree_order)):
        ordered["Accession_rename"].append(tree_order[t])
        ordered["var_order"].append(t+1)
    df_order = pd.DataFrame(ordered)

    df = pd.merge(
            df_tmp,
            df_order[["Accession_rename", "var_order"]],
            on="Accession_rename",
            how="left"
            )
    n = 66
    df[var+"_order_reverse"] = n - df["var_order"] + 1
    df_sorted = df.sort_values(by=var+"_order_reverse", ascending=True)
    return df_sorted

def read_clade_order(tree_order):
    #Read in the aves "aves.txt"
    aves = []
    with open("aves_v2.txt","r") as f:
        line = f.readline()
        while line:
            aves.append(line.strip())
            line = f.readline()
    prim = []
    #Read in the primates "primates.txt"
    with open("primates.txt","r") as f:
        line = f.readline()
        while line:
            prim.append(line.strip())
            line = f.readline()
    #Everything else is just "mammal"
    cluster_x, cluster_y, cluster_z = [], [], []
    cntr = 0
    for t in tree_order:
        if t in aves:
            cluster_z.append("A")
        elif t in prim:
            cluster_z.append("P")
        else:
            cluster_z.append("M")
        cluster_x.append(1)
        cluster_y.append(cntr)
        cntr += 1

    return cluster_x, cluster_y, cluster_z

def SetColor(z):
        if(z == "P"):
            return "#ff006e"
        elif(z == "A"):
            return "#3a86ff"
        elif(z == "M"):
            return "#ffbe0b"

def plot_interactive(anchor_name,chrom,start_coord,end_coord,step,pancounts,paircounts,genes,exons,to_highlight,bitmap):
        #genes["break"] = None
        #gene_bounds = genes.loc[:, ["start", "end", "break"]].to_numpy().flatten()
        #gene_names = genes["name"]

        tmp_lst = []
        fig = make_subplots(        
            rows=2, cols=1,
            shared_xaxes="columns",
            vertical_spacing=0.01,
            row_heights=[1,4],
            #column_widths=[2,1,5],
            #specs = [[{"rowspan":3},{},{}],
            #    [None,{},{}],
            #    [None,{},{}]]
            subplot_titles=("Pan-k-mer multiplicity","Pan-k-mer containment" )
        )

        b = bitmap#.sample(n=min(len(bitmap),50000))
        tree_shapes,tree_anno,tree_nodes,tree_order,perc_shared = create_tree(b)
        paircounts = paircounts.loc[tree_order]
        
        bin_size = int((end_coord - start_coord) / max_bins) + 1
        adjusted_bin_size = bin_size
        x = pancounts.columns * bin_size
        cntr = 0 
        hasexon = 1
        linewidth = 1 if hasexon else 15
        
        #fig.update_traces(marker_line_width = 0,
        #          selector=dict(type="bar"),row=3,col=2)
        
        #fig.add_trace(tree_nodes,row=1,col=1)
        #fig.update_layout(annotations=tree_anno,shapes=tree_shapes)
        
        y_vars, x_vars = [], []
        #This is the conserved kmer plotting section
        fig.append_trace(go.Bar(x=paircounts.columns, y=pancounts.loc[0], name=str(0),
                #legendgroup="group1", 
                #legendgrouptitle_text="Conserved K-mers",
                showlegend=False,
                marker=dict(color='grey'), 
                marker_line=dict(color='grey')
            ), row=1, col=1 )

        for i in pancounts.index[1:]:
            #print(i)
            fig.append_trace(go.Bar(x=paircounts.columns,#x=x, 
                y=pancounts.loc[i], name=str(i),
                legendgroup="group1", 
                legendgrouptitle_text="Conserved K-mers",
                marker=dict(
                    color=colors[i-1],
                    line_color=colors[i-1],
                    #colorscale="viridis",
                    #coloraxis="coloraxis",colorbar_title="Pan-Count",
                ),
                #marker_line=dict(color=colors[i-1]),
                showlegend=False,
            ), row=1, col=1 )
        #print(len(paircounts.columns))
        fig.update_layout(barmode='stack', bargap=0.0)
        fig.update_xaxes(showticklabels=False, row=2, col=1)
        sns.set_palette("viridis")
        fig.add_trace(go.Heatmap(
            z=paircounts,x=paircounts.columns, 
            y=paircounts.index, colorscale='plasma_r',#colorscale=[[0, 'rgb(0,0,255)'], [1,'rgb(255,0,0)']],
            #coloraxis="coloraxis2",showlegend=False#,colorscale=[[0, 'rgb(0,0,255)'], [1,'rgb(255,0,0)']]
        ), 
            row=2, col=1 )
        fig.update_coloraxes(showscale=False)
        fig.update_layout(font=dict(size=20))

        #Now we add the reference sequence:
        ticks = np.linspace(start_coord, end_coord+1, 10).round().astype(int)
        yvals = np.ones(len(ticks))
        

        #fig.update_yaxes(visible=False, range=[-1,2], row=1, col=3)
        title_text = "Sequence position ("+chrom+":"+str(start_coord)+"-"+str(end_coord)+")"
        fig.update_xaxes(title_text=title_text, row=2, col=1)

        #TODO don't use template, manually set background to white
        pan_cax = {
            "colorscale" : "viridis",
            "colorbar" : {"title":"Pangenome Conservation","y":0.72,"len":0.35,"yanchor":"top","title_side":"right"}
        }
        #fig.update_coloraxes(showscale=False)
        fig.update_layout(
                width =2000,
                height=2000)

        image_name=anchor_name+"."+chrom+"."+str(start_coord)+"."+str(end_coord)
        #fig.write_html(image_name+".html")
        fig.write_image(working_dir + "/" + image_name+".svg")
        fig.write_image(working_dir + "/" + image_name+".png")


genes, exons, to_highlight = [], [], [] #read_gff(chrom,start_coord,end_coord)
plot_interactive(anchor_name, chrom, start_coord, end_coord, step, pancounts, paircounts, genes, exons, to_highlight, bitmap)



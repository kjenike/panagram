from matplotlib.colors import hex2color
from plotly.subplots import make_subplots
from mycolorpy import colorlist as mcp
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

def genome_colors(index, cmap="viridis_r"):
    return mcp.gen_color(cmap=cmap,n=index.ngenomes)

def read_genome_comp(index, anchor_name):
    totals = index[anchor_name].bitfreq_chrs #.chr_occ_freq.loc[anchor_name,"total"]

    fig = make_subplots(rows=len(totals), cols=1)

    for i,(c,freqs) in enumerate(totals.iterrows()):
        perc = freqs*100
        fig.add_trace(go.Bar(x=perc.index, y=perc, marker_color=genome_colors(index), showlegend=False,), row=i+1,col=1)
    fig.update_yaxes(type="log")
    fig.update_layout(paper_bgcolor='rgba(0,0,0,0)')
    return fig

def make_all_genome_dend(index):
    dist_mat = np.zeros((index.ngenomes, index.ngenomes), np.float64)

    with open(index.genome_dist_fname) as f:
        for line in f:
            f, t, d, p, x = line.rstrip().split("\t")
            print(f,t)
            i = index.genomes[f].id
            j = index.genomes[t].id
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

def read_pangenome_comp(index):
    fig = make_subplots(rows=1, cols=1)
    colors = genome_colors(index)
    for i,g in enumerate(index.genome_names):
        fig.add_trace(go.Bar(y=index.genome_names, x=index.bitfreq_totals[i+1], name=i+1, #orientation='h',
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
    return fig

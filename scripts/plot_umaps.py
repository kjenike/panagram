import sys
import plotly.express as px
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
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import umap
import sklearn.cluster as cluster
from sklearn.cluster import DBSCAN
from collections import Counter

anchor_name = sys.argv[1]  # human
chrom_in = sys.argv[2]  # ALL
start_coord = 0  # int(sys.argv[3])#0
max_bins = int(sys.argv[3])  # 100000
eps = float(sys.argv[4])  # 1
step = 100

index = Index("/home/kjenike1/data-mschatz1/kjenike/PANAGRAM/DOGS/K51_WILD_DOGS/")

n_neighbors = int(sys.argv[5])  # 4
md = float(sys.argv[6])  # 0.0 is good
n_skips = int(sys.argv[7])  # 100
first_chrom = "chr1_RagTag"
chr_ider = "chr"
colors = mcp.gen_color(cmap="viridis_r", n=index.ngenomes)
colors_heatmap = mcp.gen_color(cmap="viridis_r", n=index.ngenomes)


def one_loc(chrom):
    end_coord = index[anchor_name].chrs.loc[chrom, "size"]
    bitmap = index.query_bitmap(anchor_name, chrom, start_coord, end_coord, n_skips)

    bin_size = max_bins  # ((end_coord - start_coord) // max_bins) + 1

    pancounts, paircounts = index.bitmap_to_bins(bitmap, bin_size)
    chrom_list = []
    starts = []
    for i in range(0, len(paircounts.T)):
        chrom_list.append(chrom)
        starts.append((i * bin_size))
    return paircounts, chrom_list, starts, bin_size


if chrom_in != "ALL":
    bins = {}
    paircounts, chrom_list, starts, bi = one_loc(chrom_in)

    bins[chrom_in] = bi
else:
    bins = {}
    chroms = list(index[anchor_name].chrs.index)
    paircounts, chrom_list, starts, bi = one_loc(first_chrom)
    paircounts = paircounts.fillna(0)

    bins[first_chrom] = bi
    for i in range(1, len(chroms)):
        c = chroms[i]

        if c.count(chr_ider) > 0:  # and int(c.split("-")[1]) > 5000000:

            paircounts_1, chrom_list_1, starts_1, bins_1 = one_loc(c)
            paircounts_1 = paircounts_1.fillna(0)

            paircounts = np.concatenate((paircounts_1, paircounts), axis=1)
            chrom_list = chrom_list + chrom_list_1
            starts = starts + starts_1
            bins[c] = bins_1
    # paircounts_2, chrom_list_2, starts_2 = one_loc("Chr2")
    # paircounts_3, chrom_list_3, starts_3 = one_loc("Chr3")
    # paircounts_4, chrom_list_4, starts_4 = one_loc("Chr4")
    # paircounts_5, chrom_list_5, starts_5 = one_loc("Chr5")
    # paircounts = np.concatenate((paircounts_1, paircounts_2, paircounts_3, paircounts_4, paircounts_5), axis=1)
    # chrom_list = chrom_list_1 + chrom_list_2 + chrom_list_3 + chrom_list_4 + chrom_list_5
    # starts = starts_1 + starts_2 + starts_3 + starts_4 + starts_5


id_colors = {"Eurasian": 1, "Non_Iberian_relict": 2, "Iberian_relict": 3, "Iberian_non_relict": 4}
# ids = []


def get_ids_from():
    with open("test.samples", "r") as f:
        line = f.readline()
        while line:
            tmp = line.strip().split("\t")
            ids.append(id_colors[tmp[0]])
            line = f.readline()
    return ids


def get_ids_from_coords():
    for i in range(0, max_bins):
        ids.append(i)
    return ids


# ids = get_ids_from_coords()


def make_pca(clusters):
    # sns.set_palette("hls", 4)
    pca = PCA(n_components=3)
    principalComponents = pca.fit_transform(paircounts.T)

    principalDf = pd.DataFrame(
        data=principalComponents,
        columns=["principal component 1", "principal component 2", "principal component 3"],
    )

    fig, ax = plt.subplots()

    scatter = ax.scatter(
        principalDf["principal component 1"],
        principalDf["principal component 2"],
        cmap="rainbow",
        c=clusters,
        alpha=0.5,
    )
    ax.set_xlabel("PC 1: " + str(pca.explained_variance_ratio_[0]))
    ax.set_ylabel("PC 2: " + str(pca.explained_variance_ratio_[1]))
    legend1 = ax.legend(*scatter.legend_elements(), loc="lower right", title="Groups")
    ax.add_artist(legend1)
    plt.legend(loc="upper right")
    plt.savefig(
        "PCAS/pca." + anchor_name + "." + chrom_in + ".bins" + str(max_bins) + ".genome.pos.png"
    )
    plt.close()


# fig_2d = px.scatter(
#    proj_2d, x=0, y=1,
#    color=df.species, labels={'color': 'species'}
# )


# Now make a umap
# n_neighbors = int(sys.argv[6])#10
# md = float(sys.argv[7])
reducer = umap.UMAP(n_neighbors=n_neighbors, min_dist=md, n_components=2, random_state=42)


# paircounts = paircounts.fillna(0)
embedding = reducer.fit_transform(paircounts.T)

# adjusted_ids = []
# end_coord = index[anchor_name].chrs.loc[chrom,"size"]
bin_size = max_bins  # int(end_coord/max_bins)
# for i in ids:
#    adjusted_ids.append(i*bin_size)


# clusters = cluster.KMeans(n_clusters=10).fit_predict(embedding)
clusters = DBSCAN(eps=eps, min_samples=1).fit_predict(embedding)


make_pca(clusters)

fig_2d = px.scatter(
    embedding,
    x=0,
    y=1,
    color=clusters,  # adjusted_ids,
    # color_discrete_sequence=px.colors.sequential.Plasma_r,
    # color_discrete_sequence=px.colors.qualitative.Spectral,
    color_continuous_scale="Spectral",  #'rainbow',
    opacity=0.7,
    # symbol = chrom_list,
    # symbol_sequence= ['circle-dot','square','diamond','cross','triangle-up','circle-dot']#'circle-open', 'circle', 'circle-open-dot', 'square', 'circle-dot']
    # auto_open=False,
    # filename='plot_result.html'#labels={'color': 'species'}
)
# fig_2d.update_traces('marker':{'color':'blue'}}, selector={'name': 'won'})
# fig_2d.update_coloraxes(showscale=False)
# fig_2d.update_traces(marker={'size': 15})
filename = (
    "UMAPS/umap.cluster."
    + anchor_name
    + "."
    + chrom_in
    + ".bins"
    + str(max_bins)
    + ".neighbors"
    + str(n_neighbors)
    + ".md"
    + str(md)
    + ".nskips."
    + str(n_skips)
    + ".genome.pos.maxbinstest.eps"
    + str(eps)
)
# fig_2d.write_html("umap.cluster."+anchor_name+"."+chrom_in+".bins"+str(max_bins)+".neighbors"+str(n_neighbors)+".md"+str(md)+".nskips."+str(n_skips)+".genome.pos.maxbinstest.eps"+str(eps)+".html")
fig_2d.write_html(filename + ".html")
fig_2d.write_image(filename + ".svg")
# fig_2d.show()
# plt.scatter(embedding[:, 0], embedding[:, 1], cmap='rainbow',c=ids,alpha=0.8)
# plt.colorbar()
# plt.savefig("umap."+anchor_name+"."+chrom+".bins"+str(max_bins)+".neighbors"+str(n_neighbors)+".genome.pos.png")
# plt.close()

# Print a bed like file


def print_stuff():
    for i in range(0, len(clusters)):
        print(
            chrom_list[i]
            + "\t"
            + str(starts[i])
            + "\t"
            + str(starts[i] + bins[chrom_list[i]] - 1)
            + "\t"
            + str(clusters[i])
            + "\t"
            + str(embedding[i, 0])
            + "\t"
            + str(embedding[i, 1])
        )


print_stuff()


def plot_interactive(
    anchor_name, chrom, start_coord, end_coord, step, pancounts, paircounts, genes
):
    # genes["break"] = None
    # gene_bounds = genes.loc[:, ["start", "end", "break"]].to_numpy().flatten()
    # gene_names = genes["name"]

    tmp_lst = []
    fig = make_subplots(
        rows=3,
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.01,
        row_heights=[1, 8, 8],
        # subplot_titles=("Ref. Sequence Position","", "",  "Conserved K-mers","" )
    )

    # We are adjusting the start and stop positions to account for the skipping.
    # The adjusted value should be the index, whereas the start_coord and end_coord are the real coordinates

    bin_size = max_bins  # int((end_coord - start_coord) / max_bins) + 1
    adjusted_bin_size = int(bin_size / step)

    x = pancounts.columns * bin_size

    cntr = 0

    # anno_types = index.genomes[anchor_name].gff_anno_types
    # hasexon = "exon" in anno_types
    hasexon = 0
    # c = np.arange(len(anno_types) + (not hasexon))
    # ann_colors = np.array(px.colors.qualitative.Prism)
    # ann_colors = ann_colors[c % len(ann_colors)]

    linewidth = 1 if hasexon else 15

    # fig.append_trace(go.Scattergl(
    #    x=gene_bounds,
    #    y=np.full(len(gene_bounds),0),
    #    line=dict(color=ann_colors[0], width=linewidth),
    #    showlegend=False,
    #    text=np.repeat(gene_names, 3),
    #    hovertemplate='<br>x:%{x}<br>m:%{text}', legendgroup="group2",
    #    name="gene"), row=2, col=1)

    # This is the conserved kmer plotting section
    fig.append_trace(
        go.Bar(
            x=x,
            y=pancounts.loc[0],
            name=str(0),
            # legendgroup="group1",
            # legendgrouptitle_text="Conserved K-mers",
            showlegend=False,
            marker=dict(color="grey"),
            marker_line=dict(color="grey"),
        ),
        row=2,
        col=1,
    )

    # t1 = time.perf_counter()

    # t0 = t1

    fig.add_trace(
        go.Scattergl(
            x=[1, 2],
            y=[1, 2],
            marker=dict(
                color=[1, index.ngenomes],
                coloraxis="coloraxis",
                colorscale="viridis",
                colorbar_title="Pan-Count",
            ),
            showlegend=False,
            opacity=0,
        ),
        row=2,
        col=1,
    )

    for i in pancounts.index[1:]:
        fig.append_trace(
            go.Bar(
                x=x,
                y=pancounts.loc[i],
                name=str(i),
                legendgroup="group1",
                legendgrouptitle_text="Conserved K-mers",
                marker=dict(
                    color=colors[i - 1],
                    line_color=colors[i - 1],
                    # colorscale="viridis",
                    # coloraxis="coloraxis",colorbar_title="Pan-Count",
                ),
                # marker_line=dict(color=colors[i-1]),
                showlegend=False,
            ),
            row=2,
            col=1,
        )

    fig.update_layout(barmode="stack", bargap=0.0)
    fig.update_xaxes(showticklabels=False, row=3, col=1)
    sns.set_palette("viridis")
    fig.add_trace(
        go.Heatmap(
            z=paircounts,
            x=paircounts.columns,
            y=paircounts.index,
            colorscale="plasma_r",  # colorscale=[[0, 'rgb(0,0,255)'], [1,'rgb(255,0,0)']],
            # coloraxis="coloraxis2",showlegend=False#,colorscale=[[0, 'rgb(0,0,255)'], [1,'rgb(255,0,0)']]
        ),
        row=3,
        col=1,
    )
    fig.update_coloraxes(showscale=False)
    # t1 = time.perf_counter()

    # t0 = t1

    # Now we add the reference sequence:
    ticks = np.linspace(start_coord, end_coord + 1, 10).round().astype(int)
    yvals = np.ones(len(ticks))
    fig.append_trace(
        go.Scattergl(
            x=ticks,
            y=yvals,
            text=ticks.astype(str),
            textposition="top center",
            showlegend=False,
            mode="lines+markers+text",
            line=dict(color="grey"),
            marker=dict(size=5, symbol="line-ns"),
        ),
        row=1,
        col=1,
    )
    # fig.update_coloraxes(showscale=False)
    # t1 = time.perf_counter()

    # t0 = t1

    fig.update_yaxes(visible=False, range=[0.9, 4], row=1, col=1)
    fig.update_xaxes(visible=False, title_text="Sequence position", row=1, col=1)
    fig.update_xaxes(title_text="Sequence position", row=3, col=1)
    # fig.update_yaxes(title_text="# of k-mers", range=[0,bin_size]  , row=3, col=1)
    fig.update_yaxes(title_text="# of k-mers", range=[0, adjusted_bin_size], row=2, col=1)

    # TODO don't use template, manually set background to white
    pan_cax = {
        "colorscale": "viridis",
        "colorbar": {
            "title": "Pangenome Conservation",
            "y": 0.72,
            "len": 0.35,
            "yanchor": "top",
            "title_side": "right",
        },
    }
    # fig.update_coloraxes(showscale=False)
    fig.update_layout(width=2000, height=800)
    fig.write_image("PLOTS/" + anchor_name + "/" + anchor_name + "." + chrom + ".png")


# anchor_name = "Col_0"
# chrom = "Chr1"
# start_coord = 0
# end_coord = 32640075
# step = 100

# index  = Index("/scratch2/katie/PANAGRAM/ARABIDOPSIS_CENT/K51_w_te_ano/")
# n_skips = 100
# bitmap = index.query_bitmap(anchor_name, chrom, start_coord, end_coord, n_skips)
# bin_size = ((end_coord - start_coord) // 200) + 1
# pancounts, paircounts = index.bitmap_to_bins(bitmap, bin_size)
# genes = []
# plot_interactive(anchor_name, chrom, start_coord, end_coord, step, pancounts, paircounts, genes)

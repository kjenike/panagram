import sys
import os.path
from mycolorpy import colorlist as mcp
from plotly.subplots import make_subplots
import seaborn as sns
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import dash
import time
import base64
from PIL import Image
from Bio import Phylo
from dash import dcc, html, Input, Output, ctx, State, no_update
from scipy.cluster import hierarchy
from io import StringIO
from scipy.cluster.hierarchy import linkage
from .index import Index
from . import figs


def view(params):
    index = Index(params.index_dir)  # Directory that contains the anchor direcotry

    anchor_name = params.genome
    if anchor_name is None:
        anchor_name = index.anchor_genomes[0]  # params.genome

    chrs = params.chrom
    if chrs is None:
        chrs = index[anchor_name].chrs.index[0]

    if params.start is None:
        params.start = 0
    if params.end is None:
        params.end = index[anchor_name].chrs.loc[chrs, "size"]

    sns.set_palette("viridis", index.ngenomes)
    colors = mcp.gen_color(cmap="viridis_r", n=index.ngenomes)

    # This will have the figure componants that we need
    layout = go.Layout(
        margin=go.layout.Margin(
            l=10, r=10, b=10, t=10  # left margin  # right margin  # bottom margin  # top margin
        )
    )

    ROOT = os.path.dirname(os.path.realpath(__file__))

    image_path = ROOT + "/assets/panagram.png"
    pil_img = Image.open(image_path)

    def b64_image(image_filename):
        with open(image_filename, "rb") as f:
            image = f.read()
        return "data:image/png;base64," + base64.b64encode(image).decode("utf-8")

    whole_genome_hists_fig = figs.read_genome_comp(index, anchor_name)
    whole_genome_hists_fig.update_layout(
        hoverlabel=dict(bgcolor="white", font_size=35, font_family="Balto"),
        title=dict(text="K-mer conservation per chromosome", x=0.5),
        font=dict(size=25, family="Balto"),
    )
    all_genomes_dend_fig = figs.make_all_genome_dend(index)
    all_genomes_dend_fig.update_layout(
        hoverlabel=dict(bgcolor="white", font_size=35, font_family="Balto"),
        font=dict(size=30, family="Balto"),
    )
    pangenome_comp_fig = figs.read_pangenome_comp(index)
    pangenome_comp_fig.update_layout(
        hoverlabel=dict(bgcolor="white", font_size=35, font_family="Balto")
    )

    chrs = index.chrs.loc[anchor_name].index[0]  # "chr1"

    config = {
        "toImageButtonOptions": {"format": "svg", "width": None, "height": None},
        "scrollZoom": False,
    }
    # tab_style = {"background-color": "lightgrey", "font-size": 36}

    # Annotation tab info
    PANGENOME_TAB = [
        html.Div(
            className="w3-half",
            children=[
                dcc.Graph(id="all_genomes", style={"font-size": 20, "height": 1500}, config=config),
            ],
        ),
        html.Div(
            className="w3-half",
            children=[
                dcc.Graph(
                    id="all_genomes_kmer_comp",
                    config=config,
                    style={"font-size": 30, "height": 1500},  # , figure = pangenome_comp_fig
                )
            ],
        ),
        html.Div(
            className="w3-third",
            children=[
                dcc.Graph(
                    id="all_genome_sizes",
                    config=config,
                    style={"font-size": 30, "height": 500},  # ,figure = pg_size_fig,
                )
            ],
        ),
        html.Div(
            className="w3-third",
            children=[
                dcc.Graph(
                    id="average_kmer_content",
                    config=config,
                    style={"font-size": 30, "height": 500},  # figure = pangenome_avg_fig,
                )
            ],
        ),
        html.Div(
            className="w3-third",
            children=[
                dcc.Graph(
                    id="pangenome_num_seqs",
                    config=config,
                    style={"font-size": 30, "height": 500},  # ,figure = pg_count_fig
                )
            ],
        ),
    ]

    # def render_anchor_tab():
    #    return
    ANCHOR_TAB = [
        # html.Div(children=[
        # html.I("Anchor genome " + anchor_name, id="anchor_labels"),
        # html.Br(),
        # html.Label('Select a chromosome: '),
        # dcc.Dropdown(chrs_list[anchor_name], chrs, style=dict(width='40%', height='110%', verticalAlign="middle", ), id="chr_select_dropdown"),
        #        html.Br(),
        #    ], style={'padding-top' : '1%', 'padding-left' : '1%', 'padding-bottom' : '1%', 'padding-right' : '1%', 'font-size':'24px', 'textAlign': 'left', "border":"2px grey solid", }),
        html.Div(
            className="w3-threequarter",
            children=[
                dcc.Graph(
                    id="genome_umap", config=config, style={"font-size": 30, "height": 900}
                )  # , figure = avg_kmer_per_chr_fig
            ],
        ),
        html.Div(
            className="w3-quarter",
            children=[
                dcc.Graph(
                    id="genome_umap_hist", config=config, style={"font-size": 30, "height": 900}
                )  # , figure = avg_kmer_per_chr_fig
            ],
        ),
        # html.Div(className="w3-quarter",children=[
        #    dcc.Graph(id="gene_content", config=config, style={"font-size": 30, "height" : 900}), #gene_content_fig)
        # ],),
        # html.Div(className="w3-quarter",children=[
        #    dcc.Graph(id="genes_per_chr",
        #    config=config)
        # ],),
        # html.Div(className="w3-quarter",children=[
        #    dcc.Graph(id="avg_kmer_chr", config=config) #, figure = avg_kmer_per_chr_fig
        # ],),
        html.Div(
            className="w3-threequarter",
            children=[
                dcc.Graph(
                    id="all_chromosomes",
                    config=config,
                    style={"height": index[anchor_name].chr_count * 250},
                )
            ],
        ),
        html.Div(
            className="w3-quarter",
            children=[
                dcc.Graph(
                    id="all_chromosomes_hists",
                    style={"height": index[anchor_name].chr_count * 250},
                    figure=whole_genome_hists_fig,
                    config=config,
                )
            ],
        ),
    ]

    # def render_chromosome_tab():
    #    return [
    CHROMOSOME_TAB = [
        html.Div(
            children=[  # summary
                html.Div(
                    className="w3-threequarter",
                    children=[
                        dcc.Graph(id="Umap", config=config, style={"height": 800, "font-size": 30}),
                    ],
                ),
                html.Div(
                    className="w3-quarter",
                    children=[
                        dcc.Graph(
                            id="Umap_genome_hist",
                            config=config,
                            style={"height": 800, "font-size": 30},
                        ),
                    ],
                ),
                html.Div(
                    className="w3-container",
                    children=[
                        # left figure
                        dcc.Graph(
                            id="chromosome", config=config, style={"font-size": 30, "height": 350}
                        )
                    ],
                ),
            ]
        ),
        html.Div(
            children=[
                html.Div(
                    className="w3-container",
                    children=[
                        dcc.Graph(
                            id="primary", config=config, style={"height": 1000, "font-size": 30}
                        )
                    ],
                )
            ]
        ),
        html.Div(
            children=[
                html.Div(
                    className="w3-container",
                    children=[
                        # html.Div(className="w3-half", children=[
                        #    dcc.Graph(id="Umap",
                        #        config=config,
                        #        style={"height": 900, "font-size": 20, "margin-right": "120px"}),
                        # ]),
                        html.Div(
                            className="w3-half",
                            children=[
                                dcc.Graph(
                                    id="Genes",
                                    # figure=gene_content_plot,
                                    config=config,
                                    style={"font-size": 30, "height": 800, "margin-left": "0px"},
                                ),
                            ],
                        ),
                        html.Div(
                            className="w3-half",
                            children=[
                                dcc.Graph(
                                    id="Third",
                                    # This is the histogram section
                                    config=config,
                                    style={"font-size": 30, "height": 800},
                                ),
                            ],
                        ),
                    ],
                )
            ]
        ),
    ]
    """
    ANNOTATION_TAB = html.Div(children=[
        html.Div(className="w3-container", children=[

            html.Div(className="w3-threequarter", children=[
                dcc.Graph(id="annotation_conservation",
                    #figure= annotation_tab_plot(annotation_tab_df, "Length", 1), #plot_anno_conserve(anchor_name),
                    config=config,
                    style={"height": 1250, "font-size": 40}
                    )
                ]),
            html.Div(className="w3-quarter", children=[
                #Control panel for the annotation tab plots 
                html.I("Color by: ", style={ "font-size": 40}),
                html.Br(),
                dcc.RadioItems(
                    options = [
                        {'label':' Log2 gene length', 'value':' Log2 gene length' },
                        {'label':' Gene length', 'value':' Gene length' },
                        {'label':' Chromosome', 'value':' Chromosome' },
                        {'label':' Bed file [FUTURE FEATURE]', 'value':' Bed file [FUTURE FEATURE]', 'disabled':True },
                        {'label':' Selected region [FUTURE FEATURE]', 'value':' Selected region [FUTURE FEATURE]', 'disabled':True }
                    ], value = " Log2 gene length",
                    #[' Log2 gene length', ' Gene length', ' Chromosome', ' Bed file [FUTURE FEATURE]'],
                    #' Log2 gene length', 
                    labelStyle={"display":"block"}, style={'font-size' : 30, }, id="radios", 
                ),
                html.Br(),
                html.I("", style={ "font-size": 30}, id="annotation_tab_clickdata_Name"),
                html.Br(),
                html.I("", style={ "font-size": 30}, id="annotation_tab_clickdata_X"),
                html.Br(),
                html.I("", style={ "font-size": 30}, id="annotation_tab_clickdata_Y"),
                html.Br(),
                html.I("", style={ "font-size": 30}, id="annotation_tab_clickdata_Color"),
                html.Br(),
                html.P("", style={ "font-size": 30, "overflow-wrap": "break-word"}, id="annotation_tab_clickdata_Attr"),
                
                ], style={'padding' : '1%', "height": 1250}),
            ])#"border":"2px grey solid"
        ])
    """
    app = dash.Dash(
        __name__,
        # external_stylesheets=[dbc.themes.CERULEAN],
        external_stylesheets=[
            "https://www.w3schools.com/w3css/4/w3.css"
        ],  # , dbc.themes.BOOTSTRAP],
        url_base_pathname=params.url_base,
        suppress_callback_exceptions=True,
    )

    app.layout = html.Div(
        [
            # html.Div(id = 'parent', children = [
            html.Div(
                style={
                    "display": "flex",  # Arrange children horizontally
                    "alignItems": "center",  # Vertically center items
                    "padding": "1%",
                    "textAlign": "left",
                },
                children=[
                    html.Img(src=pil_img, style={"height": "84px", "marginRight": "20px"}),
                    html.H1(
                        id="H1",
                        children="Panagram",
                        style={
                            "fontSize": 84,
                            "margin": 0,
                            "color": "#440154",
                            "fontWeight": "bold",
                        },
                    ),
                ],
            ),
            html.Div(
                children=[
                    html.Div(
                        className="w3-row-padding",
                        children=[
                            html.Div(
                                [
                                    html.I(f"Panagram of {index.ngenomes} genomes"),
                                    html.Br(),
                                    html.I(f"K-mer length: {index.k}"),
                                    html.Br(),
                                    html.I(f"Step size: {index.lowres_step}", id="step_size_out"),
                                ],
                                className="w3-third",
                            ),
                            html.Div(
                                [
                                    html.Label("Select a genome:"),
                                    dcc.Dropdown(
                                        options=[
                                            {"label": g, "value": g} for g in index.anchor_genomes
                                        ],
                                        value=anchor_name,
                                        id="genome_select_dropdown",
                                        style={"maxWidth": "300px"},
                                        clearable=False,
                                    ),
                                    html.Label("Select a chromosome:", style={"marginTop": "1rem"}),
                                    dcc.Dropdown(
                                        options=[
                                            {"label": c, "value": c}
                                            for c in index[anchor_name].chrs.index
                                        ],
                                        value="",
                                        id="chr_select_dropdown",
                                        style={"maxWidth": "300px"},
                                        clearable=False,
                                    ),
                                ],
                                className="w3-third",
                            ),
                            html.Div(
                                [
                                    html.I(f"{anchor_name}.{chrs}:", id="chr_name"),
                                    dcc.Input(
                                        id="Chrs_Info",
                                        debounce=True,
                                        style={"maxWidth": "300px", "marginBottom": "0.5rem"},
                                    ),
                                    html.I(id="regional_genes"),
                                    html.Br(),
                                    html.Label(
                                        "Pre-selected regions:", style={"marginTop": "1rem"}
                                    ),
                                    dcc.Dropdown(
                                        id="pre_regions_dropdown",
                                        style={"maxWidth": "300px"},
                                        clearable=False,
                                    ),
                                ],
                                className="w3-third",
                            ),
                        ],
                        style={"marginTop": "1rem", "marginBottom": "1rem"},
                    )
                ],
                style={
                    "fontSize": "30px",  # Bigger font size
                    "padding": "20px",  # Padding inside the box
                    "border": "1px solid #ccc",  # Light grey border
                    "borderRadius": "10px",  # Rounded corners
                    "boxShadow": "2px 2px 12px rgba(0,0,0,0.1)",  # Subtle shadow
                    "width": "98%",  # Optional max width for neatness
                    "boxSizing": "border-box",
                    "backgroundColor": "#eaf4fb",
                    "margin": "20px auto",  # Center horizontally with margin top/bottom
                },
            ),
            # html.H1(id = 'H1', children = 'Panagram', style = {'textAlign':'left', "font-size": 64, 'padding' : '1%'}),
            # html.Div(children=[
            #    html.Div(children = [
            #        html.I("Panagram of " + str(index.ngenomes) + " genomes"),
            #        html.Br(),
            #        html.I("K-mer length: " + str(index.k),style={"display": "inline-block",}),
            #        html.Br(),
            #        html.I("Step size: " + str(index.lowres_step), style={"display": "inline-block",}, id='step_size_out'),
            #        ], style={"display": "inline-block", 'padding-left': '1%'}),
            #    html.Div(children = [
            #        html.Label('Select a genome: '),
            #        dcc.Dropdown(index.anchor_genomes, anchor_name, style=dict(width='110%', height='100%', verticalAlign="middle"), id="genome_select_dropdown"),
            #        #html.Br(),#style=dict(width='100%', height='110%', verticalAlign="middle", )
            #        html.Label('Select a chromosome: '),
            #        dcc.Dropdown(index[anchor_name].chrs.index, "", style=dict(width='110%', height='100%', verticalAlign="middle", ), id="chr_select_dropdown"),
            #    ], style={"display": "inline-block", 'padding-left' : '10%', 'vertical-align': 'text-bottom'}),
            #    html.Div(children = [
            #        html.I(anchor_name + "." + chrs + ":", id="chr_name"),
            #        dcc.Input(id="Chrs_Info", debounce=True), #, placeholder=str(x_start_init) + "-" + str(x_stop_init)+ "          "
            #        html.Br(),
            #        html.I(id='regional_genes'), #"Genes in this region: " + str(local_genes),
            #        html.Br(),
            #        html.Label('Pre-selected regions: '),
            #        dcc.Dropdown(style=dict(width='100%', height='110%', verticalAlign="middle", ), id="pre_regions_dropdown"), #pre_bed_info[anchor_name], chrs + ":" + str(x_start_init) + "-" + str(x_stop_init) ,
            #
            #    ], style={"display": "inline-block", 'padding-left' : '10%', 'vertical-align': 'text-bottom'}),
            # ], style = {'padding-top' : '1%', 'padding-left' : '1%', 'padding-bottom' : '1%', 'padding-right' : '1%', 'font-size':'30px', "border":"2px grey solid"}),
            # html.Div(className="loader-wrapper", children=[
            dcc.Loading(
                id="loading-tabs",
                parent_className="loading_wrapper",
                type="default",
                className="dash-spinner",
                children=dcc.Tabs(
                    id="tabs",
                    value="pangenome",
                    style={"font-size": 36},
                    children=[
                        dcc.Tab(
                            value="pangenome",
                            label="Pangenome",  # style=  tab_style,
                            children=PANGENOME_TAB,
                            className="custom-tab",
                            selected_className="custom-tab--selected",
                        ),
                        dcc.Tab(
                            value="anchor",
                            label="Anchor genome",  # style= tab_style,
                            children=ANCHOR_TAB,
                            className="custom-tab",
                            selected_className="custom-tab--selected",
                        ),
                        dcc.Tab(
                            value="chromosome",
                            label="Chromosome",  # style=tab_style,
                            children=CHROMOSOME_TAB,
                            className="custom-tab",
                            selected_className="custom-tab--selected",
                        ),
                        # dcc.Tab(value="annotation", label='Annotation', style=tab_style, children=ANNOTATION_TAB),
                    ],
                    className="custom-tabs",
                ),  # style={'backgroundColor': 'transparent'}
            ),
            # ], style={'backgroundColor': 'transparent'}), #style={"maxHeight": "300vh",}),
            html.Div(id="tab-content"),
            html.Div(chrs, style={"display": "none"}, id="selected-chrom-state"),
            html.Div(anchor_name, style={"display": "none"}, id="selected-anchor-state"),
            html.Div(params.start, id="start-coord-state", style={"display": "none"}),
            html.Div(params.end, id="end-coord-state", style={"display": "none"}),
            html.Div(id="chr-genes-state", style={"display": "none"}),
            # html.Div([],id='gene-names-state',style={"display" : "none"} ),
        ]
    )  # ,])

    def annotation_tab_info(file_name):
        all_avgs = []
        all_stds = []
        all_names = []
        all_lens = []
        all_chrs = []
        all_attr = []
        all_starts = []
        all_stops = []
        with open(file_name, "r") as f:
            line = f.readline()
            while line:
                tmp = line.strip().split("\t")
                this_avg = float(tmp[4])
                all_avgs.append(this_avg)
                all_stds.append(float(tmp[5]))
                all_names.append(tmp[0])
                all_lens.append(int(tmp[3]) - int(tmp[2]))
                all_chrs.append(tmp[1])
                all_attr.append(tmp[6])
                all_starts.append(int(tmp[2]))
                all_stops.append(int(tmp[3]))
                line = f.readline()

        d = {
            "Standard_dev": all_stds,  # np.log10(all_stds),#all_stds,
            "Average": all_avgs,
            "Name": all_names,
            "Length": all_lens,
            "chr": all_chrs,
            "attr": all_attr,
            "start": all_starts,
            "end": all_stops,
        }
        df = pd.DataFrame(d)
        return df

    def annotation_tab_plot(df, color_by, log2_true):
        if log2_true == 0:
            color_me = df[color_by]
        else:
            color_me = np.log2(df[color_by])
        if color_by == "chr":
            colors = px.colors.sample_colorscale(
                "plasma",
                [
                    n / (index[anchor_name].chr_count - 1)
                    for n in range(index[anchor_name].chr_count)
                ],
            )
            fig = px.scatter(
                df,
                x="Average",
                y="Standard_dev",
                marginal_y="histogram",
                marginal_x="histogram",
                hover_name="Name",
                opacity=0.5,
                color_discrete_sequence=colors,  # px.colors.sequential.Plasma,
                color=color_me,
            )
        else:
            fig = px.scatter(
                df,
                x="Average",
                y="Standard_dev",
                marginal_y="histogram",
                marginal_x="histogram",
                hover_name="Name",
                opacity=0.5,
                color=color_me,
            )
        fig.update_layout(clickmode="event+select")
        return fig

    def get_newick(node, parent_dist, leaf_names, newick="") -> str:
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
        # Adapted from:
        # https://github.com/plotly/dash-phylogeny
        """Associates to  each clade an x-coord.
        returns dict {clade: x-coord}
        """
        xcoords = tree.depths()
        if not max(xcoords.values()):
            xcoords = tree.depths(unit_branch_lengths=True)
        return xcoords

    def get_y_coordinates(tree, dist=1.3):
        # Adapted from:
        # https://github.com/plotly/dash-phylogeny
        """
        returns  dict {clade: y-coord}
        The y-coordinates are  (float) multiple of integers (i*dist below)
        dist depends on the number of tree leafs
        """
        maxheight = tree.count_terminals()  # Counts the number of tree leafs.
        # Rows are defined by the tips/leafs
        ycoords = dict(
            (leaf, maxheight - i * dist) for i, leaf in enumerate(reversed(tree.get_terminals()))
        )

        def calc_row(clade):
            for subclade in clade:
                if subclade not in ycoords:
                    calc_row(subclade)
            ycoords[clade] = (ycoords[clade.clades[0]] + ycoords[clade.clades[-1]]) / 2

        if tree.root.clades:
            calc_row(tree.root)
        return ycoords

    def get_clade_lines(
        orientation="horizontal",
        y_curr=0,
        start_coord=0,
        x_curr=0,
        y_bot=0,
        y_top=0,
        line_color="rgb(25,25,25)",
        line_width=0.5,
    ):
        # Adapted from:
        # https://github.com/plotly/dash-phylogeny
        """define a shape of type 'line', for branch"""
        branch_line = dict(
            type="line",
            layer="above",
            xref="x",
            yref="y7",
            line=dict(color=line_color, width=line_width),
        )
        if orientation == "horizontal":
            branch_line.update(
                x0=start_coord,
                y0=y_curr,
                x1=x_curr,
                y1=y_curr,
            )
        elif orientation == "vertical":
            branch_line.update(
                x0=x_curr,
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

    def draw_clade(
        color_code,
        total_kmers,
        kmer_num,
        palette,
        clade,
        start_coord,
        line_shapes,
        line_color="grey",
        line_width=5,
        x_coords=0,
        y_coords=0,
    ):
        # Adapted from:
        # https://github.com/plotly/dash-phylogeny
        """Recursively draw the tree branches, down from the given clade"""
        x_curr = x_coords[clade]
        y_curr = y_coords[clade]
        if str(clade) != "Clade" and str(clade) != "AT":
            tmp_name = str(clade)
            line_color = color_code[
                tmp_name
            ]  # palette[int(((kmer_num[tmp_name])/total_kmers)*100)-1]
        elif clade.clades:

            line_color = palette[int(biggest_num_in_clade(clade, kmer_num)) + 10]
            # line_color = "grey"
            # Now we have to find the all of the children, and which one has the highest value
        # Draw a horizontal line from start to here
        branch_line = get_clade_lines(
            orientation="horizontal",
            y_curr=y_curr,
            start_coord=start_coord,
            x_curr=x_curr,
            line_color=line_color,
            line_width=line_width,
        )
        line_shapes.append(branch_line)
        if clade.clades:
            # Draw a vertical line connecting all children
            y_top = y_coords[clade.clades[0]]
            y_bot = y_coords[clade.clades[-1]]
            # Now we have to get the right color though. It should be the max value of any children
            line_shapes.append(
                get_clade_lines(
                    orientation="vertical",
                    x_curr=x_curr,
                    y_bot=y_bot,
                    y_top=y_top,
                    line_color=line_color,
                    line_width=line_width,
                )
            )
            # Draw descendants
            for child in clade:
                draw_clade(
                    color_code,
                    total_kmers,
                    kmer_num,
                    palette,
                    child,
                    x_curr,
                    line_shapes,
                    x_coords=x_coords,
                    y_coords=y_coords,
                )

    def create_tree(bitmap):
        # Adapted from:
        # https://github.com/plotly/dash-phylogeny

        # We need to sum up the number of times that kmers occur in each column (aka, each sample)
        kmer_num_tmp = bitmap.sum(axis=0)

        # Convert to numpy array if it's a pandas Series
        if hasattr(kmer_num_tmp, "values"):
            kmer_num_tmp = kmer_num_tmp.values

        matrix = linkage(bitmap.transpose(), method="ward", metric="euclidean")
        tree_tmp2 = hierarchy.to_tree(matrix, False)
        treedata = get_newick(tree_tmp2, tree_tmp2.dist, index.genome_names)

        palette = sns.color_palette("RdPu", 130).as_hex()
        total_kmers = max(kmer_num_tmp)  # [-1] #kmer_num_tmp["Solqui2"]
        kmer_num = {}
        color_code = {}
        kmer_num_raw = {}
        # cntr = 0
        for k in range(0, len(kmer_num_tmp)):  # kmer_num_tmp.keys():
            kmer_num[index.genome_names[k]] = float(((kmer_num_tmp[k]) / total_kmers) * 100)
            kmer_num_raw[index.genome_names[k]] = int(kmer_num_tmp[k])
            color_code[index.genome_names[k]] = palette[
                int(((kmer_num_tmp[k]) / total_kmers) * 100) + 10
            ]
            # cntr += 1
        # treedata = "(A, (B, C), (D, E))"
        handle = StringIO(treedata)
        tree = Phylo.read(handle, "newick")
        # tree = Phylo.read(tree_file, "newick")
        x_coords = get_x_coordinates(tree)
        y_coords = get_y_coordinates(tree)
        line_shapes = []
        draw_clade(
            color_code,
            total_kmers,
            kmer_num,
            palette,
            tree.root,
            0,
            line_shapes,
            line_color="blue",  #'rgb(25,25,25)',
            line_width=6,
            x_coords=x_coords,
            y_coords=y_coords,
        )
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
            # if tmp_name.count("_") > 0:
            #    tmp_name = tmp_name.split("_")[1]
            if str(cl.name) == "None":
                this_color = "grey"
                color.append(this_color)
                sizes.append(0)
            else:
                # this_color = "grey"
                this_color = color_code[tmp_name]
                color.append(this_color)
                sizes.append(0)
        # graph_title = "Pansol Phylgenetic Tree"#create_title(virus_name, nb_genome)
        intermediate_node_color = "rgb(100,100,100)"
        axis = dict(
            showline=False,
            zeroline=False,
            showgrid=False,
            showticklabels=False,
            title="",  # y title
        )
        label_legend = []
        an_x = []
        an_y = []
        for t in range(0, len(text)):
            if text[t] != None:
                label_legend.append(text[t])
                an_x.append(X[t])
                an_y.append(Y[t])
        nodes = []
        node = dict(
            type="scatter",
            x=X,
            y=Y,
            xaxis="x",
            yaxis="y7",
            mode="markers",
            marker=dict(color=color, size=sizes),
            text=text,  # vignet information of each node
            # hoverinfo='',
            showlegend=False,
        )
        nodes.append(node)

        def make_anns(x, y, text, kmer_num, i):
            tmp_txt = text  # ""
            tmp_txt += (
                " - " + str(kmer_num[tmp_txt])[:4] + "% "
            )  # (" + str(kmer_num_raw[tmp_txt]) + ")"
            return dict(
                xref="x",
                yref="y7",
                x=x,
                y=y,
                text="\t" + tmp_txt,
                showarrow=False,
                xanchor="left",
                yanchor="middle",
            )

        annotations = []
        for i in range(0, len(label_legend)):
            annotations.append(make_anns(an_x[i], an_y[i], label_legend[i], kmer_num, i))
        max_x = max(an_x)
        layout = dict(  # title=graph_title,
            paper_bgcolor="rgba(0,0,0,0)",
            # dragmode="select",
            # font=dict(family='Balto', size=22),
            # width=1000,
            # height=1500,
            autosize=True,
            showlegend=True,
            xaxis=dict(
                showline=False,  # True,
                zeroline=False,
                showgrid=False,  # True,  # To visualize the vertical lines
                ticklen=4,
                showticklabels=False,  # True,
                title="Branch Length",
                autorange=False,
                # range=[0, 0.1]
            ),
            yaxis=axis,
            # yaxis_title=None,
            hovermode="closest",
            shapes=line_shapes,
            plot_bgcolor="rgb(250,250,250)",
            legend={"x": 0, "y": 1},
            annotations=annotations,
            xaxis_range=[0, int(max_x * 1.2)],
        )
        # fig = make_subplots(
        #    rows=1, cols=1,
        #    specs=[[{"type": "scatter", 'colspan':1}]], #, None, {"type": "bar"} ]], #, {"type": "pie"} , {"type": "pie"} ]],
        #    horizontal_spacing=0.01,
        #    #subplot_titles=("This region","CDS","Exons","Genes", "Whole chromosome")
        # )
        # fig.add_trace(node, row=1, col=1)
        # fig.update_layout(layout)
        ##Sort the hist bars
        # hist_y = []#list(kmer_num_tmp.values()).sort()
        # hist_x = []

        # fig.update_yaxes(visible=False, showticklabels=False)
        # fig.update_layout(margin=dict(
        #        t=20,
        #        b=10,
        #        l=10,
        #        r=10))
        return line_shapes, annotations, node, label_legend, kmer_num  # fig,label_legend

    def get_umap_chrom_hist(anchor_name, chrom):
        umap = index.genomes[anchor_name].chrom_umaps.loc[chrom]
        color = list(umap["cluster"])
        cluster_size = umap["end"] - umap["start"]
        # Convert to numpy array if it's a pandas Series
        if hasattr(cluster_size, "values"):
            cluster_size = cluster_size.values
        x, y, z = [], [], []

        for i in range(0, max(color) + 1):
            x.append(i)
            z.append(i)
            y.append(color.count(i) * cluster_size[0])

        hist_fig = go.Figure(
            [
                go.Bar(
                    x=y,
                    y=x,
                    marker_color=z,
                    orientation="h",
                    marker={
                        "colorscale": "Portland",
                    },
                )
            ]
        )
        hist_fig.update_xaxes(title_text="Base pairs per cluster")
        hist_fig.update_xaxes(type="log")
        hist_fig.update_layout(
            plot_bgcolor="rgba(0,0,0,0)",
            font=dict(size=25),
            hoverlabel=dict(bgcolor="white", font_size=35, font_family="Balto"),
        )
        return hist_fig

    def get_umap_from_file(chrom, anchor):
        umap = index.genomes[anchor_name].chrom_umaps.loc[chrom]

        X = umap["umap1"].values if hasattr(umap["umap1"], "values") else umap["umap1"]
        Y = umap["umap2"].values if hasattr(umap["umap2"], "values") else umap["umap2"]
        color = umap["cluster"].values if hasattr(umap["cluster"], "values") else umap["cluster"]
        chrom_start = umap["start"].values if hasattr(umap["start"], "values") else umap["start"]
        chrom_stop = umap["end"].values if hasattr(umap["end"], "values") else umap["end"]

        hover_text = [
            f"Cluster: {color[i]}<br>Chromosome: {chrom}:{chrom_start[i]}-{chrom_stop[i]}"
            for i in range(len(X))
        ]  # <br>X: {X[i]}<br>Y: {Y[i]}" for i in range(len(X))]
        fig = go.Figure(
            data=go.Scatter(
                x=X,
                y=Y,
                mode="markers",
                # ids=color, #hoverinfo='text',
                hovertext=hover_text,
                marker=dict(
                    color=color,
                    # colorscale="Portland",
                    coloraxis="coloraxis",
                    showscale=True,
                    size=16,
                    opacity=0.7,
                ),
                # hovertemplate="Cluster:%{z}<br>X: %{x}<br>Y: %{y}",#'<br>x:%{x}<br>y:%{y}',
                # font=dict(size=20),
                showlegend=False,
            )
        )
        pan_cax = {
            "colorscale": "Portland",
            "colorbar": {
                "title": "DBSCAN Cluster",
                "y": 0,
                "len": 1.05,
                "yanchor": "bottom",
                "title_side": "right",
            },
        }
        # fig.update_layout(font=dict(size=20))
        fig.update_layout(
            title={
                "text": "Pangenome conservation UMAP (10kbp regions)",
                "xanchor": "center",
                "x": 0.5,
                "y": 0.85,
                "yanchor": "top",
            },
            coloraxis=pan_cax,
            plot_bgcolor="rgba(0,0,0,0)",
            # xaxis_title=sort_by[0] + " k-mer content rank",
            # yaxis_title="% of k-mers",
            # margin=dict(
            #    t=20,
            #    b=10,
            #    l=10,
            #    r=10),
            hoverlabel=dict(bgcolor="white", font_size=35, font_family="Balto"),
            font=dict(
                size=25,
            ),
        )
        fig.add_vline(x=0, line_color="black")
        fig.add_hline(y=0, line_color="black")
        fig.update_xaxes(showticklabels=False)
        fig.update_yaxes(showticklabels=False)
        # fig = go.Figure(data=[go.scatter(node)])#make_subplots()
        return fig

    def read_clusters(anchor, chrom):
        umap = index.genomes[anchor_name].chrom_umaps.loc[chrom]

        data = {"x": [], "y": [], "z": []}
        data["x"] = list(umap["start"]) + list(umap["start"])
        data["y"] = []
        for i in range(0, len(umap["start"])):
            data["y"].append(1)
        for i in range(0, len(umap["start"])):
            data["y"].append(2)
        data["z"] = list(umap["end"]) + list(umap["end"])
        return data

    def get_local_info_just_one(bar_sum_regional, anchor_name, chrs):
        x = []
        for i in range(1, index.ngenomes + 1):
            x.append(i)
        # y_whole = list(index.bitfreq_chrs.loc[anchor_name,chrs].loc[1:])
        y = [(i / sum(bar_sum_regional[1:]) * 100) for i in bar_sum_regional[1:]]
        data = go.Bar(x=x, y=y, marker_color=colors, showlegend=False)
        # fig.append_trace(go.Bar(x=x, y=y_whole, marker_color=colors, showlegend=False), row=1, col=2)
        return data

    def get_local_info(bar_sum_regional, anchor_name, chrs):
        fig = make_subplots(
            rows=1,
            cols=3,
            # specs=[[{"type": "bar", "colspan": 2}, None],
            #   [{"type": "bar"}, {"type": "bar"}]],
            subplot_titles=(
                "Whole chromosome",
                "This region",
                "Genes",
            ),
            vertical_spacing=0.1,
        )
        x = []
        for i in range(1, index.ngenomes + 1):
            x.append(i)
        # x = [1,2,3,4,5,6,7,8,9]

        # colors = ["#fde725", "#addc30", "#5ec962", "#28ae80", "#21918c", "#2c728e", "#3b528b", "#472d7b", "#440154"]
        # This region
        y = [(i / sum(bar_sum_regional[1:]) * 100) for i in bar_sum_regional[1:]]
        # y_whole=[(i/sum(bar_sum_global_tmp)*100) for i in bar_sum_global_tmp]
        y_whole = list(index.bitfreq_chrs.loc[anchor_name, chrs].loc[1:])
        # .sum(axis=1)

        # fig.add_trace(go.Bar(x=x, y=[a_i - b_i for a_i, b_i in zip(y, y_whole)], marker_color=colors, showlegend=False), row=2, col=1)
        fig.append_trace(go.Bar(x=x, y=y, marker_color=colors, showlegend=False), row=1, col=1)
        fig.append_trace(
            go.Bar(x=x, y=y_whole, marker_color=colors, showlegend=False), row=1, col=2
        )
        # Genes
        # y=[(i/sum(gene_comp[1:])*100) for i in gene_comp[1:]]
        # fig.add_trace(go.Bar(x=x, y=y_whole, marker_color="#7400b8", showlegend=False), row=1, col=4)
        totals = 0
        gene_comp = index[anchor_name].bitfreq_genes.loc[chrs] * 100

        fig.append_trace(
            go.Bar(
                x=x,
                y=[a_i - b_i for a_i, b_i in zip(gene_comp, y_whole)],
                marker_color=colors,
                showlegend=False,
            ),
            row=1,
            col=3,
        )
        # fig.update_layout(xaxis_title_text="K-mers shared in X samples", yaxis_title_text='Frequency (log)')
        fig.update_xaxes(title_text="# of genomes", row=1, col=2)
        fig.update_yaxes(title_text="Difference from whole chromosome", row=1, col=2)
        fig.update_yaxes(title_text="Difference from whole chromosome", row=1, col=3)
        fig.update_yaxes(title_text="Percent of k-mers", row=1, col=1)
        fig.update_yaxes(type="log", row=1, col=1)
        fig.update_yaxes(type="log", row=1, col=2)
        fig.update_yaxes(type="log", row=1, col=3)
        # fig.update_layout(height=1000)
        fig.update_annotations(font_size=25)
        fig.update_layout(
            plot_bgcolor="rgba(0,0,0,0)",
            hoverlabel=dict(bgcolor="white", font_size=35, font_family="Balto"),
            font=dict(size=25),
        )
        return fig

    def plot_interactive(
        bar_sum_regional,
        anchor_name,
        chrom,
        start_coord,
        end_coord,
        step,
        bitmap,
        pancounts,
        paircounts,
        genes,
    ):
        t0 = time.perf_counter()

        genes["break"] = None
        gene_bounds = genes.loc[:, ["start", "end", "break"]].to_numpy().flatten()
        gene_names = genes["name"]
        # bounds.to_numpy().flatten(),

        tmp_lst = []
        fig = make_subplots(
            rows=4,
            cols=2,
            shared_xaxes=True,  # "columns",
            # shared_yaxes=[False, False, False, True],
            shared_yaxes="rows",
            vertical_spacing=0.01,
            horizontal_spacing=0.01,
            row_heights=[2, 3, 8, 8],
            column_widths=[2, 7],
            specs=[[{}, {}], [{"rowspan": 2}, {}], [{}, {}], [{}, {}]],
            # shared_yaxes=True,
            subplot_titles=("", "", "Pan-k-mer summary in this region", "", "", "", "", ""),
        )
        coord_row = 1
        anno_row = 2
        pan_row = 3
        pair_row = 4

        # fig.update_yaxes(title_text="Difference from whole chromosome", row=1, col=2)
        # Now we are adding a histogram to the first two rows of column 1. This will require a rowspan
        hist_data = get_local_info_just_one(bar_sum_regional, anchor_name, chrom)
        fig.add_trace(hist_data, row=anno_row, col=1)
        fig.update_yaxes(type="log", row=anno_row, col=1, title="Frequency (log)")
        fig.update_xaxes(domain=[0, 0.195], row=anno_row, col=1)  # ,title="Number of samples")

        # fig.update_traces(row=anno_row, col=1,
        #        xaxis=dict(domain=[0, 0.9]))
        # fig.update_layout(
        #    yaxis=dict(
        #        type="log",
        #        domain=[0, 1],       # full height
        #        #anchor="x",
        #    ),
        #    xaxis=dict(
        #        domain=[0, 0.9],     # restrict to 90% of width
        #        #anchor="y"
        #    ),#row=anno_row,col=1
        ##)
        # start_coord = bitmap_counts

        # We are adjusting the start and stop positions to account for the skipping.
        # The adjusted value should be the index, whereas the start_coord and end_coord are the real coordinates

        # start_coord = pancnts.index[0]
        # end_coord = pancnts.index[-1] + pancnts.index.step

        b = bitmap.sample(n=min(len(bitmap), 50000))
        tree_shapes, tree_anno, tree_nodes, tree_order, perc_shared = create_tree(b)

        if params.order is not None:
            order = params.order
        else:
            order = tree_order
        paircounts = paircounts.loc[order]

        bin_size = int((end_coord - start_coord) / params.max_chr_bins) + 1
        adjusted_bin_size = int(bin_size / step)

        x = pancounts.columns * bin_size

        t1 = time.perf_counter()
        print(f"\tPivot bins {t1 - t0:0.4f} seconds")
        t0 = t1

        cntr = 0

        t1 = time.perf_counter()
        print(f"\tAnno Query {t1 - t0:0.4f} seconds")
        t0 = t1

        anno_types = index.genomes[anchor_name].gff_anno_types
        print(anno_types)
        if anno_types:
            hasexon = "exon" in anno_types

            c = np.arange(len(anno_types) + (not hasexon))
            ann_colors = np.array(px.colors.qualitative.Prism)
            ann_colors = ann_colors[c % len(ann_colors)]

            linewidth = 1 if hasexon else 15
            print("WIDHT", linewidth)

            fig.add_trace(
                go.Scattergl(
                    x=gene_bounds,
                    xaxis="x2",  # yaxis="y2",
                    y=np.full(len(gene_bounds), 0),
                    line=dict(color=ann_colors[0], width=linewidth),
                    showlegend=False,
                    text=np.repeat(gene_names, 3),
                    hovertemplate="<br>x:%{x}<br>m:%{text}",
                    legendgroup="group2",
                    name="gene",
                ),
                row=anno_row,
                col=2,
            )
            fig.update_layout(clickmode="event+select")

            t1 = time.perf_counter()
            print(f"\tGene Plot {t1 - t0:0.4f} seconds")
            t0 = t1

            # ann_colors = mcp.gen_color(cmap="turbo",n=len(index.genomes[anchor_name].gff_anno_types)-1)
            # ann_colors = np.array(["#3b528b"]+ann_colors)

            anno_names = ["gene"]

            anno = index.query_anno(
                anchor_name, chrom, start_coord, end_coord
            )  # .set_index("type_id").sort_index()
            anno["break"] = np.nan
            grp = anno.groupby("type_id")
            for t, df in anno.groupby("type_id"):
                # for i,t in enumerate(anno_types):
                #    df = anno.loc[t]
                xs = df[["start", "end", "break"]].to_numpy().flatten()
                ys = np.full(len(xs), -t * 2)
                name = df["type"].iloc[0]

                if name == "exon":
                    linewidth = 15
                else:
                    anno_names.append(name)
                    linewidth = 1

                fig.add_trace(
                    go.Scattergl(
                        x=xs,
                        y=ys,
                        xaxis="x2",  # yaxis="y2",
                        line=dict(width=linewidth, color=ann_colors[t]),
                        name=name,
                        # hoverinfo='name',
                        hovertemplate="<br>x:%{x}<br>m:%{text}",
                        text=np.repeat(df["name"], 3),
                        showlegend=False,
                        opacity=0.7,
                        marker={
                            "symbol": "line-ns",
                            "line_color": ann_colors[t],
                            "line_width": 1,
                            "size": 2,
                        },
                        mode="lines+markers",
                    ),
                    row=anno_row,
                    col=2,
                )

            # ys = np.arange(-len(ann_colors),0)+1
            # names = anno_types[::-1]+["gene"] #anno_names[::-1]
            # if len(ys) > 10:
            #    n = len(ys)//10
            #    ys = ys[::n]
            #    names = names[::n]
            # print(ys)
            # print(names)

            ticks = index.genomes[anchor_name].anno_type_ids.reset_index()
            ticks.columns = ["name", "idx"]
            ticks.loc[0, "name"] = "gene"
            ticks["y"] = -ticks["idx"]

            if len(ticks) > 5:
                n = len(ticks) // 5
                ticks = ticks.iloc[::n]

            fig.update_yaxes(
                ticktext=ticks["name"],
                tickvals=ticks["y"],
                range=[-len(ann_colors) - 0.5, 1.5],  # title="Annotation",
                showticklabels=True,
                row=anno_row,
                col=2,
            )

        t1 = time.perf_counter()
        print(f"\tAnno Plot {t1 - t0:0.4f} seconds")
        t0 = t1

        # This is the conserved kmer plotting section
        fig.add_trace(
            go.Bar(
                x=x,
                y=pancounts.loc[0],
                name=str(0),
                xaxis="x2",  # yaxis="y4",
                # legendgroup="group1",
                # legendgrouptitle_text="Conserved K-mers",
                showlegend=False,
                marker=dict(color="grey"),
                marker_line=dict(color="grey"),
            ),
            row=pan_row,
            col=2,
        )

        t1 = time.perf_counter()
        print(f"\tConserved k-mers (grey) {t1 - t0:0.4f} seconds")
        t0 = t1

        fig.add_trace(
            go.Scattergl(
                x=[1, 2],
                y=[1, 2],
                xaxis="x2",  # yaxis="y4",
                marker=dict(
                    color=[1, index.ngenomes],
                    coloraxis="coloraxis",
                    colorscale="viridis",
                    colorbar_title="Pan-Count",
                ),
                showlegend=False,
                opacity=0,
            ),
            row=pan_row,
            col=2,
        )

        for i in pancounts.index[1:]:
            fig.add_trace(
                go.Bar(
                    x=x,
                    y=pancounts.loc[i],
                    name=str(i),
                    xaxis="x2",  # yaxis="y4",
                    legendgroup="group1",
                    legendgrouptitle_text="Conserved K-mers",
                    marker=dict(
                        color=colors[i - 1],
                        line_color=colors[i - 1],
                    ),
                    showlegend=False,
                ),
                row=pan_row,
                col=2,
            )

        fig.update_layout(barmode="stack", bargap=0.0)
        fig.update_xaxes(showticklabels=False, row=pan_row, col=2)

        fig.add_trace(
            go.Heatmap(
                z=paircounts,
                x=paircounts.columns,
                y=paircounts.index,
                xaxis="x2",  # yaxis="y6",
                coloraxis="coloraxis2",
            ),
            row=pair_row,
            col=2,
        )
        # fig.update_yaxes(anchor="y6",row=4,col=1)
        # fig.update_shapes(yref="y7",row=4,col=1)

        print(tree_nodes)
        print(tree_shapes)
        # Now we add the clustering
        fig.add_trace(tree_nodes, row=pair_row, col=1)  # ,sharey = ax1)
        fig.update_layout(annotations=tree_anno, shapes=tree_shapes)
        fig.update_xaxes(showticklabels=False, row=pair_row, col=1)
        fig.update_yaxes(showticklabels=False, row=pair_row, col=1)
        # fig.add_shapes(tree_shapes,row=4,col=1)

        t1 = time.perf_counter()
        print(f"\tConserved kmers, non-grey {t1 - t0:0.4f} seconds")
        t0 = t1

        # This is where we add the UMAP
        # umap = get_umap_from_file(chrom)
        # fig.add_trace(umap, row=2,col=1)

        # Now we add the reference sequence:
        ticks = np.linspace(start_coord, end_coord + 1, 10).round().astype(int)
        yvals = np.ones(len(ticks)) * 3
        fig.add_trace(
            go.Scattergl(
                x=ticks,
                y=yvals,
                text=ticks.astype(str),
                xaxis="x2",  # yaxis="y2",
                textposition="top center",
                showlegend=False,
                mode="lines+markers+text",
                line=dict(color="grey"),
                marker=dict(size=5, symbol="line-ns"),
            ),
            row=coord_row,
            col=2,
        )

        t1 = time.perf_counter()
        print(f"\tFinishing touches {t1 - t0:0.4f} seconds")
        t0 = t1

        fig.update_yaxes(visible=False, range=[0, 5], row=coord_row, col=2)

        fig.update_xaxes(visible=False, title_text="Sequence position", row=coord_row, col=2)
        fig.update_xaxes(title_text="Sequence position", row=pair_row, col=2)
        clusters = read_clusters(anchor_name, chrom)
        fig.add_trace(
            go.Heatmap(clusters, colorscale="Portland", showscale=False), row=coord_row, col=2
        )
        # fig.update_traces(showscale=False)
        # fig.update_yaxes(title_text="# of k-mers", range=[0,bin_size]  , row=3, col=1)
        fig.update_yaxes(
            showticklabels=True,
            title_text="# of k-mers",
            range=[0, adjusted_bin_size],
            row=pan_row,
            col=2,
        )

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

        pair_cax = {
            "colorscale": "plasma_r",
            "colorbar": {
                "title": "Pairwise Conservation",
                "y": 0,
                "len": 0.35,
                "yanchor": "bottom",
                "title_side": "right",
            },
        }
        # clust_cax = {
        #    "colorscale":"Portland",
        #    "colorbar": {
        #        "title":"Clusters",
        #        "y":0,
        #        "len":0.35,
        #        "yanchor":"bottom",
        #        "title_side":"right"}
        # }
        fig.update_xaxes(range=[start_coord, end_coord], row=coord_row, col=2)
        fig.update_annotations(font_size=25)
        fig.update_layout(
            font=dict(size=20),
            plot_bgcolor="rgba(0,0,0,0)",
            coloraxis=pan_cax,
            coloraxis2=pair_cax,
            hoverlabel=dict(bgcolor="white", font_size=35, font_family="Balto"),
            #    coloraxis3=clust_cax
        )

        fig.update_layout(
            {
                "yaxis": {"matches": None},
                "yaxis2": {"matches": None},
                "yaxis3": {"matches": None},
                "yaxis4": {"matches": None},
                "yaxis5": {"matches": None},
                "yaxis6": {"matches": None},
                "yaxis7": {"matches": None},
                "yaxis8": {"matches": "y6"},
                #            "yaxis9": {"matches": "y6"},
            }
        )
        fig.update_layout(
            {
                "xaxis": {"matches": None},
                # "xaxis2": {"matches": None},
                # "xaxis3": {"matches": None},
                "xaxis2": {"matches": None},
                "xaxis3": {"matches": None},
                "xaxis4": {"matches": "x2"},
                "xaxis5": {"matches": "x2"},
                "xaxis6": {"matches": "x2"},
                "xaxis7": {"matches": None},
                "xaxis8": {"matches": "x2"},
            }
        )
        # fig.update_traces(row=4, col=2, secondary_y=False ) #selector=dict(type='bar'))
        # fig.layout['yaxis7'] = fig.layout['yaxis6']
        # for trace in fig['data']:
        #    if trace.name == 'Bottom Right Plot':
        #        trace.yaxis = 'y6'
        # fig['data'][-1]['yaxis'] = 'y6'
        # for trace in fig['data']:
        #    if trace.yaxis == 'y7':  # y7 = bottom-right
        #        trace.yaxis = 'y6'   # y6 = bottom-left
        # for trace in fig['data']:
        #    if trace.yaxis == 'y2':  # This is the bottom-right plot (based on your output)
        #        trace.yaxis = 'y6'
        # fig.layout.pop('yaxis2', None)
        # fig.update_layout(
        #        yaxis2=dict(
        #            matches='y6',          # Force y2 to match y6
        #            overlaying='y6',       # Put it on the same physical axis
        #            showticklabels=False   # Avoid duplicate ticks
        #            ))
        # for i, trace in enumerate(fig.data):
        #    #print(i)
        #    #print(trace.yaxis)
        #    if i == 16:
        #        trace.yaxis = 'y8'
        #    elif i == 17:
        #        trace.yaxis = 'y8'
        #    elif i == 18:
        #        trace.yaxis = 'y2'
        #    elif i == 19:
        #        trace.yaxis = 'y2'
        #    else:
        #        trace.yaxis = 'y6'
        for i, trace in enumerate(fig.data):
            print(
                f"Trace {i}: {trace.name if hasattr(trace, 'name') else 'No name'}, "
                f"x-axis: {trace.xaxis}, y-axis: {trace.yaxis}"
            )
        """
            for attr in dir(fig.layout):
            #if attr.startswith("yaxis") and isinstance(getattr(fig.layout, attr), go.layout.YAxis):
            if not attr.startswith("yaxis"):
                continue

            try:
                axis = getattr(fig.layout, attr)
            except AttributeError:
                continue
            if not isinstance(axis, go.layout.YAxis):
                continue
            
            if attr.startswith("yaxis"):
                #yaxis = getattr(fig.layout, attr)
                axis = getattr(fig.layout, attr)
                if not isinstance(axis, go.layout.YAxis):
                    continue
                axis_name = attr  # e.g., 'yaxis6'
                axis_id = axis_name[5:] or '1'  # 'yaxis' → y1, 'yaxis6' → y6
                axis_ref = f"y{axis_id}"
                print(f"\n🔹 {axis_name} (trace ref: {axis_ref})")
                print(f"Type: {getattr(axis, 'type', 'not specified')}")

            # Case 1: Numerical axis with defined range
            if getattr(axis, 'type', None) in (None, 'linear', 'log'):
                print(f"Range: {getattr(axis, 'range', 'None')}")
        
                # Case 2: Categorical axis
            else: #axis.type == 'category':
                categories = set()

                # Try from layout first
                if getattr(axis, 'categoryarray', None):
                    categories.update(yaxis.categoryarray)
                elif getattr(axis, 'ticktext', None):
                    categories.update(axis.ticktext)

                # Then collect from traces using this axis
                for trace in fig.data:
                    if trace.yaxis == axis_ref:
                        if isinstance(trace.y, (list, tuple)):
                            categories.update(trace.y)

                print("Categories:", sorted(categories))
        """
        # for i in range(2,8):
        #    yaxis_key = "yaxis"+str(i)
        #    yaxis = getattr(fig.layout, yaxis_key, None)#fig.layout.get("yaxis"+str(i))
        #    if yaxis and yaxis.range:
        #        print(f"{i} range: {yaxis.range}")
        #    elif trace.yaxis == 'y6' and isinstance(trace.y, (list, tuple)):
        #    else:
        #        print(f"{i} has no range explicitly set.")
        # for trace in fig.data:
        #    if trace.yaxis in ('y2', 'y7'):  # these are usually assigned to row=4,col=2
        #        trace.yaxis = 'y6'
        # for i, trace in enumerate(fig.data):
        #    print(f"Trace {i}: {trace.name}, uses y-axis: {trace.yaxis}")
        print(fig.layout["yaxis6"].domain)

        # for i, trace in enumerate(fig['data']):
        #    print(f"Trace {i}: {trace.name}, uses y-axis: {trace.yaxis}")
        # for k in fig.layout:
        #    if k.startswith("yaxis"):
        #        print(f"{k} → {fig.layout[k]}")
        t1 = time.perf_counter()
        print(f"\tTruly finishing touches {t1 - t0:0.4f} seconds")
        sys.stdout.flush()

        return fig, perc_shared

    def make_gene_whole_chr_old(x, locs):
        z_genes = [0] * (len(x) + 2)
        g = 0
        while g < len(locs):
            x_bin = int(int(locs[g]) / index.max_bin_len)
            if x_bin < len(z_genes):
                z_genes[x_bin] += 1
            g += 1
        return z_genes

    def make_gene_whole_chr(name, x, locs):
        z_genes = np.zeros(len(x) + 2, int)
        for i in x.searchsorted(locs, "right"):
            z_genes[i] += 1
        return z_genes

    def plot_chr_whole(start_coord, end_coord, anchor_name, this_chr, genes):
        df = index[anchor_name].bitsum_bins.loc[this_chr]
        div = df.sum(axis=1)

        z_1 = df[1] / div
        z_9 = df[index.ngenomes] / div
        # z_1 = index.bitsum_bins.loc[(anchor_name,this_chr), 1]
        # z_9 = index.bitsum_bins.loc[(anchor_name,this_chr), index.ngenomes]
        # y, x = [], []
        # cntr = 0
        # for xtmp in z_1:
        #    x.append(cntr)
        #    y.append(1)
        #    cntr += window_size
        x = df.index
        y = np.full(len(x), 1)

        z_genes = [0] * len(x)
        # if index.gene_tabix[anchor_name] is not None:
        if len(genes) > 0:
            bounds = genes["start"].to_numpy()  # genes.loc[:, ["start"]].to_numpy()
            z_genes = make_gene_whole_chr(anchor_name, x, bounds)

        # Umap clusters
        umap_data = read_clusters(anchor_name, this_chr)
        chr_fig = make_subplots(
            rows=4,
            cols=1,
            specs=[
                [
                    {
                        "type": "heatmap",
                    }
                ],
                [
                    {
                        "type": "heatmap",
                    }
                ],
                [
                    {
                        "type": "heatmap",
                    }
                ],
                [
                    {
                        "type": "heatmap",
                    }
                ],
            ],
            shared_xaxes=True,
            subplot_titles=(
                "Pan-k-mer clusters, k-mer density and gene density accross whole chromosome",
                "",
                "",
            ),
            vertical_spacing=0.0,
            # font=dict(
            #    size=20,
            # )
            # height=350,
        )
        chr_fig.update_annotations(font_size=30)
        chr_fig.append_trace(
            go.Heatmap(umap_data, colorscale="Portland", showscale=False, showlegend=False),
            row=1,
            col=1,
        )
        chr_fig.append_trace(
            go.Scatter(
                x=[
                    start_coord,
                    start_coord,
                    None,
                    end_coord,
                    end_coord,
                    None,
                    start_coord,
                    end_coord,
                ],
                showlegend=False,
                y=[0.5, 1.5, None, 0.5, 1.5, None, 1.45, 1.45],
                mode="lines",
                line_color="#1dd3b0",
                line_width=8,
            ),
            row=1,
            col=1,
        )

        chr_fig.append_trace(
            go.Heatmap(
                x=x,
                z=z_9,
                y=y,
                type="heatmap",
                colorscale="magma_r",
                showlegend=False,
                showscale=False,
            ),
            row=2,
            col=1,
        )
        chr_fig.append_trace(
            go.Scatter(
                x=[
                    start_coord,
                    start_coord,
                    None,
                    end_coord,
                    end_coord,
                    None,
                    start_coord,
                    end_coord,
                ],
                showlegend=False,
                y=[0.5, 1.5, None, 0.5, 1.5],  # , None, 1.45, 1.45 ],
                mode="lines",
                line_color="#1dd3b0",
                line_width=8,
            ),
            row=2,
            col=1,
        )

        chr_fig.append_trace(
            go.Heatmap(x=x, z=z_1, y=y, type="heatmap", colorscale="magma", showscale=False),
            row=3,
            col=1,
        )
        chr_fig.append_trace(
            go.Scatter(
                x=[start_coord, start_coord, None, end_coord, end_coord],
                showlegend=False,
                y=[0.5, 1.5, None, 0.5, 1.5],
                mode="lines",
                line_color="#1dd3b0",
                line_width=8,
            ),
            row=3,
            col=1,
        )

        chr_fig.append_trace(
            go.Heatmap(x=x, z=z_genes, y=y, type="heatmap", colorscale="magma_r", showscale=False),
            row=4,
            col=1,
        )
        chr_fig.append_trace(
            go.Scatter(
                x=[
                    start_coord,
                    start_coord,
                    None,
                    end_coord,
                    end_coord,
                    None,
                    start_coord,
                    end_coord,
                ],
                showlegend=False,
                y=[0.5, 1.5, None, 0.5, 1.5, None, 0.55, 0.55],
                mode="lines",
                line_color="#1dd3b0",
                line_width=8,
            ),
            row=4,
            col=1,
        )

        chr_fig.update_xaxes(
            fixedrange=True, range=[0, index.chrs.loc[anchor_name, this_chr]["size"]], row=1, col=1
        )
        chr_fig.update_yaxes(title_text="C.", range=[0.5, 1.5], showticklabels=False, row=1, col=1)

        chr_fig.update_xaxes(
            fixedrange=True, range=[0, index.chrs.loc[anchor_name, this_chr]["size"]], row=2, col=1
        )
        chr_fig.update_yaxes(title_text="U.", range=[0.5, 1.5], showticklabels=False, row=2, col=1)

        chr_fig.update_xaxes(
            fixedrange=True, range=[0, index.chrs.loc[anchor_name, this_chr]["size"]], row=3, col=1
        )
        chr_fig.update_yaxes(title_text="P.", range=[0.5, 1.5], showticklabels=False, row=3, col=1)

        chr_fig.update_xaxes(
            fixedrange=True,
            title_text="Sequence position",
            range=[0, index.chrs.loc[anchor_name, this_chr]["size"]],
            row=4,
            col=1,
        )

        chr_fig.update_yaxes(title_text="G.", range=[0.5, 1.5], showticklabels=False, row=4, col=1)
        chr_fig.update_layout(
            clickmode="event+select",
            dragmode="select",
            selectdirection="h",
            height=350,
            margin=dict(b=10, l=10, r=10),
            font=dict(size=16),
        )
        chr_fig.update_layout(
            hoverlabel=dict(bgcolor="white", font_size=35, font_family="Balto"),
        )
        return chr_fig

    def plot_whole_genome(anchor_name):
        spec = []
        sub_titles = []
        h = index[anchor_name].chr_count * 250
        chr_names = list(index[anchor_name].chrs.index)
        # print(chr_names)
        for chrom in range(index[anchor_name].chr_count):  # i in range(0, num_chrs[anchor_name]):
            spec.append(
                [
                    {
                        "type": "heatmap",
                    }
                ]
            )
            spec.append(
                [
                    {
                        "type": "heatmap",
                    }
                ]
            )
            spec.append(
                [
                    {
                        "type": "heatmap",
                    }
                ]
            )
            spec.append(
                [
                    {
                        "type": "heatmap",
                    }
                ]
            )
            sub_titles.append(chr_names[chrom])
            sub_titles.append("")
            sub_titles.append("")
            sub_titles.append("")

        wg_fig = make_subplots(
            rows=(4 * index[anchor_name].chr_count),
            cols=1,
            specs=spec,  # [[{"type": "heatmap",}], [{"type": "heatmap",}], [{"type": "heatmap"}]],
            shared_xaxes=True,
            subplot_titles=sub_titles,  # ("K-mer and gene density accross whole chromosome", "", ""),
            # y_title="TESTING",
            vertical_spacing=0.0,
        )
        umap_data, zmax = anchor_umaps_per_chr(anchor_name)
        cntr = 1
        for chrom in index.chrs.loc[anchor_name].index:
            x = list(index.bitsum_bins.loc[(anchor_name, chrom)].index)
            if len(x) != 1:
                if chrom in umap_data.keys():
                    wg_fig.append_trace(
                        go.Heatmap(
                            x=umap_data[chrom]["umap1"],
                            z=umap_data[chrom]["cluster"],
                            y=[1] * (len(umap_data[chrom]["umap1"])),
                            type="heatmap",
                            hovertemplate="Coord: %{x}<br>Cluster: %{z}",
                            colorscale="Portland",
                            zmin=0,
                            zmax=zmax,
                            showlegend=False,
                            showscale=False,
                        ),
                        # layout=go.Layout(yaxis=dict(title=dict(text="LABELME"))),
                        row=((cntr * 4) - 3),
                        col=1,
                    )
                    wg_fig.update_yaxes(title_text="Cluster", row=((cntr * 4) - 3), col=1)
                z = index.bitsum_bins.loc[(anchor_name, chrom), index.ngenomes]
                wg_fig.append_trace(
                    go.Heatmap(
                        x=x,
                        z=index.bitsum_bins.loc[(anchor_name, chrom), index.ngenomes],
                        hovertemplate="Coord: %{x}<br>Universal: %{z}",
                        y=[1] * (len(x)),
                        type="heatmap",
                        colorscale="magma_r",
                        showlegend=False,
                        showscale=False,
                    ),
                    row=((cntr * 4) - 2),
                    col=1,
                )
                wg_fig.update_yaxes(title_text="Univ.", row=((cntr * 4) - 2), col=1)
                wg_fig.append_trace(
                    go.Heatmap(
                        x=x,
                        z=index.bitsum_bins.loc[(anchor_name, chrom), 1],
                        y=[1] * (len(x)),
                        type="heatmap",
                        colorscale="magma",
                        hovertemplate="Coord: %{x}<br>Private: %{z}",
                        showscale=False,
                    ),
                    row=((cntr * 4) - 1),
                    col=1,
                )
                wg_fig.update_yaxes(title_text="Priv.", row=((cntr * 4) - 1), col=1)

            if cntr == 1:
                wg_fig.update_layout(xaxis={"side": "top"})
            cntr += 1

        wg_fig.update_layout(
            clickmode="event",
            plot_bgcolor="rgba(0,0,0,0)",
            height=h,
            hoverlabel=dict(bgcolor="white", font_size=35, font_family="Balto"),
            font=dict(size=18, family="Balto"),
        )
        wg_fig.update_annotations(font_size=25)
        wg_fig.update_yaxes(showticklabels=False)
        # wg_fig.update_layout(height=h)
        return wg_fig

    def plot_gene_content(all_genes, anchor_name, chrom, start_coord, end_coord):

        sort_by = ["Unique", "Universal"]
        f = index.bitfreq_chrs.loc[(anchor_name, chrom)]
        uniq_avg = f.loc[1]  # tmp[1]
        univ_avg = f.loc[index.ngenomes]  # tmp[-1]
        universals = all_genes[index.ngenomes]
        uniques = all_genes[1]
        sizes = all_genes["end"] - all_genes["start"]

        local_genes = all_genes[all_genes["local"]]

        colors = ["#ffd60a", "#440154"]
        # d = {
        #    "Names":all_genes["name"],
        #    "Universal": universals/sizes,
        #    "Unique":uniques/sizes
        # }
        uniq_avg = uniq_avg / 100
        univ_avg = univ_avg / 100
        # df = pd.DataFrame(d)

        x = np.arange(len(all_genes), dtype=int)  # []

        # df_sorted = df.sort_values(sort_by[-1])

        df = all_genes[["name", "start", "end", "local"]]
        df["Universal"] = all_genes[index.ngenomes] / sizes
        df["Unique"] = all_genes[1] / sizes
        df_sorted = df.sort_values(sort_by[-1])
        df_sorted["X"] = x

        cntr = 0
        fig = go.Figure(
            data=[
                go.Scattergl(
                    x=x,
                    y=df_sorted[sort_by[cntr]],
                    text=df_sorted["name"],
                    marker=dict(color=colors[cntr]),
                    name="% " + sort_by[cntr],
                    mode="markers",
                )
            ]
        )
        cntr += 1
        while cntr < len(sort_by):  # s in sort_by:
            fig.add_trace(
                go.Scattergl(
                    x=x,
                    y=df_sorted[sort_by[cntr]],
                    xaxis="x",
                    yaxis="y",
                    text=df_sorted["name"],
                    customdata=df_sorted[["start", "end"]].to_numpy(),
                    marker=dict(color=colors[cntr]),
                    name="% " + sort_by[cntr],
                    mode="markers",
                )
            )
            cntr += 1
        fig.update_layout(clickmode="event+select")
        fig.update_layout(hovermode="x unified")

        df2 = df_sorted[df_sorted["local"]]

        fig.add_trace(
            go.Scattergl(
                x=df2["X"],
                y=df2["Universal"],
                xaxis="x",
                yaxis="y",
                marker=dict(color="#FF2192", size=10),
                mode="markers",
                hoverinfo="skip",
                showlegend=False,
            )
        )
        fig.add_trace(
            go.Scattergl(
                x=df2["X"],
                y=df2["Unique"],
                xaxis="x",
                yaxis="y",
                marker=dict(color="#FF2192", size=10),
                mode="markers",
                hoverinfo="skip",
                showlegend=False,
            )
        )

        fig.add_hline(y=uniq_avg, line_dash="dash", line_color="goldenrod")
        fig.add_hline(y=univ_avg, line_dash="dash", line_color="#440154")
        fig.add_hline(y=0, line_color="black")
        fig.add_vline(x=0, line_color="black")
        # fig.update_layout(height=500)
        fig.update_layout(
            title={
                "text": "k-mer content in individual genes",
                "xanchor": "center",
                "x": 0.5,
                "y": 0.85,
                "yanchor": "top",
            },
            xaxis_title=sort_by[0] + " k-mer content rank",
            yaxis_title="% of k-mers",
            plot_bgcolor="rgba(0,0,0,0)",
            # margin=dict(
            #    t=20,
            #    b=10,
            #    l=10,
            #    r=10),
            font=dict(
                size=20,
            ),
        )
        return fig

    # Tab 1 bar chart of each genome

    # Tab 1 dendogram
    def read_genome_size_files():
        seqs = []
        sizes = []
        for l in index.genome_names:
            seqs.append([l, len(index.chrs.loc[l].index)])
            sizes.append([l, index.chrs.loc[l]["size"].sum()])

        sizes_df = pd.DataFrame(sizes, columns=["Sample", "Size"])
        seqs_df = pd.DataFrame(seqs, columns=["Sample", "Seq"])

        sizes_df = sizes_df.sort_values(by=["Size"])
        seqs_df = seqs_df.sort_values(by=["Seq"])
        return sizes_df, seqs_df

    def make_genome_count_plot(anchor_name):
        counts = index.genome_sizes["chr_count"].sort_values()
        # fig = make_subplots(rows=1, cols=1)
        fig = go.Figure(data=[go.Scattergl(x=counts.index, y=counts)])
        fig.update_yaxes(
            title_text="# of scaffolds",
        )
        fig.add_vline(x=anchor_name, line_dash="dash", line_color="darkblue")
        fig.update_layout(
            font=dict(
                size=20,
            )
        )
        return fig

    def make_genome_size_plot(anchor_name):
        lens = index.genome_sizes["length"].sort_values()
        # fig = make_subplots(rows=1, cols=1)
        fig = go.Figure(data=[go.Scattergl(x=lens.index, y=lens)])
        fig.update_yaxes(
            title_text="Size of genome",
        )
        fig.add_vline(x=anchor_name, line_dash="dash", line_color="darkblue")
        fig.update_layout(
            font=dict(
                size=20,
            )
        )

        return fig

    def plot_avgs(anchor_name):
        # fig = make_subplots(rows=1, cols=1)
        fig = go.Figure(
            data=[go.Scattergl(x=index.bitsum_totals_avg.index, y=index.bitsum_totals_avg)]
        )
        fig.add_vline(x=anchor_name, line_dash="dash", line_color="darkblue")
        fig.update_yaxes(
            title_text="Average k-mer",
        )
        fig.update_layout(
            font=dict(
                size=20,
            )
        )
        return fig

    def make_avg_kmer_fig(this_anchor):
        # fig = make_subplots(rows=1, cols=1)

        avgs = index.bitsum_chrs_avg[this_anchor]

        fig = go.Figure(data=[go.Scattergl(x=avgs.index, y=avgs)])
        fig.update_yaxes(
            title_text="Average k-mer",
        )
        fig.update_xaxes(
            title_text="Chromosome",
        )
        fig.update_layout(plot_bgcolor="rgba(0,0,0,0)", font=dict(size=20))

        fig.update_layout(
            updatemenus=[
                dict(
                    type="buttons",
                    direction="left",  # showactive=True, x=0, y=0,
                    pad={"r": 10, "t": 10},
                    showactive=True,
                    x=0,  # 0.11,
                    xanchor="left",
                    y=1.2,
                    yanchor="top",
                    buttons=list(
                        [
                            dict(
                                args=[{"yaxis.type": "linear"}], label="Linear", method="relayout"
                            ),
                            dict(args=[{"yaxis.type": "log"}], label="Log", method="relayout"),
                        ]
                    ),
                ),
            ]
        )

        return fig

    def anchor_umaps_per_chr(anchor_name):
        data = {}
        umap = index.genomes[anchor_name].genome_umap

        chroms = index.genomes[anchor_name].chrs.index
        data = {c: umap.query("chrom == @c")[["umap1", "umap2", "cluster"]] for c in chroms}
        zmax = umap["cluster"].max()

        return data, zmax

    def make_genome_umap(anchor_name):
        umap = index.genomes[anchor_name].genome_umap
        umap_x = umap["umap1"]
        umap_y = umap["umap2"]
        color = umap["cluster"]
        chroms = umap["chrom"]
        chrom_start = umap["start"]
        chrom_stop = umap["end"]

        umap_cax = {
            "colorscale": "Portland",
            "colorbar": {
                "title": "DBSCAN Cluster",
                "y": 0,
                "len": 1.0,
                "yanchor": "bottom",
                "title_side": "right",
            },
        }
        hover_text = [
            f"Cluster: {color[i]}<br>Chromosome: {chroms[i]}:{chrom_start[i]}-{chrom_stop[i]}<br>X: {umap_x[i]}<br>Y: {umap_y[i]}"
            for i in range(len(umap_x))
        ]
        umap_fig = go.Figure(
            data=[
                go.Scattergl(
                    x=umap_x,
                    y=umap_y,
                    ids=color,
                    hoverinfo="text",
                    hovertext=hover_text,  # hovertemplate = "%{ids}: <br>Popularity: %{umap_x} </br> %{umap_y}",
                    mode="markers",
                    marker=dict(
                        color=color,
                        # colorscale="Portland",
                        coloraxis="coloraxis",
                        showscale=True,
                        size=14,
                        opacity=0.7,
                    ),
                )
            ]
        )
        umap_fig.update_layout(
            plot_bgcolor="rgba(0,0,0,0)",
            font=dict(size=25),
            coloraxis=umap_cax,
        )
        umap_fig.update_layout(
            title=dict(
                text="Pangenome conservation UMAP (10kbp regions)",
                y=0.9,
                x=0.5,
                xanchor="center",
                yanchor="top",
                font=dict(size=35),
            ),
            hoverlabel=dict(bgcolor="white", font_size=35, font_family="Balto"),
        )
        umap_fig.add_vline(x=0, line_color="black")
        umap_fig.add_hline(y=0, line_color="black")

        return umap_fig

    def get_genome_umap_hist(anchor_name):
        umap = index.genomes[anchor_name].genome_umap
        color = umap["cluster"]
        cluster_size = umap["end"] - umap["start"]
        # Convert to numpy array if it's a pandas Series
        if hasattr(cluster_size, "values"):
            cluster_size = cluster_size.values
        x, y, z = [], [], []

        for i in range(0, max(color) + 1):
            x.append(i)
            z.append(i)
            y.append((color == i).sum() * cluster_size[0])

        hist_fig = go.Figure(
            [
                go.Bar(
                    x=y,
                    y=x,
                    marker_color=z,
                    orientation="h",
                    marker={
                        "colorscale": "Portland",
                    },
                )
            ]
        )
        hist_fig.update_xaxes(title_text="Base pairs per cluster")
        hist_fig.update_xaxes(type="log")
        hist_fig.update_layout(
            plot_bgcolor="rgba(0,0,0,0)",
            font=dict(size=25),
            hoverlabel=dict(bgcolor="white", font_size=35, font_family="Balto"),
        )

        return hist_fig

    def make_genes_per_chr_fig(anchor_name):
        # fig = go.Figure(data=[(rows=1, cols=1)
        # x = []
        # y = []
        # for this_chrom in index[anchor_name].chrs.index:
        #    x.append(this_chrom)
        #    try:
        #        sys.stderr.write("Quering genes 2\n")
        #        sys.stderr.flush()
        #        ngenes = len(index.query_genes(anchor_name, this_chrom, 0, index.chrs.loc[anchor_name, this_chrom]["size"]))
        #        print(ngenes)
        #        ngenes = index.chrs.loc[anchor_name, this_chrom]["gene_count"]
        #        print(ngenes)
        #        y.append(ngenes/index.chrs.loc[anchor_name, this_chrom]["size"])
        #        sys.stderr.write("Queried genes 2\n")
        #        sys.stderr.flush()
        #    except :
        #        print("failed", anchor_name + "\t" + this_chrom)
        chrs = index[anchor_name].chrs
        fig = go.Figure(data=[go.Scattergl(x=chrs.index, y=chrs["gene_count"])])
        fig.update_yaxes(
            title_text="Gene density",
        )
        fig.update_xaxes(
            title_text="Chromosome",
        )
        fig.update_layout(plot_bgcolor="rgba(0,0,0,0)", font=dict(size=20))
        fig.update_layout(
            updatemenus=[
                dict(
                    type="buttons",
                    direction="left",  # showactive=True, x=0, y=0,
                    pad={"r": 10, "t": 10},
                    showactive=True,
                    x=0,  # 0.11,
                    xanchor="left",
                    y=1.2,
                    yanchor="top",
                    buttons=list(
                        [
                            dict(
                                args=[{"yaxis.type": "linear"}], label="Linear", method="relayout"
                            ),
                            dict(args=[{"yaxis.type": "log"}], label="Log", method="relayout"),
                        ]
                    ),
                ),
            ]
        )
        return fig

    def plot_perc_shared(perc_shared, anchor):
        df = index.genomes[
            anchor
        ].total_paircounts  # pd.read_csv(perc_shared_file,index_col="name") #np.loadtxt(perc_shared_file, delimiter=",")

        # print(data.loc["Pandan_hap1"][2])
        # print(perc_shared)
        df = df.sort_values(by="frac")
        names, data, colors = [], [], []

        # for k in perc_shared.keys():
        for row in df.iterrows():
            k = row[0]
            names.append(k)
            tmp = perc_shared[k] - (df.loc[k]["frac"] * 100)
            if tmp > 0:
                colors.append("#00b4d8")
            elif tmp == 0:
                colors.append("grey")
            else:
                colors.append("#c1121f")
            data.append(tmp)
        perc_shared_fig = go.Figure([go.Bar(x=names, y=data, marker_color=colors)])
        # perc_shared_fig = go.Figure([go.Bar(x=list(perc_shared.keys()),y=list(perc_shared.values()))])
        perc_shared_fig.update_layout(
            plot_bgcolor="rgba(0,0,0,0)",
            font=dict(size=25),
            hoverlabel=dict(font_size=35, font_family="Balto"),
        )
        perc_shared_fig.update_yaxes(
            title_text="% shared globally - % shared in this region",
        )
        return perc_shared_fig

    def make_anchor_tab_legends(anchor_name):
        fig = make_subplots(rows=1, cols=1)
        fig.update_layout(plot_bgcolor="rgba(0,0,0,0)", font=dict(size=20))
        """
        a = np.array([[1, 0]])
        fig = go.Figure(
            data=[
                go.Heatmap(
                    
                    z=a,
                    colorscale="Portland",
                    showscale=True,
                    zmin=a.min(),
                    zmax=a.max(),
                    #colorbar=dict(x=0),
                ),
                #go.Heatmap(
                #    z=np.zeros_like(a), # ensure there are enough points to cover the whole image
                #    colorscale="Picnic_r",  # any colorscale that has white at 0
                #    showscale=False,
                #    zmid=0,
                #),
            ],
        )
        fig.update_layout(
            #width=75, # 75 seems to work empirically, you may have to adjust the width
            height=5,
            xaxis_showgrid=False,
            yaxis_showgrid=False,
            xaxis_zeroline=False,
            yaxis_zeroline=False,
            xaxis_visible=False,
            yaxis_visible=False,
            margin=dict(l=0, r=0, b=0, t=0), # zeroing the margins gives a tighter figure, and is also important to make the other numbers work out
        )
        fig.update_traces(colorbar_orientation='h', selector=dict(type='heatmap'))
        """
        fig.update_xaxes(visible=False)
        fig.update_yaxes(visible=False)
        return fig

    def make_gene_per_genome_fig(anchor_name):
        fig = make_subplots(rows=1, cols=1)
        colors = ["#ffd60a", "#440154"]
        gene_names = []
        universals = []
        uniques = []
        sizes = []
        genes = list()

        tic = time.perf_counter()
        genes = index.query_genes(anchor_name)
        tib = time.perf_counter()
        sys.stderr.write(f"Queried genes 3 ({tib-tic})\n")
        # genes["name"] = genes["attr"].str.extract("Name=([^;]+)")
        genes["size"] = genes["end"] - genes["start"]

        x = [i for i in range(0, len(genes))]
        genes["universal"] = genes[index.ngenomes] / genes["size"]
        genes["unique"] = genes[1] / genes["size"]
        df_sorted = genes.sort_values("universal")
        df_sorted["X"] = x
        fig.add_trace(
            go.Scattergl(
                x=x,
                y=df_sorted["unique"],
                text=df_sorted["chr"] + ":" + df_sorted["name"],
                marker=dict(color=colors[0]),
                name="Proportion " + "Unique",
                mode="markers",
            )
        )
        fig.add_trace(
            go.Scattergl(
                x=x,
                y=df_sorted["universal"],
                text=df_sorted["chr"] + ":" + df_sorted["name"],
                marker=dict(color=colors[1]),
                name="Proportion " + "Universal",
                mode="markers",
            )
        )
        fig.update_layout(plot_bgcolor="rgba(0,0,0,0)", font=dict(size=20))
        fig.update_xaxes(
            title_text="Genes",
        )
        fig.update_yaxes(
            title_text="Proportion of gene",
        )
        fig.update_layout(hovermode="x unified")

        return fig

    def get_buffer(tmp_start, tmp_stop, n_skips):
        min_dist = n_skips * params.max_chr_bins * 8
        act_dist = tmp_stop - tmp_start
        if act_dist >= min_dist:
            return 0  # The distance between start and stop are good enough
        else:  #
            return 1  # int((min_dist - act_dist)/2)+1

    def set_coords(anchor, chrom=None, start=0, end=None):
        anchor_chrs = index.chrs.loc[anchor].index
        if chrom not in anchor_chrs:
            chrom = anchor_chrs[0]
        if end is None:
            end = index.chrs.loc[(anchor, chrom), "size"]

        return anchor, chrom, start, end

    @app.callback(
        Output("selected-anchor-state", "children"),  # pangenome
        Output("selected-chrom-state", "children"),  # anchor
        Output("start-coord-state", "children"),  # start_coord div (constant?)
        Output("end-coord-state", "children"),  # x_end div (constant?)
        # Output('chr-genes-state','children'),   #x_end div (constant?)
        Output("chromosome", "relayoutData"),  # chromosome
        Output("Chrs_Info", "value"),  # chromosome
        Output("chr_name", "children"),  # chromosome
        # Output('chromosome','figure'),
        # Output('gene-names-state','children'),
        # Annotation output
        # Output('annotation_tab_clickdata_Name', 'children'),
        # Output('annotation_tab_clickdata_X', 'children'),
        # Output('annotation_tab_clickdata_Y', 'children'),
        # Output('annotation_tab_clickdata_Color', 'children'),
        # Output('annotation_tab_clickdata_Attr', 'children'),
        # Input('annotation_conservation', 'clickData'),
        Input("genome_select_dropdown", "value"),  # pangenome, anchor, chromosome
        Input("chr_select_dropdown", "value"),  # anchor, chromosome
        Input("all_chromosomes", "relayoutData"),  # anchor -> chromosome
        Input("primary", "clickData"),  # chromosome
        Input("primary", "relayoutData"),  # chromosome
        Input("chromosome", "selectedData"),  # chromosome
        Input("Genes", "clickData"),  # chromosome
        Input("Chrs_Info", "value"),  # chromosome
        Input("pre_regions_dropdown", "value"),  # chromosome
        # Input('all_chromosomes','relayoutData'), #anchor -> chromosome
        State("selected-anchor-state", "children"),
        State("selected-chrom-state", "children"),
        State("start-coord-state", "children"),
        State("end-coord-state", "children"),
        State("chr-genes-state", "children"),
        # State('chromosome','figure'),
    )
    def nav_callback(
        anchor_dropdown,
        chr_dropdown,
        anctab_chrs_relayout,
        chrtab_primary_click,
        chrtab_primary_relayout,
        chrtab_chr_select,
        chrtab_gene_click,
        user_chr_coords,
        pre_selected_region,
        anchor,
        chrom,
        start_coord,
        end_coord,
        chr_genes,
    ):
        triggered_id = ctx.triggered_id

        n_skips = 100
        anchor_init = anchor
        chrom_init = chrom
        start_init = start_coord
        end_init = end_coord

        start_coord = None
        end_coord = None
        anchor_out = None

        tmp_name, tmp_x, tmp_y, tmp_color, tmp_attr = (
            no_update,
            no_update,
            no_update,
            no_update,
            no_update,
        )  # json.dumps(clickData, indent=2)

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

        # Select start/end coordinates, triggers CHROMOSOME
        elif triggered_id == "Chrs_Info":
            start_coord = int(user_chr_coords.strip().split("-")[0])
            end_coord = int(user_chr_coords.strip().split("-")[1])

        # Select chromosome (region) in anchors tab, triggers CHROMOSOME
        elif triggered_id == "all_chromosomes":
            # Get the chromosome number
            if anctab_chrs_relayout != None and len(anctab_chrs_relayout) > 1:
                chr_num_tmp = (
                    list(anctab_chrs_relayout.keys())[0].split(".")[0].split("axis")[1]
                )  # if anctab_chrs_relayout != None and len(anctab_chrs_relayout)>1:
                if len(chr_num_tmp) == 0:
                    chr_num = 1
                else:
                    chr_num = int(int(chr_num_tmp) / 3) + 1

                chrom = index.chrs.index[chr_num - 1]
                anchor, chrom, start_coord, end_coord = set_coords(anchor, chrom)
            anchor_out = anchor if anchor != anchor_init else no_update

        # Chromosome gene plot, triggers CHROMOSOME
        elif triggered_id == "Genes":
            if chrtab_gene_click != None:
                start_coord, end_coord = chrtab_gene_click["points"][0]["customdata"]

                # print(chrtab_gene_click)
                # print(gene_coords)
                # this_gene_name = chrtab_gene_click['points'][0]['text']
                # sys.stderr.write("Quering genes 7\n")
                # sys.stderr.flush()
                # genes = index.query_genes(anchor, chrom, 0, index.chrs.loc[anchor, chrom]["size"]).set_index("name")
                # sys.stderr.write("Queried genes 7\n")
                # sys.stderr.flush()
                #
                # this_gene = genes.loc[this_gene_name]
                # start_coord = int(this_gene['start'])
                # end_coord = int(this_gene['end'])
                # print(start_coord,end_coord)

        # Chromosome top plot, triggers CHROMOSOME
        elif triggered_id == "chromosome":
            if "range" not in chrtab_chr_select:
                return (no_update,) * 12
            if "x2" in chrtab_chr_select["range"].keys():
                start_coord = int(chrtab_chr_select["range"]["x2"][0])
                end_coord = int(chrtab_chr_select["range"]["x2"][1])
            elif "x" in chrtab_chr_select["range"].keys():
                start_coord = int(chrtab_chr_select["range"]["x"][0])
                end_coord = int(chrtab_chr_select["range"]["x"][1])
            else:
                start_coord = int(chrtab_chr_select["range"]["x3"][0])
                end_coord = int(chrtab_chr_select["range"]["x3"][1])

        # Chromosome main plot, triggers CHROMOSOME
        elif triggered_id == "primary":
            print(chrtab_primary_relayout)
            if (
                chrtab_primary_relayout != None
                and "xaxis2.range[0]" in chrtab_primary_relayout.keys()
            ):
                start_coord = int(chrtab_primary_relayout["xaxis2.range[0]"])
                end_coord = int(chrtab_primary_relayout["xaxis2.range[1]"])

        # Annotation tab triggered when gene is selected in annotation plot
        """elif triggered_id == 'annotation_conservation':
            this_gene_info = annotation_tab_df.loc[annotation_tab_df['Name'] == str(anno_clickdata['points'][0]['hovertext'])]
            print(this_gene_info['chr'].iloc[0])
            print(this_gene_info['start'].iloc[0])
            print(this_gene_info['end'].iloc[0])
            tmp_name = "Name: " + str(anno_clickdata['points'][0]['hovertext'])
            tmp_x = "X: " + str(anno_clickdata['points'][0]['x'])
            tmp_y = "Y: " + str(anno_clickdata['points'][0]['y'])
            tmp_color = "Color: " + str(anno_clickdata['points'][0]['marker.color'])
            chrom = this_gene_info['chr'].iloc[0]
            start_coord = this_gene_info['start'].iloc[0]
            end_coord = this_gene_info['end'].iloc[0]
            tmp_attr = "Attributes: " + str(annotation_tab_df.loc[annotation_tab_df['Name'] == str(anno_clickdata['points'][0]['hovertext']), 'attr'].iloc[0])
        """
        if start_coord is None or end_coord is None:
            start_out = end_out = no_update
            start_coord = start_init
            end_coord = end_init

        elif (
            start_coord < 0
            or start_coord > end_coord
            or end_coord > index.chrs.loc[(anchor, chrom), "size"]
        ):
            sys.stderr.write("Error: start coord greater than end coord. Resetting\n")
            anchor = anchor_init
            chrom = chrom_init
            start_out = end_out = no_update
            start_coord = start_init
            end_coord = end_init
        else:
            start_out = start_coord
            end_out = end_coord

        if anchor_out is None:
            anchor_out = anchor  # if anchor != anchor_init else no_update

        return (
            anchor_out,
            chrom,
            start_out,
            end_out,  # chr_genes,
            chrtab_chr_select,
            update_output_div(chrom, start_coord, end_coord, anchor),
            update_out_chr(chrom, anchor),
            # tmp_name, tmp_x, tmp_y, tmp_color, tmp_attr
        )

    @app.callback(
        Output("annotation_conservation", "figure"),
        Input("radios", "value"),
    )
    def annotation_tab_callback(radio_input):
        print("Radio button")
        print(radio_input)
        log2_true = 0
        if radio_input == " Log2 gene length":
            color_by_input = "Length"
            log2_true = 1
        elif radio_input == " Gene length":
            color_by_input = "Length"
        elif radio_input == " Chromosome":
            color_by_input = "chr"
        elif radio_input == "Bed file":
            color_by_input = "bed"
        else:
            color_by_input = "selected"

        fig = go.Figure()  # annotation_tab_plot(annotation_tab_df, color_by_input, log2_true)
        return fig

    @app.callback(
        Output("all_genomes", "figure"),
        Output("all_genomes_kmer_comp", "figure"),
        Output("pangenome_num_seqs", "figure"),  # pangenome
        Output("average_kmer_content", "figure"),  # pangenome
        Output("all_genome_sizes", "figure"),  # pangenome
        Input("tabs", "value"),
        Input("selected-anchor-state", "children"),
    )
    def pangenome_callback(tab, anchor_name):
        if tab != "pangenome":
            return ({},) * 5  # + (no_update,)
        print("pangenome_callback")
        triggered_id = ctx.triggered_id
        print(triggered_id)
        if triggered_id == "tabs":
            return (
                all_genomes_dend_fig,
                pangenome_comp_fig,
                make_genome_count_plot(anchor_name),
                plot_avgs(anchor_name),
                make_genome_size_plot(anchor_name),
            )
        if triggered_id != "selected-anchor-state":
            return (no_update,) * 5  # KJ changed to 6 from 5
        return (
            all_genomes_dend_fig,
            pangenome_comp_fig,
            make_genome_count_plot(anchor_name),
            plot_avgs(anchor_name),
            make_genome_size_plot(anchor_name),
        )  # , anchor_name

    @app.callback(
        Output("all_chromosomes", "figure"),  # anchor
        # Output('gene_content', 'figure'),             #anchor
        # Output('genes_per_chr', 'figure'),            #anchor
        # Output('avg_kmer_chr', 'figure'),             #anchor
        Output("genome_umap", "figure"),  # UMAP
        Output("genome_umap_hist", "figure"),  # umap hist
        Output("chr_select_dropdown", "options"),  # anchor
        # Output('anchor_labels', 'children'),          #anchor
        # Output('selected-chrom-state', 'children'),  #anchor
        Input("tabs", "value"),
        Input("selected-anchor-state", "children"),
    )
    def anchor_callback(tab, anchor_name):
        if tab == "anchor":
            figs = (
                plot_whole_genome(
                    anchor_name
                ),  # make_anchor_tab_legends(anchor_name),#make_gene_per_genome_fig(anchor_name),
                # make_genes_per_chr_fig(anchor_name), make_avg_kmer_fig(anchor_name),
                make_genome_umap(anchor_name),
                get_genome_umap_hist(anchor_name),
            )
        else:
            figs = ({},) * 3  # Changed May 21 2025
        triggered_id = ctx.triggered_id
        return figs + (update_chromosome_list(anchor_name),)  # , anchor_name)

    @app.callback(
        Output("chromosome", "figure"),  # chromosome
        Output("primary", "figure"),  # chromosome
        Output("Third", "figure"),  # chromosome
        Output("Genes", "figure"),  # chromosome
        Output("Umap", "figure"),  # chromosome tab
        Output("Umap_genome_hist", "figure"),
        Output("regional_genes", "children"),  # chromosome
        Input("tabs", "value"),  # all
        Input("start-coord-state", "children"),  # start_coord div (constant?)
        Input("end-coord-state", "children"),  # x_end div (constant?)
        State("selected-chrom-state", "children"),
        State("selected-anchor-state", "children"),
        # State('gene-names-state','children'),
    )
    def chromosome_callback(tab, start_coord, end_coord, chrs, anchor_name):
        print("chromosome_callback\t" + str(tab))
        # tic = time.perf_counter()
        if tab != "chromosome":
            # return ({},)*5 + (no_update,)#+ (no_update,)*4
            return (no_update,) * 7  # KJ changed from 5 to 6

        triggered_id = ctx.triggered_id

        n_skips = 100
        chr_num = (
            index[anchor_name].chrs.loc[chrs, "id"] + 1
        )  # chrs_list[anchor_name].get_loc(chrs)+1
        # This should be the first time the chromosome callback is called?
        return update_all_figs(chr_num, chrs, anchor_name, 0, start_coord, end_coord, n_skips)

    def update_all_figs(chr_num, chrom, anchor_name, redo_wg, start_coord, end_coord, n_skips):
        t_start = time.perf_counter()

        sys.stderr.write("Quering genes 8\n")
        end = index.chrs.loc[anchor_name, chrom]["size"]
        sys.stderr.write(f"{anchor_name} {chrom} {end}\n")
        sys.stderr.flush()

        all_genes = index.query_genes(
            anchor_name, chrom, 0, index.chrs.loc[anchor_name, chrom]["size"]
        )
        all_genes["local"] = all_genes["start"].clip(start_coord, None) < all_genes["end"].clip(
            0, end_coord
        )

        sys.stderr.write("Queried genes 8\n")
        sys.stderr.flush()

        # bounds = all_genes.loc[:, ["start", "end"]]
        # bounds["break"] = None
        # bounds.to_numpy().flatten(),

        gene_names = all_genes["name"]  # [g.split(';')[0].split("=")[1] for g in genes['attr']]

        t1 = time.perf_counter()
        print(f"All_genes in {t1-t_start:0.4f} seconds")
        t0 = t1
        if get_buffer(start_coord, end_coord, index.lowres_step) == 1:
            n_skips = 1

        bitmap = index.query_bitmap(anchor_name, chrom, start_coord, end_coord, n_skips)

        t1 = time.perf_counter()
        print(f"query bitmap {t1 - t0:0.4f} seconds")
        t0 = t1

        t1 = time.perf_counter()
        print(f"tree plot in {t1 - t0:0.4f} seconds")
        t0 = t1

        bin_size = ((end_coord - start_coord) // params.max_chr_bins) + 1
        pan, pair = index.bitmap_to_bins(bitmap, bin_size)

        t1 = time.perf_counter()
        print(f"transform bitmap {t1 - t0:0.4f} seconds")
        toc_tmp_1 = t1
        bar_sum_regional = pan.sum(axis=1)
        fig1, perc_shared = plot_interactive(
            bar_sum_regional,
            anchor_name,
            chrom,
            start_coord,
            end_coord,
            n_skips,
            bitmap,
            pan,
            pair,
            all_genes.query("local"),
        )
        toc_tmp_2 = time.perf_counter()
        print(f"main fig in {toc_tmp_2 - toc_tmp_1:0.4f} seconds")

        fig2 = plot_perc_shared(perc_shared, anchor_name)

        # UMAP figure
        umap_fig = get_umap_from_file(chrom, anchor_name)
        umap_genome_hist = get_umap_chrom_hist(anchor_name, chrom)
        # sys.stderr.write("Quering genes 9\n")
        # sys.stderr.write(anchor_name+"\t"+chrom+"\t"+str(start_coord) + "\t"+str(end_coord)+"\n")
        # sys.stderr.flush()
        print("QUERY", anchor_name, chrom, int(start_coord), int(end_coord))
        local_genes = all_genes[all_genes["local"]]

        toc_tmp_3 = time.perf_counter()
        print(f"querying in {toc_tmp_3 - toc_tmp_2:0.4f} seconds")

        fig4 = plot_gene_content(all_genes, anchor_name, chrom, int(start_coord), int(end_coord))
        toc_tmp_31 = time.perf_counter()
        print(f"fig 4 plot in {toc_tmp_31 - toc_tmp_3:0.4f} seconds")

        chr_fig = plot_chr_whole(start_coord, end_coord, anchor_name, chrom, all_genes)
        toc_tmp_32 = time.perf_counter()
        print(f"chr fig in {toc_tmp_32 - toc_tmp_31:0.4f} seconds")

        # bar_sum_regional = pan.sum(axis=1) #bitmap_counts.value_counts().reindex(index.bitsum_index, fill_value=0)
        # fig2 = get_local_info(bar_sum_regional, anchor_name, chrom)
        toc_tmp_33 = time.perf_counter()
        print(f"fig 2 plot in {toc_tmp_33 - toc_tmp_32:0.4f} seconds")

        # if redo_wg == 1:
        big_plot = no_update  # plot_whole_genome(anchor_name)
        toc_tmp_5 = time.perf_counter()
        print(f"plots big in {toc_tmp_5 - toc_tmp_33:0.4f} seconds")

        pg_sizes_fig = make_genome_size_plot(anchor_name)
        toc_tmp_6 = time.perf_counter()
        print(f"plots sizes in {toc_tmp_6 - toc_tmp_5:0.4f} seconds")

        pg_scaffolds_fig = make_genome_count_plot(anchor_name)
        toc_tmp_7 = time.perf_counter()
        print(f"plots scfs in {toc_tmp_7 - toc_tmp_6:0.4f} seconds")

        pg_avgs = plot_avgs(anchor_name)
        toc_tmp_8 = time.perf_counter()
        print(f"plots avgs in {toc_tmp_8 - toc_tmp_7:0.4f} seconds")

        tab2_gene_density_fig = no_update  # make_genes_per_chr_fig(anchor_name)
        toc_tmp_9 = time.perf_counter()
        print(f"plots density in {toc_tmp_9 - toc_tmp_8:0.4f} seconds")

        tab2_avg_fig = no_update  # make_avg_kmer_fig(anchor_name)
        toc_tmp_10 = time.perf_counter()
        print(f"plots avg tab2 in {toc_tmp_10 - toc_tmp_9:0.4f} seconds")

        tab2_sorted_genes = no_update  # make_gene_per_genome_fig(anchor_name)
        toc_tmp_11 = time.perf_counter()
        print(f"plots sorted in {toc_tmp_11 - toc_tmp_10:0.4f} seconds")

        toc = time.perf_counter()
        print(f"Update all in {toc - t_start:0.4f} seconds")
        return (
            chr_fig,
            fig1,
            fig2,
            fig4,
            umap_fig,
            umap_genome_hist,
            update_gene_locals(local_genes, chrom, start_coord, end_coord, anchor_name),
        )

    def update_output_div(chrs, start, stop, anchor_name):
        return f"{start}-{stop}"

    def update_out_chr(chrs, anchor_name):
        return f"{anchor_name}.{chrs}:"

    def update_gene_locals(local_gene_list, chrs, start_coord, end_coord, anchor_name):
        printme = ""
        if local_gene_list is None or len(local_gene_list) == 0:
            printme = ""
        elif len(local_gene_list) == 1:
            printme += "Genes: "
            printme += local_gene_list["name"].iloc[0]  # + ": "
        elif len(local_gene_list) <= 25:
            printme += "Genes: "
            print(local_gene_list)
            for i in local_gene_list["name"]:
                printme += i + ", "
        else:
            printme = "Genes in this region: " + str(len(local_gene_list))
        return f"{printme}"

    def update_chromosome_list(anchor_name):
        return_me = [{"label": i, "value": i} for i in index.chrs.loc[anchor_name].index]
        return return_me

    app.run(host=params.host, port=params.port, debug=not params.ndebug)

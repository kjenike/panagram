from pathlib import Path
from panagram.index import Index
import plotly.express as px


def better_dir(item):
    # don't print hidden functions
    methods = dir(item)
    return [method for method in methods if not method.startswith("_")]


# def visualize(pair, output_file, inverse=False):
#     # take a look at what pair looks like after manipulation
#     # pair[pair >= 1] = 10
#     if inverse:
#         fig = px.imshow(
#             pair,
#             color_continuous_scale=px.colors.sequential.Plasma[::-1],
#             x=pair.columns,
#             y=pair.index,
#         )
#     else:
#         fig = px.imshow(pair, x=pair.columns, y=pair.index)
#     fig.write_image(output_file)
#     return

def visualize(pair, output_file, inverse=False):
    # take a look at what pair looks like after manipulation
    if inverse:
        fig = px.imshow(
            pair,
            color_continuous_scale=px.colors.sequential.Greens[
                ::-1
            ],  # px.colors.sequential.Plasma[::-1],
            x=pair.columns,
            y=pair.index,
            aspect="auto",
            zmin=0,
            # zmax=1,
        )
    else:
        fig = px.imshow(
            pair,
            x=pair.columns,
            y=pair.index,
            aspect="auto",
        )
    fig.update_layout(
        xaxis=dict(
            dtick=2000000,
        ),
    )
    fig.write_image(output_file)
    return


index_dir = "/home/nbrown62/data_mschatz1/nbrown62/panagram_data/tomato_sl4_flye"
anchor = "SL4"  # "SL5"
# chr_name = "BGV006775_MAS2.0ch11"
output_dir = Path(index_dir) / "introgression_analysis_v1"
output_dir.mkdir(parents=True, exist_ok=True)
index = Index(index_dir)
genome = index.genomes[anchor]
chrs = genome.sizes.keys()
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

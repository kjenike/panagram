import sys
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
from panagram.index import Index

#index_file = sys.argv[1]
num_bins = 200000 #this is actually closer to the bin size? 
#labels = [ "Solyc", "Solabu2", "Solaet3", "Solmac3", 
#        "Solpri1","Solcan1","Solqui2","Solmur2hap1","Solmur2hap2",
#    ]
num_samples = 2
n_skips = 100
chr_nums = 12
fai_f = sys.argv[1]
bit_file_prefix = sys.argv[2]#"Maize_testing/pan_anc_files/B73.anchor"


'''name_idx_short = {"0": 0, "1": 1, "3": 1, "5": 1, "4": 2, "6": 2, "8":2, "9":3}
name_idx_long = {}
for i in ["0","1","3","4","5","6","8","9"]:
    name_idx_long[i] = name_idx_short[i]
    for j in ["0","1","3","4","5","6","8","9"]:
        tmp_str = i+j
        tmp_total = name_idx_short[i] + name_idx_short[j]
        name_idx_long[tmp_str] = tmp_total
        for k in ["0","1","3","4","5","6","8","9"]:
            tmp_str = i+j+k
            tmp_total = name_idx_short[i] + name_idx_short[j] + name_idx_short[k]
            name_idx_long[tmp_str] = tmp_total
'''
def make_bars(cnts, Y):
    #X is the position in the chromosome 
    x = []
    y_uniq = []
    y_univ = []
    y_name_uniq = 'chr' + str(Y) + "_unique"
    y_name_univ = 'chr' + str(Y) + "_univ"
    cntr = num_bins
    while cntr < (len(cnts)*n_skips):
        x.append(cntr)
        y_uniq.append(y_name_uniq)
        y_univ.append(y_name_univ)
        cntr += num_bins
    z_unique  = [0]*(len(x)+1)
    z_univ    = [0]*(len(x)+1)

    cntr = 0
    for i in cnts:
        tmp_x = int(cntr/num_bins)
        #i = name_idx_long[i]
        if i == 1:#8:
            z_unique[tmp_x] += 1
        elif i == num_samples:
            z_univ[tmp_x] += 1
        cntr += n_skips
    return x, y_univ, y_uniq, z_unique, z_univ


anchor = Index(bit_file_prefix)
anchor_name = "Spri1"
#chr_lens = []
'''
with open(fasta_index, "r") as f:
	line = f.readline()
	while line:
		chr_lens.append(int(line.split('\t')[1]))
		line = f.readline()
'''
#chr_lens = [308452471, 243675191, 238017767, 250330460, 226353449, 
#181357234, 185808916, 182411202, 163004744, 152435371]


chr_names = []
chr_lens = []
cntr = 0
with open(fai_f, "r") as f:
    line = f.readline()
    while line and cntr < chr_nums:
    	#print(line)
    	chr_names.append(line.strip().split('\t')[0])
    	chr_lens.append(int(line.strip().split('\t')[1]))
    	cntr += 1
    	line = f.readline()

for i in range(1,len(chr_lens)+1):
	#all_x = []
	#all_y_unique = []
	#all_y_univ = []
	#all_z_unique = []
	#all_z_univ = []
	chr_start = 0
	chr_end = chr_lens[i-1]
	cntr = 0
	chr1 = anchor.query_bitmap(anchor_name, "chr" + str(i), chr_start, chr_end, 100).sum(axis=1)
	#print(chr1)
	x, y_univ, y_uniq, z_unique, z_univ = make_bars(chr1, cntr+1)

	x_string = ""
	z_univ_string = ""
	z_unique_string = ""
		
	for i in range(0, len(x)):
		x_string += str(x[i])+","
	for i in range(0, len(x)):
		z_univ_string += str(z_univ[i])+","
	for i in range(0, len(x)):
		z_unique_string += str(z_unique[i])+","
		
	print(x_string)
	print(z_univ_string)
	print(z_unique_string)

'''
for i in range(0, len(chr1)):
	x, y_univ, y_uniq, z_unique, z_univ = make_bars(cnts, cntr+1)
		
	x_string = ""
	z_univ_string = ""
	z_unique_string = ""
		
	for i in range(0, len(x)):
			x_string += str(x[i])+","
	for i in range(0, len(x)):
			z_univ_string += str(z_univ[i])+","
	for i in range(0, len(x)):
			z_unique_string += str(z_unique[i])+","
		
	print(x_string)
	print(z_univ_string)
	print(z_unique_string)

'''
'''
fig = make_subplots(rows=2, cols=1, 
	specs=[[{"type": "heatmap",}], [{"type": "heatmap",}]],
	shared_xaxes=True,
	vertical_spacing = 0.0
	)

#for k in range(0, len(all_x)):
fig.add_trace(go.Heatmap(
		z=z_unique,
		x=x,
		y=y_uniq, 
		showlegend=False, showscale=False), row=2, col=1)
fig.add_trace(go.Heatmap(
		z=z_univ,
		x=x,
		y=y_univ), row=1, col=1)
'''	

#fig.show()





























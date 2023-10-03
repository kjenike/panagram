import sys
from panagram.index import Index 
import numpy as np
import pandas as pd
import statistics
import matplotlib.pyplot as plt
import plotly.express as px

#First argument: the directory where panagram index was run
index_loc = sys.argv[1]
#Second argument: anchor name
anchor_name = sys.argv[2]
#Third argument: do you want to query the genes, annotation or bits? [gene, anno, bit]
ele = sys.argv[3]
#If you are looking at the bitmap, it's probably best to specify how large the blocks should be. By default I'll set it to 200kbp
#if ele == "bit":
    #block_size = int(sys.argv[4])
block_size = 200000
skips = 100 #You can change this as well. This is only relevant to the bitmap option. Higher number will make this go faster. At the expense of accuracy. 

index = Index(index_loc)
elements = list()

def query_something(ele, chrom):
    if ele == "gene":
        e = index.query_genes(anchor_name, chrom, 0, index.chrs.loc[anchor_name, chrom]["size"])
    elif ele == "anno":
        e = index.query_anno(anchor_name, chrom, 0, index.chrs.loc[anchor_name, chrom]["size"])
    elif ele == "bit":
        
        e = list()
        i = 0
        while i < index.chrs.loc[anchor_name, chrom]["size"]:
            e_tmp = index.query_bitmap(anchor_name, chrom, i, i+block_size, skips).sum(axis=1)
            i += block_size
            e.append(e_tmp)
    elif ele == "custom":
        #Here is where it gets really interesting! 
        # I didn't have time to make this fully runable. But I think you can figure it out? 
        # Fist, you need to get the kmers for a certain range
        start = 0
        stop = index.chrs.loc[anchor_name, chrom]["size"] #This will go through every kmer in the chromosome
        kmers = index.query_bitmap(anchor_name, chrom, start, stop, 1)
        #And now you can set the query. The below example will look for kmers that present in the first and third sample, absent in the second and fourth sample, and the remaining samples don't matter. 
        raw_locs = np.flatnonzero((kmers[:,0] == 1) & (kmers[:,1] == 0) & (kmers[:,2] == 1) & (kmers[:,3] == 0))
    else: 
        print("Please specify thing to query: gene, anno, bit")
    return e


for c in index.chrs.loc[anchor_name].index:
    e = query_something(ele, c)
    elements.append(e)

#for c in index.chrs.loc[anchor_name].index:
    #print(index.chr_bins.loc[(anchor_name,c)])
    #print(index.chr_bins.loc[(anchor_name,c)]['total'][1])
    #tmp = list(index.chr_bins.loc[(anchor_name,c)]['total'][1])
    #tmp_loc = list(index.chr_bins.loc[(anchor_name,c)].index)
    #d = {'total_occ_1':tmp, 'start':tmp_loc}
    #df = pd.DataFrame(d)
    #df_sorted = df.sort_values('total_occ_1')
    #df_sorted['chr'] = [c]*len(tmp)
    #df_all = pd.concat([df_all, df_sorted], axis=0)
    #print("****************")
    #print(c)
    #print(df_sorted.head(1))
    #print(df_sorted.tail(1))
    #tmp_loc = list(index.chr_bins.loc[(anchor_name,c)]['total_occ_30'])
    #tmp_sorted = tmp.sort_values('universal')
    #this_max = min(tmp)
    #for t in range(0,len(tmp)):
    #	if tmp[t] == this_max:
    #		this_index = t*200000
    #print(index.chrs.loc[anchor_name, c]["size"])
    #print(c)
    #this_start = 
    #this_stop = 
    #g = index.query_bitmap(anchor_name, c, )
    #print(g)
    #g = index.query_genes(anchor_name, c, 0, index.chrs.loc[anchor_name, c]["size"])
    #g = index.query_anno(anchor_name, c, 0, index.chrs.loc[anchor_name, c]["size"])
    #g = g.loc[g["type"]=="CDS"]
    #The bins 
    #print(index.chr_bins.loc[(anchor_name,c)])
    #g = []
    #for i in index.chr_bins.loc[(anchor_name,c)].T:
    #    cntr = 1
    #    tmp = []
    #    for j in list(index.chr_bins.loc[(anchor_name,c)].loc[i]):
    #        tmp.append(cntr*j)
    #        cntr += 1
    #    print(c + "\t" + str(i) + "\t" + str(sum(tmp)/sum(list(index.chr_bins.loc[(anchor_name,c)].loc[i]))))
        #g.append(sum(tmp)/sum(list(index.chr_bins.loc[(anchor_name,c)].loc[i])))


    #g = pd.DataFrame(g) 
    #g = index.chr_bins.loc[(anchor_name,c)]['total_occ_30']
    #g["bits"] = index.query_bitmap(anchor_name, chrom, start_coord, end_coord, 1).sum(axis=1)
    #print(g)
    #print(index.query_genes(anchor_name, c, 0, index.chrs.loc[anchor_name, c]["size"]))
    #print(c+"\t"+str(this_max)+"\t"+str(this_index))
    
    #genes.append(g)


#df_sorted = df_all.sort_values('total_occ_1')
#print(df_sorted.tail(40))

def find_univeral_genes(genes):
    genes = pd.concat(genes)
    #print(genes.tail(10))
    genes["name"] = genes["attr"].str.extract("Name=([^;]+)")
    genes["size"] = genes["end"] - genes["start"] #g["bits"] = index.query_bitmap(anchor_name, chrom, start_coord, end_coord, 1).sum(axis=1)
    genes['universal'] = genes['universal']/genes["size"]
    genes['unique'] = genes['unique']/genes["size"]
    df_sorted = genes.sort_values('universal')

    univ_genes = list(df_sorted.loc[df_sorted['universal'] == 1]['attr'])
    univ_chr = list(df_sorted.loc[df_sorted['universal'] == 1]['chr'])
    univ_start = list(df_sorted.loc[df_sorted['universal'] == 1]['start'])
    univ_end = list(df_sorted.loc[df_sorted['universal'] == 1]['end'])
    #print(univ_genes)
    cntr = 0
    for g in univ_genes:
        #print(univ_chr[cntr] + "\t" + str(univ_start[cntr]) + "\t" + str(univ_end[cntr]) +"\t" + g)
        cntr += 1


def find_gene_avgs(elements):
    #print(len(genes))
    #all_avgs = []
    #all_stds = []
    elements = pd.concat(elements)

    if ele != "bit":
        for i, row in elements.iterrows():
            these_kmers = index.query_bitmap(anchor_name, row['chr'], row['start'], row['end'], 1).sum(axis=1)
            if row['start'] != row['end']:
                this_std = these_kmers.std()#var()#these_kmers.std()
                this_avg = sum(these_kmers) / len(these_kmers)
                print(row['chr']+"\t"+str(row['start'])+"\t"+str(row['end'])+"\t"+str(this_avg)+"\t"+str(this_std)+'\t'+row['attr'])
    elif ele == "bit":
        for i, row in elements.iterrows():
            print(row)

    #genes["name"] = genes["attr"].str.extract("Name=([^;]+)")
    #genes["name"] = str(genes["attr"])+"_"+str(genes["start"])
    #genes["size"] = genes["end"] - genes["start"]
    #print(genes["attr"])
    
    #genes["bits"] = index.query_bitmap(anchor_name, chrom, start_coord, end_coord, 1).sum(axis=1)
    #genes['universal'] = genes['universal']/genes["size"]
    #genes['unique'] = genes['unique']/genes["size"]
    #for i, row in genes.iterrows():
        #print(row['attr'])
        #if row['attr'].count("gene_biotype=protein_coding") > 0:
        #these_kmers = index.query_bitmap(anchor_name, row['chr'], row['start'], row['end'], 1).sum(axis=1)
        #print(these_kmers)
        #print(row)
        #if row['start'] != row['end']:
        #    this_std = these_kmers.std()#var()#these_kmers.std()
        #    this_avg = sum(these_kmers) / len(these_kmers)
            #print(row)
        #    print(row["name"]+"\t"+row['chr']+"\t"+str(row['start'])+"\t"+str(row['end'])+"\t"+str(this_avg)+"\t"+str(this_std)+'\t'+row['attr'])
            #print(row)
        
            #print(row["type"]+"\t"+row['chr']+"\t"+str(row['start'])+"\t"+str(row['end'])+"\t"+str(this_avg)+"\t"+str(this_std)+"\t"+row["attr"])
        #for k in these_kmers:
        #    print(row['chr']+"\t"+ str(k) + "\t" + row['attr'])
            #this_avg = sum(these_kmers) / len(these_kmers)
            #print(row['name'] + "\t" + row['chr']+"\t"+str(row['start'])+"\t"+str(row['end'])+"\t"+str(this_avg)+"\t"+str(this_std)+'\t' + row['attr'])
        #all_avgs.append(this_avg)#print()
        #all_stds.append(this_std)
    #return #all_avgs, all_stds

def find_cpg_avgs(cpg_file):
    
    #genes = pd.concat(genes)
    all_chrs = []
    all_starts = []
    all_stops = []
    all_avgs = []
    all_stds = []
    with open(cpg_file, "r") as f:
        line = f.readline()
        while line:
            tmp = line.strip().split('\t')
            if len(tmp)>1:
                this_chr = tmp[1]
                #print(tmp)
                this_start = int(tmp[2])
                this_stop = int(tmp[3])
                this_kmers = index.query_bitmap(anchor_name, this_chr, this_start, this_stop, 1).sum(axis=1)
                this_std = this_kmers.std()#var()#these_kmers.std()
                this_avg = this_kmers.mean()
                #for k in this_kmers:
                #    print(this_chr+"\t"+str(this_start)+"\t"+str(this_stop)+"\t"+str(k))
                print(this_chr+"\t"+str(this_start)+"\t"+str(this_stop)+"\t"+str(this_avg)+"\t"+str(this_std))
                #all_chrs.append(this_chr)
                #all_starts.append(this_start)
                #all_stops.append(this_stop)
                #all_avgs.append(this_avg)
                #all_stds.append(this_std)

            line = f.readline()
def in_sweep_genes(this_name):
	if this_name in sweep_names:
		return "sweep"
	else:
		return "other"

def plot_avgs(file_name):
	all_avgs = []
	all_stds = []
	all_names = []
	all_lens = []
	is_sweep = []
	with open(file_name, "r") as f:
		line = f.readline()
		while line:
			tmp = line.strip().split('\t')
			if tmp[1] != "chrX":
				#print("hey!")
				#print(tmp)
				this_avg = float(tmp[4])
				all_avgs.append(this_avg)
				all_stds.append(float(tmp[5]))
				all_names.append(tmp[0])
				all_lens.append(int(tmp[3])-int(tmp[2]))
				is_sweep.append(in_sweep_genes(tmp[0]))
			line = f.readline()
	#plt.scatter(all_stds, all_avgs, alpha=0.4)
	#plt.ylabel('Average k-mer')
	#plt.xlabel('Standard deviation')
	#plt.hist(all_stds, bins=200, alpha=0.7)
	#plt.ylabel('frequency')
	#plt.xlabel('Standard deviation of k-mer composition in the gene')
	#plt.yscale('log')
	d = {"Standard_dev":all_stds, #np.log10(all_stds),#all_stds,
		"Average":all_avgs,
		"Name":all_names,
		"Length":all_lens,
		"Sweep":is_sweep
	}
	df = pd.DataFrame(d)
	fig = px.scatter(df,x="Average", y="Standard_dev", marginal_y='histogram', 
		marginal_x='histogram',
		hover_name="Name", 
		color="Sweep",
		opacity=0.5
		#color=np.log10(df["Length"])
		)
		#color="Length")
		#hover_data="Name")

	fig.show()
	#df_sorted = df.sort_values("Average")
	#plt.show()

find_gene_avgs(elements)


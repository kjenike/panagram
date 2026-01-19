import sys
from panagram.index import Index
from panagram import figs


anchor_name=sys.argv[1]
wd = sys.argv[2]
start_coord=0
n_skips=100

index  = Index("wd")
chroms = list(index[anchor_name].chrs.index)

def get_bits_from_chrom(c,data):
    end_coord = index[anchor_name].chrs.loc[c,"size"]
    bitmap = index.query_bitmap(anchor_name, c, start_coord, end_coord, n_skips)
    b = bitmap.sample(n=min(len(bitmap),50000))
    kmer_num_tmp = b.sum(axis=0)
    for k in data.keys():
        data[k] += int(kmer_num_tmp[k])
        #print()
        #print(k)
        #print(kmer_num_tmp[k])
    #print(kmer_num_tmp)

    #bitmap = index.query_bitmap(anchor_name, chrom, start_coord, end_coord, n_skips)
    #b = bitmap.sample(n=min(len(bitmap),50000))
    #kmer_num_tmp = b.sum(axis=0)
    '''
    colors = mcp.gen_color(cmap="viridis_r",n=index.ngenomes)
    total_kmers = max(kmer_num_tmp) #[-1] #kmer_num_tmp["Solqui2"]
    kmer_num = {}
    color_code = {}
    kmer_num_raw = {}
    #cntr = 0
    for k in range(0, len(kmer_num_tmp)):#kmer_num_tmp.keys():
        kmer_num[index.genome_names[k]] = float(((kmer_num_tmp[k])/total_kmers)*100)
        kmer_num_raw[index.genome_names[k]] = kmer_num_tmp[k]
        color_code[index.genome_names[k]] = palette[int(((kmer_num_tmp[k])/total_kmers)*100)+10]
    '''
    return data
#Initialize the dict! 
data = {}
for g in list(index.genomes):
    data[g] = 0
    #print(list(index.genomes))
for c in chroms:
    data = get_bits_from_chrom(c,data)

for k in data.keys():
    print(k + "," + str(data[k]) + "," + str(data[k]/data[anchor_name]*100))
#print(data)



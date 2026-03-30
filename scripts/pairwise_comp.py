import sys
from panagram.index import Index
from panagram import figs


anchor_name = sys.argv[1]
wd = sys.argv[2]
start_coord = 0
n_skips = 100

index = Index(wd)
chroms = list(index[anchor_name].chrs.index)


def get_bits_from_chrom(c, data):
    end_coord = index[anchor_name].chrs.loc[c, "size"]
    bitmap = index.query_bitmap(anchor_name, c, start_coord, end_coord, n_skips)
    b = bitmap.sample(n=min(len(bitmap), 50000))
    kmer_num_tmp = b.sum(axis=0)
    for k in data.keys():
        data[k] += int(kmer_num_tmp[k])
    return data


# Initialize the dict!
data = {}
for g in list(index.genomes):
    data[g] = 0
    # print(list(index.genomes))
for c in chroms:
    data = get_bits_from_chrom(c, data)

for k in data.keys():
    print(k + "," + str(data[k]) + "," + str(data[k] / data[anchor_name] * 100))

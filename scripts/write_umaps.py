import sys
import pandas as pd
import numpy as np
import time
from panagram.index import Index
import umap
import sklearn.cluster as cluster
from sklearn.cluster import DBSCAN

anchor_name = sys.argv[1]#human
chrom_in = sys.argv[2]#ALL
bin_size = int(sys.argv[3])#100000
eps = float(sys.argv[4])#1

n_neighbors = int(sys.argv[5])#4
md = float(sys.argv[6])#0.0 is good
stepsize = int(sys.argv[7])#100 

index  = Index(".")

paircounts = {}
for chrom in index[anchor_name].chrs.index:
    bitmap = index.query_bitmap(anchor_name, chrom, step=stepsize)
    df = index.bitmap_to_paircount_bins(bitmap,bin_size).T.fillna(0)
    paircounts[chrom] = df
paircounts = pd.concat(paircounts,names=["chrom","start"])

reducer = umap.UMAP(n_neighbors=n_neighbors, min_dist=md, n_components=2, random_state=42)

embedding = reducer.fit_transform(paircounts.to_numpy())

clusters = DBSCAN(eps = eps, min_samples = 1).fit_predict(embedding)

out = pd.DataFrame(embedding, index=paircounts.index, columns=["umap1","umap2"]).reset_index()
out["end"] = out["start"] + bin_size
out["cluster"] = clusters
print(out[["chrom","start","end","umap1","umap2","cluster"]].to_csv(sep="\t",index=False))

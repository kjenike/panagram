import sys                          
import pysam                        
import py_kmc_api as kmc            
import bgzip
import numpy as np
from time import time

print("init")
                                    
ngenomes = sys.argv[1]              
                                    
t = time()
db = kmc.KMCFile()                 
db.OpenForRA(sys.argv[2])           
print("db", time()-t)
t = time()
                                    
fasta = pysam.FastaFile(sys.argv[3])

out_prefix = sys.argv[4]
pac_out = open(f"{out_prefix}.pank", "wb")
ann_out = open(f"{out_prefix}.pana", "w")

pac_bgz = bgzip.BGZipWriter(pac_out)

vec = kmc.CountVec()

offs = 0

nbits = 32#max(8, int(8 * 2**np.ceil(np.log2(ngenomes/8))))
ann_out.write(f"{ngenomes}\t{nbits}\n")
                                    
for name in fasta.references:       
    seq = fasta.fetch(name, 0)#, 10000000) 
    print("fetch", time()-t)
    t = time()

    db.GetCountersForRead(seq, vec)
    print("count", time()-t)
    t = time()

    arr = np.array(vec, dtype="uint16")
    print("numpy", time()-t)
    t = time()

    pac_bgz.write(arr.tobytes())
    print("write", time()-t)
    t = time()

    ann_out.write(f"{offs}\t{name}\n")

    offs += len(arr)

ann_out.close()
pac_bgz.close()
pac_out.close()

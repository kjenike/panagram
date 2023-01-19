import sys
import os
import subprocess
import numpy as np
import pandas as pd
from .kmc import py_kmc_api
import bgzip
import gzip
import csv
import glob
from collections import defaultdict
from time import time
from Bio import bgzf, SeqIO

ROOT_DIR = os.path.dirname(os.path.realpath(__file__))
KMC_DIR = os.path.join(ROOT_DIR, "kmc")

BGZ_SUFFIX = "pank"
IDX_SUFFIX = "panx"
#BGZ_IDX_SUFFIX = "/anchor.panx"
#BGZ_100NT_SUFFIX = "/anchor.100nt"
ANN_SUFFIX = "pana"
    

class KmerBitmap:
    def __init__(self, prefix, mode="r"):
        self.prefix = prefix

        self.genomes = list()
        self.genome_sizes = dict()
        self.offsets = dict()
        self.bitmap_lens = defaultdict(int)

        if mode == "r":
            self._init_read()
        elif mode == "w":
            self._init_write()
        else:
            raise ValueError("Mode must be 'r' or 'w'")

    def bgz_fname(self, step): 
        return f"{self.prefix}.{step}.{BGZ_SUFFIX}"

    def idx_fname(self, step): 
        return f"{self.prefix}.{step}.{IDX_SUFFIX}"

    @property
    def ann_fname(self): 
        return f"{self.prefix}.{ANN_SUFFIX}"

    def _init_read(self):
        self.write_mode = False

        self.steps = list()
        for fname in glob.glob(f"{self.prefix}.*.{BGZ_SUFFIX}"):
            step = int(fname.split(".")[-2])
            self.steps.append(step)

        with open(self.ann_fname) as self.ann:
            self._set_genomes(self.ann.readline().strip().split(("\t")))

            offs = [0 for s in self.steps]
            for line in self.ann:
                gi, size, name = line.strip().split("\t")
                gi, size = map(int, (gi,size))


                for i,step in enumerate(self.steps):
                    self.offsets[(step,gi,name)] = offs[i]
                    offs[i] += int(np.ceil(size / step))

            for i,step in enumerate(self.steps):
                self.bitmap_lens[step] = offs[i]

        t = time()
        #self.bgz_blocks = self.load_bgz_blocks(self.idx_fname(1))
        #self.bgz = bgzf.BgzfReader(self.bgz_fname(1), "rb")

        #self.bgz_100nt_blocks = self.load_bgz_blocks(self.idx_fname(100))
        #self.bgz_100nt = bgzf.BgzfReader(self.bgz_fname(100), "rb")
        self.blocks = {s : self.load_bgz_blocks(self.idx_fname(s)) for s in self.steps}
        self.bitmaps = {s : bgzf.BgzfReader(self.bgz_fname(s), "rb") for s in self.steps}

    def load_bgz_blocks(self, fname):
        idx_in = open(fname, "rb")
        nblocks = np.fromfile(idx_in, "uint64", 1)[0]
        dtype = [("rstart", "uint64"), ("dstart", "uint64")]
        blocks = np.zeros(int(nblocks)+1, dtype=dtype)
        blocks[1:] = np.fromfile(idx_in, dtype, nblocks)
        return blocks.astype([("rstart", int), ("dstart", int)])

    def _set_genomes(self, genomes):
        self.genomes = genomes
        self.genome_ids = {g:i for i,g in enumerate(genomes)}
        self.ngenomes = len(genomes)
        self.nbytes = int(np.ceil(self.ngenomes / 8))

    def query(self, genome, name, start, end, step=1):
        print(self.genome_ids)
        gi = self.genome_ids[genome]
        seqid = (gi,name)

        bstep = 1
        for s in self.steps:
            if step % s == 0:
                bstep = max(bstep, s)

        pac = self._query(seqid, start, end, step, bstep)

        ret = np.unpackbits(pac, bitorder="little", axis=1)[:,:self.ngenomes]

        return ret

    def _query(self, name, start, end, step, bstep):
        byte_start = self.nbytes * (self.offsets[(bstep,)+name] + (start//bstep))
        length  = (end - start) // bstep

        step = step // bstep

        blk = np.searchsorted(self.blocks[bstep]["dstart"], byte_start, side="right")-1
        blk_offs = byte_start - self.blocks[bstep]["dstart"][blk]
        blk_start = self.blocks[bstep]["rstart"][blk]

        self.bitmaps[bstep].seek(bgzf.make_virtual_offset(blk_start, blk_offs))
        buf = self.bitmaps[bstep].read(length * self.nbytes)

        pac = np.frombuffer(buf, "uint8").reshape((length, self.nbytes))

        if step > 1:
            return pac[::step]
        else:
            return pac
    
    @property
    def genome_names(self):
        return list(self.genomes)

    def _init_write(self):
        self.write_mode = True
        #self.bgz = bgzip.BGZipWriter(open(self.bgz_fname, "wb"))
        #self.bgz = bgzf.BgzfWriter(self.bgz_fname(1), "wb")
        #self.bgz_100nt = bgzf.BgzfWriter(self.bgz_fname(100), "wb")

        self.ann = open(self.ann_fname, "w")

    def write(self, kmc_files, genome_tsv, steps=[1,100], anchors=None):
        self.bitmaps = {s : bgzf.BgzfWriter(self.bgz_fname(s), "wb") for s in steps}

        self.kmc_dbs = list()
        for fname in kmc_files:
            self.kmc_dbs.append(py_kmc_api.KMCFile())
            self.kmc_dbs[-1].OpenForRA(fname)

        self.steps = steps

        self.fastas = dict()
        genomes = list()
        with open(genome_tsv) as genome_file:
            genomes_in = csv.reader(genome_file, delimiter="\t")
            for name, fasta in genomes_in:
                genomes.append(name)
                self.fastas[name] = fasta

        self._set_genomes(genomes)

        self.ann.write("\t".join(self.genomes) + "\n")

        for i, (name, fasta) in enumerate(self.fastas.items()):
            if anchors is None or name in anchors:
                self._load_fasta(i, name, fasta)


    def _load_fasta(self, gi, name, fname):

        with gzip.open(fname, "rt") as fasta:
            vec = py_kmc_api.CountVec()
            t = time()
            #for seq_name in fasta.references: 
            #    seq = fasta.fetch(seq_name)#, 0, 10000000) 
            for rec in SeqIO.parse(fasta, "fasta"):
                seq_name = rec.id
                seq = str(rec.seq)

                #arrs = list()
                #arrs_100nt = list()
                byte_arrs = defaultdict(list)

                def nbytes(i):
                    if self.nbytes <= 4:
                        return self.nbytes
                    elif i == len(self.kmc_dbs)-1:
                        return self.nbytes % 4
                    else:
                        return 4
                        
                for ki,db in enumerate(self.kmc_dbs): 

                    db.GetCountersForRead(seq, vec)

                    arr = np.array(vec, dtype="uint32")
                    arr = arr.view("uint8").reshape((len(arr),4))
                    
                    arr = arr[:,:nbytes(ki)]

                    for s in self.steps:
                        byte_arrs[s].append(arr[::s])
                        
                    #arr_100nt = arr[::100]
                    #arrs.append(arr)
                    #arrs_100nt.append(arr_100nt)

                size = None
                for step,arr in byte_arrs.items():
                    arr = np.concatenate(arr, axis=1)
                    self.bitmaps[step].write(arr.tobytes())
                    self.bitmap_lens[step] += len(arr)
                    if step == 1:
                        size = len(arr)

                #arr = np.concatenate(arrs, axis=1)
                #arr_100nt = np.concatenate(arrs_100nt, axis=1)
                #self.bgz.write(arr.tobytes())
                #self.bgz_100nt.write(arr_100nt.tobytes())


                #self.ann.write(f"{self.ref_len}\t{self.ref_len_100nt}\t{name}.{seq_name}\n")
                self.ann.write(f"{gi}\t{size}\t{seq_name}\n")
                #self.bitmap_lens[1] += len(arr)
                #self.bitmap_lens[100] += len(arr_100nt)

                print(f"Anchored {name} {seq_name}", time()-t)
                t = time()

    def close(self):
        self.ann.close()
        for f in self.bitmaps.values():
            f.close()
        #self.bgz.close()
        #self.bgz_100nt.close()

        #if self.write_mode:
        #    self.kmc_db.close()

#def PangenomeIndex:
#    def __init__(self, genomes):
#        if isinstance(genomes, str):
#            self.genomes = pd.read_csv(fname, names=["name", "fasta"])
#        elif isinstance(genomes, pd.DataFrame):
#            self.genomes = genomes
#        else:
#            raise ValueError("genomes must be filename or DataFrame")
#        
#        for i,row in self.genomes.iterrows():
#            if not os.path.exists(row["fasta"]+".fai"):
#                subprocess.check_call(["samtools", "faidx", row["fasta"]])
#
#
#    def _load_genome_tsv(self, fname):
#        with open(genome_tsv) as genome_file:
#            genomes_in = csv.reader(genome_file, delimiter="\t")
#            for name, fasta in genomes_in:
#                genomes.append(name)
#                self.fastas[name] = fasta

def index(genomes, out_dir, k):

    def init_dir(name):
        d = os.path.join(out_dir, name)
        os.makedirs(d, exist_ok=True)

        return d

    count_dir = init_dir("kmc_count")
    onehot_dir = init_dir("kmc_onehot")
    bitvec_dir = init_dir("kmc_bitvec")
    tmp_dir = init_dir("tmp")

    i = 0
    db_count = 1
    samp_count = 0

    genome_dbs = list()
    
    with open(genomes) as genome_file:
        genomes_in = csv.reader(genome_file, delimiter="\t")
        for name, fasta in genomes_in:
            if i >= 32:
                i = 0
                db_count += 1


            onehot_id = str(2**i)

            print(f"Counting {name}")
            
            count_db = os.path.join(count_dir, name)
            onehot_db = os.path.join(onehot_dir, name)

            genome_dbs.append((name, os.path.abspath(onehot_db)))

            #subprocess.check_call([
            #    f"{KMC_DIR}/kmc", f"-k{k}", 
            #    "-t4", "-m8", "-ci1", "-cs10000000", "-fm",
            #    fasta, count_db, tmp_dir
            #])

            #subprocess.check_call([
            #    f"{KMC_DIR}/kmc_tools", "-t4", "transform",
            #    count_db, "set_counts", onehot_id, onehot_db
            #])

            i += 1
            samp_count += 1

    bitvec_dbs = list()

    for i in range(db_count):
        h = (i+1)*32

        if db_count == 1 or i < db_count:
            t=32
        else:
            t = samp_count-32

        opdef_fname = os.path.join(bitvec_dir, f"{i}.opdef.txt")
        bitvec_fname = os.path.abspath(os.path.join(bitvec_dir, f"{i}"))

        with open(opdef_fname, "w") as opdefs:
            opdefs.write("INPUT:\n")
            for name, db in genome_dbs:
                opdefs.write(f"{name} = {db}\n")
            opdefs.write(f"OUTPUT:\n{bitvec_fname} = {genome_dbs[0][0]}")
            for name,_ in genome_dbs[1:]:
                opdefs.write(f" + {name}")
            opdefs.write("\n-ocsum\n")

        opdefs.close()

        #subprocess.check_call([
        #    f"{KMC_DIR}/kmc_tools", "complex", opdef_fname
        #])
        
        bitvec_dbs.append(bitvec_fname)

    bits = KmerBitmap(f"{out_dir}/anchor", "w")
    bits.write(bitvec_dbs, genomes, steps=[1, 100, 1000])
    bits.close()

    for step in bits.steps:
        subprocess.check_call([
            "bgzip", "-r", f"{out_dir}/anchor.{step}.pank", 
                     "-I", f"{out_dir}/anchor.{step}.panx"])

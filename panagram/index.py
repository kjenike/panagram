import sys
import os
import subprocess
import numpy as np
import pandas as pd
import bgzip
import gzip
import csv
import glob
import pysam
from collections import defaultdict
from time import time
from Bio import bgzf, SeqIO
import toml
import multiprocessing as mp
ROOT_DIR = os.path.dirname(os.path.realpath(__file__))
KMC_DIR = os.path.join(ROOT_DIR, "kmc")

BGZ_SUFFIX = "bgz"
IDX_SUFFIX = "gzi"
ANN_SUFFIX = "pana"

MODES = {"r","w"}

#TODO
#mash
#genomes: k-mer composition for genome
#bins: precomputed kmer types for whole genome (200kb)
#chromosome: k-mer types for each chromosome
#genes: for each gene %universal+%unique
#exons+repeats: track, type, id, metadata for random access display

class Index:
    def __init__(self, prefix=None, conf=None):
        if conf is not None:
            self._load_conf(conf)
            self._init_write()
        elif prefix is not None:
            self.prefix = prefix
            self._init_read()
        else:
            raise ValueError(f"Must specify a prefix (read mode) or config file (write mode)")
        #if mode == "r":
        #    self._init_read()
        #elif mode == "w":
        #    self._init_write()
        #if mode not in MODES:
        #    raise ValueError(f"'mode' must be in {MODES}")

    PARAMS = [
        "k", 
        "anchor_only",
        "bitmap_resolutions",
        ("out_dir", "prefix"), 
        ("kmc", "kmc_args"),
        ("genomes", "genome_fastas")
    ]
    def _load_conf(self, conf):
        if conf is None: return
        elif isinstance(conf, str):
            conf = toml.load(conf)

        for p in self.PARAMS:
            if isinstance(p, tuple):
                name,tgt = p
            else:
                name = tgt = p
            setattr(self, tgt, conf[name])

        self.genome_fastas = pd.Series(self.genome_fastas, name="fasta")


    def _init_read(self):
        self.write_mode = False
        self._init_dirs()
        self.chrs = pd.read_csv(f"{self.prefix}/chrs.csv").set_index(["genome","chr"])
        self.genomes = self.chrs.index.unique("genome")

        genome_names = self.chrs.index.unique("genome")
        self.bitmaps = {genome : KmerBitmap(self, genome) for genome in genome_names}

    def _init_write(self):
        self.write_mode = True
        self._init_dirs()
    
    def _load_kmc(self, files):
        from .kmc import py_kmc_api
        self.kmc = py_kmc_api
        self.kmc_dbs = list()
        for fname in files:
            self.kmc_dbs.append(self.kmc.KMCFile())
            self.kmc_dbs[-1].OpenForRA(fname)

    def write(self):
        genomes = list()
        self.genome_ids = dict()
        for i,(name,fasta) in enumerate(self.genome_fastas.items()):
            self.genome_ids[name] = i

            if not os.path.exists(fasta+".fai"):
                cmd = ["samtools", "faidx", fasta]
                subprocess.check_call(cmd)
            fa = pysam.FastaFile(fasta)

            genomes.append(pd.DataFrame(
                [(name, chrom, i, fa.get_reference_length(chrom)) 
                 for chrom in fa.references], 
                columns=["genome", "chr", "id", "size"]))
        self.chrs = pd.concat(genomes).set_index(["genome", "chr"])
        self.genomes = self.chrs.index.unique("genome")

        if self.anchor_only:
            pre = f""
            suf = ".kmc_pre"
            bitvec_dbs = [f[:-len(suf)] 
                for f in glob.glob(f"{self.bitvec_dir}/*{suf}")]
        else:
            bitvec_dbs = self._run_kmc()

        self._load_kmc(bitvec_dbs)

        ngenomes = len(self.genome_fastas)
        self.bitmaps = dict()
        for name,fasta in self.genome_fastas.items():
            self.bitmaps[name] = KmerBitmap(self, name, "w")

        #Calculate per-chromosome k-mer occurence counts
        for i in range(len(self.genomes)):
            self.chrs[f"total_occ_{i+1}"] = 0

        for (genome,chrom),size in self.chrs["size"].items():
            popcnts = self.query_bitmap(genome, chrom, 0, size, 100).sum(axis=1)
            occs, counts = np.unique(popcnts, return_counts=True)
            for occ,count in zip(occs, counts):
                self.chrs.loc[(genome,chrom), f"total_occ_{occ}"] = count

        self.chrs.to_csv(f"{self.prefix}/chrs.csv")


    def close(self):
        for b in self.bitmaps.values():
            b.close()

    def init_dir(self, path):
        d = os.path.join(self.prefix, path)
        if self.write_mode:
            os.makedirs(d, exist_ok=True)
        return d

    def _init_dirs(self):
        self.count_dir =  self.init_dir("kmc_count")
        self.onehot_dir = self.init_dir("kmc_onehot")
        self.bitvec_dir = self.init_dir("kmc_bitvec")
        self.tmp_dir =    self.init_dir("tmp")
        self.anchor_dir = self.init_dir("anchors")

    def query_bitmap(self, genome, chrom, start, end, step):
        return self.bitmaps[genome].query(chrom, start, end, step)

    @staticmethod
    def _run_kmc_genome(args):
        args, k, i, name, fasta, count_db, onehot_db, tmp_dir = args
        memory = args["memory"]
        threads = args["threads"]

        onehot_id = str(2**i)

        subprocess.check_call([
            f"{KMC_DIR}/kmc", f"-k{k}", 
            f"-t{threads}", f"-m{memory}", "-ci1", "-cs10000000", "-fm",
            fasta, count_db, tmp_dir
        ])

        subprocess.check_call([
            f"{KMC_DIR}/kmc_tools", "-t4", "transform",
            count_db, "set_counts", onehot_id, onehot_db
        ])

        return name, onehot_db


    def _iter_kmc_genome_args(self):
        i = 0
        db_count = 0
        for name,fasta in self.genome_fastas.items():
            if i >= 32:
                i = 0
                db_count += 1

            count_db = os.path.join(self.count_dir, name)
            onehot_db = os.path.join(self.onehot_dir, name)
            tmp_dir = self.init_dir(f"tmp/{name}")
            yield (self.kmc_args, self.k, i, name, fasta, count_db, onehot_db, tmp_dir)

            i += 1

    def _run_kmc(self):

        i = 0
        db_count = 1
        samp_count = 0

        genome_dbs = list()
        
        if self.kmc_args["processes"] == 1:
            for args in self._iter_kmc_genome_args():
                genome_dbs.append(self._run_kmc_genome(args))
        else:
            with mp.Pool(processes=self.kmc_args["processes"]) as pool:
                for db in pool.imap(self._run_kmc_genome, self._iter_kmc_genome_args(), chunksize=1):
                    genome_dbs.append(db)

        bitvec_dbs = list()

        for i in range(db_count):
            h = (i+1)*32

            if db_count == 1 or i < db_count:
                t=32
            else:
                t = samp_count-32

            opdef_fname = os.path.join(self.bitvec_dir, f"{i}.opdef.txt")
            bitvec_fname = os.path.abspath(os.path.join(self.bitvec_dir, f"{i}"))

            with open(opdef_fname, "w") as opdefs:
                opdefs.write("INPUT:\n")
                for name, db in genome_dbs:
                    opdefs.write(f"{name} = {db}\n")
                opdefs.write(f"OUTPUT:\n{bitvec_fname} = {genome_dbs[0][0]}")
                for name,_ in genome_dbs[1:]:
                    opdefs.write(f" + {name}")
                opdefs.write("\n-ocsum\n")

            opdefs.close()

            subprocess.check_call([
                f"{KMC_DIR}/kmc_tools", "complex", opdef_fname
            ])
            
            bitvec_dbs.append(bitvec_fname)
        return bitvec_dbs

class KmerBitmap:
    def __init__(self, index, anchor, mode="r"):
        self.prefix = os.path.join(index.anchor_dir, anchor)
        self.anchor_name = anchor

        self.anchor_id = index.chrs.loc[anchor]["id"].iloc[0]
        self.sizes = index.chrs.loc[anchor]["size"]
        self.genomes = index.genomes
        self.ngenomes = len(self.genomes)
        self.nbytes = int(np.ceil(self.ngenomes / 8))

        self._init_steps(index)
        step_sizes = pd.DataFrame({step : np.ceil(self.sizes / step) for step in self.steps}, dtype=int)
        self.offsets = step_sizes.cumsum().shift(fill_value=0) 

        self.seq_lens = defaultdict(dict)
        self.bitmap_lens = defaultdict(int)

        if mode == "w":
            self._init_write(index, [1,100])
        elif mode != "r":
            raise ValueError("Mode must be 'r' or 'w'")

        self._init_read(index)

    def _init_steps(self, index):
        if hasattr(index, "bitmap_resolutions"):
            self.steps = index.bitmap_resolutions
        else:
            self.steps = list()
            for fname in glob.glob(f"{self.prefix}.*.{BGZ_SUFFIX}"):
                step = int(fname.split(".")[-2])
                self.steps.append(step)
        

    def bgz_fname(self, step): 
        return f"{self.prefix}.{step}.{BGZ_SUFFIX}"

    def idx_fname(self, step): 
        return f"{self.prefix}.{step}.{IDX_SUFFIX}"

    @property
    def ann_fname(self): 
        return f"{self.prefix}.{ANN_SUFFIX}"

    def _init_read(self, index):
        self.write_mode = False

        self.blocks = {s : self.load_bgz_blocks(self.idx_fname(s)) for s in self.steps}
        self.bitmaps = {s : bgzf.BgzfReader(self.bgz_fname(s), "rb") for s in self.steps}

    @property
    def chrs(self):
        return self.sizes.index

    def seq_len(self, seq_name):
        return self.sizes.loc[seq_name]

    def load_bgz_blocks(self, fname):
        idx_in = open(fname, "rb")
        nblocks = np.fromfile(idx_in, "uint64", 1)[0]
        dtype = [("rstart", "uint64"), ("dstart", "uint64")]
        blocks = np.zeros(int(nblocks)+1, dtype=dtype)
        blocks[1:] = np.fromfile(idx_in, dtype, nblocks)
        return blocks.astype([("rstart", int), ("dstart", int)])

    def query(self, name, start, end, step=1):
        bstep = 1
        for s in self.steps:
            if step % s == 0:
                bstep = max(bstep, s)

        pac = self._query(name, start, end, step, bstep)

        ret = np.unpackbits(pac, bitorder="little", axis=1)[:,:self.ngenomes]

        return ret

    def _query(self, name, start, end, step, bstep):
        byte_start = self.nbytes * (self.offsets.loc[name,bstep] + (start//bstep))
        length  = (end - start) // bstep

        step = step // bstep

        blk = np.searchsorted(self.blocks[bstep]["dstart"], byte_start, side="right")-1
        blk_offs = byte_start - self.blocks[bstep]["dstart"][blk]
        blk_start = self.blocks[bstep]["rstart"][blk]

        self.bitmaps[bstep].seek(bgzf.make_virtual_offset(blk_start, blk_offs))
        buf = self.bitmaps[bstep].read(length * self.nbytes)

        pac = np.frombuffer(buf, "uint8").reshape((len(buf)//self.nbytes, self.nbytes))

        if step > 1:
            return pac[::step]
        else:
            return pac
    
    @property
    def genome_names(self):
        return list(self.genomes)

    def _get_kmc_counts(self, db, seq):
        vec = self.kmc.CountVec()
        db.GetCountersForRead(seq, vec)

        pac = np.array(vec, dtype="uint32")
        return pac.view("uint8").reshape((len(pac),4))

    def _init_write(self, index, steps):
        self.kmc_dbs = index.kmc_dbs
        self.kmc = index.kmc
        self.anchor_only = index.anchor_only

        self.bitmaps = {s : bgzf.BgzfWriter(self.bgz_fname(s), "wb") for s in steps}

        gi = self.anchor_id
        name = self.anchor_name
        fasta_fname = index.genome_fastas.loc[self.anchor_name]

        #with gzip.open(fname, "rt") as fasta:
        with open(fasta_fname, "r") as fasta:
            t = time()
            #for seq_name in fasta.references: 
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
                    pacbytes = self._get_kmc_counts(db, seq)
                    pacbytes = pacbytes[:,:nbytes(ki)]

                    for s in self.steps:
                        byte_arrs[s].append(pacbytes[::s])

                size = None
                for step,arr in byte_arrs.items():
                    arr = np.concatenate(arr, axis=1)
                    self.bitmaps[step].write(arr.tobytes())
                    self.bitmap_lens[step] += len(arr)
                    if step == 1:
                        size = len(arr)

                self.seq_lens[gi][seq_name] = size

                print(f"Anchored {name} {seq_name}", time()-t)
                t = time()

        self.close()

        for step in self.steps:
            subprocess.check_call([
                "bgzip", "-r", self.bgz_fname(step), "-I", self.idx_fname(step)])

    def close(self):
        for f in self.bitmaps.values():
            f.close()


def index(conf): #genomes, out_dir, k):

    idx = Index(conf=conf)
    idx.write()#args)
    idx.close()


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
        self.genomes = pd.read_csv(f"{self.anchor_dir}/index.csv").set_index(["genome","chr"])

        genome_names = self.genomes.index.unique("genome")
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
        #self.genome_fastas = pd.read_csv(args.genomes, names=["name", "fasta"], sep="\t").set_index("name")["fasta"]
        #self.k = args.k
        
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
        self.genomes = pd.concat(genomes).set_index(["genome", "chr"])

        self.genomes.to_csv(f"{self.anchor_dir}/index.csv")

        #if args.anchor_only:
        #    pre = f""
        #    suf = ".kmc_pre"
        #    bitvec_dbs = [f[:-len(suf)] 
        #        for f in glob.glob(f"{args.out_dir}/kmc_bitvec/*{suf}")]
        #else:
        bitvec_dbs = self.run_kmc_bitvec()

        self._load_kmc(bitvec_dbs)

        self.anchors = dict()

        ngenomes = len(self.genome_fastas)
        self.bitmaps = dict()
        for name,fasta in self.genome_fastas.items():
            self.bitmaps[name] = KmerBitmap(self, name, "w")

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

    def run_kmc_bitvec(self):

        i = 0
        db_count = 1
        samp_count = 0

        genome_dbs = list()
        
        for name,fasta in self.genome_fastas.items():
            if i >= 32:
                i = 0
                db_count += 1


            onehot_id = str(2**i)

            count_db = os.path.join(self.count_dir, name)
            onehot_db = os.path.join(self.onehot_dir, name)

            genome_dbs.append((name, os.path.abspath(onehot_db)))

            subprocess.check_call([
                f"{KMC_DIR}/kmc", f"-k{self.k}", 
                "-t4", "-m8", "-ci1", "-cs10000000", "-fm",
                fasta, count_db, self.tmp_dir
            ])

            subprocess.check_call([
                f"{KMC_DIR}/kmc_tools", "-t4", "transform",
                count_db, "set_counts", onehot_id, onehot_db
            ])

            i += 1
            samp_count += 1

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
        self.anchor_id = index.genomes.loc[anchor]["id"].iloc[0]
        self.sizes = index.genomes.loc[anchor]["size"]
        self.genomes = index.genomes.index.unique("genome")
        self.ngenomes = len(self.genomes)
        self.nbytes = int(np.ceil(self.ngenomes / 8))

        self.seq_lens = defaultdict(dict)
        self.offsets = dict()
        self.bitmap_lens = defaultdict(int)

        

        if mode == "r":
            self._init_read(index)
        elif mode == "w":
            self._init_write(index, [1,100])
        else:
            raise ValueError("Mode must be 'r' or 'w'")

    def bgz_fname(self, step): 
        return f"{self.prefix}.{step}.{BGZ_SUFFIX}"

    def idx_fname(self, step): 
        return f"{self.prefix}.{step}.{IDX_SUFFIX}"

    @property
    def ann_fname(self): 
        return f"{self.prefix}.{ANN_SUFFIX}"

    def _init_read(self, index):
        self.write_mode = False

        self.steps = list()
        for fname in glob.glob(f"{self.prefix}.*.{BGZ_SUFFIX}"):
            step = int(fname.split(".")[-2])
            self.steps.append(step)

        offs = [0 for s in self.steps]
        for name,size in self.sizes.items():
            for i,step in enumerate(self.steps):
                self.offsets[(step,name)] = offs[i]
                offs[i] += int(np.ceil(size / step))
            
        for i,step in enumerate(self.steps):
            self.bitmap_lens[step] = offs[i]

        #with open(self.ann_fname) as self.ann:
        #    self._set_genomes(self.ann.readline().strip().split(("\t")))
        #    for line in self.ann:
        #        gi, size, name = line.strip().split("\t")
        #        gi, size = map(int, (gi,size))

        #        self.seq_lens[gi][name] = size

        #        for i,step in enumerate(self.steps):
        #            self.offsets[(step,gi,name)] = offs[i]
        #            offs[i] += int(np.ceil(size / step))

        t = time()
        self.blocks = {s : self.load_bgz_blocks(self.idx_fname(s)) for s in self.steps}
        self.bitmaps = {s : bgzf.BgzfReader(self.bgz_fname(s), "rb") for s in self.steps}

    def genome_seqs(self, genome):
        gi = self.genome_ids[genome]
        return list(self.seq_lens[gi].keys())


    def seq_len(self, genome, seq_name):
        gi = self.genome_ids[genome]
        return self.seq_lens[gi][seq_name]

    def genome_len(self, genome):
        gi = self.genome_ids[genome]
        return sum(self.seq_lens[gi].values())

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
        byte_start = self.nbytes * (self.offsets[(bstep,name)] + (start//bstep))
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

    def _get_kmc_counts(self, db, seq):
        vec = self.kmc.CountVec()
        db.GetCountersForRead(seq, vec)

        pac = np.array(vec, dtype="uint32")
        return pac.view("uint8").reshape((len(pac),4))

    def _init_write(self, index, steps):
        self.bitmaps = {s : bgzf.BgzfWriter(self.bgz_fname(s), "wb") for s in steps}
        self.kmc_dbs = index.kmc_dbs
        self.kmc = index.kmc

        self.steps = steps

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

                #arr = np.concatenate(arrs, axis=1)
                #arr_100nt = np.concatenate(arrs_100nt, axis=1)
                #self.bgz.write(arr.tobytes())
                #self.bgz_100nt.write(arr_100nt.tobytes())


                #self.ann.write(f"{self.ref_len}\t{self.ref_len_100nt}\t{name}.{seq_name}\n")
                #self.ann.write(f"{gi}\t{size}\t{seq_name}\n")
                #self.bitmap_lens[1] += len(arr)
                #self.bitmap_lens[100] += len(arr_100nt)

                print(f"Anchored {name} {seq_name}", time()-t)
                t = time()

    def close(self):
        #self.ann.close()
        for f in self.bitmaps.values():
            f.close()

        for step in self.steps:
            subprocess.check_call([
                "bgzip", "-r", self.bgz_fname(step), "-I", self.idx_fname(step)])
        #self.bgz.close()
        #self.bgz_100nt.close()

        #if self.write_mode:
        #    self.kmc_db.close()


def index(conf): #genomes, out_dir, k):

    idx = Index(conf=conf)
    idx.write()#args)
    idx.close()


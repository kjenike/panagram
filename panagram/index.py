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
from types import SimpleNamespace
ROOT_DIR = os.path.dirname(os.path.realpath(__file__))
KMC_DIR = os.path.join(ROOT_DIR, "kmc")

BGZ_SUFFIX = "bgz"
IDX_SUFFIX = "gzi"
ANCHOR_DIR = "anchors"

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
        self.conf = conf
        self.prefix = prefix
        if conf is not None:
            if prefix is None:
                self.prefix = conf.prefix
            self._init_write()
        elif prefix is not None:
            self.conf = SimpleNamespace()
            self.conf.prefix = prefix
            self._init_read()
        else:
            raise ValueError(f"Must specify a prefix (read mode) or config file (write mode)")


    def _init_read(self):
        self.write_mode = False
        self._init_dirs()
        self.chrs = pd.read_csv(f"{self.prefix}/chrs.csv").set_index(["genome","chr"])
        self.genomes = self.chrs.index.unique("genome")

        genome_names = self.chrs.index.unique("genome")
        self.bitmaps = {genome : KmerBitmap(self.conf, genome, "r", self.chrs) for genome in genome_names}

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

    @staticmethod
    def run_anchor(args):
        return KmerBitmap(*args).anchor_name

    def _init_genomes(self):
        self.genomes = self.chrs.index.unique("genome")
        self.ngenomes = len(self.genomes)
        self._total_occ_idx = pd.Index([f"total_occ_{i+1}" for i in range(self.ngenomes)])
        self._gene_occ_idx = pd.Index([f"gene_occ_{i+1}" for i in range(self.ngenomes)])

    def write(self):
        genomes = list()
        self.genome_ids = dict()
        self.genome_files = pd.DataFrame({
            "fasta" : pd.Series(self.conf.fasta),
            #"gene" : pd.Series(self.conf.gene)
        })

        if hasattr(self.conf, "gene"):
            self.genome_files["gene"] = self.conf.gene


        for i,(name,fasta) in enumerate(self.genome_files["fasta"].items()):
            self.genome_ids[name] = i

            if not os.path.exists(fasta+".fai"):
                cmd = ["samtools", "faidx", fasta]
                subprocess.check_call(cmd)
            fa = pysam.FastaFile(fasta)

            genomes.append(pd.DataFrame(
                [(name, chrom, i, fa.get_reference_length(chrom)-self.conf.k+1) 
                 for chrom in fa.references], 
                columns=["genome", "chr", "id", "size"]))

        self.chrs = pd.concat(genomes).set_index(["genome", "chr"])
        self._init_genomes()
        self.chrs.to_csv(f"{self.prefix}/chrs.csv")

        if self.conf.anchor_only or self.conf.anno_only:
            pre = f""
            suf = ".kmc_pre"
            bitvec_dbs = [f[:-len(suf)] 
                for f in glob.glob(f"{self.bitvec_dir}/*{suf}")]

        else:
            bitvec_dbs = self._run_kmc()
        
        def iter_anchor_args():
            for name,fasta in self.genome_files["fasta"].items():
                yield (self.conf, name, "w", self.chrs, fasta, bitvec_dbs)

        if not self.conf.anno_only:
            if self.conf.processes == 1:
                for args in iter_anchor_args():
                    print("Anchored", self.run_anchor(args))
            else:
                with mp.Pool(processes=self.conf.processes) as pool:
                    for name in pool.imap_unordered(self.run_anchor, iter_anchor_args(), chunksize=1):
                        print(f"Anchored {name}")
                        sys.stdout.flush()

        self.bitmaps = {
            name : KmerBitmap(self.conf, name, "r", self.chrs) for name in self.genomes
        }

        #Calculate per-chromosome k-mer occurence counts
        self.chrs[self._total_occ_idx] = 0
        self.chrs[self._gene_occ_idx] = 0

        for (genome,chrom),size in self.chrs["size"].items():
            counts = self.query_occ_counts(genome,chrom,0,size,100,self._total_occ_idx)
            self.chrs.loc[(genome,chrom), counts.index] = counts

        self.chrs.to_csv(f"{self.prefix}/chrs.csv")

        if "gene" in self.genome_files.columns:
            for g in self.genomes:
                self._load_genes(g)
        self.chrs.to_csv(f"{self.prefix}/chrs.csv")

    def query_occ_counts(self, genome, chrom, start, end, step, idx):
        popcnts = self.query_bitmap(genome, chrom, start, end, step).sum(axis=1)
        occs, counts = np.unique(popcnts, return_counts=True)
        return pd.Series(index=idx[occs-1], data=counts).reindex(idx, fill_value=0)

    TABIX_COLS = ["chr","start","end","type","attr"]
    def _iter_gff(self, fname):
        print(fname)
        for df in pd.read_csv(
            fname, 
            sep="\t", comment="#", chunksize=10000,
            names = ["chr", "source", "type", "start", "end", "score", "strand", "phase", "attr"],
            usecols = self.TABIX_COLS): yield df[self.TABIX_COLS]

    def _write_tabix(self, df, genome, typ):
        prefix = os.path.join(self.anno_dir, f"{genome}.{typ}")
        bed = prefix+".bed"
        tbx = bed+".bgz"

        df.to_csv(bed, sep="\t", header=None, index=False)
        pysam.tabix_compress(bed, tbx, True)
        pysam.tabix_index(tbx, True, 0,1,2)
        #os.remove(bed)

    def _load_genes(self, genome):
        genes = list()
        exons = list()


        for df in self._iter_gff(self.genome_files.loc[genome, "gene"]):
            genes.append(df[df["type"] == "gene"])
            exons.append(df[df["type"] == "exon"])

        def _merge_dfs(dfs):
            return pd.concat(dfs).sort_values(["chr","start"]).reset_index(drop=True)

        self._write_tabix(_merge_dfs(exons), genome, "exon")

        genes = _merge_dfs(genes)
        genes["gene_occ_1"] = 0
        univ = f"gene_occ_{self.ngenomes}"
        genes[univ] = 0
        for i,g in genes.iterrows():
            counts = self.query_occ_counts(genome, g["chr"], g["start"], g["end"], 1, self._gene_occ_idx)
            genes.loc[i, "gene_occ_1"] += counts.loc["gene_occ_1"]
            genes.loc[i, f"gene_occ_{self.ngenomes}"] += counts.loc[f"gene_occ_{self.ngenomes}"]
            self.chrs.loc[(genome,g["chr"]), counts.index] += counts
            #for occ,count in counts.items():
            #    self.chrs.loc[(genome,g["chr"]), occ] += count

        self._write_tabix(genes, genome, "gene")


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
        self.anchor_dir = self.init_dir(ANCHOR_DIR)
        self.anno_dir = self.init_dir("anno")

    def query_bitmap(self, genome, chrom, start, end, step=1):
        return self.bitmaps[genome].query(chrom, start, end, step)

    @staticmethod
    def _run_kmc_genome(args):
        conf, i, name, fasta, count_db, onehot_db, tmp_dir = args

        onehot_id = str(2**i)

        subprocess.check_call([
            f"{KMC_DIR}/kmc", f"-k{conf.k}", 
            f"-t{conf.kmc.threads}", f"-m{conf.kmc.memory}", "-ci1", "-cs10000000", "-fm",
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
        for name,fasta in self.genome_files["fasta"].items():
            if i >= 32:
                i = 0
                db_count += 1

            count_db = os.path.join(self.count_dir, name)
            onehot_db = os.path.join(self.onehot_dir, name)
            tmp_dir = self.init_dir(f"tmp/{name}")
            yield (self.conf, i, name, fasta, count_db, onehot_db, tmp_dir)

            i += 1

    def _run_kmc(self):

        i = 0
        db_count = 1
        samp_count = 0

        genome_dbs = list()
        if self.conf.kmc.processes == 1:
            for args in self._iter_kmc_genome_args():
                genome_dbs.append(self._run_kmc_genome(args))
        else:
            with mp.Pool(processes=self.conf.kmc.processes) as pool:
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
    def __init__(self, conf, anchor, mode="r", chrs=None, fasta=None, kmc_dbs=None):
        self.conf = conf
        self.prefix = os.path.join(conf.prefix, ANCHOR_DIR, anchor)
        self.anchor_name = anchor
        self.fasta = fasta

        self.anchor_id = chrs.loc[anchor]["id"].iloc[0]
        self.sizes = chrs.loc[anchor]["size"]
        self.genomes = chrs.index.unique(0)
        self.ngenomes = len(self.genomes)
        self.nbytes = int(np.ceil(self.ngenomes / 8))

        self._init_steps(conf)
        step_sizes = pd.DataFrame({step : np.ceil(self.sizes / step) for step in self.steps}, dtype=int)
        self.offsets = step_sizes.cumsum().shift(fill_value=0) 

        self.seq_lens = defaultdict(dict)
        self.bitmap_lens = defaultdict(int)

        if mode == "w":
            self._init_write(kmc_dbs)
        elif mode == "r":
            self._init_read(index)
        else:
            raise ValueError("Mode must be 'r' or 'w'")


    def _init_steps(self, conf):
        if hasattr(conf, "bitmap_resolutions"):
            self.steps = conf.bitmap_resolutions
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

    def _get_kmc_counts(self, db, seq):
        vec = self.kmc.CountVec()
        db.GetCountersForRead(seq, vec)

        pac = np.array(vec, dtype="uint32")
        return pac.view("uint8").reshape((len(pac),4))
    
    def _load_kmc(self, kmc_dbs):
        from .kmc import py_kmc_api
        self.kmc = py_kmc_api
        self.kmc_dbs = list()
        for db in kmc_dbs:
            if isinstance(db, str):
                self.kmc_dbs.append(self.kmc.KMCFile())
                self.kmc_dbs[-1].OpenForRA(db)
            else:
                self.kmc_dbs.append(db)

    def _init_write(self, kmc_dbs=None):
        self._load_kmc(kmc_dbs)

        #self.bitmaps = {s : bgzf.BgzfWriter(self.bgz_fname(s), "wb") for s in self.steps}
        self.bitmaps = {s : bgzip.BGZipWriter(open(self.bgz_fname(s), "wb"))for s in self.steps}

        gi = self.anchor_id
        name = self.anchor_name

        if self.fasta.endswith(".gz") or self.fasta.endswith(".bgz"):
            opn = lambda f: gzip.open(f, "rt")
        else:
            opn = lambda f: open(f, "r")

        with opn(self.fasta) as fasta:
        #with open(self.fasta, "r") as fasta:
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
                sys.stdout.write(f"Anchored {seq_name}\n")

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


import sys
import os
import os.path
from os import path
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

TABIX_COLS = ["chr","start","end","type","attr"]
TABIX_TYPES = {"start" : int, "end" : int}
GENE_TABIX_COLS = TABIX_COLS + ["unique","universal"]
GENE_TABIX_TYPES = {"start" : int, "end" : int, "unique" : int, "universal" : int}
TABIX_SUFFIX = ".bgz"

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
        #self.chrs = pd.read_csv(f"{self.prefix}/chrs.csv").set_index(["genome","chr"])
        #self.genomes = self.chrs.index.unique("genome")
        self._load_chrs()

        self.bitmaps = dict()
        self.gene_tabix = dict()
        self.anno_tabix = dict()

        for g in self.genomes:
            self.bitmaps[g] = KmerBitmap(self.conf, g, "r", self.chrs)
            self.gene_tabix[g] = self._load_tabix(g, "gene")
            self.anno_tabix[g] = self._load_tabix(g, "anno")

        if not hasattr(self.conf, "chr_bin_kbp"):
            fnames = glob.glob(f"{self.prefix}/bins_*kbp.pkl")
            if len(fnames) != 1:
                raise RuntimeError(f"Exactly one chromosome bin file must be present, found {len(fnames)}: {fnames}")
            self.conf.chr_bin_kbp = int(fnames[0].split("_")[-1][:-7])

        self._load_chr_bins()

    def _load_chr_bins(self):
        arr = np.fromfile(self.chr_bins_fname, "uint32")
        arr = arr.reshape((len(arr)//self.ngenomes, self.ngenomes))

        idx = pd.MultiIndex.from_frame(pd.concat([self.chr_bin_coords(*ch).drop(columns=["end"]) for ch in self.chrs.index]))
        self.chr_bins = pd.DataFrame(arr, columns=self._total_occ_idx, index=idx)
        #self.chr_bins = pd.read_pickle(self.chr_bins_fname)

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
        bitmap = KmerBitmap(*args)
        bitmap.close()
        return bitmap.anchor_name

    def _load_chrs(self):
        self.chrs = pd.read_csv(f"{self.prefix}/chrs.csv").set_index(["genome","chr"])
        self._init_genomes()

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
        })

        if hasattr(self.conf, "gff"):
            self.genome_files["gff"] = pd.Series(self.conf.gff)

        if self.conf.anno_only:
            self._load_chrs()
        else:
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

        print("Computing chromosome summaries")
        self.chrs[self._total_occ_idx] = 0
        self.chrs[self._gene_occ_idx] = 0

        chr_bins_out = open(self.chr_bins_fname, "wb")
        #bin_coords = list()
        #chr_bin_occs = list()

        prev_genome = None
        for (genome,chrom),size in self.chrs["size"].items():
            if genome != prev_genome:
                print(genome)
                prev_genome = genome

            occs = self.query_bitmap(genome,chrom,0,size,100).sum(axis=1)
            chr_counts = self.count_occs(occs)
            self.chrs.loc[(genome,chrom), self._total_occ_idx] = chr_counts

            coords = self.chr_bin_coords(genome,chrom)
            #bin_coords.append(coords.drop(columns="end"))

            for i,c in coords.iterrows():
                st = c["start"] // self.conf.lowres_step
                en = c["end"] // self.conf.lowres_step
                bin_counts = self.count_occs(occs[st:en])
                chr_bins_out.write(bin_counts.tobytes())

        chr_bins_out.close()

        self.chrs.to_csv(f"{self.prefix}/chrs.csv")

        print("Computing gene summaries")
        if "gff" in self.genome_files.columns:
            for g in self.genomes:
                self._load_gffs(g)

        self.chrs.to_csv(f"{self.prefix}/chrs.csv")

    @property
    def chr_bins_fname(self):
        return f"{self.prefix}/bins_{self.conf.chr_bin_kbp}kbp.bin"

    def chr_bin_coords(self, genome, chrom):
        binlen = self.conf.chr_bin_kbp*1000
        size = self.chrs.loc[(genome,chrom),"size"]
        starts = np.arange(0,size,binlen)
        ends = np.clip(starts+binlen, 0, size)
        ret = pd.DataFrame({"genome" : genome, "chr" : chrom, "start" : starts, "end" : ends})
        return ret


    def query_occ_counts(self, genome, chrom, start, end, step=1):
        occs = self.query_bitmap(genome, chrom, start, end, step).sum(axis=1)
        return self.count_occs(occs)
    
    def count_occs(self, occs):
        ret = np.zeros(self.ngenomes, "uint32")
        occs, counts = np.unique(occs, return_counts=True)
        ret[occs-1] = counts
        return ret
        #return pd.Series(index=idx[occs-1], data=counts).reindex(idx, fill_value=0)

    def _iter_gff(self, fname):
        print(fname)
        for df in pd.read_csv(
            fname, 
            sep="\t", comment="#", chunksize=10000,
            names = ["chr","source","type","start","end","score","strand","phase","attr"],
            usecols = TABIX_COLS): yield df[TABIX_COLS]

    def _load_tabix(self, genome, type_):
        fname = self.tabix_fname(genome, type_)
        if not os.path.exists(fname):
            return None
        return pysam.TabixFile(fname, parser=pysam.asTuple())

    def tabix_fname(self, genome, typ):
        return os.path.join(self.anno_dir, f"{genome}.{typ}.bed{TABIX_SUFFIX}") 

    def _write_tabix(self, df, genome, typ):
        tbx = self.tabix_fname(genome, typ)
        bed = tbx[:-len(TABIX_SUFFIX)]

        df.to_csv(bed, sep="\t", header=None, index=False)
        pysam.tabix_compress(bed, tbx, True)
        pysam.tabix_index(tbx, True, 0,1,2)
        #os.remove(bed)

    def _load_gffs(self, genome):
        genes = list()
        annos = list()
        fname = self.genome_files.loc[genome, "gff"]
        if pd.isna(fname): return

        for df in self._iter_gff(fname):
            genes.append(df[df["type"].isin(self.conf.gff_gene_types)])
            annos.append(df[df["type"].isin(self.conf.gff_anno_types)])

        def _merge_dfs(dfs):
            return pd.concat(dfs).sort_values(["chr","start"]).reset_index(drop=True)

        self._write_tabix(_merge_dfs(annos), genome, "anno")

        genes = _merge_dfs(genes)
        genes["unique"] = 0
        genes["universal"] = 0
        for i,g in genes.iterrows():
            counts = self.query_occ_counts(genome, g["chr"], g["start"], g["end"])
            genes.loc[i, "unique"] += counts[0]
            genes.loc[i, "universal"] += counts[-1]
            self.chrs.loc[(genome,g["chr"]), self._gene_occ_idx] += counts
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

    def query_genes(self, genome, chrom, start, end):
        if not genome in self.gene_tabix:
            return pd.DataFrame(columns=GENE_TABIX_COLS)
        rows = self.gene_tabix[genome].fetch(chrom, start, end)
        return pd.DataFrame(rows, columns=GENE_TABIX_COLS).astype(GENE_TABIX_TYPES)

    def query_anno(self, genome, chrom, start, end):
        if not genome in self.anno_tabix:
            return pd.DataFrame(columns=TABIX_COLS)
        rows = self.anno_tabix[genome].fetch(chrom, start, end)
        return pd.DataFrame(rows, columns=TABIX_COLS).astype(TABIX_TYPES)
    
    #def _query_tabix(self, tabix genome, chrom, start, end):

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
        if hasattr(conf, "lowres_step"):
            self.steps = [1, conf.lowres_step]
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


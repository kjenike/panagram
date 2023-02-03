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

import dataclasses
from simple_parsing import field
from typing import Any, List, Tuple, Type, Union
import argparse

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

REQUIRED_FILES = ["panagram.toml", "chrs.csv", "bins_*.bin", "anchors/*bgz", "anchors/*gzi"]

MODES = {"r","w"}

#TODO
#store mash distances {prefix}/dists.csv

@dataclasses.dataclass
class KMC:
    """Parameters for KMC kmer counting"""
    memory: int = field(default=8, dest="main.kmc.memory")
    processes: int = field(default=4, dest="main.kmc.processes")
    threads: int = field(default=1, dest="main.kmc.threads")

    #Use existing KMC count and onehot genome databases if present
    use_existing: bool = field(action="store_true", default=False)

@dataclasses.dataclass
class Index:
    """Anchor KMC bitvectors to reference FASTA files to create pan-kmer bitmap"""

    #configuration file (toml)
    input: str = field(positional=True, metavar="config_file")

    prefix: str = field(alias=["-o"], default=None)

    #K-mer length
    k: int = field(alias=["-k"], default=21)

    #Number of processes
    processes: int = field(alias=["-p"], default=1)

    #Step size for low-resolution pan-kmer bitmap (in nucleotides, larger step = lower resolution)
    lowres_step: int = 100

    #Size of chromosome-scale occurence count bins in kilobases 
    chr_bin_kbp: int = 200

    gff_gene_types: List[str] = field(default_factory=lambda: ["gene"])
    gff_anno_types: List[str] = field(default=None)

    #Only perform anchoring and annotation
    anchor_only: bool = False

    #Only perform annotaion
    anno_only: bool = False

    kmc: KMC = KMC()

    fasta: dict = field(default_factory=dict,help=argparse.SUPPRESS)
    gff: dict = field(default_factory=dict,help=argparse.SUPPRESS)

    #Dummy parameters to force KMC params to be in "kmc.*" format
    use_existing: int = field(default=1,help=argparse.SUPPRESS)
    threads: int = field(default=1,help=argparse.SUPPRESS)
    memory: int = field(default=1,help=argparse.SUPPRESS)

    def _load_dict(self, root, vals):
        for key,val in vals.items():
            dest = getattr(root, key, None)
            if dataclasses.is_dataclass(dest):
                if isinstance(val, dict):
                    self._load_dict(dest, val)
                elif dataclasses.is_dataclass(val):
                    setattr(root, key, val)
                else:
                    raise ValueError(f"{key} must be dict or dataclass, found {key} = {val}")
            else:
                setattr(root, key, val)

    def run(self):
        self.write()#args)
        self.close()

    @property
    def params(self):
        return dataclasses.asdict(self)

    def __post_init__(self):

        self.write_mode = False

        if self.input is not None:
            if os.path.isdir(self.input):
                self.prefix = self.input

            elif os.path.exists(self.input):
                d = toml.load(self.input)
                self._load_dict(self, d)
                self.write_mode = True

            self.input = None

        if self.prefix is None:
            raise ValueError("Must specify index prefix or valid configuration file")

        #self.write_mode = hasattr(self, "fasta")
        #for pattern in REQUIRED_FILES:
        #    if len(glob.glob(pattern)) == 0:
        #        print(pattern)
        #        self.write_mode = True
        #        break

        self.count_dir =  self.init_dir("kmc_count")
        self.onehot_dir = self.init_dir("kmc_onehot")
        self.bitvec_dir = self.init_dir("kmc_bitvec")
        self.tmp_dir =    self.init_dir("tmp")
        self.anchor_dir = self.init_dir(ANCHOR_DIR)
        self.anno_dir = self.init_dir("anno")

        if not self.write_mode:
            self._init_read()

    def _init_read(self):
        self._load_chrs()

        self._load_dict(self, toml.load(self.index_config_file))

        self.bitmaps = dict()
        self.gene_tabix = dict()
        self.anno_tabix = dict()

        for g in self.genomes:
            self.bitmaps[g] = KmerBitmap(self.params, g, "r", self.chrs)
            self.gene_tabix[g] = self._load_tabix(g, "gene")
            self.anno_tabix[g] = self._load_tabix(g, "anno")

        if self.chr_bin_kbp is None:
            fnames = glob.glob(f"{self.prefix}/bins_*kbp.bin")
            if len(fnames) != 1:
                raise RuntimeError(f"Exactly one chromosome bin file must be present, found {len(fnames)}: {fnames}")
            self.chr_bin_kbp = int(fnames[0].split("_")[-1][:-7])

        self._load_chr_bins()

    def _load_chr_bins(self):
        arr = np.fromfile(self.chr_bins_fname, "uint32")
        arr = arr.reshape((len(arr)//self.ngenomes, self.ngenomes))

        idx = pd.MultiIndex.from_frame(pd.concat([self.chr_bin_coords(*ch).drop(columns=["end"]) for ch in self.chrs.index]))
        self.chr_bins = pd.DataFrame(arr, columns=self._total_occ_idx, index=idx)
        #self.chr_bins = pd.read_pickle(self.chr_bins_fname)

    def init_dir(self, path):
        d = os.path.join(self.prefix, path)
        if self.write_mode:
            os.makedirs(d, exist_ok=True)
        return d


    
    #def _load_kmc(self, files):
    #    from .kmc import py_kmc_api
    #    self.kmc = py_kmc_api
    #    self.kmc_dbs = list()
    #    for fname in files:
    #        self.kmc_dbs.append(self.kmc.KMCFile())
    #        self.kmc_dbs[-1].OpenForRA(fname)

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
            "fasta" : pd.Series(self.fasta),
        })

        if hasattr(self, "gff"):
            self.genome_files["gff"] = pd.Series(self.gff)

        if self.anno_only:
            self._load_chrs()
        else:
            for i,(name,fasta) in enumerate(self.genome_files["fasta"].items()):
                self.genome_ids[name] = i

                if not os.path.exists(fasta+".fai"):
                    cmd = ["samtools", "faidx", fasta]
                    subprocess.check_call(cmd)
                fa = pysam.FastaFile(fasta)

                genomes.append(pd.DataFrame(
                    [(name, chrom, i, fa.get_reference_length(chrom)-self.k+1) 
                     for chrom in fa.references], 
                    columns=["genome", "chr", "id", "size"]))

            self.chrs = pd.concat(genomes).set_index(["genome", "chr"])
            self._init_genomes()
            self.chrs.to_csv(f"{self.prefix}/chrs.csv")

        if self.anchor_only or self.anno_only:
            pre = f""
            suf = ".kmc_pre"
            bitvec_dbs = [f[:-len(suf)] 
                for f in glob.glob(f"{self.bitvec_dir}/*{suf}")]

        else:
            bitvec_dbs = self._run_kmc()
        
        def iter_anchor_args():
            for name,fasta in self.genome_files["fasta"].items():
                yield (self.params, name, "w", self.chrs, fasta, bitvec_dbs)

        if not self.anno_only:
            if self.processes == 1:
                for args in iter_anchor_args():
                    print("Anchored", self.run_anchor(args))
            else:
                with mp.Pool(processes=self.processes) as pool:
                    for name in pool.imap_unordered(self.run_anchor, iter_anchor_args(), chunksize=1):
                        print(f"Anchored {name}")
                        sys.stdout.flush()

        self.bitmaps = {
            name : KmerBitmap(self.params, name, "r", self.chrs) for name in self.genomes
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
                prev_genome = genome
                print(genome)

            occs = self.query_bitmap(genome,chrom,0,size,100).sum(axis=1)
            chr_counts = self.count_occs(occs)
            self.chrs.loc[(genome,chrom), self._total_occ_idx] = chr_counts

            coords = self.chr_bin_coords(genome,chrom)
            #bin_coords.append(coords.drop(columns="end"))

            for i,c in coords.iterrows():
                st = c["start"] // self.lowres_step
                en = c["end"] // self.lowres_step
                bin_counts = self.count_occs(occs[st:en])
                chr_bins_out.write(bin_counts.tobytes())

        chr_bins_out.close()

        self.chrs.to_csv(f"{self.prefix}/chrs.csv")

        if "gff" in self.genome_files.columns:

            if self.gff_anno_types is None:
                self.all_anno_types = pd.Index([])

            print("Computing gene summaries")
            for g in self.genomes:
                print(g)
                self._load_gffs(g)

            if self.gff_anno_types is None:
                self.gff_anno_types = list(self.all_anno_types)

        self.chrs.to_csv(f"{self.prefix}/chrs.csv")

        with open(self.index_config_file, "w") as conf_out:
            toml.dump(self.params, conf_out)

    @property
    def index_config_file(self):
        return f"{self.prefix}/panagram.toml"

    @property
    def chr_bins_fname(self):
        return f"{self.prefix}/bins_{self.chr_bin_kbp}kbp.bin"

    def chr_bin_coords(self, genome, chrom):
        binlen = self.chr_bin_kbp*1000
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
            gmask = df["type"].isin(self.gff_gene_types)
            genes.append(df[gmask])

            if self.gff_anno_types is not None:
                annos.append(df[df["type"].isin(self.gff_anno_types)])
            else:
                annos.append(df[~gmask])

        def _merge_dfs(dfs):
            return pd.concat(dfs).sort_values(["chr","start"]).reset_index(drop=True)

        annos = _merge_dfs(annos)
        self._write_tabix(annos, genome, "anno")

        if self.gff_anno_types is None:
            self.all_anno_types = self.all_anno_types.union(annos["type"].unique())

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


    def query_bitmap(self, genome, chrom, start=None, end=None, step=1):
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
        conf, db_i, i, name, fasta, count_db, onehot_db, tmp_dir = args

        onehot_id = str(2**i)

        def should_build(db):
            if not (conf["kmc"]["use_existing"] and
                    os.path.exists(db+".kmc_pre") and
                    os.path.exists(db+".kmc_suf")):
                return True
            print(f"Using exisitng KMC db: {db}")
            return False

        if should_build(count_db):
            subprocess.check_call([
                f"{KMC_DIR}/kmc", f"-k{conf['k']}", 
                f"-t{conf['kmc']['threads']}", 
                f"-m{conf['kmc']['memory']}", 
                "-ci1", "-cs10000000", "-fm",
                fasta, count_db, tmp_dir
            ])

        if should_build(onehot_db):
            subprocess.check_call([
                f"{KMC_DIR}/kmc_tools", "-t4", "transform",
                count_db, "set_counts", onehot_id, onehot_db
            ])

        return db_i, name, onehot_db


    def _iter_kmc_genome_args(self):
        i = 0
        db_i = 0
        for name,fasta in self.genome_files["fasta"].items():
            if i >= 32:
                i = 0
                db_i += 1

            count_db = os.path.join(self.count_dir, name)
            onehot_db = os.path.join(self.onehot_dir, name)
            tmp_dir = self.init_dir(f"tmp/{name}")
            yield (self.params, db_i, i, name, fasta, count_db, onehot_db, tmp_dir)

            i += 1

    def _run_kmc(self):

        i = 0
        samp_count = len(self.genome_files)
        db_count = int(np.ceil(samp_count / 32))

        genome_dbs = [list() for i in range(db_count)]
        if self.kmc.processes == 1:
            for args in self._iter_kmc_genome_args():
                i,name,db = self._run_kmc_genome(args)
                genome_dbs[i].append((name,db))
        else:
            with mp.Pool(processes=self.kmc.processes) as pool:
                for i,name,db in pool.imap(self._run_kmc_genome, self._iter_kmc_genome_args(), chunksize=1):
                    genome_dbs[i].append((name,db))

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
                for name, db in genome_dbs[i]:
                    opdefs.write(f"{name} = {db}\n")
                opdefs.write(f"OUTPUT:\n{bitvec_fname} = {genome_dbs[i][0][0]}")
                for name,_ in genome_dbs[i][1:]:
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
        self.prefix = os.path.join(conf["prefix"], ANCHOR_DIR, anchor)
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
            self._init_read()
        else:
            raise ValueError("Mode must be 'r' or 'w'")


    def _init_steps(self, conf):
        if "lowres_step" in conf:
            self.steps = [1, conf["lowres_step"]]
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

    def _init_read(self):
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

    def query(self, name, start=None, end=None, step=1):
        bstep = 1
        for s in self.steps:
            if step % s == 0:
                bstep = max(bstep, s)

        if start is None:
            start = 0

        if end is None:
            end = self.seq_len(name)

        pac = self._query(name, start, end, step, bstep)

        ret = np.unpackbits(pac, bitorder="little", axis=1)[:,:self.ngenomes]

        return ret

    def _query(self, name, start, end, step, bstep):
        byte_start = self.nbytes * (self.offsets.loc[name,bstep] + (start//bstep))
        length  = int((end - start) // bstep)

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
        try:
            from .kmc import py_kmc_api
        except ModuleNotFoundError:
            raise ModuleNotFoundError("py_kmc_api failed to install. See https://github.com/kjenike/panagram#readme for more information")

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

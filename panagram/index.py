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
from collections import defaultdict, Counter
from time import time
from Bio import bgzf, SeqIO
import yaml
import multiprocessing as mp
from types import SimpleNamespace
import shutil
import snakemake
import re
import logging

import dataclasses
from simple_parsing import field
from simple_parsing.helpers import Serializable
from typing import Any, List, Tuple, Type, Union
import argparse

logger = logging.getLogger(__name__)
def init_logger(logfile):
    logging.basicConfig(
            filename=logfile, level=logging.INFO, 
            format='[%(asctime)s %(levelname)s] %(message)s',
            datefmt="%Y-%m-%d %H:%M:%S"
    )

    #def handler(typ, val, tb):
    #    logger.exception(repr(val))

    #sys.excepthook = handler

SRC_DIR = os.path.dirname(os.path.realpath(__file__))
EXTRA_DIR = os.path.join(SRC_DIR, "extra")
SNAKEFILE = os.path.join(SRC_DIR, "workflow", "Snakefile")
NAME_REGEX = "[A-Za-z0-9_-]+"

BGZ_SUFFIX = "gz"
IDX_SUFFIX = "gzi"
ANCHOR_DIR = "anchor"

GFF_COLS = ["chr","start","end","type","attr"]
TABIX_COLS = ["chr","start","end","type","name"]
TABIX_TYPES = {"start" : int, "end" : int}
GENE_COLS = ["chr","start","end","name"]#,"unique","universal"]
GENE_TABIX_TYPES = {"start" : int, "end" : int, "unique" : int, "universal" : int}
TABIX_SUFFIX = ".gz"

MODES = {"r","w"}

@dataclasses.dataclass
class KMC(Serializable):
    """Parameters for KMC kmer counting"""
    memory: int = field(default=8, dest="main.kmc.memory")
    threads: int = field(default=1, dest="main.kmc.threads")

    #Use existing KMC count and onehot genome databases if present
    use_existing: bool = field(action="store_true", default=False)

@dataclasses.dataclass
class Index(Serializable):
    """Anchor KMC bitvectors to reference FASTA files to create pan-kmer bitmap"""

    #configuration file (yaml)
    input: str = field(positional=True, metavar="config_file")

    mode: str = field(default=None, help=argparse.SUPPRESS)

    prefix: str = field(alias=["-o"], default=None)

    #K-mer length
    k: int = field(alias=["-k"], default=21)

    #Number of processes
    cores: int = field(alias=["-c"], default=1)

    #Step size for low-resolution pan-kmer bitmap (in nucleotides, larger step = lower resolution)
    lowres_step: int = 100

    #Size of chromosome-scale occurence count bins in kilobases
    max_bin_kbp: int = 200
    min_bin_count: int = 100

    max_view_chrs: int = 50

    gff_gene_types: List[str] = field(default_factory=lambda: ["gene"], help="GFF features to store locations and conservation scores")
    gff_anno_types: List[str] = field(default=None, help="GFF features of which to store locations, but not conservation scores")
    gff_name: str = field(default="Name", help="GFF attributes to store in annotation indexes")

    #Subset of genome IDs to generate anchor genomes for. Will use all genomes as anchors if not specified
    anchor_genomes: List[str] = field(default=None)

    prepare: bool = field(alias=["-p"], default=False, action="store_true")

    kmc: KMC = field(default_factory=lambda:KMC())

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

    @property
    def params(self):
        return dataclasses.asdict(self)
    
    @property
    def bitsum_index(self):
        return pd.RangeIndex(0,self.ngenomes+1)

    def kmc_prefix(self, *names):
        return os.path.join(self.get_subdir("kmc"), ".".join(names))

    def get_subdir(self,name):
        return os.path.join(self.prefix, name)

    @property
    def snakefile(self):
        return os.path.join(self.prefix,"Snakefile")

    #panagram index command
    def run(self):
        self.init_config()
        print('Wrote config.yaml and samples.tsv')

        if not os.path.exists(self.snakefile):
            print(f'Wrote {self.snakefile}')
            shutil.copy(SNAKEFILE, self.snakefile)
        else:
            print(f'Using existing {self.snakefile}')

        args = ["--cores", f"{self.cores}", "all"]
        argstr = " ".join(args)

        if self.prepare:
            print(f"Prepared. Run 'snakemake {argstr}' to build index")
        else:
            print(f"Running 'snakemake {argstr}'")
            snakemake.main(args)

        self.close()

    def __getitem__(self, k):
        return self.genomes[k]

    def __post_init__(self):
        if not (self.mode is None or self.mode in MODES):
            raise ValueError("Invalid mode '{self.mode}', must be 'r' or 'w'")

        if self.mode is None:
            self.write_mode = os.path.isfile(self.input)
        else:
            self.write_mode = self.mode == "w"

        if self.write_mode:
            if self.prefix is None:
                self.prefix = os.path.dirname(self.input)
            if len(self.prefix) == 0:
                self.prefix = "."

            if os.path.isdir(self.input):
                self.prefix = self.input
                if not (os.path.isfile(self.config_fname) and os.path.isfile(self.samples_fname)):
                    raise ValueError("Index write directory not initialized")
                self.input = self.samples_fname
                self.load_config()

            elif os.path.isfile(self.input):
                self.init_config()

            else:
                raise ValueError("Index input must be sample TSV or initialized directory")

        else:
            if not os.path.isdir(self.input):
                raise ValueError("Index input must be directory mode='r'")
            self.prefix = self.input

        os.chdir(self.prefix)
        self.prefix = ""

        self.samples = pd.read_table(self.samples_fname).set_index("name")

        self.genomes = dict()
        p = self.params
        for name,row in self.samples.iterrows():
            self.genomes[name] = Genome(self, row["id"], name, row["fasta"], row["gff"], row["anchor"], write=self.write_mode)

        self.chrs = None

        if not self.write_mode:
            self._init_read()

    @property
    def genome_names(self):
        return self.samples.index

    def init_config(self):

        samples = pd.read_table(self.input)#[["name","fasta","gff"]].set_index("name")
        if not "name" in samples.columns or not "fasta" in samples.columns:
            raise ValueError("Input samples must contain 'name' and 'fasta' column headers")
        if not "gff" in samples:
            samples["gff"] = pd.NA

        invalid = ~samples["name"].str.fullmatch(NAME_REGEX)
        if np.any(invalid):
            bad = "', '".join(samples["name"][invalid])
            raise ValueError(f"Invalid genome names: '{bad}'\nMust match r'{NAME_REGEX}'.")

        samples = samples[["name","fasta","gff"]].set_index("name").dropna(how="all")
        samples["id"] = np.arange(len(samples), dtype=int)

        if self.anchor_genomes is None:
            if "anchor" in samples:
                self.anchor_genomes = list(samples.query("anchor").index)
            else:
                self.anchor_genomes = list(samples["fasta"].dropna().index)

        samples["anchor"] = samples.index.isin(self.anchor_genomes)

        samples.to_csv(self.samples_fname, sep="\t")

        self.write_config()


    def _init_read(self):
        self.load_config()

        self.ngenomes = len(self.samples)

        self.chrs = pd.concat({
            genome : self.genomes[genome].chrs
            for genome in self.anchor_genomes
        }, names=["genome","chr"])#.sort_index()

        self.bitsum_bins = pd.concat({
            genome : self.genomes[genome].bitsum_bins #_read_bitsum_bins().droplevel("end")
            for genome in self.anchor_genomes
        }, names=["genome","chr","start"]).sort_index()

        self.bitsum_chrs = pd.concat({
            genome : self.genomes[genome].bitsum_chrs #_read_bitsum_bins().droplevel("end")
            for genome in self.anchor_genomes
        }, axis=0)

        self.bitfreq_chrs = pd.concat({
            genome : self.genomes[genome].bitfreq_chrs #_read_bitsum_bins().droplevel("end")
            for genome in self.anchor_genomes
        }, axis=0)
        self.bitsum_totals = pd.concat({
            genome : self.genomes[genome].bitsum_total #_read_bitsum_bins().droplevel("end")
            for genome in self.anchor_genomes
        }, axis=1).T
        self.bitfreq_totals = self.bitsum_totals.divide(self.bitsum_totals.sum(axis=1),axis=0)

        n = np.array(self.bitfreq_totals.columns)
        self.bitsum_totals_avg = (self.bitfreq_totals*n).sum(axis=1).sort_values()
        self.bitsum_chrs_avg = (self.bitfreq_chrs*n).sum(axis=1).sort_values()

        g = self.chrs["size"].groupby("genome")
        self.genome_sizes = pd.DataFrame({
            "length" : g.sum(),
            "chr_count" : g.count()
        })

    def __getitem__(self, genome):
        return self.genomes[genome]

    def write_config(self, exclude=["prefix"]):
        prms = self.params
        for p in exclude:
            del prms[p]

        with open(self.config_fname, "w") as conf_out:
            yaml.dump(prms, conf_out)

    def load_config(self):
        with open(self.config_fname) as f:
            self._load_dict(self, yaml.load(f,yaml.Loader))

    @property
    def config_fname(self):
        return os.path.join(self.prefix, "config.yaml")

    @property
    def samples_fname(self):
        return os.path.join(self.prefix, "samples.tsv")

    @property
    def genome_dist_fname(self):
        return os.path.join(self.prefix, "genome_dist.tsv")

    def bitsum_count(self, occs):
        ret = np.zeros(self.ngenomes, "uint32")
        occs, counts = np.unique(occs, return_counts=True)
        ret[occs-1] = counts
        return ret
        #return pd.Series(index=idx[occs-1], data=counts).reindex(idx, fill_value=0)

    def close(self):
        for b in self.genomes.values():
            b.close()

    def query_bitmap(self, genome, chrom, start=None, end=None, step=1):
        return self.genomes[genome].query(chrom, start, end, step)

    def query_genes(self, genome, chrom=None, start=None, end=None):
        return self.genomes[genome].query_genes(chrom, start, end)

    def query_anno(self, genome, chrom, start, end):
        return self.genomes[genome].query_anno(chrom, start, end)

    @property
    def kmc_bitvec_count(self):
        return int(np.ceil(len(self.samples) / 32.0))

    @property
    def opdef_filenames(self):
        return [self.kmc_prefix(f"opdef{i}.txt") for i in range(self.kmc_bitvec_count)]

    @property
    def bitvec_prefixes(self):
        return [self.kmc_prefix(f"bitvec{i}") for i in range(self.kmc_bitvec_count)]

    @property
    def steps(self):
        return (1, self.lowres_step)

    def init_opdefs(self):
        genome_dbs = [list()]
        i = 0
        for name,fasta in self.samples["fasta"].items():
            if i == 32:
                genome_dbs.append(list())
                i = 0
            genome_dbs[-1].append((name,self.kmc_prefix(name,"onehot")))
            i += 1


        for i,fname in enumerate(self.opdef_filenames):
            bitvec_fname = self.kmc_prefix(f"bitvec{i}")
            with open(fname, "w") as opdefs:
                opdefs.write("INPUT:\n")
                for name, db in genome_dbs[i]:
                    opdefs.write(f"{name} = {db}\n")
                opdefs.write(f"OUTPUT:\n{bitvec_fname} = {genome_dbs[i][0][0]}")
                for name,_ in genome_dbs[i][1:]:
                    opdefs.write(f" + {name}")
                opdefs.write("\n-ocsum\n")

    def get_anchor_filenames(self, name):
        return self.genomes[name].anchor_filenames

    @property
    def anchor_filenames(self):
        ret = list()
        for g in self.genomes.values():
            ret += g.anchor_filenames
        return ret
    
    def bitmap_to_bins(self, bitmap, binlen):
        df = bitmap.set_index(bitmap.index // binlen)#.groupby(level=1)

        pancounts = df.sum(axis=1).reset_index().value_counts()
        pancount_bins = pancounts.unstack(level=1,fill_value=0).T
        pancount_bins = pancount_bins.reindex(self.bitsum_index, fill_value=0)

        paircount_bins = df.groupby(level=0).sum()
        paircount_bins = paircount_bins.set_index(paircount_bins.index*binlen).T
        paircount_bins = paircount_bins.div(paircount_bins.max(axis=0),axis=1)

        return pancount_bins, paircount_bins
    
    def bitmap_to_pancount(self, bitmap):
        return pd.Series(bitmap.to_numpy().sum(axis=1),index=bitmap.index)
    
    def bitmap_to_paircount(self, bitmap):
        return pd.Series(bitmap.to_numpy().sum(axis=1),index=bitmap.index)
    
    def pancount_to_bins(self, pancnts, binlen):
        bin_counts = pd.DataFrame({
            "count" : pancnts.to_numpy(),
            "bin" : pancnts.index // binlen
        }).value_counts(sort=False)
        return bin_counts.unstack(level=1,fill_value=0).reindex(self.bitsum_index, fill_value=0)
    

class Genome:
    def __init__(self, idx, id, name, fasta=None, gff=None, anchor=None, write=False):
        self.index = idx
        self.samples = idx.samples
        self.params = idx.params
        self.prefix = os.path.join(self.params["prefix"], ANCHOR_DIR, name)
        self.id = id#chrs.loc[name]["id"].iloc[0]
        self.name = name
        self.fasta = fasta
        self.gff = gff
        self.write_mode = write
        self.anchored = anchor if (anchor is not None) else (fasta is not None)
        self.annotated = not pd.isna(gff)

        self.genome_names = idx.genome_names
        self.ngenomes = len(self.samples)
        self.nbytes = int(np.ceil(self.ngenomes / 8))
        self.bitmaps = None
        self.chrs = None

        if not self.anchored:
            return

        self._init_steps()

        self.seq_lens = defaultdict(dict)
        self.bitmap_lens = defaultdict(int)

        self._init_anno_types()

        if os.path.exists(self.chrs_fname):
            self.load_chrs()
        elif self.fasta is not None:
            self.init_chrs(self.fasta)
        else:
            sys.stderr.write(f"Warning: failed to initialze '{self.name}' chromosomes")
            self.chrs = None

        if not self.write_mode:
            self.init_read()

    @property
    def anchor_filenames(self):
        if not self.anchored:
            return []
        ret = [self.chrs_fname, self.bins_fname]
        for s in self.steps:
            ret += [self.bitmap_gz_fname(s), self.bitmap_gzi_fname(s)]

        if self.annotated:
            ret.append(self.chr_genes_fname)
            for t in ["gene","anno"]:
                ret += [self.tabix_fname(t), self.tabix_idx_fname(t)]
        return ret

    @property
    def chrs_fname(self):
        return os.path.join(self.prefix, "chrs.tsv")

    @property
    def bins_fname(self):
        return os.path.join(self.prefix, f"bitsum.bins.tsv")

    @property
    def chr_genes_fname(self):
        return os.path.join(self.prefix, f"bitsum.genes.tsv")

    @property
    def anno_types_fname(self):
        return os.path.join(self.prefix, "anno_types.txt")

    def bitmap_gz_fname(self, step):
        return os.path.join(self.prefix, f"bitmap.{step}.{BGZ_SUFFIX}")

    def bitmap_gzi_fname(self, step):
        return os.path.join(self.prefix, f"bitmap.{step}.{IDX_SUFFIX}")

    def bed_tmp_fname(self, typ):
        return os.path.join(self.prefix, f"{typ}.bed")

    def tabix_fname(self, typ):
        return self.bed_tmp_fname(typ)+".gz" #os.path.join(self.prefix, f"{typ}.bed{TABIX_SUFFIX}")

    def tabix_idx_fname(self, typ):
        return self.tabix_fname(typ)+".csi"
    
    @property
    def bitsum_index(self):
        return pd.RangeIndex(0,self.ngenomes+1)

    @property
    def gene_tabix_cols(self):
        return GENE_COLS + [1,self.ngenomes]#list(self.bitsum_index)
        #return TABIX_COLS + list(self.bitsum_index)

    @property
    def gene_tabix_types(self):
        r = {"start" : int, "end" : int}
        for i in [1,self.ngenomes]:#self.bitsum_index:
            r[i] = int
        return r#TABIX_COLS + list(self.bitsum_index)
    
    @property
    def chr_count(self):
        return len(self.chrs)

    def init_chrs(self, fasta):
        fa = pysam.FastaFile(fasta)
        chrs = pd.DataFrame(
            [(i, name, fa.get_reference_length(name)-self.params['k']+1)
             for i,name in enumerate(fa.references)],
            columns=["id", "name", "size"]).set_index("name")
        self.set_chrs(chrs)

        return chrs

    def write_chrs(self):
        self.chrs.to_csv(self.chrs_fname, sep="\t")

    def load_chrs(self):
        self.set_chrs(pd.read_table(self.chrs_fname,index_col="name"))

    def set_chrs(self, chrs):
        self.chrs = chrs
        if not "gene_count" in self.chrs.columns:
            self.chrs["gene_count"] = 0
        self.sizes = chrs["size"]

        step_sizes = pd.DataFrame({step : np.ceil(self.sizes / step) for step in self.steps}, dtype=int)
        self.offsets = step_sizes.cumsum().shift(fill_value=0)

    def _init_steps(self):
        if "lowres_step" in self.params:
            self.steps = [1, self.params["lowres_step"]]
        else:
            self.steps = list()
            for fname in glob.glob(f"{self.prefix}.*.{BGZ_SUFFIX}"):
                step = int(fname.split(".")[-2])
                self.steps.append(step)

    def init_read(self):
        self.blocks = {s : self.load_bgz_blocks(self.bitmap_gzi_fname(s)) for s in self.steps}
        self.bitmaps = {s : bgzf.BgzfReader(self.bitmap_gz_fname(s), "rb") for s in self.steps}

        self.bitsum_bins = self._read_bitsum_bins()
        self.bitsum_chrs = self.bitsum_bins.groupby("chr").sum()
        self.bitsum_total = self.bitsum_bins.sum()

        sum2freq = lambda df: df.divide(df.sum(axis=1), axis=0)
        self.bitfreq_bins = sum2freq(self.bitsum_bins)
        self.bitfreq_chrs = sum2freq(self.bitsum_chrs)

        self.gene_tabix = self._load_tabix("gene")
        self.anno_tabix = self._load_tabix("anno")
        self.annotated = self.gene_tabix is not None or self.anno_tabix is not None
        
        if self.annotated:
            self.bitsum_genes = pd.read_table(self.chr_genes_fname).set_index("chr")
            self.bitfreq_genes = sum2freq(self.bitsum_genes)
        else:
            self.bitfreq_genes = self.bitsum_genes = pd.DataFrame(0, index=self.chrs.index, columns=self.gene_tabix_cols)

    def _load_tabix(self, type_):
        fname = self.tabix_fname(type_)
        if not os.path.exists(fname):
            return None

        index_fname = self.tabix_idx_fname(type_)
        return pysam.TabixFile(fname, parser=pysam.asTuple(), index=index_fname)


    def _read_bitsum_bins(self):
        df = pd.read_table(self.bins_fname)
        df["chr"] = self.chrs.index[df["chr"]]
        df.set_index(["chr","start"],inplace=True)
        df.columns = df.columns.astype(int)
        return df

    def seq_len(self, seq_name):
        return self.sizes.loc[seq_name]
    
    def _gffattr(self,df,name):
        return df["attr"].str.extract(f"{name}=([^;]+)", re.IGNORECASE)[0]

    def _iter_gff(self):
        for df in pd.read_csv(
                self.gff,
                sep="\t", comment="#", chunksize=10000,
                names = ["chr","source","type","start","end","score","strand","phase","attr"],
                usecols = GFF_COLS):
            
            df["id"] = self._gffattr(df, "ID")

            #df["id"] = gffattr(df, "ID")
            #df["name"] = pd.NA
            #for attr in self.params["gff_name"]:
            #    isna = df.index[df["name"].isna()]
            #    if len(isna) > 0:
            #        names = self._gffattr(df.loc[isna], attr)
            #        df.loc[isna,"name"] = names
            #    else:
            #        break

            yield df#[TABIX_COLS]

    def _init_anno_types(self):
        if self.params["gff_anno_types"] is not None:
            anno_types = [self.params["gff_anno_types"]]

        if os.path.exists(self.anno_types_fname) and not self.write_mode:
            with open(self.anno_types_fname) as f:
                anno_types = [l.strip() for l in f]
        else:
            self.gff_anno_types = None
            self.anno_type_ids = None
            return
        
        if "exon" in anno_types:
            if anno_types[0] != "exon":
                anno_types = ["exon"] + [a for a in anno_types if a != "exon"]
            id0 = 0
        else:
            id0 = 1

        self.gff_anno_types = set(anno_types)
        self.anno_type_ids = pd.Series({a : id0+i for i,a in enumerate(anno_types)})
        

    def _write_anno_types(self):
        with open(self.anno_types_fname, "w") as f:
            for t in self.gff_anno_types:
                f.write(f"{t}\n")

    def init_gff(self, filename=None):
        if filename is None:
            filename = self.gff
        if pd.isna(filename): return

        if self.params["gff_anno_types"] is not None:
            gff_anno_types = set(self.params["gff_anno_types"])
        else:
            gff_anno_types = None

        genes = list()
        annos = list()

        for df in self._iter_gff():
            gmask = df["type"].isin(self.params["gff_gene_types"])
            genes.append(df[gmask])

            if gff_anno_types is not None:
                annos.append(df[df["type"].isin(gff_anno_types)])
            else:
                annos.append(df[~gmask])

        def _merge_dfs(dfs):
            return pd.concat(dfs).sort_values(["chr","start"]).reset_index(drop=True)

        annos = _merge_dfs(annos)#.set_index("type").sort_index()
        genes = _merge_dfs(genes)
        genes["name"] = self._gffattr(genes, self.params["gff_name"]).fillna(genes["id"])

        parents = self._gffattr(annos, "Parent")#.set_index(annos["id"]).dropna()
        anno_ids = annos.reset_index().dropna().set_index("id")["index"]
        gene_names = genes[["id","name"]].dropna().set_index("id")["name"]

        p = parents.isin(anno_ids.index)

        n = 0
        while p.any():
            parents[p] = parents.loc[anno_ids.loc[parents[p]]].to_numpy()
            p = parents.isin(anno_ids.index)
            n += 1

        annos["name"] = gene_names.loc[parents].to_numpy()

        annos = annos[annos["type"] != "transcript"][TABIX_COLS].drop_duplicates()

        self._write_tabix(annos, "anno")

        if gff_anno_types is None:
            self.gff_anno_types = set(annos["type"].unique())
        else:
            self.gff_anno_types = gff_anno_types.intersection(annos["type"])
        self._write_anno_types()

        for i in self.bitsum_index:
            genes[i] = 0
        return genes.set_index(["chr","start","end"]).sort_index()

    def _write_tabix(self, df, typ):
        tbx = self.tabix_fname(typ)
        bed = self.bed_tmp_fname(typ) #tbx[:-len(TABIX_SUFFIX)]

        df.to_csv(bed, sep="\t", header=None, index=False)
        pysam.tabix_compress(bed, tbx, True)
        pysam.tabix_index(tbx, True, 0,1,2, csi=True)

    def load_bgz_blocks(self, fname):
        with open(fname, "rb") as idx_in:
            nblocks = np.fromfile(idx_in, "uint64", 1)[0]
            dtype = [("rstart", "uint64"), ("dstart", "uint64")]
            blocks = np.zeros(int(nblocks)+1, dtype=dtype)
            blocks[1:] = np.fromfile(idx_in, dtype, nblocks)
        return blocks.astype([("rstart", int), ("dstart", int)])

    def bitsum_count(self, occs):
        return pd.Series(occs).value_counts()

    def query(self, name, start=None, end=None, step=1):
        t = time()
        bstep = 1
        for s in self.steps:
            if step % s == 0:
                bstep = max(bstep, s)

        if start is None:
            start = 0

        if end is None:
            end = self.seq_len(name)

        pac = self._query_bytes(name, start, end-1, step, bstep)
        bits = self._bytes_to_bits(pac)
        #idx = np.arange(start, end, step, dtype=int)#[:len(bits)-1]
        idx = pd.RangeIndex(start, end, step)
        df = pd.DataFrame(bits, index=idx, columns=self.genome_names)
        return df

    def _bytes_to_bits(self, pac):
        return np.unpackbits(pac, bitorder="little", axis=1)[:,:self.ngenomes]

    def _query_bytes(self, name, start, end, step, bstep):
        byte_start = self.nbytes * (self.offsets.loc[name,bstep] + (start//bstep))
        length  = int((end - start) // bstep) + 1#bstep

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


    def _load_kmc(self, kmc_dbs):
        try:
            from .extra import py_kmc_api
        except ModuleNotFoundError:
            raise ModuleNotFoundError("py_kmc_api failed to install. See https://github.com/kjenike/panagram#readme for more information")

        self.kmc = py_kmc_api
        dbs = list()
        for db in kmc_dbs:
            if isinstance(db, str):
                dbs.append(self.kmc.KMCFile())
                dbs[-1].OpenForRA(db)
            else:
                dbs.append(db)
        return dbs

    def query_genes(self, chrom=None, start=None, end=None):#, attrs=["name","id"]):
        t = time()
        if self.gene_tabix is None:
            #return pd.DataFrame(columns=list(self.gene_tabix_cols) + ["attr"])
            rows = []
        else:
            try:
                rows = self.gene_tabix.fetch(chrom, start, end, multiple_iterators=True)
            except ValueError:
                rows = []

        rows = list(rows)
        ret = pd.DataFrame(rows, columns=self.gene_tabix_cols).astype(self.gene_tabix_types)

        #if "attr" in ret.columns:
        #    attr = lambda a: ret["attr"].str.extract(f"{a}=([^;]+)", re.IGNORECASE)
        #    names = attr("Name")
        #    ids = attr("ID")
        #    names[names.isna()] = ids[names.isna()]
        #    ret["name"] = names
        #else:
        #    ret["name"] = ""

        return ret
    
    def query_anno(self, chrom, start, end):
        t0 = time()
        if self.anno_tabix is None:
            return pd.DataFrame(columns=TABIX_COLS)
        try:
            rows = self.anno_tabix.fetch(chrom, start, end)
        except ValueError:
            rows = []
        #t1 = time()
        #print("fetch",t1-t0)
        #t0 = t1

        rows = list(rows)
        #t1 = time()
        #print("list", t1-t0)
        #t0 = t1

        df =  pd.DataFrame(rows, columns=TABIX_COLS)#.drop_duplicates()

        #t1 = time()
        #print("df", t1-t0)
        #t0 = t1

        df = df.astype(TABIX_TYPES)
        df["type_id"] = self.anno_type_ids[df["type"]].to_numpy()

        #t1 = time()
        #print("type", t1-t0)
        #t0 = t1

        return df

    def iter_fasta(self):
        if self.fasta.endswith(".gz") or self.fasta.endswith(".bgz"):
            opn = lambda f: gzip.open(f, "rt")
        else:
            opn = lambda f: open(f, "r")

        with opn(self.fasta) as fasta:
            for rec in SeqIO.parse(fasta, "fasta"):
                yield rec

    def _query_kmc_bytes(self, db_i, seq):
        db = self.kmc_dbs[db_i]
        vec = self.kmc.CountVec()
        db.GetCountersForRead(seq, vec)

        pac32 = np.array(vec, dtype="uint32")
        pac8 = pac32.view("uint8").reshape((len(pac32),4))

        if self.nbytes <= 4:
            n = self.nbytes
        elif db_i == len(self.kmc_dbs) and self.nbytes % 4 > 0:
            n = self.nbytes % 4
        else:
            n = 4

        return pac8[:,:n]

    def _write_bitmap(self, name, seq):
        byte_arrs = defaultdict(list)

        for db in range(len(self.kmc_dbs)):
            pacbytes = self._query_kmc_bytes(db, seq)
            for s in self.steps:
                a = pacbytes[::s]
                byte_arrs[s].append(a)

        size = None
        arrs = dict()
        for step,arr in byte_arrs.items():
            arrs[step] = np.concatenate(arr, axis=1)
            self.bitmaps[step].write(arrs[step].tobytes())
            self.bitmap_lens[step] += len(arrs[step])
            if step == 1:
                size = len(arr)

        self.seq_lens[self.id][name] = size

        return self._bytes_to_bits(arrs[1])

    def run_annotate(self, gff_file=None, logfile=None, nogene=False):
        logging.basicConfig(
                filename=logfile, level=logging.INFO, 
                format='[ %(asctime)s %(levelname)7s ] %(message)s',
                datefmt="%Y-%m-%d %H:%M:%S"
        )

        gene_df = self.init_gff(gff_file)

        if nogene: return

        #for chrom,df in gene_df.groupby(level="chr"):
        for chrom in gene_df.index.unique("chr"):
            if not chrom in self.sizes.index:
                logger.warning(f"Skipping gene at {chrom}:{start}-{end}, chromosome not found")
                continue
            df = gene_df.loc[chrom]
            st = df.index.get_level_values("start").min()
            en = min(self.sizes[chrom],df.index.get_level_values("end").max())

            bitsum = self.query(chrom, st, en).sum(axis=1)

            for start,end in df.index:
                if end <= start or start < 0 or end-st > len(bitsum):
                    logger.warning(f"Skipping gene at {chrom}:{start}-{end}, coordinates out-of-bounds")
                    continue
                occ, counts = np.unique(bitsum[start-st:end-st], return_counts=True)
                gene_df.loc[(chrom,start,end),occ] += counts

            self.bitsum_genes = gene_df.groupby("chr",sort=False)[self.bitsum_index].sum()#.sort_index()
            self.bitsum_genes.to_csv(self.chr_genes_fname, sep="\t")

            gene_tabix = gene_df.reset_index()[self.gene_tabix_cols]
            self._write_tabix(gene_tabix, "gene")



    def run_anchor(self, bitvecs, logfile=None):
        logging.basicConfig(
                filename=logfile, level=logging.INFO, 
                format='[ %(asctime)s %(levelname)7s ] %(message)s',
                datefmt="%Y-%m-%d %H:%M:%S"
        )

        if not self.anchored:
            logger.info(f"Skipping non-anchor genome '{self.name}'")
            return

        self.kmc_dbs = self._load_kmc(bitvecs)
        logger.info(f"KMC Database Loaded")

        if self.annotated:
            gene_df = self.init_gff()#.groupby("chr")
            chr_genes = gene_df.index.get_level_values(0).value_counts()
            logger.info(f"Annotation pre-processed")
        else:
            chr_genes = pd.Series([0])
        self.chrs["gene_count"] = chr_genes.reindex(self.chrs.index,fill_value=0)

        self.bitmaps = {s : bgzip.BGZipWriter(open(self.bitmap_gz_fname(s), "wb"))for s in self.steps}
        bin_occs = dict()
        logger.info(f"Anchoring Started")

        for i,rec in enumerate(self.iter_fasta()):
            chrom = rec.id
            bitmap = self._write_bitmap(chrom, str(rec.seq))

            bitsum = bitmap.sum(axis=1)
            bin_occs[i] = self.bin_bitsum(bitsum) 

            logger.info(f"Anchored {chrom}")

            if self.annotated and chrom in chr_genes:
                for _,start,end in gene_df.loc[[chrom]].index:
                    if end <= start or start < 0 or end > len(bitsum):
                        logger.warning(f"Skipping gene at {chrom}:{start}-{end}, coordinates out-of-bounds")
                        continue
                    occ, counts = np.unique(bitsum[start:end], return_counts=True)
                    gene_df.loc[(chrom,start,end),occ] += counts
                logger.info(f"Annotated {chrom}")

            t = time()

        if self.annotated:
            gene_tabix = gene_df.reset_index()[self.gene_tabix_cols]
            self._write_tabix(gene_tabix, "gene")
            self.bitsum_genes = gene_df.groupby("chr",sort=False)[self.bitsum_index].sum()#.sort_index()
            self.bitsum_genes.to_csv(self.chr_genes_fname, sep="\t")

        bin_occs = pd.concat(bin_occs,names=["chr","start","end"]).droplevel("end")
        bin_occs.to_csv(self.bins_fname, sep="\t")

        self.write_chrs()

        self.close()

        for step in self.steps:
            subprocess.check_call([
                "bgzip", "-rI", self.bitmap_gzi_fname(step), self.bitmap_gz_fname(step)])

    def bin_bitsum(self, bitsum):
        binlen = self.params["max_bin_kbp"]*1000
        if len(bitsum) / binlen < self.params["min_bin_count"]:
            binlen = len(bitsum) // self.params["min_bin_count"]

        starts = np.arange(0,len(bitsum),binlen)
        ends = np.clip(starts+binlen, 0, len(bitsum))
        coords = pd.MultiIndex.from_arrays([starts, ends])
        return pd.DataFrame([
            pd.value_counts(bitsum[st:en]).reindex(self.bitsum_index, fill_value=0) for st,en in coords
        ], index=coords).astype(int)

    def close(self):
        if self.bitmaps is not None:
            for f in self.bitmaps.values():
                f.close()

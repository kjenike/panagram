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
import yaml
import multiprocessing as mp
from types import SimpleNamespace
import shutil
import snakemake

import dataclasses
from simple_parsing import field
from simple_parsing.helpers import Serializable
from typing import Any, List, Tuple, Type, Union
import argparse

SRC_DIR = os.path.dirname(os.path.realpath(__file__))
EXTRA_DIR = os.path.join(SRC_DIR, "extra")
SNAKEFILE = os.path.join(SRC_DIR, "workflow", "Snakefile")

BGZ_SUFFIX = "bgz"
IDX_SUFFIX = "gzi"
ANCHOR_DIR = "anchor"

TABIX_COLS = ["chr","start","end","type","attr"]
TABIX_TYPES = {"start" : int, "end" : int}
GENE_TABIX_COLS = TABIX_COLS + ["unique","universal"]
GENE_TABIX_TYPES = {"start" : int, "end" : int, "unique" : int, "universal" : int}
TABIX_SUFFIX = ".bgz"

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
    chr_bin_kbp: int = 200

    gff_gene_types: List[str] = field(default_factory=lambda: ["gene"])
    gff_anno_types: List[str] = field(default=None)

    #Subset of genome IDs to generate anchor genomes for. Will use all genomes as anchors if not specified
    anchor_genomes: List[str] = field(default=None)

    prepare: bool = field(alias=["-p"], default=False, action="store_true")

    #Only perform anchoring and annotation
    anchor_only: bool = field(default=False, action="store_true")

    #Only perform annotaion
    anno_only: bool = field(default=False, action="store_true")

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

    def kmc_prefix(self, *names):
        return os.path.join(self.get_subdir("kmc"), ".".join(names))

    def get_subdir(self,name):
        return os.path.join(self.prefix, name)

    @property
    def tmp_dir(self):
        return self.get_subdir("tmp")

    @property
    def anchor_dir(self):
        return self.get_subdir("anchor")

    @property
    def anno_dir(self):
        return self.get_subdir("anno")

    @property
    def mash_dir(self):
        return self.get_subdir("mash")

    @property
    def chrs_file(self):
        return os.path.join(self.prefix, "chrs.tsv")

    @property
    def snakefile(self):
        return os.path.join(self.prefix,"Snakefile")

    #panagram index command
    def run(self):
        #self.write()
        #self.close()
        #return
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
    
    def __post_init__(self):

        if not (self.mode is None or self.mode in MODES):
            raise ValueError("Invalid mode '{self.mode}', must be 'r' or 'w'")

        if self.mode is None:
            self.write_mode = os.path.isfile(self.input)
        else:
            self.write_mode = self.mode == "w"

        if self.write_mode:
            if os.path.isdir(self.input):
                self.prefix = self.input
                if not (os.path.isfile(self.index_config_file) and os.path.isfile(self.index_samples_file)):
                    raise ValueError("Index write directory not initialized")

            elif not os.path.isfile(self.input):
                raise ValueError("Index input must be sample CSV/TSV or initialized directory")
            
            else:
                if self.prefix is None:
                    self.prefix = os.path.dirname(self.input)
                if len(self.prefix) == 0:
                    self.prefix = "."
                samples = pd.read_table(self.input)[["name","fasta","gff"]]
                samples.index.name = "id"
                samples.to_csv(self.index_samples_file, sep="\t")

            os.chdir(self.prefix)
            self.prefix = ""

            self.write_config()

        else:
            if not os.path.isdir(self.input):
                raise ValueError("Index input must be directory mode='r'")
            self.prefix = self.input
        self.samples = pd.read_table(self.index_samples_file).set_index("name")

        self.genomes = dict()
        p = self.params
        for name,row in self.samples.iterrows():
            self.genomes[name] = Genome(self, row["id"], name, row["fasta"], row["gff"], write=self.write_mode)

        self.chrs = None

        if not self.write_mode:
            self._init_read()


    def _init_read(self):
        self._load_chrs()

        with open(self.index_config_file) as f:
            self._load_dict(self, yaml.load(f,yaml.Loader))

        self.gene_tabix = dict()
        self.anno_tabix = dict()

        for g in self.anchor_genomes:
            self.genomes[g].init_read(self.chrs)# = Genome(self.params, g, "r", self.chrs)
            self.gene_tabix[g] = self._load_tabix(g, "gene")
            self.anno_tabix[g] = self._load_tabix(g, "anno")

        self.chr_bins = pd.concat({
            genome : self.genomes[genome].get_bin_counts().droplevel("end")
            for genome in self.anchor_genomes
        }, names=["genome","chr","start"]).sort_index()

    def init_dir(self, path):
        d = os.path.join(self.prefix, path)
        #if self.write_mode:
        #    os.makedirs(d, exist_ok=True)
        return d

    def __getitem__(self, genome):
        return self.genomes[genome]

    @staticmethod
    def run_anchor(args):
        print(args)
        bitmap = Genome(*args)
        bitmap.init_write(["kmc/bitvec0"])
        bitmap.close()
        return bitmap.name

    #TODO
    def _load_chrs(self):
        self.chrs = pd.read_csv(self.chrs_file, sep="\t").set_index(["genome","chr"])
        names = self.chrs.columns.str

        total_cols = names.startswith("total_occ_")
        gene_cols = names.startswith("gene_occ_")
        occ_cols = total_cols | gene_cols
        self.chr_occ = self.chrs.loc[:,total_cols | gene_cols].copy()
        if np.any(self.chr_occ):
            cols = self.chr_occ.columns.str.split("_occ_", expand=True)
            self.chr_occ.columns = cols.set_levels(cols.levels[1].astype(int), level=1)

            self.genome_occ = self.chr_occ.groupby(level="genome", sort=False).sum()
            normalize = lambda df: df.divide(df.sum(axis=1), axis=0)
            
            self.genome_occ_freq = pd.concat({
                "total" : normalize(self.genome_occ["total"]),
                "gene" : normalize(self.genome_occ["gene"]),
            }, axis=1)
            
            self.chr_occ_freq = pd.concat({
                "total" : normalize(self.chr_occ["total"]),
                "gene" : normalize(self.chr_occ["gene"]),
            }, axis=1)

            self._init_genomes()

            self.genome_occ_avg = (self.genome_occ_freq["total"]*self._occ_idx).sum(axis=1).sort_values()
            self.chr_occ_avg = (self.chr_occ_freq["total"]*self._occ_idx).sum(axis=1).sort_values()
        self.chrs = self.chrs[self.chrs.columns[~(total_cols | gene_cols)]]

    @property
    def genome_occs():
        return self.genome_occs.divide(self.genome_occs.sum(axis=1))

    def _init_genomes(self):
        if self.anchor_genomes is None:
            self.anchor_genomes = self.chrs.query("size > 0").index.unique("genome")
        self.ngenomes = len(self.samples)

        g = self.chrs["size"].groupby("genome")
        self.genome_sizes = pd.DataFrame({
            "length" : g.sum(),
            "chr_count" : g.count()
        })

        self._occ_idx = pd.RangeIndex(1, self.ngenomes+1)
        self._total_occ_idx = pd.MultiIndex.from_product([["total"], self._occ_idx]) 
        self._gene_occ_idx = pd.MultiIndex.from_product([["gene"], self._occ_idx])  

    def _run_mash(self):
        mash_files = list()
        for i,(name,fasta) in enumerate(self.samples["fasta"].dropna().items()):
            cmd =[f"{EXTRA_DIR}/mash", "sketch", "-C", name, "-o", f"{self.mash_dir}/{name}", "-r", "-s", "10000", fasta]
            subprocess.check_call(cmd)
            mash_files.append(f"{self.mash_dir}/{name}.msh")

        cmd = [f"{EXTRA_DIR}/mash", "triangle", "-C", "-E"] + mash_files
        triangle_fname = f"{self.mash_dir}/triangle.txt"
        with open(self.genome_dist_fname, "w") as mash_out:
            subprocess.check_call(cmd, stdout=mash_out)

    @property
    def genome_dist_fname(self):
        return os.path.join(self.prefix, "genome_dist.tsv")

    def write(self):
        if self.anno_only:
            self._load_chrs()
        else:
            self.init_chrs()

        if self.anchor_only or self.anno_only:
            pre = f""
            suf = ".kmc_pre"
            #bitvec_dbs = [f"{self.bitvec_dir}/{i}" for i in range(self.kmc_bitvec_count)]
            bitvec_dbs = [self.kmc_prefix(f"{i}") for i in range(self.kmc_bitvec_count)]
            #[f[:-len(suf)] 
            #    for f in glob.glob(f"{self.bitvec_dir}/*{suf}")]

        else:
            bitvec_dbs = self._run_kmc()
        
        def iter_anchor_args():
            #for name,fasta in self.samples["fasta"].dropna().items():
            for name,row in self.samples.iterrows():
                if self.anchor_genomes is None or name in self.anchor_genomes:
                    #print("INTER", row["gff"])
                    yield (self, row["id"], name, row["fasta"], row["gff"], True)

        if not self.anno_only:
            if self.cores == 1:
                for args in iter_anchor_args():
                    print("Anchored", self.run_anchor(args))
            else:
                with mp.Pool(processes=self.cores) as pool:
                    for name in pool.imap_unordered(self.run_anchor, iter_anchor_args(), chunksize=1):
                        #print(f"Anchored {name}")
                        sys.stdout.flush()

        p = self.params
        self.genomes = {
            name : Genome(self, i, name, chrs=self.chrs) for i,name in enumerate(self.anchor_genomes)
        }

        #print("Computing chromosome summaries")
        self.chrs[self._total_occ_idx] = 0
        self.chrs[self._gene_occ_idx] = 0

        #chr_bins_out = open(self.chr_bins_fname_old, "wb")
        #bin_coords = list()
        #chr_bin_occs = list()
        rows = list()

        prev_genome = None
        for (genome,chrom),size in self.chrs["size"].items():
            if size == 0: 
                continue

            if genome != prev_genome:
                prev_genome = genome

            occs = self.query_bitmap(genome,chrom,0,size,1).sum(axis=1)
            chr_counts = self.count_occs(occs)
            self.chrs.loc[(genome,chrom), self._total_occ_idx] = chr_counts

        self._write_chrs()

        self._run_mash()

        if "gff" in self.samples.columns:

            if self.gff_anno_types is None:
                self.all_anno_types = pd.Index([])

            for g in self.samples.index:
                self.genomes[g].init_gff()
                self._load_gffs(g)

            if self.gff_anno_types is None:
                self.gff_anno_types = list(self.all_anno_types)

        #self.chrs.to_csv(f"{self.prefix}/chrs.csv")
        self._write_chrs()

        self.write_config()
        #with open(self.index_config_file, "w") as conf_out:
        #    yaml.dump(self.params, conf_out)

    def write_config(self, exclude=["prefix"]):
        prms = self.params
        for p in exclude:
            del prms[p]
        with open(self.index_config_file, "w") as conf_out:
            yaml.dump(prms, conf_out)

    @property
    def index_config_file(self):
        return os.path.join(self.prefix, "panagram.yaml")

    @property
    def index_samples_file(self):
        return os.path.join(self.prefix, "samples.tsv")

    @property
    def chr_bins_fname_old(self):
        return os.path.join(self.prefix, f"bins_{self.chr_bin_kbp}kbp.bin")

    @property
    def chr_bins_fname(self):
        return os.path.join(self.prefix, f"bins_{self.chr_bin_kbp}kbp.npy")

    def query_occ_counts(self, genome, chrom, start, end, step=1):
        return self.genomes[genome].query_occ_counts(chrom, start, end, step)
        #occs = self.query_bitmap(genome, chrom, start, end, step).sum(axis=1)
        #return self.count_occs(occs)
    
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

        index_fname = fname+".csi"
        if not os.path.exists(index_fname):
            index_fname = fname+".tbi"
            if not os.path.exists(index_fname):
                raise FileNotFoundError("Index file not found: '{fname}.csi' or '{fname.tbi}' must be preset")

        return pysam.TabixFile(fname, parser=pysam.asTuple(), index=index_fname)
    
    def get_chrs(self):
        if self.chrs is None:
            self._load_chrs()
        return self.chrs
    
    def init_chrs(self):
        genomes = list()
        genome_id = 0
        for name,fasta in self.samples["fasta"].items():#, "fastq"]]
            #if pd.isnull(fasta) == pd.isnull(fastq):
            #    raise ValueError(f"Must specify either FASTA or FASTQ for {name} (not both)")

            genome_id += 1

            is_anchor = (self.anchor_genomes is None or name in self.anchor_genomes) #pd.isnull(fastq) and 

            if is_anchor:
                if not os.path.exists(fasta+".fai"):
                    cmd = ["samtools", "faidx", fasta]
                    subprocess.check_call(cmd)
                fa = pysam.FastaFile(fasta)

                genomes.append(pd.DataFrame(
                    [(name, chrom, genome_id, fa.get_reference_length(chrom)-self.k+1) 
                     for chrom in fa.references], 
                    columns=["genome", "chr", "id", "size"]))

            else:
                genomes.append(pd.DataFrame(
                    {"genome" : [name], "id" : genome_id, "chr" : None, "size" : 0}
                ))

        self.chrs = pd.concat(genomes).set_index(["genome", "chr"])
        self._init_genomes()
        #self.chrs.to_csv(f"{self.prefix}/chrs.csv")
        self._write_chrs()
    
    def _write_chrs(self):
        out = self.chrs.copy()
        cols_out = ["_occ_".join(map(str, c)) if isinstance(c, tuple) else c
                    for c in out.columns]

        out = out.set_axis(cols_out, axis="columns")
        out.to_csv(self.chrs_file, sep="\t")

    def tabix_fname(self, genome, typ):
        #return os.path.join("anchor", genome, "{typ}.bed{TABIX_SUFFIX}") 
        return os.path.join("anchor", genome, "{typ}.bed{TABIX_SUFFIX}") 

    def _load_gffs(self, genome):
        genes = list()
        annos = list()
        fname = self.samples.loc[genome, "gff"]
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
        #self._write_tabix(annos, genome, "anno")
        self.genomes[genome]._write_tabix(annos, "anno")

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

        self.genomes[genome]._write_tabix(genes, "gene")


    def close(self):
        for b in self.genomes.values():
            b.close()


    def query_bitmap(self, genome, chrom, start=None, end=None, step=1):
        return self.genomes[genome].query(chrom, start, end, step)

    def query_genes(self, genome, chrom, start, end, attrs=["name"]):
        if self.gene_tabix.get(genome, None) is None:
            return pd.DataFrame(columns=GENE_TABIX_COLS+attrs)
        try:
            rows = self.gene_tabix[genome].fetch(chrom, start, end)
        except ValueError:
            rows = []

        ret = pd.DataFrame(rows, columns=GENE_TABIX_COLS).astype(GENE_TABIX_TYPES)
        for a in attrs:
            ret[a.lower()] = ret["attr"].str.extract(f"{a}=([^;]+)")
        
        return ret

    def query_anno(self, genome, chrom, start, end):
        if self.gene_tabix.get(genome, None) is None:
            return pd.DataFrame(columns=TABIX_COLS)
        try:
            rows = self.anno_tabix[genome].fetch(chrom, start, end)
        except ValueError:
            rows = []
        return pd.DataFrame(rows, columns=TABIX_COLS).astype(TABIX_TYPES)
    
    
    #def _query_tabix(self, tabix genome, chrom, start, end):

    @staticmethod
    def _run_kmc_genome(args):
        conf, db_i, i, name, fasta, count_db, onehot_db, tmp_dir, fasta_in = args

        onehot_id = str(2**i)

        def should_build(db):
            if not (conf["kmc"]["use_existing"] and
                    os.path.exists(db+".kmc_pre") and
                    os.path.exists(db+".kmc_suf")):
                return True
            #print(f"Using exisitng KMC db: {db}")
            return False

        in_arg = "-fm" if fasta_in else "-fq"

        if should_build(count_db):
            cmd = [
                f"{EXTRA_DIR}/kmc", f"-k{conf['k']}", 
                f"-t{conf['kmc']['threads']}", 
                f"-m{conf['kmc']['memory']}", 
                "-ci1", "-cs10000000", in_arg,
                fasta, count_db, tmp_dir
            ]
            subprocess.check_call(cmd)

        if should_build(onehot_db):
            cmd = [
                f"{EXTRA_DIR}/kmc_tools", "-t4", "transform",
                count_db, "set_counts", onehot_id, onehot_db
            ]
            subprocess.check_call(cmd)

        return db_i, name, onehot_db


    def _iter_kmc_genome_args(self):
        i = 0
        db_i = 0
        #for name,fasta in self.samples["fasta"].items():
        for name,fasta in self.samples["fasta"].items():
            if i >= 32:
                i = 0
                db_i += 1

            count_db = self.kmc_prefix(name, "count")
            onehot_db = self.kmc_prefix(name, "onehot")
            tmp_dir = self.init_dir(f"tmp/{name}")
            yield (self.params, db_i, i, name, fasta, count_db, onehot_db, tmp_dir, fasta)

            i += 1
    
    @property
    def kmc_bitvec_count(self):
        return int(np.ceil(len(self.samples) / 32.0))
    
    @property
    def opdef_filenames(self):
        return [self.kmc_prefix(f"opdef{i}.txt") for i in range(self.kmc_bitvec_count)]
    
    @property
    def bitvec_filenames(self):
        return [self.kmc_prefix(f"bitvec{i}.txt") for i in range(self.kmc_bitvec_count)]
    
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

    def _run_kmc(self):
        i = 0
        samp_count = len(self.samples)
        kmc_bitvec_count = int(np.ceil(samp_count / 32))

        genome_dbs = [list() for i in range(kmc_bitvec_count)]
        #print("PROCESS", self.kmc.processes)
        if self.cores == 1:
            for args in self._iter_kmc_genome_args():
                i,name,db = self._run_kmc_genome(args)
                genome_dbs[i].append((name,db))
        else:
            with mp.Pool(processes=self.cores) as pool:
                for i,name,db in pool.imap(self._run_kmc_genome, self._iter_kmc_genome_args(), chunksize=1):
                    genome_dbs[i].append((name,db))

        bitvec_dbs = list()

        for i in range(kmc_bitvec_count):
            h = (i+1)*32

            if kmc_bitvec_count == 1 or i < kmc_bitvec_count:
                t=32
            else:
                t = samp_count-32

            opdef_fname = self.kmc_prefix(f"opdef{i}.txt")
            bitvec_fname = self.kmc_prefix(f"bitvec{i}")

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
                f"{EXTRA_DIR}/kmc_tools", "complex", opdef_fname
            ])
            
            bitvec_dbs.append(bitvec_fname)
        return bitvec_dbs

class Genome:
    def __init__(self, idx, id, name, fasta=None, gff=None, write=False, chrs=None):
        self.samples = idx.samples
        self.params = idx.params
        self.prefix = os.path.join(self.params["prefix"], ANCHOR_DIR, name)
        self.id = id#chrs.loc[name]["id"].iloc[0]
        self.name = name
        self.fasta = fasta
        self.gff = gff
        self.write_mode = write
        
        self.ngenomes = len(self.samples)
        self.nbytes = int(np.ceil(self.ngenomes / 8))
        self.genomes = None
        self.chrs = None

        self._init_steps()

        self.seq_lens = defaultdict(dict)
        self.bitmap_lens = defaultdict(int)
        if chrs is not None:
            self.init_chrs(chrs)

        #if write_mode:
        #    self.init_write(kmc_dbs)
        #else:
        #    self.init_read()

    @property
    def chrs_file(self):
        return os.path.join(self.prefix, "chrs.tsv")
    
    def init_chrs(self):
        fa = pysam.FastaFile(self.fasta)
        chrs = pd.DataFrame(
            [(self.name, chrom, self.id, fa.get_reference_length(chrom)-self.params['k']+1) 
             for chrom in fa.references], 
            columns=["genome", "chr", "id", "size"])
        chrs.to_csv(self.chrs_file, sep="\t", index=False)
        self.set_chrs(chrs)
        return chrs

    def load_chrs(self):
        self.set_chrs(pd.read_table(self.chrs_file))

    def set_chrs(self, chrs):
        self.chrs = chrs
        print(chrs)
        self.sizes = chrs["size"]

        step_sizes = pd.DataFrame({step : np.ceil(self.sizes / step) for step in self.steps}, dtype=int)
        self.offsets = step_sizes.cumsum().shift(fill_value=0) 

        if not self.write_mode:
            self.init_read()

    def _init_steps(self):
        if "lowres_step" in self.params:
            self.steps = [1, self.params["lowres_step"]]
        else:
            self.steps = list()
            for fname in glob.glob(f"{self.prefix}.*.{BGZ_SUFFIX}"):
                step = int(fname.split(".")[-2])
                self.steps.append(step)

    @property
    def bins_fname(self): 
        kb = self.params["chr_bin_kbp"]
        return os.path.join(self.prefix, f"bitmap.{kb}kb.tsv")

    def bgz_fname(self, step): 
        return os.path.join(self.prefix, f"bitmap.{step}.{BGZ_SUFFIX}")

    def idx_fname(self, step): 
        return os.path.join(self.prefix, f"bitmap.{step}.{IDX_SUFFIX}")

    def init_read(self):
        self.blocks = {s : self.load_bgz_blocks(self.idx_fname(s)) for s in self.steps}
        self.genomes = {s : bgzf.BgzfReader(self.bgz_fname(s), "rb") for s in self.steps}

    def get_bin_counts(self):
        df = pd.read_table(self.bins_fname)
        df["chr"] = self.chrs[df["chr"]]
        df.set_index(["chr","start","end"],inplace=True)
        df.columns = df.columns.astype(int)
        return df

    def seq_len(self, seq_name):
        return self.sizes.loc[seq_name]

    def _iter_gff(self):
        for df in pd.read_csv(
            self.gff, 
            sep="\t", comment="#", chunksize=10000,
            names = ["chr","source","type","start","end","score","strand","phase","attr"],
            usecols = TABIX_COLS): yield df[TABIX_COLS]

    def init_gff(self):
        genes = list()
        annos = list()
        #fname = self.samples.loc[genome, "gff"]
        #print("GEEF", self.gff)
        if pd.isna(self.gff): return

        for df in self._iter_gff():
            gmask = df["type"].isin(self.params["gff_gene_types"])
            genes.append(df[gmask])

            if self.params["gff_anno_types"] is not None:
                annos.append(df[df["type"].isin(self.params["gff_anno_types"])])
            else:
                annos.append(df[~gmask])

        def _merge_dfs(dfs):
            return pd.concat(dfs).sort_values(["chr","start"]).reset_index(drop=True)

        annos = _merge_dfs(annos)
        self._write_tabix(annos, "anno")

        if self.params["gff_anno_types"] is None:
            self.params["gff_anno_types"] = self.params["gff_anno_types"].union(annos["type"].unique())

        genes = _merge_dfs(genes)
        genes["unique"] = 0
        genes["universal"] = 0
        for i,g in genes.iterrows():
            counts = self.query_occ_counts(g["chr"], g["start"], g["end"])
            genes.loc[i, "unique"] += counts[0]
            genes.loc[i, "universal"] += counts[-1]

        self._write_tabix(genes, "gene")

    def tabix_fname(self, typ):
        #return os.path.join("anchor", self.name, "{typ}.bed{TABIX_SUFFIX}") 
        return os.path.join("anno", f"{self.name}.{typ}.bed{TABIX_SUFFIX}") 
    
    def _write_tabix(self, df, typ):
        tbx = self.tabix_fname(typ)
        bed = tbx[:-len(TABIX_SUFFIX)]

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

    def count_occs(self, occs):
        ret = np.zeros(self.ngenomes, "uint32")
        occs, counts = np.unique(occs, return_counts=True)
        ret[occs-1] = counts
        return ret
    
    def query_occ_counts(self, name, start=None, end=None, step=1):
        occs = self.query(name,start,end,step).sum(axis=1)
        return self.count_occs(occs)

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

        self.genomes[bstep].seek(bgzf.make_virtual_offset(blk_start, blk_offs))
        buf = self.genomes[bstep].read(length * self.nbytes)

        pac = np.frombuffer(buf, "uint8").reshape((len(buf)//self.nbytes, self.nbytes))

        #pac = pac[:,::-1]

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
            from .extra import py_kmc_api
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

    def init_write(self, kmc_dbs=None):
        self._load_kmc(kmc_dbs)

        #self.genomes = {s : bgzf.BgzfWriter(self.bgz_fname(s), "wb") for s in self.steps}
        self.genomes = {s : bgzip.BGZipWriter(open(self.bgz_fname(s), "wb"))for s in self.steps}
        bin_counts = dict()

        gi = self.id
        name = self.name

        if self.fasta.endswith(".gz") or self.fasta.endswith(".bgz"):
            opn = lambda f: gzip.open(f, "rt")
        else:
            opn = lambda f: open(f, "r")

        with opn(self.fasta) as fasta:
        #with open(self.fasta, "r") as fasta:
            t = time()
            #for seq_name in fasta.references: 
            for i,rec in enumerate(SeqIO.parse(fasta, "fasta")):
                seq_name = rec.id
                seq = str(rec.seq)

                #arrs = list()
                #arrs_100nt = list()
                byte_arrs = defaultdict(list)

                def nbytes(k):
                    if self.nbytes <= 4:
                        return self.nbytes
                    elif k == len(self.kmc_dbs)-1 and self.nbytes % 4 > 0:
                        return self.nbytes % 4
                    else:
                        return 4

                sizes = defaultdict(list)
                for ki,db in enumerate(self.kmc_dbs): 
                    sys.stdout.flush()
                    pacbytes = self._get_kmc_counts(db, seq)
                    pacbytes = pacbytes[:,:nbytes(ki)]

                    for s in self.steps:
                        a = pacbytes[::s]
                        byte_arrs[s].append(a)
                        sizes[s].append(len(a))

                size = None
                arrs = dict()
                for step,arr in byte_arrs.items():
                    arrs[step] = np.concatenate(arr, axis=1)
                    self.genomes[step].write(arrs[step].tobytes())
                    self.bitmap_lens[step] += len(arrs[step])
                    if step == 1:
                        size = len(arr)

                #pacbytes
                occs = np.unpackbits(arrs[1], bitorder="little", axis=1)[:,:self.ngenomes].sum(axis=1)

                binlen = self.params["chr_bin_kbp"]*1000
                starts = np.arange(0,len(occs),binlen)
                ends = np.clip(starts+binlen, 0, len(occs))
                coords = pd.MultiIndex.from_arrays([starts, ends])
                cols = np.arange(self.ngenomes)+1
                bin_counts[i] = pd.DataFrame([
                    pd.value_counts(occs[st:en]).reindex(cols) for st,en in coords
                ], index=coords)

                self.seq_lens[gi][seq_name] = size
                sys.stdout.write(f"Anchored {seq_name}\n")
                sys.stdout.flush()

                t = time()

        self.close()
        
        bin_counts = pd.concat(bin_counts,names=["chr","start","end"])
        bin_counts.to_csv(self.bins_fname, sep="\t")

        for step in self.steps:
            subprocess.check_call([
                "bgzip", "-rI", self.idx_fname(step), self.bgz_fname(step)])
        #self._init_read()
        #self.init_gff()

    def close(self):
        if self.genomes is not None:
            for f in self.genomes.values():
                f.close()

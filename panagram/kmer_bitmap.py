import sys
import csv
import numpy as np
import pysam                        
from . import py_kmc_api as kmc            
import bgzip
import gzip
from time import time
from Bio import bgzf, SeqIO

BGZ_SUFFIX = ".pank"
BGZ_IDX_SUFFIX = ".panx"
BGZ_100NT_SUFFIX = ".100nt"
ANN_SUFFIX = ".pana"

class KmerBitmapBgz:
    def __init__(self, prefix, kmc_dbs=None, genome_tsv=None, anchors=None):
        self.prefix = prefix
        self.bgz_fname = prefix + BGZ_SUFFIX
        self.bgz_100nt_fname = prefix + BGZ_100NT_SUFFIX + BGZ_SUFFIX
        self.ann_fname = prefix + ANN_SUFFIX

        self.offsets = dict()
        self.ref_len = 0
        self.ref_len_100nt = 0

        if kmc_dbs is None and genome_tsv is None:
            self._init_read()
        else:
            self._init_write(kmc_dbs, genome_tsv, anchors)

    def _init_read(self):
        self.write_mode = False

        with open(self.ann_fname) as self.ann:
            self.ngenomes, self.nbytes = map(int, self.ann.readline().split())
            #self._set_ngenomes(ngenomes)
            #self.ngenomes = len(self.fastas)

            self.genomes = list()

            nt_offs = 0
            for line in self.ann:
                nt_offs, nt_offs_100nt, seq_name = line.split()[:4]
                self.offsets[seq_name] = (int(nt_offs), int(nt_offs_100nt))

                genome = seq_name[:seq_name.find(".")]
                if len(self.genomes) == 0 or genome != self.genomes[-1]:
                    self.genomes.append(genome)

            self.ref_len = nt_offs
            self.ref_len_100nt = nt_offs_100nt

        t = time()
        self.bgz_blocks = self.load_bgz_blocks(self.prefix+BGZ_IDX_SUFFIX)
        self.bgz = bgzf.BgzfReader(self.bgz_fname, "rb")

        self.bgz_100nt_blocks = self.load_bgz_blocks(self.prefix + BGZ_100NT_SUFFIX + BGZ_IDX_SUFFIX)
        self.bgz_100nt = bgzf.BgzfReader(self.bgz_100nt_fname, "rb")

    def load_bgz_blocks(self, fname):
        idx_in = open(fname, "rb")
        nblocks = np.fromfile(idx_in, "uint64", 1)[0]
        dtype = [("rstart", "uint64"), ("dstart", "uint64")]
        blocks = np.zeros(int(nblocks)+1, dtype=dtype)
        blocks[1:] = np.fromfile(idx_in, dtype, nblocks)
        return blocks.astype([("rstart", int), ("dstart", int)])

    def query(self, name, start, end, step=1):
        #byte_start = self.nbytes * (self.offsets[name] + start)
        #length  = end - start

        if step > 1 and step % 100 == 0:
            pac = self._query_100nt(name, start, end, step)
        else:
            pac = self._query_1nt(name, start, end, step)
        #    pac = self._query_100nt(byte_start, length, step)
        #else:
        #    pac = self._query_1nt(byte_start, length, step)

        ret = np.unpackbits(pac, bitorder="little", axis=1)[:,:self.ngenomes]

        return ret

    #def _query_1nt(self, byte_start, length, step):
    def _query_1nt(self, name, start, end, step):
        byte_start = self.nbytes * (self.offsets[name][0] + start)
        length  = end - start

        blk = np.searchsorted(self.bgz_blocks["dstart"], byte_start, side="right")-1
        blk_offs = byte_start - self.bgz_blocks["dstart"][blk]
        blk_start = self.bgz_blocks["rstart"][blk]

        self.bgz.seek(bgzf.make_virtual_offset(blk_start, blk_offs))
        buf = self.bgz.read(length * self.nbytes)

        pac = np.frombuffer(buf, "uint8").reshape((length, self.nbytes))

        if step > 1:
            return pac[::step]
        else:
            return pac

    #def _query_100nt(self, byte_start, length, step):
    def _query_100nt(self, name, start, end, step):
        byte_start = self.nbytes * (self.offsets[name][1] + (start//100))
        length  = (end - start) // 100

        #byte_start = byte_start // 100
        #length = length // 100
        step = step // 100

        blk = np.searchsorted(self.bgz_100nt_blocks["dstart"], byte_start, side="right")-1
        blk_offs = byte_start - self.bgz_100nt_blocks["dstart"][blk]
        blk_start = self.bgz_100nt_blocks["rstart"][blk]

        self.bgz_100nt.seek(bgzf.make_virtual_offset(blk_start, blk_offs))
        buf = self.bgz_100nt.read(length * self.nbytes)

        pac = np.frombuffer(buf, "uint8").reshape((length, self.nbytes))

        if step > 1:
            return pac[::step]
        else:
            return pac
    
    @property
    def genome_names(self):
        return list(self.genomes)

    def _init_write(self, kmc_files, genome_tsv, anchors):
        self.write_mode = True
        #self.bgz = bgzip.BGZipWriter(open(self.bgz_fname, "wb"))
        self.bgz = bgzf.BgzfWriter(self.bgz_fname, "wb")
        self.bgz_100nt = bgzf.BgzfWriter(self.bgz_100nt_fname, "wb")

        self.ann = open(self.ann_fname, "w")

        self.kmc_dbs = list()
        for fname in kmc_files:
            self.kmc_dbs.append(kmc.KMCFile())
            self.kmc_dbs[-1].OpenForRA(fname)

        self.genomes = list()
        self.fastas = dict()
        with open(genome_tsv) as genome_file:
            genomes_in = csv.reader(genome_file, delimiter="\t")
            for name, fasta in genomes_in:
                self.genomes.append(name)
                self.fastas[name] = fasta

        self.ngenomes = len(self.fastas)
        self.nbytes = int(np.ceil(self.ngenomes / 8))

        self.ann.write(f"{self.ngenomes}\t{self.nbytes}\n")

        for name, fasta in self.fastas.items():
            if anchors is None or name in anchors:
                self._load_fasta(name, fasta)

    def _load_fasta(self, name, fname):

        #with pysam.FastaFile(fname) as fasta:
        with gzip.open(fname, "rt") as fasta:
            vec = kmc.CountVec()
            t = time()
            #for seq_name in fasta.references: 
            #    seq = fasta.fetch(seq_name)#, 0, 10000000) 
            for rec in SeqIO.parse(fasta, "fasta"):
                seq_name = rec.id
                seq = str(rec.seq)

                arrs = list()
                arrs_100nt = list()

                def nbytes(i):
                    if self.nbytes <= 4:
                        return self.nbytes
                    elif i == len(self.kmc_dbs)-1:
                        return self.nbytes % 4
                    else:
                        return 4
                        

                for i,db in enumerate(self.kmc_dbs): 

                    db.GetCountersForRead(seq, vec)

                    arr = np.array(vec, dtype="uint32")
                    arr = arr.view("uint8").reshape((len(arr),4))
                    
                    arr = arr[:,:nbytes(i)]

                    arr_100nt = arr[::100]

                    arrs.append(arr)
                    arrs_100nt.append(arr_100nt)

                arr = np.concatenate(arrs, axis=1)
                arr_100nt = np.concatenate(arrs_100nt, axis=1)

                self.bgz.write(arr.tobytes())
                self.bgz_100nt.write(arr_100nt.tobytes())

                self.ann.write(f"{self.ref_len}\t{self.ref_len_100nt}\t{name}.{seq_name}\n")
                self.ref_len += len(arr)
                self.ref_len_100nt += len(arr_100nt)

                print(f"Anchored {name}.{seq_name}", time()-t)
                t = time()

    def close(self):
        self.ann.close()
        self.bgz.close()
        self.bgz_100nt.close()

        #if self.write_mode:
        #    self.kmc_db.close()

class KmerBitmapRaw:
    def __init__(self, prefix, kmc=None, genome_tsv=None, anchors=None):
        self.prefix = prefix
        self.bitmap_fname = prefix + BGZ_SUFFIX
        self.ann_fname = prefix + ANN_SUFFIX

        self.offsets = dict()
        self.ref_len = 0

        if kmc is None and genome_tsv is None:
            self._init_read()
        else:
            self._init_write(kmc, genome_tsv, anchors)

    def _init_read(self):
        self.write_mode = False

        with open(self.ann_fname) as self.ann:
            self.ngenomes, self.nbytes = map(int, self.ann.readline().split())

            self.genomes = list()

            nt_offs = 0
            for line in self.ann:
                nt_offs, genome, seq_name = line.split()[:3]
                self.offsets[seq_name] = int(nt_offs)
                if len(self.genomes) == 0 or genome != self.genomes[-1]:
                    self.genomes.append(genome)
            self.ref_len = nt_offs


        self.bitmap = open(self.bitmap_fname, "rb")

    def query(self, name, start, end, step=1):
        byte_start = self.nbytes * (self.offsets[name] + start)
        length  = end - start
        self.bitmap.seek(byte_start)
        buf = self.bitmap.read(length * self.nbytes)

        pac = np.frombuffer(buf, "uint8").reshape((length, self.nbytes))
        if step > 1:
            pac = pac[::step]
        ret = np.unpackbits(pac, bitorder="little", axis=1)[:,:self.ngenomes]

        return ret
    
    @property
    def genome_names(self):
        return list(self.genomes)


    def _init_write(self, kmc_file, genome_tsv, anchors):
        self.write_mode = True
        #self.bitmap = bitmapip.BGZipWriter(open(self.bitmap_fname, "wb"))
        self.bitmap = open(self.bitmap_fname, "wb")
        self.ann = open(self.ann_fname, "w")

        self.kmc_db = kmc.KMCFile()        
        self.kmc_db.OpenForRA(kmc_file)

        self.genomes = list()
        self.fastas = dict()
        with open(genome_tsv) as genome_file:
            genomes_in = csv.reader(genome_file, delimiter="\t")
            for name, fasta in genomes_in:
                self.genomes.append(name)
                self.fastas[name] = fasta

        self.ngenomes = len(self.fastas)
        self.nbytes = int(np.ceil(self.ngenomes / 8))

        self.ann.write(f"{self.ngenomes}\t{self.nbytes}\n")

        for name, fasta in self.fastas.items():
            if anchors is None or name in anchors:
                self._load_fasta(name, fasta)

    def _load_fasta(self, name, fname):

        #with pysam.FastaFile(fname) as fasta:
        with gzip.open(fname, "rt") as fasta:
            vec = kmc.CountVec()
            t = time()
            for rec in SeqIO.parse(fasta, "fasta"):
                seq_name = rec.id
                seq = str(rec.seq)

                self.kmc_db.GetCountersForRead(seq, vec)

                arr = np.array(vec.value, dtype="uint32")
                arr = arr.view("uint8").reshape((len(arr),4))
                if self.nbytes < 4:
                    arr = arr[:,:self.nbytes]

                self.bitmap.write(arr.tobytes())

                self.ann.write(f"{self.ref_len}\t{name}\t{seq_name}\n")
                self.ref_len += len(arr)

                print(f"Anchored {seq_name}", time()-t)
                t = time()

    def close(self):
        self.ann.close()
        self.bitmap.close()

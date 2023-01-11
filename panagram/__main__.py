import sys
import argparse
from dataclasses import dataclass
from simple_parsing import ArgumentParser

from .kmer_bitmap import KmerBitmapBgz as KmerBitmap

from .view import view

def parse_coords(coords):
    name, coords = coords.split(":")
    start, end = map(int, coords.split("-"))
    return name, start, end

def comma_split(s):
    return s.split(",")

#@dataclass
#class Main:
#    index_dir: str
#    k: int = 21
#
#@dataclass
#class Index:
#    genome_tsv: str
#    kmc_files: str


def anchor_opts(sp):
    p = sp.add_parser("anchor")
    p.add_argument("genome_tsv", help="Anchor genome fasta files")
    p.add_argument("kmc_prefixes", help="KMC Database Prefix", nargs="+")
    p.add_argument("-o", "--out-prefix", required=True, help="Anchor output prefix")
    p.add_argument("-g", "--genome-count", type=int, default=None)
    p.add_argument("-a", "--anchor-genomes", default=None, type=comma_split)

def bitdump_opts(sp):
    p = sp.add_parser("bitdump")
    p.add_argument("in_prefix", help="Anchor genome input prefix")
    p.add_argument("coords", default=None, type=parse_coords, help="Coordinates to query")
    p.add_argument("step", nargs="?", default=1, type=int)
    p.add_argument("-v", "--verbose", action="store_true")

def view_opts(sp):
    p = sp.add_parser("view")
    p.add_argument("config", help="Config file")

parser = argparse.ArgumentParser("Panagram")
sp = parser.add_subparsers(dest="subcmd")
anchor_opts(sp)
bitdump_opts(sp)
view_opts(sp)

def anchor(args): #prefix, kmc, *fastas):
    bits = KmerBitmap(args.out_prefix, args.kmc_prefixes, args.genome_tsv, args.anchor_genomes)
    bits.close()

def bitdump(args):
    bitmap = KmerBitmap(args.in_prefix)
    name, start, end = args.coords
    bits = bitmap.query(name, start, end, args.step)

    if args.verbose:
        print(" ".join(bitmap.genome_names))
        for i in range(len(bits)):
            print(" ".join(bits[i].astype(str)))
    else:
        print(bits)

    bitmap.close()

cmds = {"view", "anchor", "bitdump"}

def main():
    args = parser.parse_args()
    cmd = args.subcmd
    if cmd not in cmds:
        raise ValueError(f"Unknown command: {cmd}")
    globals()[cmd](args)

#if __name__ == "__main__":
#    main()

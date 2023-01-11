import sys
import argparse

from simple_parsing import ArgumentParser, field
from dataclasses import dataclass
from typing import Any, List, Tuple, Type, Union

@dataclass
class View:
    """Display panagram viewer"""
    config: str = field(positional=True)

    def run(self):
        from .view import view
        view(self.config)

@dataclass
class Index:
    """Create anchor genomes"""
    genomes: str = field(positional=True)
    kmc_files: List[str] = field(positional=True)
    out_dir: str = field(alias=["-o"], required=True)
    k: int = field(alias=["-k"], default=21)

    def run(self):
        from .kmer_bitmap import KmerBitmapBgz as KmerBitmap
        bits = KmerBitmap(
            self.out_dir, 
            self.kmc_files, 
            self.genomes)
        bits.close()

@dataclass
class Bitdump:
    """Query pan-kmer bitmap"""
    index_dir: str = field(positional=True)
    coords: str = field(positional=True)
    step: str = field(positional=True, nargs="?", default=1)
    verbose: bool = field(alias=["-v"], default=False)

    def run(self):
        from .kmer_bitmap import KmerBitmapBgz as KmerBitmap
        bitmap = KmerBitmap(self.index_dir)
        name, start, end = parse_coords(self.coords)
        bits = bitmap.query(name, start, end, self.step)

        if self.verbose:
            print(" ".join(bitmap.genome_names))
            for i in range(len(bits)):
                print(" ".join(bits[i].astype(str)))
        else:
            print(bits)

        bitmap.close()

@dataclass
class Main:
    """Panagram!!"""
    cmd: Union[View, Index, Bitdump]

    def run(self):
        return self.cmd.run()

def parse_coords(coords):
    name, coords = coords.split(":")
    start, end = map(int, coords.split("-"))
    return name, start, end

def comma_split(s):
    return s.split(",")

def main():
    parser = ArgumentParser()
    parser.add_arguments(Main, dest="main")
    parser.parse_args().main.run()

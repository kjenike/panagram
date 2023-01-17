import sys
import os
import csv
import subprocess
from .kmer_bitmap import KmerBitmapBgz as KmerBitmap

ROOT_DIR = os.path.dirname(os.path.realpath(__file__))
KMC_DIR = os.path.join(ROOT_DIR, "kmc")


def index(genomes, out_dir, k):
    from .kmer_bitmap import KmerBitmapBgz as KmerBitmap

    def init_dir(name):
        d = os.path.join(out_dir, name)
        os.makedirs(d, exist_ok=True)

        return d

    count_dir = init_dir("kmc_count")
    onehot_dir = init_dir("kmc_onehot")
    bitvec_dir = init_dir("kmc_bitvec")
    tmp_dir = init_dir("tmp")

    i = 0
    db_count = 1
    samp_count = 0

    genome_dbs = list()
    
    with open(genomes) as genome_file:
        genomes_in = csv.reader(genome_file, delimiter="\t")
        for name, fasta in genomes_in:
            if i >= 32:
                i = 0
                db_count += 1


            onehot_id = str(2**i)

            print(f"Counting {name}")
            
            count_db = os.path.join(count_dir, name)
            onehot_db = os.path.join(onehot_dir, name)

            genome_dbs.append((name, os.path.abspath(onehot_db)))

            subprocess.check_call([
                f"{KMC_DIR}/kmc", f"-k{k}", 
                "-t4", "-m8", "-ci1", "-cs10000000", "-fm",
                fasta, count_db, tmp_dir
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

        opdef_fname = os.path.join(bitvec_dir, f"{i}.opdef.txt")
        bitvec_fname = os.path.abspath(os.path.join(bitvec_dir, f"{i}"))

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

    bits = KmerBitmap(f"{out_dir}/", bitvec_dbs, genomes)
    bits.close()

    subprocess.check_call([
        "bgzip", "-r", f"{out_dir}/anchor.pank", 
                 "-I", f"{out_dir}/anchor.panx"])

    subprocess.check_call([
        "bgzip", "-r", f"{out_dir}/anchor.100nt.pank", 
                 "-I", f"{out_dir}/anchor.100nt.panx"])

from setuptools import setup
from setuptools.command.build_py import build_py
import subprocess
import os
import glob
import sys

import platform
from urllib.request import urlretrieve
import shutil
import tarfile

ROOT_DIR = os.getcwd()
KMC_DIR = os.path.join(ROOT_DIR, "KMC")

SUBMODS = ["KMC"]

def call_cmd(*cmd):
    s = " ".join(cmd)
    sys.stdout.write(f'Calling "{s}"\n')
    ret = subprocess.run(cmd)
    return ret.returncode == 0

#kmc_api = Pybind11Extension(
#    "panagram.kmc",
#
#    sources = [ 
#        "KMC/py_kmc_api.cpp",
#        "KMC/kmc_api/kmc_file.cpp",
#        "KMC/kmc_api/kmer_api.cpp",
#        "KMC/kmc_api/mmer.cpp"
#    ],
#
#    include_dirs = [
#        "./submods",
#        "./submods/hdf5/include", 
#        "./submods/fast5/include",
#        "./submods/pdqsort",
#        "./submods/toml11",
#        get_pybind_include()
#    ],
#
#
#    libraries = ["bwa", "z", "dl", "m"],
#
#    extra_compile_args = ["-O3", "-Wall", "-m64", "-std=c++14"],
#
#    define_macros = [("PYBIND", None)]#, ("PYDEBUG", None)]
#)

KMC_BINS = ["kmc", "kmc_tools", "kmc_dump"]
KMC_URLS = {
    "Linu" : "https://github.com/refresh-bio/KMC/releases/download/v3.2.1/KMC3.2.1.linux.tar.gz",
    "Darwin" : "https://github.com/refresh-bio/KMC/releases/download/v3.2.1/KMC3.2.1.mac.tar.gz",
    "Windows" : "https://github.com/refresh-bio/KMC/releases/download/v3.2.1/KMC3.2.1.windows.tar.gz"
}

class pre_build(build_py):
    def run(self):
        make_cmd = ["make", "-j", "4", "-C", KMC_DIR, "py_kmc_api"]

        plat = platform.system()
        kmc_url = KMC_URLS.get(plat, None)
        if kmc_url is None:
            print(f"Warning: no KMC binaries found for system {plat}, attempting to compile")
            make_cmd.append("all")

        if call_cmd(*make_cmd):
            kmc_dest = f"{ROOT_DIR}/panagram/kmc"
            shutil.rmtree(kmc_dest)
            call_cmd("mv", "-f", f"{KMC_DIR}/bin", kmc_dest)

            if kmc_url is not None:
                tar_fname = f"{kmc_dest}/kmc.tar.gz"
                urlretrieve(kmc_url, tar_fname)
                kmc_tar = tarfile.open(tar_fname)
                kmc_tar.extractall(kmc_dest)
                kmc_tar.close()
                os.remove(tar_fname)

                for b in KMC_BINS:
                    os.rename(f"{kmc_dest}/bin/{b}", f"{kmc_dest}/{b}")
                shutil.rmtree(f"{kmc_dest}/bin")
                shutil.rmtree(f"{kmc_dest}/include")
        else:
            print("Warning: KMC failed to install. 'panagram index' will not be functional, but 'panagram view' will work. See https://github.com/kjenike/panagram#readme for more information")


        build_py.run(self)

if __name__ == "__main__":
    setup(
        cmdclass={'build_py': pre_build},
        packages=["panagram", "panagram.kmc"],
    )

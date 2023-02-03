from setuptools import setup
from setuptools.command.build_py import build_py
from setuptools.command.build import SubCommand
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

KMC_BINS = ["kmc", "kmc_tools", "kmc_dump"]
KMC_URLS = {
    "Linux" : "https://github.com/refresh-bio/KMC/releases/download/v3.2.1/KMC3.2.1.linux.tar.gz",
    "Darwin" : "https://github.com/refresh-bio/KMC/releases/download/v3.2.1/KMC3.2.1.mac.tar.gz",
    #"Windows" : "https://github.com/refresh-bio/KMC/releases/download/v3.2.1/KMC3.2.1.windows.tar.gz"
}

MASH_URLS = { 
    "Linux" : "https://github.com/marbl/Mash/releases/download/v2.3/mash-Linux64-v2.3.tar",
    "Darwin" : "https://github.com/marbl/Mash/releases/download/v2.3/mash-OSX64-v2.3.tar",
}

class pre_build(build_py):
    bulid_lib = ROOT_DIR
    
    def run(self):
        call_cmd("make", "-C", KMC_DIR, "clean")
        make_cmd = ["make", "-C", KMC_DIR, "py_kmc_api"]

        plat = platform.system()
        kmc_url = KMC_URLS.get(plat, None)
        if kmc_url is None:
            print(f"Warning: no KMC binaries found for system {plat}, attempting to compile")
            make_cmd.append("all")

        extra_dir = f"{ROOT_DIR}/panagram/extra"

        if call_cmd(*make_cmd):
            shutil.rmtree(extra_dir)
            call_cmd("mv", "-f", f"{KMC_DIR}/bin", extra_dir)

            open(f"{extra_dir}/__init__.py","w").close()

            if kmc_url is not None:
                tar_fname = f"{extra_dir}/kmc.tar.gz"
                urlretrieve(kmc_url, tar_fname)
                kmc_tar = tarfile.open(tar_fname)
                kmc_tar.extractall(extra_dir)
                kmc_tar.close()
                os.remove(tar_fname)

                for b in KMC_BINS:
                    os.rename(f"{extra_dir}/bin/{b}", f"{extra_dir}/{b}")
                shutil.rmtree(f"{extra_dir}/bin")
                shutil.rmtree(f"{extra_dir}/include")
        else:
            print("Warning: KMC failed to install. 'panagram index' will not be functional, but 'panagram view' will work. See https://github.com/kjenike/panagram#readme for more information")

        mash_url = MASH_URLS.get(plat, None)
        if mash_url is None:
            print(f"Warning: no MASH libraries found for system {plat}. You will be unable to run panagram index")
        else:
            tar_fname = f"{extra_dir}/mash.tar.gz"
            urlretrieve(mash_url, tar_fname)
            tar = tarfile.open(tar_fname)
            tar.extractall(extra_dir)
            tar.close()
            os.remove(tar_fname)

            dirname = os.path.basename(mash_url)[:-4]

            os.rename(f"{extra_dir}/{dirname}/mash", f"{extra_dir}/mash")
            shutil.rmtree(f"{extra_dir}/{dirname}")

        build_py.run(self)

    #def get_output_mapping(self):
    #    print("GET OUTPUT MAPPING")
    #    return build_py.get_output_mapping()
    #
    #def get_outputs(self):
    #    print("GET OUTPUTS")
    #    return build_py.get_outputs(self)

    #def get_source_files(self):
    #    print("SOURCES OUT", build_py.get_source_files(self))
    #    return build_py.get_source_files(self) + ["KMC/py_kmc_api/py_kmc_api.cpp"]

    #def initialize_options(self):
    #    return build_py.initialize_options(self)

if __name__ == "__main__":
    setup(
        cmdclass={'build_py': pre_build},
        packages=["panagram", "panagram.extra"],
    )

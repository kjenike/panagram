from setuptools import setup
from setuptools.command.build_py import build_py
import subprocess
import os
import glob

ROOT_DIR = os.getcwd()
KMC_DIR = os.path.join(ROOT_DIR, "KMC")

SUBMODS = ["KMC"]

class pre_build(build_py):
    def run(self):
        subprocess.check_call([
            "make", "-j", "4", "-C", KMC_DIR, "py_kmc_api"#, "all" 
        ])

        #kmc_libs = glob.glob(KMC_DIR+"/bin/py_kmc_api*.so")
        #assert(len(kmc_libs) == 1)
        #mv = ["mv"] + kmc_libs + [ROOT_DIR+"/panagram"] 

        mv = ["mv", f"{KMC_DIR}/bin", f"{ROOT_DIR}/panagram/kmc"]
        subprocess.check_call(mv)
        build_py.run(self)

if __name__ == "__main__":
    setup(
        cmdclass={'build_py': pre_build},
        packages=["panagram"],
    )

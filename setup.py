from setuptools import setup
from setuptools.command.build_py import build_py
import subprocess
import os
import glob
import sys

ROOT_DIR = os.getcwd()
KMC_DIR = os.path.join(ROOT_DIR, "KMC")

SUBMODS = ["KMC"]

def call_cmd(*cmd):
    s = " ".join(cmd)
    sys.stdout.write(f'Calling "{s}"\n')
    subprocess.check_call(cmd)

class pre_build(build_py):
    def run(self):
        call_cmd("make", "-j", "4", "-C", KMC_DIR, "py_kmc_api", "all")
        call_cmd("mv", f"{KMC_DIR}/bin", f"{ROOT_DIR}/panagram/kmc")
        build_py.run(self)

if __name__ == "__main__":
    setup(
        cmdclass={'build_py': pre_build},
        packages=["panagram"],
    )

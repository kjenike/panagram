import os

# snakemake sets OMP_NUM_THREADS to the rule's thread count, and rule umaps declares 0
# (it holds no core); libgomp rejects 0 with a warning, so drop it before numpy loads.
if os.environ.get("OMP_NUM_THREADS") == "0":
	del os.environ["OMP_NUM_THREADS"]

from .__main__ import main
from .__about__ import (
    __title__,
    __version__,
    __summary__,
    __author__,
    __uri__,
    __license__,
    __copyright__,
)

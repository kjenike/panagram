# panagram test suite

Run from anywhere with a Python (3.11-3.13) that has panagram's dependencies and pytest:

```bash
PYTHONPATH=<clone> OMP_NUM_THREADS=4 python -m pytest -q <clone>/tests
```

`PYTHONPATH` must point at the clone so the code under test is the clone's; the suite exits
if `panagram` is imported from anywhere else. All genomes are synthetic and seeded
(`tests/synthetic.py`); no external data is used.

## KMC, mash and `EXTRA_DIR`

`panagram/extra/` in a clone has no binaries. The anchoring fixtures need `kmc`, `kmc_tools`,
`mash` and `py_kmc_api` plus `samtools` and `bgzip` on `PATH`. They look in
`$PANAGRAM_EXTRA_DIR`, then in an installed `panagram` distribution, and skip (never error) if
neither has them. The `panagram_overlay` fixture copies the clone's package into a temporary
directory, links the binaries into its `extra/`, and runs `panagram` and `snakemake` in
subprocesses with `PYTHONPATH` set to that copy, so `panagram.index.EXTRA_DIR` is correct in the
Snakemake rules. If `py_kmc_api` cannot load because a conda `libstdc++` is older than the system
one, the fixture retries with the system `libstdc++` in `LD_PRELOAD`.

## Fixtures (`tests/conftest.py`)

| fixture | scope | gives |
| --- | --- | --- |
| `synthetic_genomes` | session | `Genomes`: FASTA, GFF3 and `samples.tsv` for five genomes (alpha..epsilon, epsilon without GFF), `k = 15`; ground truth in `.data`, `.owners()`, `.expected_bits()`, `.oracle_bits()` |
| `prepared_index` | session | `Path` of an index directory after `panagram index samples.tsv -k 15 --prepare` |
| `fresh_prepared_index` | function | the same, in a private directory a test may modify |
| `anchored_index` | session | `Path` of a complete index (`snakemake --cores 2 all`), built once, about 20 s; read-only |
| `anchored_reader` | function | `panagram.index.Index` in read mode over `anchored_index`; cwd is the index directory for the test |
| `run_panagram` | session | `run_panagram(*args, cwd=None)` runs the `panagram` CLI in a subprocess, raising on failure |
| `panagram_overlay`, `extra_dir` | session | subprocess environment and the binaries directory |

`open_index(path)` opens an index and restores the working directory (`Index` calls `os.chdir`).
An autouse fixture restores the working directory after every test.

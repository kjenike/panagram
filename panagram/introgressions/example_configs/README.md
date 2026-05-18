# Additional Config File Information

This folder contains the 2-way and 3-way introgression parameter config files used for Panagram's
paper. The command used was `panagram intros <path to config.yaml> --sweep`. For 2-way calling,
the `--sweep` flag runs the introgression caller for the following set of thresholds:
0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95.
For 3-way calling, these thresholds are used instead:
0.0, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.36, 0.4, 0.44, 0.48, 0.52, 0.56, 0.6, 0.64,
0.68.

The introgression caller was tested on 3 pangenomes: a simulated *Arabidopsis thaliana*
pangenome, a *Solanum lycopersicum* pangenome, and a *Solanum aethiopicum* pangenome.

The simulated pangenome was generated with the following command:

```bash
panagram intros simulate \
  --ref ./Athal_edited.fasta \
  --out-folder ./arabidopsis \
  --num-introgressions 10 \
  --introgression-size-min 2000000 \
  --introgression-size-max 6000000 \
  --rel-ins-size-min 1 \
  --rel-ins-size-max 1000 \
  --rel-del-size-min 1 \
  --rel-del-size-max 500 \
  --rel-sub-rate 3.3e-3 \
  --rel-ins-rate 3.3e-3 \
  --rel-del-rate 3.3e-3 \
  --mut-ins-size-min 1 \
  --mut-ins-size-max 1000 \
  --mut-del-size-min 1 \
  --mut-del-size-max 500 \
  --mut-rate-start 1e-3 \
  --mut-sub-rate 1e-2 \
  --mut-ins-rate 1e-2 \
  --mut-del-rate 1e-2 \
  --seed 42
```

`Athal_edited.fasta` is the TAIR10 reference for *Arabidopsis thaliana*, edited so that the
chromosomes are named with the convention `chr*` and so chromosomes from the
mitochondria and choloroplasts are excluded.
Once generated and indexed, the following config files were used for introgression calling:

2-way calling: <br>
simulator_2way_config.yaml

3-way calling: <br>
simulator_Gen1_3way_config.yaml <br>
simulator_Gen2_3way_config.yaml <br>
simulator_Gen3_3way_config.yaml <br>
simulator_Gen4_3way_config.yaml <br>
simulator_Gen5_3way_config.yaml <br>
simulator_Gen6_3way_config.yaml

Usually, you can use one config file for all accessions for 2-way and 3-way calling. For the
simulator, 3-way calling was done separately for each accession in order to use the ground
truth file corresponding to the coordinate space of each individual accession. The parameters
in the 3-way config
files are otherwise the same. For the other pangenomes, liftover was used in order to convert
3-way results to reference space because ground truth files were only available in that space.

*Solanum lycopersicum's* pangenome used the following config files:

2-way calling: <br>
lycopersicum_2way_config.yaml

3-way calling: <br>
lycopersicum_3way_config.yaml

*Solanum aethiopicum's* pangenome used the following config files:

2-way calling: <br>
aethiopicum_2way_config.yaml

3-way calling: <br>
aethiopicum_3way_config.yaml

See the [introgressions README](../README.md) for a breakdown of what parameters mean
and how to use them.

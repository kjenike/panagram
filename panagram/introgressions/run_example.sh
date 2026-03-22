#!/bin/bash

set -e

# create a simulated wild relative and offspring using chr1 from Arabidopsis thaliana
echo "Simulating introgressions..."
mkdir -p ./example/simulated_data
python simulate_introgressions.py \
  --ref ./example/FASTAS/athaliana_chr1.fasta \
  --out-folder ./example/simulated_data \
  --num-introgressions 2 \
  --introgression-size-min 3000000 \
  --introgression-size-max 7000000 \
  --rel-ins-size-min 1 \
  --rel-ins-size-max 1000 \
  --rel-del-size-min 1 \
  --rel-del-size-max 500 \
  --mut-ins-size-min 1 \
  --mut-ins-size-max 1000 \
  --mut-del-size-min 1 \
  --mut-del-size-max 500 \
  --seed 7 \
  > ./example/simulated_data/simulator.log 2>&1

# run panagram on the simulated data
mv ./example/simulated_data/*.fasta ./example/FASTAS/
cd ./example
panagram index samples.tsv -k 31 --prepare
snakemake --cores 1 all

# create heatmaps for Reference and Offspring after omitting fixed kmers
anchors=("Reference" "OffspringGen1" "OffspringGen2" "OffspringGen3" "OffspringGen4" "OffspringGen5" "OffspringGen6")

for anchor in "${anchors[@]}"; do
    echo "Creating heatmap for $anchor..."
    python ../create_heatmap.py \
        --index-dir "$(pwd)" \
        --anchor "$anchor" \
        --groups group.tsv
done

cd ..

python bed_to_txt.py \
--gt_bed_file ./example/simulated_data/athaliana_chr1_0_introgressions.bed \
--index_dir ./example \
--ref Reference \
--wild_type WildRelative \
--wild_type_group WT

python introgression_runner.py ./example/2way_example_config.yaml
python introgression_runner.py ./example/3way_example_config.yaml

echo "Done! Check out the kmer simularity plots in example/panagram_visuals/ and the introgression outputs in example/introgressions/!"

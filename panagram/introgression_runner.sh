#!/bin/bash
# Runs the introgression caller once with parameters set by the user

set -e

base_output_dir="/home/nbrown62/data_mschatz1/nbrown62/panagram_data/tomato_sl4/introgressions"
index_dir="/home/nbrown62/data_mschatz1/nbrown62/panagram_data/tomato_sl4"
tsv="/home/nbrown62/data_mschatz1/nbrown62/panagram_data/tomato_sl4/group.tsv"
call_flags="--vis --rmf --bin 500000 -a M82 -c chr4 -p REF --urf SL4 --gnm -1 --sft mean --ssz 4"
threshold=0.7
call_dir_stem=test
postprocess_flags="-a fgap rmbn --bin 500000 --min 2 --gap 1 --ref SL4 --paf /home/nbrown62/data_mschatz1/nbrown62/panagram_data/tomato_sl4/alignments"
gdt_dir="/home/nbrown62/data_mschatz1/nbrown62/CallIntrogressions_data/tomato_sl4_paper_subset_500kb"
score_flags="--vis -g SP SLC -a fgap rmbn --bin 500000 --min 2 --gap 1 --ref SL4"

call_dir="${base_output_dir}/${call_dir_stem}/${call_dir_stem}_${threshold}"
postprocess_dir="${call_dir}/postprocessed"
score_dir="${call_dir}/scored_bins"
mkdir -p "$call_dir"

python call_introgressions.py --thr "$threshold" $call_flags \
    --idx "$index_dir" \
    --tsv "$tsv" \
    --out "$call_dir"

python postprocess_introgressions.py $postprocess_flags \
    --bed "$call_dir" \
    --idx "$index_dir" \
    --out "$postprocess_dir"

python score_introgressions.py $score_flags \
    --how bins \
    --pre "$postprocess_dir" \
    --gdt "$gdt_dir" \
    --idx "$index_dir" \
    --out "$score_dir"

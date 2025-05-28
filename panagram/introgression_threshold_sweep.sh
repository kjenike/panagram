#!/bin/bash
# Sweep through thresholds 0.1-0.9 to create a ROC/PR curve for introgressions

set -e

base_output_dir="/home/nbrown62/data_mschatz1/nbrown62/panagram_data/tomato_sl4"
index_dir="/home/nbrown62/data_mschatz1/nbrown62/panagram_data/tomato_sl4"
tsv="/home/nbrown62/data_mschatz1/nbrown62/panagram_data/tomato_sl4/group.tsv"
gdt_dir="/home/nbrown62/data_mschatz1/nbrown62/CallIntrogressions_data/tomato_sl4_paper_subset"
call_flags="--rmf -g SLL -p SLC SP REF"
postprocess_flags="-a lift fgap rmbn --ref SL4 --paf /home/nbrown62/data_mschatz1/nbrown62/panagram_data/tomato_sl4/alignments"
score_flags="--how bins -g SP SLC -a fgap --vis"

for i in $(seq 1 9); do
    threshold=$(echo "scale=1; $i / 10" | bc)
    call_dir="${base_output_dir}/rmf_fgap_rmbn_bins_0${threshold}"
    postprocess_dir="${call_dir}/postprocessed"
    score_dir="${call_dir}/scored"
    echo "Score dir: $score_dir"

    python call_introgressions.py --thr "$threshold" $call_flags \
    --idx "$index_dir" \
    --tsv "$tsv" \
    --out "$call_dir"

    python postprocess_introgressions.py $postprocess_flags \
    --bed "$call_dir" \
    --idx "$index_dir" \
    --out "$postprocess_dir"

    python score_introgressions.py $score_flags \
    --pre "$postprocess_dir" \
    --gdt "$gdt_dir" \
    --idx "$index_dir" \
    --out "$score_dir"

done

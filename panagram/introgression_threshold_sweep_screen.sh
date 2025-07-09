#!/bin/bash
# Sweep through thresholds 0.1-0.9 to create a ROC/PR curve for introgressions
# For slower running jobs

set -e

base_output_dir="/home/nbrown62/data_mschatz1/nbrown62/panagram_data/tomato_sl4/introgressions"
index_dir="/home/nbrown62/data_mschatz1/nbrown62/panagram_data/tomato_sl4"
tsv="/home/nbrown62/data_mschatz1/nbrown62/panagram_data/tomato_sl4/group.tsv"
call_flags="--vis --rmf --bin 500000 -g SLL -p REF --urf SL4 --gnm -1 --sft mean --ssz 4"
thresholds=(0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95)
call_dir_stem=pre_500kb_500kbgt
postprocess_flags="-a fgap rmbn --bin 500000 --min 2 --gap 1 --ref SL4 --paf /home/nbrown62/data_mschatz1/nbrown62/panagram_data/tomato_sl4/alignments"
gdt_dir="/home/nbrown62/data_mschatz1/nbrown62/CallIntrogressions_data/tomato_sl4_paper_subset_20"
score_flags="--vis -g SP SLC -a fgap rmbn --bin 500000 --min 2 --gap 1 --ref SL4 --thr 0"

rm -rf "${base_output_dir:?}/${call_dir_stem}"
for threshold in "${thresholds[@]}"; do
    call_dir="${base_output_dir}/${call_dir_stem}/${call_dir_stem}_${threshold}"
    postprocess_dir="${call_dir}/postprocessed"
    score_dir="${call_dir}/scored_bins"
    logfile="${call_dir}/sweep.log"
    screen_name="thr_${threshold}"
    mkdir -p "$call_dir"

    screen -d -m -L -Logfile "$logfile" -S "$screen_name" bash -c "
        python call_introgressions.py --thr \"$threshold\" $call_flags \
            --idx \"$index_dir\" \
            --tsv \"$tsv\" \
            --out \"$call_dir\"

        python postprocess_introgressions.py $postprocess_flags \
            --bed \"$call_dir\" \
            --idx \"$index_dir\" \
            --out \"$postprocess_dir\"

        python score_introgressions.py $score_flags \
            --how bins \
            --pre \"$postprocess_dir\" \
            --gdt \"$gdt_dir\" \
            --idx \"$index_dir\" \
            --out \"$score_dir\"

        touch \"${call_dir}/done.marker\"
    "
done

# Wait for all done.marker files to exist
for threshold in "${thresholds[@]}"; do
    call_dir="${base_output_dir}/${call_dir_stem}/${call_dir_stem}_${threshold}"
    while [ ! -f "${call_dir}/done.marker" ]; do
        echo "Waiting for $threshold..."
        sleep 15
    done
done
echo "All screens completed."

call_dir="${base_output_dir}/${call_dir_stem}"

python visualize_introgressions.py \
-v prc prcc prca shtmp \
--dir "$call_dir" \
--how bins

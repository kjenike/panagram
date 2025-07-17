import sys
from pathlib import Path
import shutil
import subprocess
import time
import yaml

def parse_config(config_path):
    with config_path.open() as f:
        config = yaml.safe_load(f)

    general = config.get("general")
    calling = config.get("calling")
    postprocessing = config.get("postprocessing")
    scoring = config.get("scoring")

    output_dir = Path(general["output_dir_base"]) / general["output_dir_stem"]
    output_dir.mkdir(parents=True, exist_ok=True)
    # Copy the config file into the output directory
    config_copy_path = output_dir / "intro_config.yaml"
    shutil.copy(config_path, config_copy_path)

    index_dir = Path(general["index_dir"])
    tsv = Path(general["tsv"])
    bin_size = general["bin"]
    ref = general["ref"]

    call_run = calling["run"]
    call_rmf = calling["rmf"]
    call_thr = calling["thr"]
    call_stp = calling["stp"]
    call_grp = calling["grp"]
    call_anc = calling["anc"]
    call_chr = calling["chr"]
    call_cmp = calling["cmp"]
    call_urf = calling["urf"]
    call_gnm = calling["gnm"]
    call_trm = calling["trm"]
    call_sft = calling["sft"]
    call_ssz = calling["ssz"]
    call_cst = calling["cst"]
    call_isc = calling["isc"]
    call_vis = calling["vis"]

    post_run = postprocessing["run"]
    post_act = postprocessing["act"]
    post_min = postprocessing["min"]
    post_gap = postprocessing["gap"]
    post_map = postprocessing["map"]
    post_paf = postprocessing["paf"]

    score_gdt_dir = Path(scoring.get("gdt_dir", ""))
    score_run = scoring["run"]
    score_act = scoring["act"]
    score_min = scoring["min"]
    score_gap = scoring["gap"]
    score_how = scoring["how"]
    score_thr = scoring["thr"]
    score_cmp = scoring["cmp"]
    score_vis = scoring["vis"]
    score_othr = scoring["othr"]

    # Build flags for each step
    call_flags = None
    if call_run:
        call_flags = [
            f"--idx {index_dir}",
            f"--tsv {tsv}",
            f"--grp {call_grp}" if call_grp else None,
            f"--anc {call_anc}" if call_anc else None,
            f"--chr {call_chr}" if call_chr else None,
            f"--cmp {call_cmp}",
            f"--bin {bin_size}",
            f"--stp {call_stp}",
            f"--gnm {call_gnm}" if call_gnm else None,
            f"--trm {call_trm}" if call_trm else None,
            f"--sft {call_sft}" if call_sft else None,
            f"--ssz {call_ssz}" if call_ssz else None,
            f"--urf {ref}" if call_urf else None,
            "--rmf" if call_rmf else None,
            "--cst" if call_cst else None,
            "--isc" if call_isc else None,
            "--vis" if call_vis else None,
        ]
        call_flags = [f for f in call_flags if f]

    postprocess_flags = None
    if post_run:
        postprocess_flags = [
            f"--idx {index_dir}",
            f"--act {' '.join(post_act)}",
            f"--bin {bin_size}",
            f"--min {post_min}" if post_min else None,
            f"--gap {post_gap}" if post_gap else None,
            f"--map {post_map}" if post_map else None,
            f"--paf {post_paf}" if post_paf else None,
            f"--ref {ref}",
        ]
        postprocess_flags = [f for f in postprocess_flags if f]

    score_flags = None
    if score_run:
        score_flags = [
            f"--idx {index_dir}",
            f"--gdt {score_gdt_dir}",
            f"--how {score_how}",
            f"--thr {score_thr}",
            f"--othr {score_othr}" if score_othr else None,
            f"--cmp {' '.join(score_cmp)}",
            f"--act {' '.join(score_act)}",
            f"--bin {bin_size}",
            f"--min {score_min}" if score_min else None,
            f"--gap {score_gap}" if score_gap else None,
            f"--ref {ref}",
            "--vis" if score_vis else None,
        ]
        score_flags = [f for f in score_flags if f]
    return call_flags, postprocess_flags, score_flags, output_dir, call_thr

def run_introgression_pipeline(call_flags, postprocess_flags, score_flags, output_dir, call_thr):
    # Create output directories
    call_dir = output_dir / f"{output_dir.name}_{call_thr}"
    call_dir.mkdir(parents=True, exist_ok=True)

    # Run pipeline
    if call_flags:
        call_flags += [f"--out {call_dir}", f"--thr {call_thr}"]
        subprocess.run(f"python call_introgressions.py {' '.join(call_flags)}", shell=True, check=True)
    if postprocess_flags:
        postprocess_dir = call_dir / "postprocessed"
        postprocess_flags += [f"--bed {call_dir}", f"--out {postprocess_dir}"]
        subprocess.run(f"python postprocess_introgressions.py {' '.join(postprocess_flags)}", shell=True, check=True)
    if score_flags:
        how_to_score = None
        for flag in score_flags:
            if flag.startswith("--how "):
                how_to_score = flag.split(" ")[1]
                break
        score_dir = call_dir / f"scored_{how_to_score}"
        score_flags += [f"--pre {postprocess_dir}", f"--out {score_dir}"]
        subprocess.run(f"python score_introgressions.py {' '.join(score_flags)}", shell=True, check=True)
    print("Introgressions analysis complete.")
    return


def run_introgression_sweep(call_flags, postprocess_flags, score_flags, output_dir):
    # Warn user if output directory already exists; ask permission to delete
    if output_dir.exists():
        response = input(f"Output directory {output_dir} already exists. Do you want to delete it? (y/n): ")
        if response.lower() == 'y':
            shutil.rmtree(output_dir)
            print(f"Deleted existing output directory: {output_dir}")
        else:
            print("Exiting without running the sweep.")
            return

    thresholds = [
        0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5,
        0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
    ]

    for thr in thresholds:
        thr_call_flags = call_flags.copy()
        thr_postprocess_flags = postprocess_flags.copy()
        thr_score_flags = score_flags.copy()

        call_dir = output_dir / f"{output_dir.name}_{thr}"
        call_dir.mkdir(parents=True, exist_ok=True)

        # Build the command to run the pipeline for this threshold
        cmd = []
        if thr_call_flags:
            thr_call_flags += [f"--out {call_dir}", f"--thr {thr}"]
            cmd += [f"python call_introgressions.py {' '.join(thr_call_flags)}"]
        if thr_postprocess_flags:
            postprocess_dir = call_dir / "postprocessed"
            thr_postprocess_flags += [f"--bed {call_dir}", f"--out {postprocess_dir}"]
            cmd += [f"python postprocess_introgressions.py {' '.join(thr_postprocess_flags)}"]
        if thr_score_flags:
            how_to_score = None
            for flag in thr_score_flags:
                if flag.startswith("--how "):
                    how_to_score = flag.split(" ")[1]
                    break
            score_dir = call_dir / f"scored_{how_to_score}"
            thr_score_flags += [f"--pre {postprocess_dir}", f"--out {score_dir}"]
            cmd += [f"python score_introgressions.py {' '.join(thr_score_flags)}"]

        cmd = " && ".join(cmd)
        cmd += f"; touch {call_dir / 'done.marker'}"

        log_file = call_dir / "sweep.log"
        screen_name = f"thr_{thr}"

        # Add a screen wrapper for the command
        # This allows the command to run in a detached screen session in parallel
        # Wrap cmd in single quotes so the full command is executed in the screen session
        screen_cmd = f"screen -d -m -L -Logfile {log_file} -S {screen_name} bash -c '{cmd}'"
        subprocess.run(screen_cmd, shell=True, check=False)

    # Wait for all done.marker files to exist
    for thr in thresholds:
        call_dir = output_dir / f"{output_dir.name}_{thr}"
        done_marker = call_dir / "done.marker"
        while not done_marker.exists():
            print(f"Waiting for {thr}...")
            time.sleep(15)
    print("All screens completed.")

    # Run visualization step after all thresholds are done
    vis_cmd = f"python visualize_introgressions.py -v prc prcc prca shtmp --dir {output_dir} --how bins"
    subprocess.run(vis_cmd, shell=True, check=True)
    print("Sweep complete.")
    return

if __name__ == "__main__":
    # Get config file path from command line argument
    if len(sys.argv) < 2:
        print("Usage: python introgression_runner.py <config.yaml> [--sweep]")
        sys.exit(1)

    config_path = Path(sys.argv[1])
    if not config_path.is_file():
        print(f"Config file {config_path} does not exist.")
        sys.exit(1)
    call_flags, postprocess_flags, score_flags, output_dir, call_thr = parse_config(config_path)

    if len(sys.argv) > 2:
        if sys.argv[2] == "--sweep":
            run_introgression_sweep(call_flags, postprocess_flags, score_flags, output_dir)
        else:
            print("Unknown command. Use '--sweep' to run pipeline for various thresholds or no argument to run the pipeline with one threshold.")
    else:
        run_introgression_pipeline(call_flags, postprocess_flags, score_flags, output_dir, call_thr)

import sys
from pathlib import Path
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import yaml


def parse_config(config_path):
    """Parse the introgression configuration file.

    Args:
        config_path (Path): Path to the YAML configuration file

    Returns:
        dict: parsed configuration as a dictionary
    """

    with config_path.open() as f:
        config = yaml.safe_load(f)

    general = config.get("general")
    calling = config.get("calling")
    postprocessing = config.get("postprocessing")
    scoring = config.get("scoring")

    output_dir = Path(general["output_dir"]).resolve()
    index_dir = Path(general["index_dir"]).resolve()
    tsv = Path(general["tsv"]).resolve()
    bin_size = general["bin"]
    ref = general["ref"]
    threads = general["threads"]

    call_run = calling["run"]
    call_rmf = calling["rmf"]
    call_rmu = calling["rmu"]
    call_ogrp = calling["ogrp"]
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
    call_edg = calling["edg"]
    call_vis = calling["vis"]

    post_run = postprocessing["run"]
    post_act = postprocessing["act"]
    post_min = postprocessing["min"]
    post_gap = postprocessing["gap"]
    post_map = postprocessing["map"]
    post_paf = postprocessing["paf"]

    score_gdt_dir = Path(scoring["gdt"]).resolve() if scoring["gdt"] else None
    score_run = scoring["run"]
    score_act = scoring["act"]
    score_min = scoring["min"]
    score_gap = scoring["gap"]
    score_thr = scoring["thr"]
    score_cmp = scoring["cmp"]
    score_vis = scoring["vis"]

    # Build flags for each step
    call_flags = None
    if call_run:
        if call_rmu is True:
            call_rmu = "--rmu true"
        else:
            call_rmu = f"--rmu {' '.join(call_rmu)}" if call_rmu else None

        call_flags = [
            f"--out {output_dir}",
            f"--idx {index_dir}",
            f"--tsv {tsv}",
            f"--grp {call_grp}" if call_grp else None,
            f"--anc {' '.join(call_anc)}" if call_anc else None,
            f"--chr {' '.join(call_chr)}" if call_chr else None,
            f"--cmp {' '.join(call_cmp)}",
            f"--bin {bin_size}",
            f"--stp {call_stp}",
            f"--gnm {call_gnm}" if call_gnm is not None else None,
            f"--trm {call_trm}" if call_trm is not None else None,
            f"--sft {call_sft}" if call_sft is not None else None,
            f"--ssz {call_ssz}" if call_ssz is not None else None,
            "--urf" if call_urf else None,
            call_rmu,
            f"--ogrp {' '.join(call_ogrp)}" if call_ogrp else None,
            f"--ref {ref}",
            "--rmf" if call_rmf else None,
            "--edg" if call_edg else None,
            "--vis" if call_vis else None,
            f"--threads {threads}",
        ]
        call_flags = [f for f in call_flags if f]

    postprocess_flags = None
    if post_run:
        postprocess_flags = [
            f"--idx {index_dir}",
            f"--act {' '.join(post_act)}" if post_act else None,
            f"--bin {bin_size}",
            f"--min {post_min}" if post_min is not None else None,
            f"--gap {post_gap}" if post_gap is not None else None,
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
            f"--thr {score_thr}" if score_thr is not None else None,
            f"--cmp {' '.join(score_cmp)}" if score_cmp else None,
            f"--act {' '.join(score_act)}" if score_act else None,
            f"--bin {bin_size}",
            f"--min {score_min}" if score_min is not None else None,
            f"--gap {score_gap}" if score_gap is not None else None,
            f"--ref {ref}",
            "--vis" if score_vis else None,
            f"--grp {tsv}" if score_vis else None,
        ]
        score_flags = [f for f in score_flags if f]

    return call_flags, postprocess_flags, score_flags, output_dir, call_thr, call_cmp, threads


def run_introgression_pipeline(
    call_flags, postprocess_flags, score_flags, output_dir, call_thr, call_cmp, threads, sweep
):
    """Run the introgression pipeline.

    Args:
        call_flags (list): flags for the calling step
        postprocess_flags (list): flags for the postprocessing step
        score_flags (list): flags for the scoring step
        output_dir (Path): directory for output files
        call_thr (list): list of thresholds to run for the calling step
        call_cmp (list): comparison groups for the calling step
        threads (int): number of threads to use
        sweep (bool): whether to run a sweep of thresholds
    """
    # decide what set if thresholds to use
    if sweep:
        if call_cmp == ["REF"]:
            print("Running sweep for 2-way comparison. Thresholds will be between 0 and 1.")
            call_thr = [
                0.1,
                0.15,
                0.2,
                0.25,
                0.3,
                0.35,
                0.4,
                0.45,
                0.5,
                0.55,
                0.6,
                0.65,
                0.7,
                0.75,
                0.8,
                0.85,
                0.9,
                0.95,
            ]

        else:
            print("Running sweep for 3-way comparison. Thresholds will be between 0 and 0.7.")
            call_thr = [
                0.0,
                0.04,
                0.08,
                0.12,
                0.16,
                0.2,
                0.24,
                0.28,
                0.32,
                0.36,
                0.4,
                0.44,
                0.48,
                0.52,
                0.56,
                0.6,
                0.64,
                0.68,
            ]

    # Call already runs chromosomes in parallel, so no additional parallelization is needed
    # Call also internally handles multiple thresholds
    call_thr = [float(thr) for thr in call_thr]  # ensure thresholds are in their float form
    if call_flags:
        call_flags += [f"--thr {' '.join(str(thr) for thr in call_thr)}"]
        subprocess.run(
            f"python call_introgressions.py {' '.join(call_flags)}", shell=True, check=True
        )

    def run_postprocess_and_score(thr):
        thr_postprocess_flags = postprocess_flags.copy() if postprocess_flags else None
        thr_score_flags = score_flags.copy() if score_flags else None

        call_dir = output_dir / f"{output_dir.name}_{thr}"
        if not call_dir.exists():
            raise ValueError(
                f"Expected call output directory {call_dir} does not exist. Check if the calling step completed successfully before running postprocessing."
            )

        postprocess_dir = call_dir / "postprocessed"
        log_file = call_dir / "caller.log"

        with log_file.open("w") as log_handle:
            if thr_postprocess_flags:
                thr_postprocess_flags += [f"--bed {call_dir / 'raw'}", f"--out {postprocess_dir}"]
                result = subprocess.run(
                    f"python postprocess_introgressions.py {' '.join(thr_postprocess_flags)}",
                    shell=True,
                    check=False,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                )
                log_handle.write(result.stdout)
                log_handle.write("\n")
                log_handle.flush()

                if result.returncode != 0:
                    raise subprocess.CalledProcessError(
                        result.returncode, result.args, result.stdout, result.stderr
                    )

            if thr_score_flags:
                score_dir = call_dir / "scored"
                thr_score_flags += [f"--pre {postprocess_dir}", f"--out {score_dir}"]
                result = subprocess.run(
                    f"python score_introgressions.py {' '.join(thr_score_flags)}",
                    shell=True,
                    check=False,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                )

                log_handle.write(result.stdout)
                log_handle.write("\n")
                log_handle.flush()

                if result.returncode != 0:
                    raise subprocess.CalledProcessError(
                        result.returncode, result.args, result.stdout, result.stderr
                    )

    # Run postprocessing and scoring steps for each threshold in parallel
    if postprocess_flags or score_flags:
        print("Running postprocessing and/or scoring steps for each threshold...")
        threads = min(threads, len(call_thr))

        with ThreadPoolExecutor(max_workers=threads) as executor:
            futures = {executor.submit(run_postprocess_and_score, thr): thr for thr in call_thr}

            for future in as_completed(futures):
                thr = futures[future]
                try:
                    future.result()  # raise any exception that occurred
                except subprocess.CalledProcessError as e:
                    print(f"An error occurred while processing threshold {thr}: {e}")
                    print(
                        "Check the caller.log files in each threshold's output directory for more details."
                    )
                    # Cancel remaining futures
                    for f in futures:
                        f.cancel()
                    sys.exit(1)

    # Run visualization step after all thresholds are done
    if score_flags and any(flag == "--vis" for flag in score_flags) and sweep:
        vis_cmd = f"python visualize_introgressions.py -v prc prcc prca mcc shtmp --dir {output_dir} --thresholds {' '.join(str(thr) for thr in call_thr)}"
        subprocess.run(vis_cmd, shell=True, check=True)

    print("Introgressions analysis complete.")
    return


if __name__ == "__main__":
    # Get config file path from command line argument
    if len(sys.argv) < 2:
        print("Usage: python introgression_runner.py <config.yaml> [--sweep]")
        sys.exit(0)

    sweep = False
    if len(sys.argv) > 2:
        if sys.argv[2] == "--sweep":
            sweep = True
        else:
            print(
                "Unknown command. Use '--sweep' to run pipeline for a preset list of thresholds or no argument to run the pipeline with your list of thresholds."
            )

    config_path = Path(sys.argv[1])
    if not config_path.is_file():
        print(f"Config file {config_path} does not exist.")
        sys.exit(1)
    call_flags, postprocess_flags, score_flags, output_dir, call_thr, call_cmp, threads = (
        parse_config(config_path)
    )

    # Check if the output directory already exists
    if output_dir.exists() and len(list(output_dir.iterdir())) > 1:
        response = input(
            f"Output directory {output_dir} already exists. Do you want to delete it? (y/n): "
        )
        if response.lower() == "y":
            shutil.rmtree(output_dir)
            print(f"Deleted existing output directory: {output_dir}")
        elif response.lower() == "n":
            print("Running pipeline with existing output directory. Results may be overwritten.")
        else:
            print("Invalid response. Exiting.")
            sys.exit(1)

    # Copy the config file into the output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    config_copy_path = output_dir / "intro_config.yaml"
    shutil.copy(config_path, config_copy_path)

    if (
        postprocess_flags
        and (len(call_thr) > 1 or sweep)
        and any(flag.startswith("--act") and "lift" in flag for flag in postprocess_flags)
    ):
        if not any(flag.startswith("--paf") for flag in postprocess_flags):
            raise ValueError(
                "PAF folder must be provided when using liftover with multiple thresholds. Try running the pipeline with a single threshold first to generate the PAF folder, or run alignment manually and specify the PAF folder."
            )

    run_introgression_pipeline(
        call_flags, postprocess_flags, score_flags, output_dir, call_thr, call_cmp, threads, sweep
    )

#!/usr/bin/env bash

# Usage: bash run_pipeline.sh [--check] [target.toml ...]
#
# Runs the whole pipeline: builds the shared CisBP reference, then calls TFs and
# assigns PWMs for each target organism. Targets are config files, given as a path or
# as a file name in config/; with none given, every config/*.toml except
# reference.toml is run. --check also runs check_sr_vs_cisbp.py, which compares the
# method with CisBP's own inferences.
#
# Run inside the conda environment (conda activate igem-pwm-assignment). Downloads
# happen once, into data/raw/; later runs only check their hashes. Each step logs to
# logs/, and a failing step stops the run.

set -euo pipefail

# reading the options; target paths are made absolute before moving to the repo root
run_check=false
targets=()
for arg in "$@"; do
    case "$arg" in
        --check) run_check=true ;;
        -h|--help) sed -n '3,13p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        -*) echo "unknown option: $arg (see --help)" >&2; exit 2 ;;
        *) if [ -f "$arg" ]; then targets+=("$(realpath "$arg")"); else targets+=("$arg"); fi ;;
    esac
done

cd "$(dirname "$(realpath "$0")")"

# every target config, if none were given
if [ ${#targets[@]} -eq 0 ]; then
    for config in config/*.toml; do
        if [ "$(basename "$config")" != "reference.toml" ]; then
            targets+=("$(basename "$config")")
        fi
    done
fi

# checking the tools before starting anything long
for tool in python hmmscan hmmalign hmmfetch hmmpress gffread curl git; do
    if ! command -v "$tool" > /dev/null; then
        echo "$tool not found; activate the environment first: conda activate igem-pwm-assignment" >&2
        exit 1
    fi
done
if ! python -c "import Bio, numpy, pandas" 2> /dev/null; then
    echo "python lacks biopython, numpy or pandas; activate the environment first: conda activate igem-pwm-assignment" >&2
    exit 1
fi

# the shared reference
python scripts/build_reference.py

# each target organism
for target in "${targets[@]}"; do
    python scripts/call_tfs.py --target "$target"
    python scripts/assign_motifs.py --target "$target"
done

# the optional method check
if [ "$run_check" = true ]; then
    python scripts/check_sr_vs_cisbp.py
fi

echo "done: results/ holds one folder per target; logs/ holds the logs"

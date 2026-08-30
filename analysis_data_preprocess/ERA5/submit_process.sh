#!/bin/bash
# Build one ERA5 time series per array task.
#
#   conda activate e3sm_diags_dev_py313
#   export ERA5_DIR=$HOME/e3sm_diags/analysis_data_preprocess/ERA5
#   cd $SCRATCH/analysis_data_e3sm_diags/ERA5_v2/slurm   # logs land here
#   sbatch --array=0-43%12 $ERA5_DIR/submit_process.sh
#
# The array index is a position in the `variables` table of era5_variables.yml,
# so the task list has to be renumbered if that table gains or loses an entry.
# List the current mapping with:
#
#   python -c "import yaml; print(*enumerate(yaml.safe_load(open('era5_variables.yml'))['variables']), sep='\n')"
#
# Splitting by variable rather than running the serial loop matters because the
# eight pressure-level fields (ta, ua, va, wap, hus, hur, zg, tro3) are ~212 GB
# of the 237 GB of raw input; back to back they dominate the wall clock, and one
# failure part way through costs the whole run. Per-variable tasks also make
# `--overwrite` unnecessary on a retry: `process` skips outputs already on disk.
#
# Unlike `submit_download.sh` this is real I/O and compute rather than waiting on
# ECMWF, so it runs in `shared` (charged, fractional nodes) rather than `xfer`.
# The raw files are chunked at ~8 MB, so `process` streams them through dask and
# stays well under the memory below; the request is sized for the 3D fields,
# which read ~23 GB and write about as much.

#SBATCH --job-name=era5_process
#SBATCH --qos=shared
#SBATCH --constraint=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=%x_%A_%a.log

set -euo pipefail

index=${SLURM_ARRAY_TASK_ID:?submit with --array=<index> or --array=<first>-<last>}

# Slurm runs a copy of this script from its spool directory, so `$0` cannot
# locate the pipeline. Submit from the ERA5 directory, or point ERA5_DIR at it.
era5_dir="${ERA5_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
if [[ ! -f "${era5_dir}/era5_pipeline.py" ]]; then
    echo "era5_pipeline.py not in ${era5_dir}: submit from the ERA5 directory or set ERA5_DIR" >&2
    exit 1
fi
cd "${era5_dir}"

start_year=${ERA5_START_YEAR:-1979}
end_year=${ERA5_END_YEAR:-2025}

variable=$(python - "${index}" <<'PY'
import sys
import yaml

names = list(yaml.safe_load(open("era5_variables.yml"))["variables"])
index = int(sys.argv[1])
if not 0 <= index < len(names):
    sys.exit(f"array index {index} is outside 0-{len(names) - 1}")
print(names[index])
PY
)

echo "host=$(hostname) index=${index} variable=${variable} start=$(date -Is)"
python -u era5_pipeline.py process \
    --start-year "${start_year}" --end-year "${end_year}" \
    --variables "${variable}"
echo "done=$(date -Is)"

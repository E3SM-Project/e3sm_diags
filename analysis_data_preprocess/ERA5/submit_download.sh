#!/bin/bash
# Retrieve one ERA5 year per array task.
#
#   conda activate e3sm_diags_dev_py313
#   export ERA5_DIR=$HOME/e3sm_diags/analysis_data_preprocess/ERA5
#   cd $SCRATCH/analysis_data_e3sm_diags/ERA5_v2/slurm   # logs land here
#   sbatch --array=1979-2025%3 $ERA5_DIR/submit_download.sh
#
# The array index is the year, so extending the record later is
# `sbatch --array=2026 submit_download.sh`. Retrievals already on disk are
# skipped, so a task that hits the wall clock can simply be resubmitted.
#
# This runs in the `xfer` QOS, which executes on a login node, is free of
# charge and allows 12 hours -- the right fit for a job that spends most of its
# wall time queued at ECMWF rather than computing. Note that xfer rejects
# `-N/--nodes`, defaults to 2 GB of memory, and is only for data movement: run
# `process` and `climo` under `--qos=shared` instead.
#
# A year is two retrievals, one per CDS dataset, but the CDS queues each one for
# as long as it likes -- waits of ~1-2 hours were seen in August 2026 -- so the
# wall clock is set generously.
#
# Concurrency is capped (%3) because the CDS applies per-user request limits,
# which bind well before the 15-job xfer limit does.

#SBATCH --job-name=era5_download
#SBATCH --qos=xfer
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%A_%a.log

set -euo pipefail

year=${SLURM_ARRAY_TASK_ID:?submit with --array=<year> or --array=<first>-<last>}

# Slurm exports the submitting environment, so the dev env must be active when
# `sbatch` runs rather than activated here.
if ! python -c "import cdsapi" 2>/dev/null; then
    echo "cdsapi not importable: activate the e3sm_diags dev env before submitting" >&2
    exit 1
fi

# Only compute nodes need the NERSC proxy to reach the CDS; a login node
# reaches it directly. Set ERA5_USE_PROXY=1 if this is ever moved to `shared`.
if [[ "${ERA5_USE_PROXY:-0}" == "1" ]]; then
    export https_proxy=http://proxy.nersc.gov:3128
    export http_proxy=http://proxy.nersc.gov:3128
fi

# Slurm runs a copy of this script from its spool directory, so `$0` cannot
# locate the pipeline. Submit from the ERA5 directory, or point ERA5_DIR at it.
era5_dir="${ERA5_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
if [[ ! -f "${era5_dir}/era5_pipeline.py" ]]; then
    echo "era5_pipeline.py not in ${era5_dir}: submit from the ERA5 directory or set ERA5_DIR" >&2
    exit 1
fi
cd "${era5_dir}"

echo "host=$(hostname) year=${year} start=$(date -Is)"
python -u era5_pipeline.py download --start-year "${year}" --end-year "${year}"
echo "done=$(date -Is)"

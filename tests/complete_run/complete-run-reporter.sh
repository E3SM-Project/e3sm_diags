#!/bin/bash
# NERSC scrontab reporter. It only inspects completed Slurm jobs and artifacts.
set -euo pipefail

if [[ $# -ne 1 ]]; then
    printf '%s\n' 'Usage: complete-run-reporter.sh /absolute/path/to/config.env' >&2
    exit 2
fi

source "$1"
: "${REPOSITORY:?}"
: "${CONDA_BASE:?}"
: "${CONTROLLER_ENV_PREFIX:?}"
: "${RESULTS_ROOT:?}"
: "${E3SM_DIAGS_REPOSITORY_ID:?}"
: "${E3SM_DIAGS_CATEGORY_ID:?}"
: "${E3SM_DIAGS_TOKEN_FILE:?}"

mkdir -p "$RESULTS_ROOT/automation"
exec 9>"$RESULTS_ROOT/automation/reporter.lock"
if ! flock -n 9; then
    printf '%s\n' 'A complete-run reporter is already active; exiting.'
    exit 0
fi

export CARTOPY_DATA_DIR="${CARTOPY_DATA_DIR:-}"
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$CONTROLLER_ENV_PREFIX"
cd "$REPOSITORY"
python -m tests.complete_run.reporter \
    --results-root "$RESULTS_ROOT" \
    --repository-id "$E3SM_DIAGS_REPOSITORY_ID" \
    --category-id "$E3SM_DIAGS_CATEGORY_ID" \
    --token-file "$E3SM_DIAGS_TOKEN_FILE"

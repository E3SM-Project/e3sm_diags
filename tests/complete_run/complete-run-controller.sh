#!/bin/bash
# NERSC scrontab controller. Its configuration file is maintained outside Git.
set -euo pipefail

if [[ $# -ne 1 ]]; then
    printf '%s\n' 'Usage: complete-run-controller.sh /absolute/path/to/config.env' >&2
    exit 2
fi

source "$1"
: "${REPOSITORY:?}"
: "${CONDA_BASE:?}"
: "${CONTROLLER_ENV_PREFIX:?}"
: "${RESULTS_ROOT:?}"
: "${PSCRATCH:?PSCRATCH is required for transient complete-run files}"
: "${SLURM_ACCOUNT:?}"

# Conda environments and detached worktrees contain many small files. Keep
# these disposable artifacts out of the constrained home filesystem by default.
WORKTREE_ROOT="${WORKTREE_ROOT:-$PSCRATCH/e3sm_diags/complete-run/worktrees}"
ENVIRONMENT_ROOT="${ENVIRONMENT_ROOT:-$PSCRATCH/e3sm_diags/complete-run/environments}"

# Standard cron cannot express every second Monday across month boundaries.
ISO_WEEK=$((10#$(TZ=America/Los_Angeles date +%V)))
if (( ISO_WEEK % 2 != 0 )); then
    printf '%s\n' 'Skipping odd ISO week; complete runs are biweekly.'
    exit 0
fi
if [[ "$(TZ=America/Los_Angeles date +%H)" != "06" ]]; then
    printf '%s\n' 'Skipping UTC schedule entry outside 06:00 Pacific.'
    exit 0
fi

# scrontab guidance requires controller jobs to discard inherited job settings
# before submitting a separate Perlmutter diagnostics allocation. Preserve the
# allocation policy explicitly supplied by the external controller config.
unset SLURM_MEM_PER_CPU SLURM_OPEN_MODE
for variable in "${!SLURM_@}"; do
    case "$variable" in
        SLURM_ACCOUNT|SLURM_CONSTRAINT|SLURM_NODES|SLURM_QOS|SLURM_WALLTIME) ;;
        *) unset "$variable" ;;
    esac
done

mkdir -p "$RESULTS_ROOT/automation"
exec 8>"$RESULTS_ROOT/automation/controller-environment.lock"
if ! flock -n 8; then
    printf '%s\n' 'A complete-run environment update is active; exiting.'
    exit 0
fi
exec 9>"$RESULTS_ROOT/automation/controller.lock"
if ! flock -n 9; then
    printf '%s\n' 'A complete-run controller is already active; exiting.'
    exit 0
fi

# The Cartopy activation hook appends to this variable. Initialize it before
# activating Conda because this controller intentionally uses ``set -u``.
export CARTOPY_DATA_DIR="${CARTOPY_DATA_DIR:-}"
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$CONTROLLER_ENV_PREFIX"
cd "$REPOSITORY"

python -m tests.complete_run.automation \
    --repo "$REPOSITORY" \
    --results-root "$RESULTS_ROOT" \
    --worktree-root "$WORKTREE_ROOT" \
    --environment-root "$ENVIRONMENT_ROOT" \
    --account "$SLURM_ACCOUNT" \
    --qos "${SLURM_QOS:-regular}" \
    --nodes "${SLURM_NODES:-1}" \
    --walltime "${SLURM_WALLTIME:-02:00:00}" \
    --constraint "${SLURM_CONSTRAINT:-cpu}"

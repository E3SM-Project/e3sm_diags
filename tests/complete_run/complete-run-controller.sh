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
: "${E3SM_DIAGS_REPOSITORY_ID:?}"
: "${E3SM_DIAGS_CATEGORY_ID:?}"
: "${E3SM_DIAGS_TOKEN_FILE:?}"

# Conda environments and detached worktrees contain many small files. Keep
# these disposable artifacts out of the constrained home filesystem by default.
WORKTREE_ROOT="${WORKTREE_ROOT:-$PSCRATCH/e3sm_diags/complete-run/worktrees}"
ENVIRONMENT_ROOT="${ENVIRONMENT_ROOT:-$PSCRATCH/e3sm_diags/complete-run/environments}"

# Standard cron cannot express every second Monday across month boundaries.
ISO_WEEK=$((10#$(date +%V)))
if (( ISO_WEEK % 2 != 0 )); then
    printf '%s\n' 'Skipping odd ISO week; complete runs are biweekly.'
    exit 0
fi

# scrontab guidance requires controller jobs to discard inherited job settings
# before submitting a separate Perlmutter diagnostics allocation.
unset SLURM_MEM_PER_CPU SLURM_OPEN_MODE
for variable in "${!SLURM_@}"; do
    unset "$variable"
done

mkdir -p "$RESULTS_ROOT/automation"
exec 9>"$RESULTS_ROOT/automation/controller.lock"
if ! flock -n 9; then
    printf '%s\n' 'A complete-run controller is already active; exiting.'
    exit 0
fi

source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$CONTROLLER_ENV_PREFIX"
cd "$REPOSITORY"

COMPLETION_FILE="$RESULTS_ROOT/automation/controller-completion.json"
rm -f "$COMPLETION_FILE"
AUTOMATION_EXIT=0
python -m tests.complete_run.automation \
    --repo "$REPOSITORY" \
    --results-root "$RESULTS_ROOT" \
    --worktree-root "$WORKTREE_ROOT" \
    --environment-root "$ENVIRONMENT_ROOT" \
    --account "$SLURM_ACCOUNT" \
    --qos "${SLURM_QOS:-regular}" \
    --walltime "${SLURM_WALLTIME:-01:00:00}" \
    --completion-file "$COMPLETION_FILE" || AUTOMATION_EXIT=$?

RUN_ROOT=$(python -c 'import json, sys; print(json.load(open(sys.argv[1]))["run_root"])' "$COMPLETION_FILE")
REPORTS=("$RUN_ROOT"/comparison/*/comparison-report.json)
if [[ -f "${REPORTS[0]}" ]]; then
    COMPARISON_REPORT="${REPORTS[0]}"
    RECEIPT="$(dirname "$COMPARISON_REPORT")/publication-receipt.json"
    PUBLISH_EXIT=0
    SHOULD_PUBLISH=$(python -c 'import json, sys; print(int(json.load(open(sys.argv[1]))["summary"]["failure_count"] > 0))' "$COMPARISON_REPORT")
    if [[ "$SHOULD_PUBLISH" -eq 1 ]]; then
        python -m tests.complete_run.report publish \
            --markdown "$RUN_ROOT/automation-report.md" \
            --receipt "$RECEIPT" \
            --repository-id "$E3SM_DIAGS_REPOSITORY_ID" \
            --category-id "$E3SM_DIAGS_CATEGORY_ID" \
            --token-file "$E3SM_DIAGS_TOKEN_FILE" || PUBLISH_EXIT=$?
    fi
    python -m tests.complete_run.report render \
        --status "$RUN_ROOT/status.json" \
        --comparison-report "$COMPARISON_REPORT" \
        --output-dir "$RUN_ROOT"
    if [[ "$AUTOMATION_EXIT" -eq 0 && "$PUBLISH_EXIT" -ne 0 ]]; then
        exit "$PUBLISH_EXIT"
    fi
fi

exit "$AUTOMATION_EXIT"

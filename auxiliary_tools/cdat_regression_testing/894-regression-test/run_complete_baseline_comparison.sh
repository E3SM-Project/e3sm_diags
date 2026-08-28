#!/usr/bin/env bash
# One-off pre-merge complete-run utility for HPC use.
#
# On NERSC, request an interactive Slurm allocation first:
#   salloc --nodes 1 --qos interactive --time 04:00:00 --constraint cpu --account=e3sm
#
# Then run:
#   bash auxiliary_tools/cdat_regression_testing/894-regression-test/run_complete_baseline_comparison.sh
#
# It creates separate Conda environments for main and the current branch when
# needed, then installs each checkout before running an apples-to-apples
# comparison. Existing environments are reused unchanged. Main package code
# uses this branch's complete-run workflow scripts.

set -euo pipefail

DEFAULT_RESULTS_ROOT="/global/cfs/cdirs/e3sm/www/e3sm_diags/complete-run-test"
main_ref="main"
results_root="$DEFAULT_RESULTS_ROOT"
num_workers=24
keep_main_worktree=false
main_worktree=""
worktree_created=false

_usage() {
    cat <<'EOF'
Usage: run_complete_baseline_comparison.sh [OPTIONS]

Run a complete diagnostics comparison between this non-main branch and main
package code using this branch's complete-run workflow scripts.

This command requires a Slurm compute allocation and Conda in PATH. It creates
revision-specific E3SM Diags Conda environments from each checkout's
conda-env/dev.yml when missing, then runs make install. Existing environments
are only activated and reused.

Options:
  --main-ref REF          Main revision to test (default: main)
  --results-root DIR      Existing root for CFS results
                          (default: /global/cfs/cdirs/e3sm/www/e3sm_diags/complete-run-test)
  --num-workers N         Diagnostics worker count (default: 24)
  --keep-main-worktree    Retain the temporary detached main worktree
  --help                  Show this help and exit
EOF
}

_die() {
    printf 'Error: %s\n' "$*" >&2
    exit 1
}

_activate_conda() {
    local conda_hook
    local restore_nounset=false

    conda_hook="$("$conda_exe" shell.bash hook)"
    eval "$conda_hook"
    if [[ $- == *u* ]]; then
        set +u
        restore_nounset=true
    fi
    if ! conda activate "$1"; then
        if [[ "$restore_nounset" == true ]]; then
            set -u
        fi
        return 1
    fi
    if [[ "$restore_nounset" == true ]]; then
        set -u
    fi
}

_prepare_environment() {
    local checkout=$1
    local environment=$2

    if "$conda_exe" env list | awk -v environment="$environment" \
        '$1 == environment || ($1 == "*" && $2 == environment) { found = 1 } END { exit !found }'; then
        printf 'Reusing Conda environment: %s\n' "$environment"
        _activate_conda "$environment"
        return
    fi
    (
        cd "$checkout"
        make env NAME="$environment"
    )
    _activate_conda "$environment"
    (
        cd "$checkout"
        make install
    )
}

_cleanup() {
    local status=$?

    if [[ "$worktree_created" == true && "$keep_main_worktree" == false ]]; then
        git -C "$source_root" worktree remove --force "$main_worktree" || \
            printf 'Warning: could not remove temporary worktree: %s\n' "$main_worktree" >&2
    elif [[ "$worktree_created" == true ]]; then
        printf 'Retained temporary main worktree: %s\n' "$main_worktree"
    elif [[ -n "$main_worktree" && -d "$main_worktree" ]]; then
        rmdir "$main_worktree" || \
            printf 'Warning: could not remove temporary directory: %s\n' "$main_worktree" >&2
    fi
    exit "$status"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --main-ref)
            [[ $# -ge 2 && -n "$2" ]] || _die "--main-ref requires a revision."
            main_ref="$2"
            shift 2
            ;;
        --results-root)
            [[ $# -ge 2 && -n "$2" ]] || _die "--results-root requires a directory."
            results_root="$2"
            shift 2
            ;;
        --num-workers)
            [[ $# -ge 2 && "$2" =~ ^[1-9][0-9]*$ ]] || \
                _die "--num-workers requires a positive integer."
            num_workers="$2"
            shift 2
            ;;
        --keep-main-worktree)
            keep_main_worktree=true
            shift
            ;;
        --help)
            _usage
            exit 0
            ;;
        *)
            _die "Unknown or malformed option: $1"
            ;;
    esac
done

[[ -n "${SLURM_JOB_ID:-}" ]] || _die "A Slurm compute allocation is required (SLURM_JOB_ID is unset)."
conda_exe="$(type -P conda || true)"
[[ -n "$conda_exe" && -x "$conda_exe" ]] || \
    _die "An executable Conda installation is required in PATH."

source_root="$(git rev-parse --show-toplevel 2>/dev/null)" || \
    _die "Run this script from inside a Git worktree."
if ! worktree_status="$(git -C "$source_root" status --porcelain)"; then
    _die "Unable to determine whether the current worktree is clean."
fi
[[ -z "$worktree_status" ]] || _die "The current worktree must be clean."

branch="$(git -C "$source_root" branch --show-current)"
[[ -n "$branch" ]] || _die "The current checkout is detached; a non-main branch is required."
[[ "$branch" != "main" ]] || _die "Run this utility from a non-main branch."
branch_sha="$(git -C "$source_root" rev-parse HEAD)"
branch_short_sha="$(git -C "$source_root" rev-parse --short HEAD)"

[[ -d "$results_root" && -w "$results_root" ]] || \
    _die "Results root does not exist, is not a directory, or is not writable: $results_root"
main_sha="$(git -C "$source_root" rev-parse --verify "${main_ref}^{commit}")" || \
    _die "Main ref does not resolve to a commit: $main_ref"
[[ "$main_sha" != "$branch_sha" ]] || \
    _die "Main ref resolves to the current branch commit; refusing self-comparison."

main_worktree="$(mktemp -d "${TMPDIR:-/tmp}/e3sm-diags-main-XXXXXX")"
trap _cleanup EXIT
# Older Git versions require the destination not to exist, even when it is empty.
rmdir "$main_worktree" || _die "Could not prepare temporary worktree path: $main_worktree"
git -C "$source_root" worktree add --detach "$main_worktree" "$main_sha"
worktree_created=true

# Overlay only workflow utilities. The e3sm_diags package remains at main.
git -C "$main_worktree" checkout "$branch_sha" -- tests/complete_run/
main_short_sha="$(git -C "$main_worktree" rev-parse --short HEAD)"
timestamp="$(date -u +%Y%m%d-%H%M%S)"
branch_label="$(printf '%s' "$branch" | tr -cs 'A-Za-z0-9._-' '-')"
main_conda_env="e3sm_diags_main_${main_short_sha}"
branch_conda_env="e3sm_diags_${branch_label}_${branch_short_sha}"
main_results="$results_root/${timestamp}-main-${main_short_sha}"
branch_results="$results_root/${timestamp}-${branch_label}-${branch_short_sha}"

for result_dir in "$main_results" "$branch_results"; do
    [[ ! -e "$result_dir" && ! -L "$result_dir" ]] || \
        _die "Refusing to overwrite existing results directory: $result_dir"
done

(
    cd "$main_worktree"
    _prepare_environment "$main_worktree" "$main_conda_env"
    python -m tests.complete_run.run \
        --results-dir "$main_results" \
        --num-workers "$num_workers" \
        --workflow-revision "$branch_sha"
)
(
    cd "$source_root"
    _prepare_environment "$source_root" "$branch_conda_env"
    python -m tests.complete_run.run \
        --results-dir "$branch_results" \
        --num-workers "$num_workers" \
        --workflow-revision "$branch_sha"
)
(
    cd "$source_root"
    _activate_conda "$branch_conda_env"
    python -m tests.complete_run.compare \
        --dev-dir "$branch_results" \
        --baseline-dir "$main_results" \
        --write-diff-pngs
)

printf '\nComparison succeeded.\n'
printf 'Main result:   %s\n' "$main_results"
printf 'Branch result: %s\n' "$branch_results"
printf '\nOptional reviewed-baseline promotion command (not run automatically):\n'
printf 'python -m tests.complete_run.baseline promote --run-dir %q --channel main\n' "$main_results"

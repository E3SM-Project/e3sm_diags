#!/usr/bin/env bash
# One-off login-node submitter for the complete-run pre-merge comparison.
#
# It submits main and branch diagnostics runs in parallel, then a dependent
# netCDF comparison. Activate no interactive allocation: Slurm owns all work.
#
# Example Usage
# -------------
# bash auxiliary_tools/cdat_regression_testing/894-regression-test/submit_complete_baseline_comparison.sh --conda-env ed_dev_894_baseline

set -euo pipefail

DEFAULT_RESULTS_ROOT="/global/cfs/cdirs/e3sm/www/e3sm_diags/complete-run-test"
main_ref="main"
results_root="$DEFAULT_RESULTS_ROOT"
num_workers=24
account="e3sm"
qos="regular"
walltime="04:00:00"
constraint="cpu"
nodes=1
conda_env="ed_dev_894_baseline"

_usage() {
    cat <<'EOF'
Usage: submit_complete_baseline_comparison.sh [OPTIONS]

Submit a one-off pre-merge complete-run comparison from a login node. Main and
current-branch runs are submitted in parallel; comparison runs after both pass.

Options:
  --main-ref REF          Main revision to test (default: main)
  --results-root DIR      Existing writable CFS results root
                          (default: /global/cfs/cdirs/e3sm/www/e3sm_diags/complete-run-test)
  --num-workers N         Diagnostics worker count (default: 24)
  --account ACCOUNT       Slurm account (default: e3sm)
  --qos QOS               Slurm QoS (default: regular)
  --time HH:MM:SS         Slurm walltime (default: 04:00:00)
  --constraint VALUE      Slurm node constraint (default: cpu)
  --nodes N               Slurm node count (default: 1)
  --conda-env NAME        Conda environment (default: ed_dev_894_baseline)
  --help                  Show this help and exit
EOF
}

_die() {
    printf 'Error: %s\n' "$*" >&2
    exit 1
}

_write_common_header() {
    cat <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

checkout=""
checkout_created=false
_activate_conda() {
    local conda_hook
    conda_hook="$("$conda_exe" shell.bash hook)"
    eval "$conda_hook"
    conda activate "$conda_env"
}

_cleanup() {
    local status=$?
    if [[ "$checkout_created" == true && -n "$checkout" ]]; then
        rm -rf -- "$checkout" || \
            printf 'Warning: could not remove temporary checkout: %s\n' "$checkout" >&2
    fi
    exit "$status"
}

_prepare_checkout() {
    local revision=$1
    checkout="$(mktemp -d "${TMPDIR:-/tmp}/e3sm-diags-complete-XXXXXX")"
    checkout_created=true
    trap _cleanup EXIT
    if ! git clone --no-checkout --shared "$source_root" "$checkout"; then
        printf 'Error: Git refused the empty temporary checkout directory: %s\n' \
            "$checkout" >&2
        return 1
    fi
    if ! git -C "$checkout" checkout --detach "$revision"; then
        printf 'Error: could not detached-checkout revision %s in %s\n' \
            "$revision" "$checkout" >&2
        return 1
    fi
}
EOF
}

_write_assignments() {
    printf 'source_root=%q\n' "$source_root"
    printf 'branch_sha=%q\n' "$branch_sha"
    printf 'main_sha=%q\n' "$main_sha"
    printf 'conda_exe=%q\n' "$conda_exe"
    printf 'conda_env=%q\n' "$conda_env"
    printf 'num_workers=%q\n' "$num_workers"
    printf 'main_results=%q\n' "$main_results"
    printf 'branch_results=%q\n' "$branch_results"
}

submitted_job_ids=()
submissions_complete=false

_cleanup_submissions() {
    local status=${1:-$?}

    if [[ "$submissions_complete" == false && ${#submitted_job_ids[@]} -gt 0 ]]; then
        printf 'Cancelling partially submitted jobs: %s\n' "${submitted_job_ids[*]}" >&2
        scancel "${submitted_job_ids[@]}" || \
            printf 'Warning: could not cancel all partially submitted jobs.\n' >&2
    fi
    trap - EXIT INT TERM
    exit "$status"
}

_submit_job() {
    local output
    local variable_name=$1
    local job_id

    shift
    if ! output="$(sbatch --parsable "$@")"; then
        _die "sbatch submission failed for ${variable_name}."
    fi
    if [[ ! "$output" =~ ^([0-9]+)(\;[A-Za-z0-9._-]+)?$ ]]; then
        _die "Unexpected sbatch --parsable output for ${variable_name}: ${output@Q}"
    fi
    job_id="${BASH_REMATCH[1]}"
    printf -v "$variable_name" '%s' "$job_id"
    submitted_job_ids+=("$job_id")
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
        --account|--qos|--constraint|--conda-env)
            [[ $# -ge 2 && -n "$2" ]] || _die "$1 requires a value."
            case "$1" in
                --account) account="$2" ;;
                --qos) qos="$2" ;;
                --constraint) constraint="$2" ;;
                --conda-env) conda_env="$2" ;;
            esac
            shift 2
            ;;
        --time)
            [[ $# -ge 2 && "$2" =~ ^[0-9]{2}:[0-5][0-9]:[0-5][0-9]$ ]] || \
                _die "--time requires HH:MM:SS."
            walltime="$2"
            shift 2
            ;;
        --nodes)
            [[ $# -ge 2 && "$2" =~ ^[1-9][0-9]*$ ]] || \
                _die "--nodes requires a positive integer."
            nodes="$2"
            shift 2
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

command -v sbatch >/dev/null 2>&1 || _die "sbatch is required on the login node."
command -v scancel >/dev/null 2>&1 || _die "scancel is required for transactional submission cleanup."
conda_exe="$(type -P conda || true)"
[[ -n "$conda_exe" && -x "$conda_exe" ]] || \
    _die "An executable Conda installation is required in PATH."
if ! "$conda_exe" run -n "$conda_env" python --version >/dev/null; then
    _die "Selected Conda environment is not usable via ${conda_exe}: ${conda_env}"
fi
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
main_sha="$(git -C "$source_root" rev-parse --verify "${main_ref}^{commit}")" || \
    _die "Main ref does not resolve to a commit: $main_ref"
[[ "$main_sha" != "$branch_sha" ]] || \
    _die "Main ref resolves to the current branch commit; refusing self-comparison."
[[ -d "$results_root" && -w "$results_root" ]] || \
    _die "Results root does not exist, is not a directory, or is not writable: $results_root"

timestamp="$(date -u +%Y%m%d-%H%M%S)"
branch_label="$(printf '%s' "$branch" | tr -cs 'A-Za-z0-9._-' '-')"
main_short_sha="${main_sha:0:12}"
job_dir="$(mktemp -d "$results_root/complete-baseline-jobs-${timestamp}-${branch_short_sha}-XXXXXX")"
run_token="${job_dir##*/}"
main_results="$results_root/${timestamp}-main-${main_short_sha}-${run_token}"
branch_results="$results_root/${timestamp}-${branch_label}-${branch_short_sha}-${run_token}"
main_script="$job_dir/run-main.sh"
branch_script="$job_dir/run-branch.sh"
compare_script="$job_dir/compare.sh"

{
    _write_common_header
    _write_assignments
    cat <<'EOF'
_prepare_checkout "$main_sha"
git -C "$checkout" checkout "$branch_sha" -- tests/complete_run/
cd "$checkout"
mkdir "$main_results"
_activate_conda
python -m tests.complete_run.run \
    --results-dir "$main_results" \
    --num-workers "$num_workers" \
    --workflow-revision "$branch_sha"
EOF
} >"$main_script"

{
    _write_common_header
    _write_assignments
    cat <<'EOF'
_prepare_checkout "$branch_sha"
cd "$checkout"
mkdir "$branch_results"
_activate_conda
python -m tests.complete_run.run \
    --results-dir "$branch_results" \
    --num-workers "$num_workers" \
    --workflow-revision "$branch_sha"
EOF
} >"$branch_script"

{
    _write_common_header
    _write_assignments
    cat <<'EOF'
_prepare_checkout "$branch_sha"
cd "$checkout"
_activate_conda
python -m tests.complete_run.compare \
    --dev-dir "$branch_results" \
    --baseline-dir "$main_results" \
    --write-diff-pngs
printf '\nComparison succeeded.\n'
printf 'Main result:   %s\n' "$main_results"
printf 'Branch result: %s\n' "$branch_results"
printf '\nOptional reviewed-baseline promotion command (not run automatically):\n'
printf 'python -m tests.complete_run.baseline promote --run-dir %q --channel main\n' "$main_results"
EOF
} >"$compare_script"

chmod +x "$main_script" "$branch_script" "$compare_script"
slurm_args=(
    --account="$account"
    --qos="$qos"
    --time="$walltime"
    --constraint="$constraint"
    --nodes="$nodes"
)
trap _cleanup_submissions EXIT
trap '_cleanup_submissions 130' INT TERM
_submit_job main_job_id "${slurm_args[@]}" --output="$job_dir/main-%j.out" "$main_script"
_submit_job branch_job_id "${slurm_args[@]}" --output="$job_dir/branch-%j.out" "$branch_script"
_submit_job compare_job_id "${slurm_args[@]}" \
    --dependency="afterok:${main_job_id}:${branch_job_id}" \
    --kill-on-invalid-dep=yes \
    --output="$job_dir/compare-%j.out" "$compare_script"
submissions_complete=true

printf 'Submitted main job:       %s\n' "$main_job_id"
printf 'Submitted branch job:     %s\n' "$branch_job_id"
printf 'Submitted comparison job: %s\n' "$compare_job_id"
printf 'Main result directory:    %s\n' "$main_results"
printf 'Branch result directory:  %s\n' "$branch_results"
printf 'Job scripts and logs:     %s\n' "$job_dir"
printf '\nMonitor jobs:\n'
printf 'squeue -j %s,%s,%s\n' "$main_job_id" "$branch_job_id" "$compare_job_id"
printf 'sacct -j %s,%s,%s --format=JobID,JobName,State,Elapsed,ExitCode\n' \
    "$main_job_id" "$branch_job_id" "$compare_job_id"

#!/usr/bin/env bash
# Submit one full E3SM Diags run for each Python environment created by
# 1_create_python_release_envs_xr2026070.sh. Run from the repository root:
#
#   bash auxiliary_tools/debug/1048-py314-stall-cont/investigations/260720-xarray-lock-stall-across-py-versions/3_run_python_release_full_lcrc_xr2026070.sh
#
# To render and validate the five batch scripts without submitting them:
#
#   bash auxiliary_tools/debug/1048-py314-stall-cont/investigations/260720-xarray-lock-stall-across-py-versions/3_run_python_release_full_lcrc_xr2026070.sh --dry-run
#
# Each generated script uses its matching Conda environment and a distinct
# results directory under /lcrc/group/e3sm/public_html/ac.tvo/. The NetCDF3
# file-lock workaround is disabled in every generated script.

set -euo pipefail

readonly TEMPLATE_ENV="ed_1048_xr_latest_2026070_py3143"
readonly RESULTS_ROOT="/lcrc/group/e3sm/public_html/ac.tvo"
readonly TARGET_ENVS=(
  "ed_1048_xr_2026070_py31314"
  "ed_1048_xr_2026070_py3143"
  "ed_1048_xr_2026070_py3144"
  "ed_1048_xr_2026070_py3145"
  "ed_1048_xr_2026070_py3146"
)

DRY_RUN=0
SCRIPT_DIR=""
REPO_ROOT=""
TEMPLATE_SCRIPT=""
RUN_DIR=""
SUMMARY_FILE=""


usage() {
  cat <<'EOF'
Usage: 3_run_python_release_full_lcrc_xr2026070.sh [--dry-run]

Render and submit the existing full E3SM Diags LCRC batch script for the five
Python environments using xarray=2026.7.0. Each generated script substitutes:
  - The Conda environment used by `conda activate`.
  - The environment-specific results directory.

The generated scripts retain:
  export E3SM_DIAGS_DISABLE_CLIMO_LOCK_WORKAROUND=1

Options:
  --dry-run     Render and validate scripts without submitting Slurm jobs.
  -h, --help    Show this help text.
EOF
}


log() {
  printf '[run_python_release_full_lcrc_xr2026070] %s\n' "$*"
}


die() {
  printf '[run_python_release_full_lcrc_xr2026070] Error: %s\n' "$*" >&2
  exit 1
}


require_cmd() {
  local cmd=$1
  command -v "${cmd}" >/dev/null 2>&1 || die "Required command not found: ${cmd}"
}


parse_args() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --dry-run)
        DRY_RUN=1
        shift
        ;;
      -h|--help)
        usage
        exit 0
        ;;
      *)
        die "Unknown argument: $1"
        ;;
    esac
  done
}


init_paths() {
  local timestamp

  SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
  REPO_ROOT=$(git -C "${SCRIPT_DIR}" rev-parse --show-toplevel)
  TEMPLATE_SCRIPT="${REPO_ROOT}/auxiliary_tools/debug/1048-py314-stall-cont/parallel-lcrc/bash-scripts/e3sm_diags_atm_monthly_180x360_aave_model_vs_obs_1985-2014.py314_ed_1048_xr_latest_2026070_py3143.bash"
  [[ -f "${TEMPLATE_SCRIPT}" ]] || die "Missing template script: ${TEMPLATE_SCRIPT}"

  timestamp=$(date +%Y%m%dT%H%M%S)
  RUN_DIR="${SCRIPT_DIR}/full-runs/python-release-xr2026070-lcrc-${timestamp}"
  SUMMARY_FILE="${RUN_DIR}/summary.tsv"
  mkdir -p "${RUN_DIR}/scripts"
}


env_exists() {
  local env_name=$1
  conda env list | awk 'NR > 2 {gsub(/\*/, "", $1); if ($1 != "") print $1}' | grep -Fxq "${env_name}"
}


validate_template() {
  local activation_count
  local results_count

  activation_count=$(grep -Fxc \
    "source \"\${HOME}/miniforge3/etc/profile.d/conda.sh\"; conda activate ${TEMPLATE_ENV}" \
    "${TEMPLATE_SCRIPT}")
  results_count=$(grep -Fxc \
    "results_dir=${RESULTS_ROOT}/${TEMPLATE_ENV}/\${tag}_\${Y1}-\${Y2}" \
    "${TEMPLATE_SCRIPT}")

  [[ ${activation_count} -eq 1 ]] \
    || die "Expected one template Conda activation, found ${activation_count}"
  [[ ${results_count} -eq 1 ]] \
    || die "Expected one template results_dir assignment, found ${results_count}"
  grep -Fqx "export E3SM_DIAGS_DISABLE_CLIMO_LOCK_WORKAROUND=1" \
    "${TEMPLATE_SCRIPT}" \
    || die "Template does not disable the file-lock workaround"
}


preflight() {
  local env_name

  require_cmd conda
  require_cmd git
  require_cmd grep
  require_cmd sed

  if [[ ${DRY_RUN} -eq 0 ]]; then
    require_cmd sbatch
  fi

  for env_name in "${TARGET_ENVS[@]}"; do
    env_exists "${env_name}" || die "Conda environment not found: ${env_name}"
  done

  validate_template
}


render_script() {
  local env_name=$1
  local rendered_script=$2
  local results_dir="${RESULTS_ROOT}/${env_name}/\${tag}_\${Y1}-\${Y2}"

  sed \
    -e "s#conda activate ${TEMPLATE_ENV}#conda activate ${env_name}#" \
    -e "s#results_dir=${RESULTS_ROOT}/${TEMPLATE_ENV}/\\\${tag}_\\\${Y1}-\\\${Y2}#results_dir=${results_dir}#" \
    "${TEMPLATE_SCRIPT}" > "${rendered_script}"

  grep -Fq "conda activate ${env_name}" "${rendered_script}" \
    || die "Failed to set Conda environment in ${rendered_script}"
  grep -Fqx "results_dir=${results_dir}" "${rendered_script}" \
    || die "Failed to set results_dir in ${rendered_script}"
  grep -Fqx "export E3SM_DIAGS_DISABLE_CLIMO_LOCK_WORKAROUND=1" \
    "${rendered_script}" \
    || die "File-lock workaround is not disabled in ${rendered_script}"
  if grep -Fq "${TEMPLATE_ENV}" "${rendered_script}" \
    && [[ "${env_name}" != "${TEMPLATE_ENV}" ]]; then
    die "Template environment remains in ${rendered_script}"
  fi

  chmod +x "${rendered_script}"
}


submit_runs() {
  local env_name
  local rendered_script
  local results_dir
  local provenance_log
  local job_id
  local result

  printf 'env\tresults_dir\tprovenance_log\tjob_id\tresult\ttotal_run_time\tresult_detail\tscript\n' > "${SUMMARY_FILE}"

  for env_name in "${TARGET_ENVS[@]}"; do
    rendered_script="${RUN_DIR}/scripts/e3sm_diags_full_${env_name}.bash"
    results_dir="${RESULTS_ROOT}/${env_name}/model_vs_obs_1985-2014_units"
    provenance_log="${results_dir}/prov/e3sm_diags_run.log"
    render_script "${env_name}" "${rendered_script}"

    if [[ ${DRY_RUN} -eq 1 ]]; then
      job_id="DRY_RUN"
      result="dry_run"
      log "Rendered ${env_name}"
    else
      job_id=$(sbatch --parsable "${rendered_script}")
      result="submitted"
      log "Submitted ${env_name} as job ${job_id}"
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${env_name}" "${results_dir}" "${provenance_log}" "${job_id}" \
      "${result}" "NA" "NA" "${rendered_script}" \
      >> "${SUMMARY_FILE}"
  done
}


main() {
  parse_args "$@"
  init_paths
  preflight

  log "Template: ${TEMPLATE_SCRIPT}"
  log "Run directory: ${RUN_DIR}"
  if [[ ${DRY_RUN} -eq 1 ]]; then
    log "Dry run: no jobs will be submitted"
  fi

  submit_runs
  printf '\nSummary file: %s\n' "${SUMMARY_FILE}"
}


main "$@"
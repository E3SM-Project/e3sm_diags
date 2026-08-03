#!/usr/bin/env bash
# Run the minimum QA reproduction on LCRC for each Python environment created
# by 1_create_python_release_envs_xr2026070.sh. Run from the repository root:
#
#   bash auxiliary_tools/debug/1048-py314-stall-cont/investigations/260720-xarray-lock-stall-across-py-versions/2_run_python_release_qa_lcrc_xr2026070.sh
#
# Defaults: five QA iterations, a 75-minute timeout, and five concurrent jobs.
# To override them, for example:
#
#   bash auxiliary_tools/debug/1048-py314-stall-cont/investigations/260720-xarray-lock-stall-across-py-versions/2_run_python_release_qa_lcrc_xr2026070.sh --repro-runs 3 --timeout 60m --max-active-jobs 3
#
# The runner can enable or disable the NetCDF3 file-lock workaround. It is
# disabled by default to preserve the original regression-test behavior.

set -euo pipefail

readonly XARRAY_VERSION="2026.7.0"
readonly DEFAULT_SRUN_TIME="01:10:00"
readonly DEFAULT_TIMEOUT="75m"
readonly DEFAULT_REPRO_RUNS="5"
readonly DEFAULT_MAX_ACTIVE_JOBS="10"
readonly DEFAULT_CLIMO_LOCK_WORKAROUND="disabled"
readonly TARGET_ENVS=(
  "ed_1048_xr_2026070_py31312"
  "ed_1048_xr_2026070_py31313"
  "ed_1048_xr_2026070_py31314"
  "ed_1048_xr_2026070_py3140"
  "ed_1048_xr_2026070_py3141"
  "ed_1048_xr_2026070_py3142"
  "ed_1048_xr_2026070_py3143"
  "ed_1048_xr_2026070_py3144"
  "ed_1048_xr_2026070_py3145"
  "ed_1048_xr_2026070_py3146"
)

ACCOUNT=""
SRUN_TIME="${DEFAULT_SRUN_TIME}"
RUN_TIMEOUT="${DEFAULT_TIMEOUT}"
REPRO_RUNS="${DEFAULT_REPRO_RUNS}"
MAX_ACTIVE_JOBS="${DEFAULT_MAX_ACTIVE_JOBS}"
CLIMO_LOCK_WORKAROUND="${DEFAULT_CLIMO_LOCK_WORKAROUND}"

SCRIPT_DIR=""
REPO_ROOT=""
QA_SCRIPT=""
CONDA_BASE=""
RUN_DIR=""
SUMMARY_FILE=""

declare -A STATUS_FILES=()
declare -A RUN_LOGS=()
declare -A SRUN_STDERR_LOGS=()
declare -A COMMAND_SCRIPTS=()
declare -A RUN_PIDS=()
declare -A FINAL_JOB_ID=()
declare -A FINAL_STATUS=()
declare -A FINAL_EXIT_CODE=()
declare -A FINAL_ERROR_NOTE=()


usage() {
  cat <<'EOF'
Usage: 2_run_python_release_qa_lcrc_xr2026070.sh [options]

Run the #1048 minimum QA reproduction with xarray=2026.7.0 in these Python
environments:
  - ed_1048_xr_2026070_py31312 (python=3.13.12)
  - ed_1048_xr_2026070_py31313 (python=3.13.13)
  - ed_1048_xr_2026070_py31314 (python=3.13.14)
  - ed_1048_xr_2026070_py3140  (python=3.14.0)
  - ed_1048_xr_2026070_py3141  (python=3.14.1)
  - ed_1048_xr_2026070_py3142  (python=3.14.2)
  - ed_1048_xr_2026070_py3143  (python=3.14.3)
  - ed_1048_xr_2026070_py3144  (python=3.14.4)
  - ed_1048_xr_2026070_py3145  (python=3.14.5)
  - ed_1048_xr_2026070_py3146  (python=3.14.6)

The environments run concurrently through `srun --pty`. Each generated command
sets E3SM_DIAGS_DISABLE_CLIMO_LOCK_WORKAROUND from the selected workaround mode.

Options:
  --account NAME          Slurm account to charge. Omit for site defaults.
  --time HH:MM:SS         Slurm walltime per allocation (default: 01:10:00)
  --timeout DURATION      QA timeout inside each allocation (default: 75m)
  --repro-runs N          Number of qa.py iterations per env (default: 5)
  --max-active-jobs N     Maximum concurrent allocations (default: 10)
  --climo-lock-workaround MODE
                          enabled or disabled (default: disabled)
  -h, --help              Show this help text.
EOF
}


log() {
  printf '[run_python_release_qa_lcrc_xr2026070] %s\n' "$*"
}


die() {
  printf '[run_python_release_qa_lcrc_xr2026070] Error: %s\n' "$*" >&2
  exit 1
}


require_cmd() {
  local cmd=$1
  command -v "${cmd}" >/dev/null 2>&1 || die "Required command not found: ${cmd}"
}


parse_args() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --account)
        [[ $# -ge 2 ]] || die "--account requires a value"
        ACCOUNT=$2
        shift 2
        ;;
      --time)
        [[ $# -ge 2 ]] || die "--time requires a value"
        SRUN_TIME=$2
        shift 2
        ;;
      --timeout)
        [[ $# -ge 2 ]] || die "--timeout requires a value"
        RUN_TIMEOUT=$2
        shift 2
        ;;
      --repro-runs)
        [[ $# -ge 2 ]] || die "--repro-runs requires a value"
        REPRO_RUNS=$2
        shift 2
        ;;
      --max-active-jobs)
        [[ $# -ge 2 ]] || die "--max-active-jobs requires a value"
        MAX_ACTIVE_JOBS=$2
        shift 2
        ;;
      --climo-lock-workaround)
        [[ $# -ge 2 ]] || die "--climo-lock-workaround requires a value"
        CLIMO_LOCK_WORKAROUND=$2
        shift 2
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


normalize_integer_arg() {
  local name=$1
  local value=$2

  [[ "${value}" =~ ^[0-9]+$ ]] || die "${name} must be an integer"
  [[ "${value}" -ge 1 ]] || die "${name} must be at least 1"
}


validate_climo_lock_workaround() {
  case "${CLIMO_LOCK_WORKAROUND}" in
    enabled|disabled) ;;
    *) die "--climo-lock-workaround must be enabled or disabled" ;;
  esac
}


climo_lock_workaround_override() {
  if [[ "${CLIMO_LOCK_WORKAROUND}" == "disabled" ]]; then
    printf '%s\n' "1"
  else
    printf '%s\n' "0"
  fi
}


init_paths() {
  local timestamp

  SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
  REPO_ROOT=$(git -C "${SCRIPT_DIR}" rev-parse --show-toplevel)
  QA_SCRIPT="${REPO_ROOT}/auxiliary_tools/debug/1048-py314-stall-cont/parallel-lcrc/min-scripts/qa.py"
  [[ -f "${QA_SCRIPT}" ]] || die "Missing QA script: ${QA_SCRIPT}"

  CONDA_BASE=$(conda info --base)
  [[ -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]] \
    || die "Missing conda activation script under ${CONDA_BASE}"

  timestamp=$(date +%Y%m%dT%H%M%S)
  RUN_DIR="${SCRIPT_DIR}/runs/python-release-xr2026070-lcrc-${timestamp}"
  SUMMARY_FILE="${RUN_DIR}/summary.tsv"
  mkdir -p "${RUN_DIR}/commands" "${RUN_DIR}/logs" \
    "${RUN_DIR}/srun" "${RUN_DIR}/status"
}


env_exists() {
  local env_name=$1
  conda env list | awk 'NR > 2 {gsub(/\*/, "", $1); if ($1 != "") print $1}' | grep -Fxq "${env_name}"
}


python_version_for_env() {
  local env_name=$1

  case "${env_name}" in
    ed_1048_xr_2026070_py31312) printf '%s\n' "3.13.12" ;;
    ed_1048_xr_2026070_py31313) printf '%s\n' "3.13.13" ;;
    ed_1048_xr_2026070_py31314) printf '%s\n' "3.13.14" ;;
    ed_1048_xr_2026070_py3140) printf '%s\n' "3.14.0" ;;
    ed_1048_xr_2026070_py3141) printf '%s\n' "3.14.1" ;;
    ed_1048_xr_2026070_py3142) printf '%s\n' "3.14.2" ;;
    ed_1048_xr_2026070_py3143) printf '%s\n' "3.14.3" ;;
    ed_1048_xr_2026070_py3144) printf '%s\n' "3.14.4" ;;
    ed_1048_xr_2026070_py3145) printf '%s\n' "3.14.5" ;;
    ed_1048_xr_2026070_py3146) printf '%s\n' "3.14.6" ;;
    *) die "Unknown environment: ${env_name}" ;;
  esac
}


preflight() {
  local env_name
  local host_name

  require_cmd conda
  require_cmd git
  require_cmd srun
  require_cmd timeout
  normalize_integer_arg "--repro-runs" "${REPRO_RUNS}"
  normalize_integer_arg "--max-active-jobs" "${MAX_ACTIVE_JOBS}"
  validate_climo_lock_workaround

  host_name=$(hostname -s 2>/dev/null || hostname)
  case "${host_name}" in
    chrlogin*|nid*|chrysalis*) ;;
    *)
      log "Host ${host_name} does not look like a Chrysalis/LCRC node."
      log "The srun --pty launch shape may not work on this host."
      ;;
  esac

  for env_name in "${TARGET_ENVS[@]}"; do
    env_exists "${env_name}" || die "Conda environment not found: ${env_name}"
  done
}


prepare_env_run() {
  local env_name=$1
  local python_version
  local climo_lock_workaround_override
  local command_script="${RUN_DIR}/commands/${env_name}.command.sh"
  local run_log="${RUN_DIR}/logs/${env_name}.run.log"
  local srun_stderr="${RUN_DIR}/srun/${env_name}.stderr.log"
  local status_file="${RUN_DIR}/status/${env_name}.env_status"

  python_version=$(python_version_for_env "${env_name}")
  climo_lock_workaround_override=$(climo_lock_workaround_override)
  COMMAND_SCRIPTS["${env_name}"]="${command_script}"
  RUN_LOGS["${env_name}"]="${run_log}"
  SRUN_STDERR_LOGS["${env_name}"]="${srun_stderr}"
  STATUS_FILES["${env_name}"]="${status_file}"

  cat > "${command_script}" <<EOF
#!/usr/bin/env bash

set -euo pipefail

readonly ENV_NAME="${env_name}"
readonly EXPECTED_PYTHON="${python_version}"
readonly EXPECTED_XARRAY="${XARRAY_VERSION}"
readonly REPO_ROOT="${REPO_ROOT}"
readonly QA_SCRIPT="${QA_SCRIPT}"
readonly RUN_LOG="${run_log}"
readonly STATUS_FILE="${status_file}"
readonly RUN_TIMEOUT="${RUN_TIMEOUT}"
readonly REPRO_RUNS="${REPRO_RUNS}"
readonly CLIMO_LOCK_WORKAROUND="${CLIMO_LOCK_WORKAROUND}"
readonly CLIMO_LOCK_WORKAROUND_OVERRIDE="${climo_lock_workaround_override}"

exec > >(tee -a "\${RUN_LOG}") 2>&1

source "${CONDA_BASE}/etc/profile.d/conda.sh"
set +u
conda activate "\${ENV_NAME}"
set -u
cd "\${REPO_ROOT}"

export CARTOPY_DATA_DIR="\${CARTOPY_DATA_DIR:-}"
export E3SM_DIAGS_DISABLE_CLIMO_LOCK_WORKAROUND="\${CLIMO_LOCK_WORKAROUND_OVERRIDE}"
export E3SM_DIAGS_REPRO_RUNS="\${REPRO_RUNS}"

echo "Started: \$(date --iso-8601=seconds)"
echo "Host: \$(hostname)"
echo "SLURM_JOB_ID: \${SLURM_JOB_ID:-unknown}"
echo "Conda env: \${CONDA_DEFAULT_ENV:-unknown}"
echo "QA script: \${QA_SCRIPT}"
echo "Timeout: \${RUN_TIMEOUT}"
echo "Climo lock workaround: \${CLIMO_LOCK_WORKAROUND}"
echo "E3SM_DIAGS_DISABLE_CLIMO_LOCK_WORKAROUND=\${E3SM_DIAGS_DISABLE_CLIMO_LOCK_WORKAROUND}"
echo "E3SM_DIAGS_REPRO_RUNS=\${E3SM_DIAGS_REPRO_RUNS}"

if ! python - <<'PY'
import platform
import sys

import xarray as xr

expected_python = "${python_version}"
expected_xarray = "${XARRAY_VERSION}"
actual_python = platform.python_version()
actual_xarray = xr.__version__

print(f"Python version: {actual_python}")
print(f"Python executable: {sys.executable}")
print(f"xarray version: {actual_xarray}")
print(f"xarray module path: {xr.__file__}")

if actual_python != expected_python:
    raise SystemExit(f"Expected Python {expected_python}, got {actual_python}")
if actual_xarray != expected_xarray:
    raise SystemExit(f"Expected xarray {expected_xarray}, got {actual_xarray}")
PY
then
  cat > "\${STATUS_FILE}" <<STATUS_EOF
env_name=\${ENV_NAME}
python_version=\${EXPECTED_PYTHON}
job_id=\${SLURM_JOB_ID:-}
status=failed
exit_code=2
run_log=\${RUN_LOG}
error_note=runtime version validation failed
ended_at=\$(date --iso-8601=seconds)
STATUS_EOF
  exit 2
fi

exit_code=0
status="completed"
timeout --signal=TERM --kill-after=2m "\${RUN_TIMEOUT}" \
  python "\${QA_SCRIPT}" || exit_code=\$?

case "\${exit_code}" in
  0) status="completed" ;;
  124|137) status="timeout_stall" ;;
  *) status="failed" ;;
esac

cat > "\${STATUS_FILE}" <<STATUS_EOF
env_name=\${ENV_NAME}
python_version=\${EXPECTED_PYTHON}
job_id=\${SLURM_JOB_ID:-}
status=\${status}
exit_code=\${exit_code}
run_log=\${RUN_LOG}
error_note=NA
ended_at=\$(date --iso-8601=seconds)
STATUS_EOF

echo "Final status: \${status}"
echo "Exit code: \${exit_code}"
exit "\${exit_code}"
EOF

  chmod +x "${command_script}"
}


extract_job_id() {
  local env_name=$1
  local run_log=${RUN_LOGS["${env_name}"]}
  local job_id=""

  if [[ -f "${run_log}" ]]; then
    job_id=$(grep -Eo 'SLURM_JOB_ID:[[:space:]]*[0-9]+' "${run_log}" \
      | tail -n 1 | grep -Eo '[0-9]+' | tail -n 1) || true
  fi
  printf '%s\n' "${job_id:-NA}"
}


srun_error_note() {
  local env_name=$1
  local stderr_log=${SRUN_STDERR_LOGS["${env_name}"]}
  local note=""

  if [[ -f "${stderr_log}" ]]; then
    note=$(tail -n 5 "${stderr_log}" \
      | tr '\n' ' ' | sed 's/[[:space:]]\+/ /g; s/^ //; s/ $//')
  fi
  printf '%s\n' "${note:-NA}"
}


run_env() {
  local env_name=$1
  local command_script=${COMMAND_SCRIPTS["${env_name}"]}
  local srun_stderr=${SRUN_STDERR_LOGS["${env_name}"]}
  local srun_args=(--pty --nodes=1 --time="${SRUN_TIME}")

  if [[ -n "${ACCOUNT}" ]]; then
    srun_args+=(--account="${ACCOUNT}")
  fi

  log "Starting ${env_name}"
  srun "${srun_args[@]}" /bin/bash "${command_script}" 2>"${srun_stderr}" || true
}


read_final_status() {
  local env_name=$1
  local status_file=${STATUS_FILES["${env_name}"]}
  local key
  local value

  FINAL_JOB_ID["${env_name}"]="NA"
  FINAL_STATUS["${env_name}"]="failed"
  FINAL_EXIT_CODE["${env_name}"]="NA"
  FINAL_ERROR_NOTE["${env_name}"]="NA"

  if [[ ! -f "${status_file}" ]]; then
    FINAL_JOB_ID["${env_name}"]=$(extract_job_id "${env_name}")
    FINAL_ERROR_NOTE["${env_name}"]=$(srun_error_note "${env_name}")
    if [[ ! -f "${RUN_LOGS[${env_name}]}" && -s "${SRUN_STDERR_LOGS[${env_name}]}" ]]; then
      FINAL_STATUS["${env_name}"]="submit_failed"
    fi
    return
  fi

  while IFS='=' read -r key value; do
    case "${key}" in
      job_id) FINAL_JOB_ID["${env_name}"]=${value:-NA} ;;
      status) FINAL_STATUS["${env_name}"]=${value:-failed} ;;
      exit_code) FINAL_EXIT_CODE["${env_name}"]=${value:-NA} ;;
      error_note) FINAL_ERROR_NOTE["${env_name}"]=${value:-NA} ;;
    esac
  done < "${status_file}"
}


run_batch() {
  local env_name
  local batch_envs=("$@")

  for env_name in "${batch_envs[@]}"; do
    prepare_env_run "${env_name}"
    run_env "${env_name}" &
    RUN_PIDS["${env_name}"]=$!
  done

  for env_name in "${batch_envs[@]}"; do
    wait "${RUN_PIDS[${env_name}]}" || true
    read_final_status "${env_name}"
    log "Finalized ${env_name}: ${FINAL_STATUS[${env_name}]}"
  done
}


run_all_envs() {
  local env_name
  local batch=()

  for env_name in "${TARGET_ENVS[@]}"; do
    batch+=("${env_name}")
    if [[ ${#batch[@]} -ge ${MAX_ACTIVE_JOBS} ]]; then
      run_batch "${batch[@]}"
      batch=()
    fi
  done

  if [[ ${#batch[@]} -gt 0 ]]; then
    run_batch "${batch[@]}"
  fi
}


write_summary() {
  local env_name

  {
    printf 'env\tpython\txarray\tclimo_lock_workaround\tjob_id\tstatus\texit_code\tlog_path\terror_note\n'
    for env_name in "${TARGET_ENVS[@]}"; do
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${env_name}" \
        "$(python_version_for_env "${env_name}")" \
        "${XARRAY_VERSION}" \
        "${CLIMO_LOCK_WORKAROUND}" \
        "${FINAL_JOB_ID[${env_name}]:-NA}" \
        "${FINAL_STATUS[${env_name}]:-unknown}" \
        "${FINAL_EXIT_CODE[${env_name}]:-NA}" \
        "${RUN_LOGS[${env_name}]:-NA}" \
        "${FINAL_ERROR_NOTE[${env_name}]:-NA}"
    done
  } > "${SUMMARY_FILE}"
}


print_summary() {
  local env_name

  printf '\n%-10s  %-12s  %-10s  %-14s  %-6s  %s\n' \
    "Python" "Xarray" "Job ID" "Status" "Exit" "Log"
  printf '%-10s  %-12s  %-10s  %-14s  %-6s  %s\n' \
    "----------" "------------" "----------" "--------------" "------" "---"

  for env_name in "${TARGET_ENVS[@]}"; do
    printf '%-10s  %-12s  %-10s  %-14s  %-6s  %s\n' \
      "$(python_version_for_env "${env_name}")" \
      "${XARRAY_VERSION}" \
      "${FINAL_JOB_ID[${env_name}]:-NA}" \
      "${FINAL_STATUS[${env_name}]:-unknown}" \
      "${FINAL_EXIT_CODE[${env_name}]:-NA}" \
      "${RUN_LOGS[${env_name}]:-NA}"
  done

  printf '\nSummary file: %s\n' "${SUMMARY_FILE}"
}


main() {
  local env_name
  local overall_exit=0

  parse_args "$@"
  require_cmd conda
  require_cmd git
  init_paths
  preflight

  log "Run directory: ${RUN_DIR}"
  log "xarray=${XARRAY_VERSION}, timeout=${RUN_TIMEOUT}, reproductions=${REPRO_RUNS}"
  log "Climo lock workaround: ${CLIMO_LOCK_WORKAROUND}"
  log "Submitting up to ${MAX_ACTIVE_JOBS} active allocation(s)"

  run_all_envs
  write_summary
  print_summary

  for env_name in "${TARGET_ENVS[@]}"; do
    if [[ "${FINAL_STATUS[${env_name}]}" != "completed" ]]; then
      overall_exit=1
      break
    fi
  done

  exit "${overall_exit}"
}


main "$@"

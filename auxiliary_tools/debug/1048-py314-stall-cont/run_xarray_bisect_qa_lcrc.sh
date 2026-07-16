#!/usr/bin/env bash
# bash auxiliary_tools/debug/1048-py314-stall-cont/run_xarray_bisect_qa_lcrc.sh

set -euo pipefail

readonly DEFAULT_SRUN_TIME="01:00:00"
readonly DEFAULT_TIMEOUT="30m"
readonly DEFAULT_REPRO_RUNS="10"

readonly TARGET_ENVS=(
  "ed_1048_xr_before_018ad08b"
  "ed_1048_xr_after_018ad08b"
  "ed_1048_xr_before_0a2d81c7"
  "ed_1048_xr_after_0a2d81c7"
)

ACCOUNT=""
SRUN_TIME="${DEFAULT_SRUN_TIME}"
RUN_TIMEOUT="${DEFAULT_TIMEOUT}"
REPRO_RUNS="${DEFAULT_REPRO_RUNS}"

SCRIPT_DIR=""
REPO_ROOT=""
QA_SCRIPT=""
RUNS_ROOT=""
RUN_DIR=""
SUMMARY_FILE=""

declare -A STATUS_FILES=()
declare -A RUN_LOGS=()
declare -A SRUN_STDERR_LOGS=()
declare -A COMMAND_SCRIPTS=()
declare -A FINAL_JOB_ID=()
declare -A FINAL_STATUS=()
declare -A FINAL_EXIT_CODE=()
declare -A FINAL_STALL_SIGNAL=()
declare -A FINAL_ERROR_NOTE=()


usage() {
  cat <<'EOF'
Usage: run_xarray_bisect_qa_lcrc.sh [options]

Run the 4 xarray bisect environments for the #1048 Chrysalis/LCRC minimum
repro using `srun --pty --nodes=1 --time=01:00:00 /bin/bash`.

Default launch shape per environment:
  srun --pty --nodes=1 --time=01:00:00 /bin/bash
  conda activate <env>
  python auxiliary_tools/debug/1048-py314-stall-cont/parallel-lcrc/min-scripts/qa.py

The runs execute sequentially so each `srun --pty` session has a clean
interactive allocation.

Options:
  --account NAME          Slurm account to charge. Omit to use site defaults.
  --time HH:MM:SS         Slurm walltime per `srun` allocation (default: 01:00:00)
  --timeout DURATION      Command timeout inside each allocation (default: 30m)
  --repro-runs N          E3SM_DIAGS_REPRO_RUNS value (default: 10)
  -h, --help              Show this help text.
EOF
}


log() {
  printf '[run_xarray_bisect_qa_lcrc] %s\n' "$*"
}


die() {
  printf '[run_xarray_bisect_qa_lcrc] Error: %s\n' "$*" >&2
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
  QA_SCRIPT="${REPO_ROOT}/auxiliary_tools/debug/1048-py314-stall-cont/parallel-lcrc/min-scripts/qa.py"
  [[ -f "${QA_SCRIPT}" ]] || die "Missing QA script: ${QA_SCRIPT}"

  RUNS_ROOT="${SCRIPT_DIR}/runs"
  timestamp=$(date +%Y%m%dT%H%M%S)
  RUN_DIR="${RUNS_ROOT}/xarray-bisect-lcrc-${timestamp}"
  SUMMARY_FILE="${RUN_DIR}/summary.tsv"

  mkdir -p "${RUN_DIR}/commands" "${RUN_DIR}/logs" "${RUN_DIR}/srun" "${RUN_DIR}/status"
}


normalize_integer_arg() {
  local name=$1
  local value=$2

  [[ "${value}" =~ ^[0-9]+$ ]] || die "${name} must be an integer"
  if [[ "${value}" -lt 1 ]]; then
    die "${name} must be at least 1"
  fi
}


env_exists() {
  local env_name=$1

  conda env list | awk 'NR > 2 {gsub(/\*/, "", $1); if ($1 != "") print $1}' | grep -Fxq "${env_name}"
}


preflight() {
  local env_name
  local host_name

  require_cmd bash
  require_cmd git
  require_cmd conda
  require_cmd srun
  require_cmd timeout

  normalize_integer_arg "--repro-runs" "${REPRO_RUNS}"

  host_name=$(hostname -s 2>/dev/null || hostname)
  case "${host_name}" in
    chrlogin*|nid*|chrysalis*)
      ;;
    *)
      log "Host ${host_name} does not look like a Chrysalis/LCRC login or compute node."
      log "If this is not Chrysalis, the default `srun --pty` flow may be invalid."
      ;;
  esac

  for env_name in "${TARGET_ENVS[@]}"; do
    env_exists "${env_name}" || die "Conda environment not found: ${env_name}"
  done
}


status_file_for_env() {
  local env_name=$1
  printf '%s/status/%s.env_status\n' "${RUN_DIR}" "${env_name}"
}


run_log_for_env() {
  local env_name=$1
  printf '%s/logs/%s.run.log\n' "${RUN_DIR}" "${env_name}"
}


srun_stderr_for_env() {
  local env_name=$1
  printf '%s/srun/%s.stderr.log\n' "${RUN_DIR}" "${env_name}"
}


command_script_for_env() {
  local env_name=$1
  printf '%s/commands/%s.command.sh\n' "${RUN_DIR}" "${env_name}"
}


render_command_script() {
  local env_name=$1
  local command_script=$2
  local run_log=$3
  local status_file=$4

  cat > "${command_script}" <<EOF
#!/usr/bin/env bash

set -euo pipefail
set -o pipefail

readonly ENV_NAME="${env_name}"
readonly REPO_ROOT="${REPO_ROOT}"
readonly QA_SCRIPT="${QA_SCRIPT}"
readonly RUN_LOG="${run_log}"
readonly STATUS_FILE="${status_file}"
readonly RUN_TIMEOUT="${RUN_TIMEOUT}"
readonly REPRO_RUNS="${REPRO_RUNS}"

mkdir -p "\$(dirname "\${RUN_LOG}")" "\$(dirname "\${STATUS_FILE}")"

source "\${HOME}/miniforge3/etc/profile.d/conda.sh"
conda activate "${env_name}"
cd "\${REPO_ROOT}"

export E3SM_DIAGS_DISABLE_CLIMO_LOCK_WORKAROUND=1
export E3SM_DIAGS_REPRO_RUNS="\${E3SM_DIAGS_REPRO_RUNS:-\${REPRO_RUNS}}"

exec > >(tee -a "\${RUN_LOG}") 2>&1

echo "Started: \$(date --iso-8601=seconds)"
echo "Host: \$(hostname)"
echo "SLURM_JOB_ID: \${SLURM_JOB_ID:-unknown}"
echo "Conda env: \${CONDA_DEFAULT_ENV:-unknown}"
echo "QA script: \${QA_SCRIPT}"
echo "Timeout: \${RUN_TIMEOUT}"
echo "E3SM_DIAGS_DISABLE_CLIMO_LOCK_WORKAROUND=\${E3SM_DIAGS_DISABLE_CLIMO_LOCK_WORKAROUND}"
echo "E3SM_DIAGS_REPRO_RUNS=\${E3SM_DIAGS_REPRO_RUNS}"

exit_code=0
status="completed"

timeout --signal=TERM --kill-after=2m "\${RUN_TIMEOUT}" python "\${QA_SCRIPT}" || exit_code=\$?

case "\${exit_code}" in
  0)
    status="completed"
    ;;
  124|137)
    status="timeout_stall"
    ;;
  *)
    status="failed"
    ;;
esac

cat > "\${STATUS_FILE}" <<STATUS_EOF
env_name=\${ENV_NAME}
job_id=\${SLURM_JOB_ID:-}
status=\${status}
exit_code=\${exit_code}
run_log=\${RUN_LOG}
ended_at=\$(date --iso-8601=seconds)
STATUS_EOF

echo "Final status: \${status}"
echo "Exit code: \${exit_code}"
exit "\${exit_code}"
EOF

  chmod +x "${command_script}"
}


read_status_file() {
  local env_name=$1
  local status_file=${STATUS_FILES["${env_name}"]}
  local key
  local value
  local job_id=""
  local status=""
  local exit_code=""

  [[ -f "${status_file}" ]] || return 1

  while IFS='=' read -r key value; do
    case "${key}" in
      job_id)
        job_id=${value}
        ;;
      status)
        status=${value}
        ;;
      exit_code)
        exit_code=${value}
        ;;
    esac
  done < "${status_file}"

  [[ -n "${status}" ]] || return 1

  FINAL_JOB_ID["${env_name}"]="${job_id:-NA}"
  FINAL_STATUS["${env_name}"]="${status}"
  FINAL_EXIT_CODE["${env_name}"]="${exit_code:-NA}"
  FINAL_ERROR_NOTE["${env_name}"]="NA"
  return 0
}


extract_job_id_from_log() {
  local env_name=$1
  local log_path=${RUN_LOGS["${env_name}"]}
  local job_id

  if [[ -f "${log_path}" ]]; then
    job_id=$(
      grep -Eo 'SLURM_JOB_ID:[[:space:]]*[0-9]+' "${log_path}" \
        | tail -n 1 \
        | grep -Eo '[0-9]+' \
        | tail -n 1
    ) || true
    if [[ -n "${job_id:-}" ]]; then
      printf '%s\n' "${job_id}"
      return
    fi
  fi

  printf 'NA\n'
}


extract_srun_error_note() {
  local env_name=$1
  local stderr_log=${SRUN_STDERR_LOGS["${env_name}"]}
  local note=""

  if [[ -f "${stderr_log}" ]]; then
    note=$(tail -n 5 "${stderr_log}" | tr '\n' ' ' | sed 's/[[:space:]]\+/ /g; s/^ //; s/ $//')
  fi

  printf '%s\n' "${note:-NA}"
}


detect_stall_signal() {
  local log_path=$1
  local last_start
  local last_done

  [[ -f "${log_path}" ]] || {
    printf 'no\n'
    return
  }

  last_start=$(
    grep -nE 'Climo backend (open_mfdataset|open_dataset) start' "${log_path}" \
      | tail -n 1 \
      | cut -d: -f1
  ) || true

  last_done=$(
    grep -nE 'Climo backend (open_mfdataset|open_dataset)( retry drop_time)? done' "${log_path}" \
      | tail -n 1 \
      | cut -d: -f1
  ) || true

  if [[ -n "${last_start}" && ( -z "${last_done}" || "${last_start}" -gt "${last_done}" ) ]]; then
    printf 'yes\n'
  else
    printf 'no\n'
  fi
}


run_env() {
  local env_name=$1
  local command_script
  local run_log
  local status_file
  local srun_stderr
  local srun_exit=0
  local srun_args=(--pty --nodes=1 --time="${SRUN_TIME}")

  if [[ -n "${ACCOUNT}" ]]; then
    srun_args+=(--account="${ACCOUNT}")
  fi

  command_script=$(command_script_for_env "${env_name}")
  run_log=$(run_log_for_env "${env_name}")
  status_file=$(status_file_for_env "${env_name}")
  srun_stderr=$(srun_stderr_for_env "${env_name}")

  STATUS_FILES["${env_name}"]="${status_file}"
  RUN_LOGS["${env_name}"]="${run_log}"
  SRUN_STDERR_LOGS["${env_name}"]="${srun_stderr}"
  COMMAND_SCRIPTS["${env_name}"]="${command_script}"

  render_command_script "${env_name}" "${command_script}" "${run_log}" "${status_file}"

  log "Starting ${env_name}"
  if srun "${srun_args[@]}" /bin/bash "${command_script}" 2>"${srun_stderr}"; then
    srun_exit=0
  else
    srun_exit=$?
  fi

  if ! read_status_file "${env_name}"; then
    FINAL_JOB_ID["${env_name}"]="$(extract_job_id_from_log "${env_name}")"
    FINAL_EXIT_CODE["${env_name}"]="${srun_exit}"
    FINAL_ERROR_NOTE["${env_name}"]="$(extract_srun_error_note "${env_name}")"
    if [[ "${srun_exit}" == "124" || "${srun_exit}" == "137" ]]; then
      FINAL_STATUS["${env_name}"]="timeout_stall"
    elif [[ ! -f "${run_log}" && -s "${srun_stderr}" ]]; then
      FINAL_STATUS["${env_name}"]="submit_failed"
    else
      FINAL_STATUS["${env_name}"]="failed"
    fi
  fi

  if [[ "${FINAL_STATUS[${env_name}]}" == "timeout_stall" ]]; then
    FINAL_STALL_SIGNAL["${env_name}"]="$(detect_stall_signal "${run_log}")"
  else
    FINAL_STALL_SIGNAL["${env_name}"]="no"
  fi

  log "Finalized ${env_name}: ${FINAL_STATUS[${env_name}]}"
}


write_summary() {
  local env_name

  {
    printf 'env\tjob_id\tstatus\texit_code\tstall_signal\tlog_path\terror_note\n'
    for env_name in "${TARGET_ENVS[@]}"; do
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${env_name}" \
        "${FINAL_JOB_ID[${env_name}]:-NA}" \
        "${FINAL_STATUS[${env_name}]:-unknown}" \
        "${FINAL_EXIT_CODE[${env_name}]:-NA}" \
        "${FINAL_STALL_SIGNAL[${env_name}]:-no}" \
        "${RUN_LOGS[${env_name}]:-NA}" \
        "${FINAL_ERROR_NOTE[${env_name}]:-NA}"
    done
  } > "${SUMMARY_FILE}"
}


print_summary() {
  local env_name

  printf '\n'
  printf '%-32s  %-10s  %-14s  %-9s  %-12s  %s\n' \
    "Environment" "Job ID" "Status" "Exit" "Stall Signal" "Log"
  printf '%-32s  %-10s  %-14s  %-9s  %-12s  %s\n' \
    "--------------------------------" "----------" "--------------" "---------" "------------" "---"

  for env_name in "${TARGET_ENVS[@]}"; do
    printf '%-32s  %-10s  %-14s  %-9s  %-12s  %s\n' \
      "${env_name}" \
      "${FINAL_JOB_ID[${env_name}]:-NA}" \
      "${FINAL_STATUS[${env_name}]:-unknown}" \
      "${FINAL_EXIT_CODE[${env_name}]:-NA}" \
      "${FINAL_STALL_SIGNAL[${env_name}]:-no}" \
      "${RUN_LOGS[${env_name}]:-NA}"
  done

  printf '\nSummary file: %s\n' "${SUMMARY_FILE}"
  printf 'If any run is `submit_failed`, check the per-env files under %s/srun/.\n' "${RUN_DIR}"
}


main() {
  local env_name
  local overall_exit=0

  parse_args "$@"
  init_paths
  preflight

  log "Run directory: ${RUN_DIR}"
  if [[ -n "${ACCOUNT}" ]]; then
    log "Using srun --pty with account=${ACCOUNT}, time=${SRUN_TIME}"
  else
    log "Using srun --pty with site-default account, time=${SRUN_TIME}"
  fi
  log "Per-run timeout=${RUN_TIMEOUT}, E3SM_DIAGS_REPRO_RUNS=${REPRO_RUNS}"

  for env_name in "${TARGET_ENVS[@]}"; do
    run_env "${env_name}"
  done

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

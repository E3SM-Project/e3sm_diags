#!/usr/bin/env bash
# bash auxiliary_tools/debug/1048-py314-stall-cont/run_xarray_bisect_qa_nersc.sh

set -euo pipefail

readonly DEFAULT_ACCOUNT="e3sm"
readonly DEFAULT_CONSTRAINT="cpu"
readonly DEFAULT_SLURM_TIME="02:00:00"
readonly DEFAULT_TIMEOUT="30m"
readonly DEFAULT_MAX_ACTIVE_JOBS=3
readonly POLL_INTERVAL_SECONDS=15

readonly TARGET_ENVS=(
  "ed_1048_xr_before_018ad08b"
  "ed_1048_xr_after_018ad08b"
  "ed_1048_xr_before_0a2d81c7"
  "ed_1048_xr_after_0a2d81c7"
)

ACCOUNT="${DEFAULT_ACCOUNT}"
CONSTRAINT="${DEFAULT_CONSTRAINT}"
SLURM_TIME="${DEFAULT_SLURM_TIME}"
RUN_TIMEOUT="${DEFAULT_TIMEOUT}"
MAX_ACTIVE_JOBS="${DEFAULT_MAX_ACTIVE_JOBS}"

SCRIPT_DIR=""
REPO_ROOT=""
QA_SCRIPT=""
RUNS_ROOT=""
RUN_DIR=""
SUMMARY_FILE=""

declare -A JOB_IDS=()
declare -A JOB_ACTIVE=()
declare -A JOB_FINALIZED=()
declare -A STATUS_FILES=()
declare -A LOG_PATHS=()
declare -A SLURM_STDOUT_PATHS=()
declare -A SLURM_STDERR_PATHS=()
declare -A BATCH_SCRIPTS=()
declare -A FINAL_STATUS=()
declare -A FINAL_EXIT_CODE=()
declare -A FINAL_STALL_SIGNAL=()


usage() {
  cat <<'EOF'
Usage: run_xarray_bisect_qa_nersc.sh [options]

Submit and monitor the 4 xarray bisect environments for the #1048 NERSC
parallel minimum repro while respecting the 3-job concurrency limit.

Options:
  --account NAME          Slurm account to charge (default: e3sm)
  --constraint VALUE      Slurm constraint value (default: cpu)
  --time HH:MM:SS         Slurm walltime per job (default: 02:00:00)
  --timeout DURATION      Command timeout inside each job (default: 30m)
  --max-active-jobs N     Max concurrent jobs; capped at 3 (default: 3)
  -h, --help              Show this help text.
EOF
}


log() {
  printf '[run_xarray_bisect_qa_nersc] %s\n' "$*"
}


die() {
  printf '[run_xarray_bisect_qa_nersc] Error: %s\n' "$*" >&2
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
      --constraint)
        [[ $# -ge 2 ]] || die "--constraint requires a value"
        CONSTRAINT=$2
        shift 2
        ;;
      --time)
        [[ $# -ge 2 ]] || die "--time requires a value"
        SLURM_TIME=$2
        shift 2
        ;;
      --timeout)
        [[ $# -ge 2 ]] || die "--timeout requires a value"
        RUN_TIMEOUT=$2
        shift 2
        ;;
      --max-active-jobs)
        [[ $# -ge 2 ]] || die "--max-active-jobs requires a value"
        MAX_ACTIVE_JOBS=$2
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
  QA_SCRIPT="${REPO_ROOT}/auxiliary_tools/debug/1048-py314-stall-cont/parallel-nersc/min-scripts/qa.py"
  [[ -f "${QA_SCRIPT}" ]] || die "Missing QA script: ${QA_SCRIPT}"

  RUNS_ROOT="${SCRIPT_DIR}/runs"
  timestamp=$(date +%Y%m%dT%H%M%S)
  RUN_DIR="${RUNS_ROOT}/xarray-bisect-${timestamp}"
  SUMMARY_FILE="${RUN_DIR}/summary.tsv"

  mkdir -p "${RUN_DIR}/batch" "${RUN_DIR}/logs" "${RUN_DIR}/slurm" "${RUN_DIR}/status"
}


normalize_max_active_jobs() {
  [[ "${MAX_ACTIVE_JOBS}" =~ ^[0-9]+$ ]] || die "--max-active-jobs must be an integer"
  if [[ "${MAX_ACTIVE_JOBS}" -lt 1 ]]; then
    die "--max-active-jobs must be at least 1"
  fi
  if [[ "${MAX_ACTIVE_JOBS}" -gt 3 ]]; then
    log "Capping --max-active-jobs from ${MAX_ACTIVE_JOBS} to 3 to respect NERSC limits"
    MAX_ACTIVE_JOBS=3
  fi
}


env_exists() {
  local env_name=$1

  conda env list | awk 'NR > 2 {gsub(/\*/, "", $1); if ($1 != "") print $1}' | grep -Fxq "${env_name}"
}


preflight() {
  local env_name

  require_cmd bash
  require_cmd git
  require_cmd conda
  require_cmd sbatch
  require_cmd squeue
  require_cmd sacct
  require_cmd timeout

  normalize_max_active_jobs

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


slurm_stdout_for_env() {
  local env_name=$1
  printf '%s/slurm/%s.stdout.log\n' "${RUN_DIR}" "${env_name}"
}


slurm_stderr_for_env() {
  local env_name=$1
  printf '%s/slurm/%s.stderr.log\n' "${RUN_DIR}" "${env_name}"
}


batch_script_for_env() {
  local env_name=$1
  printf '%s/batch/%s.sbatch.sh\n' "${RUN_DIR}" "${env_name}"
}


render_batch_script() {
  local env_name=$1
  local batch_script=$2
  local run_log=$3
  local status_file=$4

  cat > "${batch_script}" <<EOF
#!/usr/bin/env bash

set -euo pipefail
set -o pipefail

readonly ENV_NAME="${env_name}"
readonly REPO_ROOT="${REPO_ROOT}"
readonly QA_SCRIPT="${QA_SCRIPT}"
readonly RUN_LOG="${run_log}"
readonly STATUS_FILE="${status_file}"
readonly RUN_TIMEOUT="${RUN_TIMEOUT}"

mkdir -p "\$(dirname "\${RUN_LOG}")" "\$(dirname "\${STATUS_FILE}")"

source "\${HOME}/miniforge3/etc/profile.d/conda.sh"
conda activate "${env_name}"
cd "\${REPO_ROOT}"

export E3SM_DIAGS_DISABLE_CLIMO_LOCK_WORKAROUND=1
export E3SM_DIAGS_REPRO_RUNS="\${E3SM_DIAGS_REPRO_RUNS:-10}"

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

  chmod +x "${batch_script}"
}


submit_env_job() {
  local env_name=$1
  local batch_script
  local status_file
  local run_log
  local slurm_stdout
  local slurm_stderr
  local job_id
  local job_name

  batch_script=$(batch_script_for_env "${env_name}")
  status_file=$(status_file_for_env "${env_name}")
  run_log=$(run_log_for_env "${env_name}")
  slurm_stdout=$(slurm_stdout_for_env "${env_name}")
  slurm_stderr=$(slurm_stderr_for_env "${env_name}")
  job_name="xr1048_${env_name#ed_1048_xr_}"

  STATUS_FILES["${env_name}"]="${status_file}"
  LOG_PATHS["${env_name}"]="${run_log}"
  SLURM_STDOUT_PATHS["${env_name}"]="${slurm_stdout}"
  SLURM_STDERR_PATHS["${env_name}"]="${slurm_stderr}"
  BATCH_SCRIPTS["${env_name}"]="${batch_script}"
  JOB_FINALIZED["${env_name}"]=0

  render_batch_script "${env_name}" "${batch_script}" "${run_log}" "${status_file}"

  if job_id=$(sbatch \
    --parsable \
    --job-name "${job_name}" \
    --account "${ACCOUNT}" \
    --constraint "${CONSTRAINT}" \
    --time "${SLURM_TIME}" \
    --output "${slurm_stdout}" \
    --error "${slurm_stderr}" \
    "${batch_script}"); then
    JOB_IDS["${env_name}"]="${job_id}"
    JOB_ACTIVE["${env_name}"]=1
    log "Submitted ${env_name} as job ${job_id}"
  else
    JOB_IDS["${env_name}"]=""
    JOB_ACTIVE["${env_name}"]=0
    JOB_FINALIZED["${env_name}"]=1
    FINAL_STATUS["${env_name}"]="submit_failed"
    FINAL_EXIT_CODE["${env_name}"]="NA"
    FINAL_STALL_SIGNAL["${env_name}"]="no"
    log "Submission failed for ${env_name}"
  fi
}


active_job_count() {
  local env_name
  local count=0

  for env_name in "${TARGET_ENVS[@]}"; do
    if [[ "${JOB_ACTIVE[${env_name}]:-0}" == "1" ]]; then
      ((count+=1))
    fi
  done

  printf '%s\n' "${count}"
}


read_status_file() {
  local env_name=$1
  local status_file=${STATUS_FILES["${env_name}"]}
  local key
  local value
  local status=""
  local exit_code=""

  [[ -f "${status_file}" ]] || return 1

  while IFS='=' read -r key value; do
    case "${key}" in
      status)
        status=${value}
        ;;
      exit_code)
        exit_code=${value}
        ;;
    esac
  done < "${status_file}"

  [[ -n "${status}" ]] || return 1

  FINAL_STATUS["${env_name}"]="${status}"
  FINAL_EXIT_CODE["${env_name}"]="${exit_code:-NA}"
  return 0
}


map_sacct_state_to_status() {
  local state=$1

  case "${state}" in
    COMPLETED)
      printf 'completed\n'
      ;;
    TIMEOUT)
      printf 'timeout_stall\n'
      ;;
    *)
      printf 'failed\n'
      ;;
  esac
}


read_sacct_fallback() {
  local env_name=$1
  local job_id=${JOB_IDS["${env_name}"]}
  local sacct_line
  local state
  local exit_code

  [[ -n "${job_id}" ]] || return 1

  sacct_line=$(sacct -n -P -X -j "${job_id}" -o State,ExitCode | head -n 1 | tr -d ' ')
  [[ -n "${sacct_line}" ]] || return 1

  state=${sacct_line%%|*}
  exit_code=${sacct_line##*|}
  [[ -n "${state}" ]] || return 1

  state=${state%%+*}
  state=${state%% *}

  FINAL_STATUS["${env_name}"]="$(map_sacct_state_to_status "${state}")"
  FINAL_EXIT_CODE["${env_name}"]="${exit_code%%:*}"
  return 0
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


finalize_env_if_ready() {
  local env_name=$1

  if [[ "${JOB_FINALIZED[${env_name}]:-0}" == "1" ]]; then
    return
  fi

  if read_status_file "${env_name}" || read_sacct_fallback "${env_name}"; then
    if [[ "${FINAL_STATUS[${env_name}]}" == "timeout_stall" ]]; then
      FINAL_STALL_SIGNAL["${env_name}"]="$(detect_stall_signal "${LOG_PATHS[${env_name}]}")"
    else
      FINAL_STALL_SIGNAL["${env_name}"]="no"
    fi
    JOB_FINALIZED["${env_name}"]=1
    JOB_ACTIVE["${env_name}"]=0
    log "Finalized ${env_name}: ${FINAL_STATUS[${env_name}]}"
  fi
}


refresh_job_states() {
  local env_name
  local job_id
  local squeue_output

  for env_name in "${TARGET_ENVS[@]}"; do
    if [[ -z "${JOB_IDS[${env_name}]:-}" || "${JOB_FINALIZED[${env_name}]:-0}" == "1" ]]; then
      continue
    fi

    job_id=${JOB_IDS["${env_name}"]}
    squeue_output=$(squeue -h -j "${job_id}" -o '%A')

    if [[ -n "${squeue_output}" ]]; then
      JOB_ACTIVE["${env_name}"]=1
    else
      JOB_ACTIVE["${env_name}"]=0
      finalize_env_if_ready "${env_name}"
    fi
  done
}


all_jobs_finalized() {
  local env_name

  for env_name in "${TARGET_ENVS[@]}"; do
    if [[ "${JOB_FINALIZED[${env_name}]:-0}" != "1" ]]; then
      return 1
    fi
  done

  return 0
}


write_summary() {
  local env_name

  {
    printf 'env\tjob_id\tstatus\texit_code\tstall_signal\tlog_path\n'
    for env_name in "${TARGET_ENVS[@]}"; do
      printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${env_name}" \
        "${JOB_IDS[${env_name}]:-NA}" \
        "${FINAL_STATUS[${env_name}]:-unknown}" \
        "${FINAL_EXIT_CODE[${env_name}]:-NA}" \
        "${FINAL_STALL_SIGNAL[${env_name}]:-no}" \
        "${LOG_PATHS[${env_name}]:-NA}"
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
      "${JOB_IDS[${env_name}]:-NA}" \
      "${FINAL_STATUS[${env_name}]:-unknown}" \
      "${FINAL_EXIT_CODE[${env_name}]:-NA}" \
      "${FINAL_STALL_SIGNAL[${env_name}]:-no}" \
      "${LOG_PATHS[${env_name}]:-NA}"
  done

  printf '\nSummary file: %s\n' "${SUMMARY_FILE}"
}


submit_and_monitor_jobs() {
  local next_index=0
  local total_envs=${#TARGET_ENVS[@]}

  while (( next_index < total_envs )) || ! all_jobs_finalized; do
    refresh_job_states

    while (( next_index < total_envs )) && (( $(active_job_count) < MAX_ACTIVE_JOBS )); do
      submit_env_job "${TARGET_ENVS[$next_index]}"
      ((next_index+=1))
      refresh_job_states
    done

    if (( next_index < total_envs )) || ! all_jobs_finalized; then
      sleep "${POLL_INTERVAL_SECONDS}"
    fi
  done
}


main() {
  local env_name
  local overall_exit=0

  parse_args "$@"
  init_paths
  preflight

  log "Run directory: ${RUN_DIR}"
  log "Submitting up to ${MAX_ACTIVE_JOBS} active job(s) at a time"
  log "Per-job Slurm time=${SLURM_TIME}, command timeout=${RUN_TIMEOUT}"

  for env_name in "${TARGET_ENVS[@]}"; do
    JOB_ACTIVE["${env_name}"]=0
    JOB_FINALIZED["${env_name}"]=0
  done

  submit_and_monitor_jobs
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

#!/usr/bin/env bash
# Run the minimum #1048 reproduction in one Python environment with temporary
# faulthandler, RSS, and GC telemetry injected through sitecustomize.py.
#
# Default target:
#   ed_1048_xr_2026070_py3144
#
# Run from the repository root:
#   bash auxiliary_tools/debug/1048-py314-stall-cont/investigations/260720-xarray-lock-stall-across-py-versions/4_run_py3144_gc_trace_min_lcrc.sh
#
# Useful variants:
#   bash auxiliary_tools/debug/1048-py314-stall-cont/investigations/260720-xarray-lock-stall-across-py-versions/4_run_py3144_gc_trace_min_lcrc.sh --env ed_1048_xr_2026070_py3143
#   bash auxiliary_tools/debug/1048-py314-stall-cont/investigations/260720-xarray-lock-stall-across-py-versions/4_run_py3144_gc_trace_min_lcrc.sh --timeout 75m --trace-interval 60 --rss-interval 15

set -euo pipefail

readonly XARRAY_VERSION="2026.7.0"
readonly DEFAULT_ENV_NAME="ed_1048_xr_2026070_py3144"
readonly DEFAULT_SRUN_TIME="01:10:00"
readonly DEFAULT_TIMEOUT="60m"
readonly DEFAULT_REPRO_RUNS="3"
readonly DEFAULT_TRACE_INTERVAL_SECONDS="120"
readonly DEFAULT_RSS_INTERVAL_SECONDS="30"

ACCOUNT=""
ENV_NAME="${DEFAULT_ENV_NAME}"
SRUN_TIME="${DEFAULT_SRUN_TIME}"
RUN_TIMEOUT="${DEFAULT_TIMEOUT}"
REPRO_RUNS="${DEFAULT_REPRO_RUNS}"
TRACE_INTERVAL_SECONDS="${DEFAULT_TRACE_INTERVAL_SECONDS}"
RSS_INTERVAL_SECONDS="${DEFAULT_RSS_INTERVAL_SECONDS}"

SCRIPT_DIR=""
REPO_ROOT=""
QA_SCRIPT=""
CONDA_BASE=""
RUN_DIR=""
COMMAND_SCRIPT=""
RUN_LOG=""
SRUN_STDERR_LOG=""
SUMMARY_FILE=""
STATUS_FILE=""
PROBE_DIR=""


usage() {
  cat <<'EOF'
Usage: 4_run_py3144_gc_trace_min_lcrc.sh [options]

Run the minimum E3SM Diags #1048 QA reproduction with Python-level telemetry:
  - parent faulthandler traceback dumps every N seconds
  - parent VmRSS/VmHWM samples from /proc/self/status
  - gc.get_count(), gc.get_threshold(), and compact gc.get_stats() samples
  - forked-worker PID breadcrumbs and inherited SIGUSR1 traceback dumps

The script defaults to ed_1048_xr_2026070_py3144 because Python 3.14.4 stalls
in the July 20 matrix. Use --env to compare another pre-created environment.

Options:
  --env NAME               Conda environment to run (default: py3144 env)
  --account NAME           Slurm account to charge. Omit for site defaults.
  --time HH:MM:SS          Slurm walltime for the allocation (default: 01:10:00)
  --timeout DURATION       QA timeout inside the allocation (default: 60m)
  --repro-runs N           Number of qa.py iterations (default: 1)
  --trace-interval SEC     faulthandler dump interval; 0 disables (default: 120)
  --rss-interval SEC       RSS/GC sample interval; 0 disables (default: 30)
  -h, --help               Show this help text.
EOF
}


log() {
  printf '[run_py3144_gc_trace_min_lcrc] %s\n' "$*"
}


die() {
  printf '[run_py3144_gc_trace_min_lcrc] Error: %s\n' "$*" >&2
  exit 1
}


require_cmd() {
  local cmd=$1
  command -v "${cmd}" >/dev/null 2>&1 || die "Required command not found: ${cmd}"
}


parse_args() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --env)
        [[ $# -ge 2 ]] || die "--env requires a value"
        ENV_NAME=$2
        shift 2
        ;;
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
      --trace-interval)
        [[ $# -ge 2 ]] || die "--trace-interval requires a value"
        TRACE_INTERVAL_SECONDS=$2
        shift 2
        ;;
      --rss-interval)
        [[ $# -ge 2 ]] || die "--rss-interval requires a value"
        RSS_INTERVAL_SECONDS=$2
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
}


python_version_for_env() {
  case "$1" in
    ed_1048_xr_2026070_py31314) printf '%s\n' "3.13.14" ;;
    ed_1048_xr_2026070_py3143) printf '%s\n' "3.14.3" ;;
    ed_1048_xr_2026070_py3144) printf '%s\n' "3.14.4" ;;
    ed_1048_xr_2026070_py3145) printf '%s\n' "3.14.5" ;;
    ed_1048_xr_2026070_py3146) printf '%s\n' "3.14.6" ;;
    *) die "Unknown Python-release test environment: $1" ;;
  esac
}


env_exists() {
  local env_name=$1
  conda env list | awk 'NR > 2 {gsub(/\*/, "", $1); if ($1 != "") print $1}' | grep -Fxq "${env_name}"
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
  RUN_DIR="${SCRIPT_DIR}/confirmation-runs/gc-trace-min-lcrc-${timestamp}-${ENV_NAME}"
  COMMAND_SCRIPT="${RUN_DIR}/commands/${ENV_NAME}.command.sh"
  RUN_LOG="${RUN_DIR}/logs/${ENV_NAME}.run.log"
  SRUN_STDERR_LOG="${RUN_DIR}/srun/${ENV_NAME}.stderr.log"
  SUMMARY_FILE="${RUN_DIR}/summary.tsv"
  STATUS_FILE="${RUN_DIR}/status/${ENV_NAME}.env_status"
  PROBE_DIR="${RUN_DIR}/sitecustomize"

  mkdir -p "${RUN_DIR}/commands" "${RUN_DIR}/logs" "${RUN_DIR}/srun" \
    "${RUN_DIR}/status" "${PROBE_DIR}"
}


preflight() {
  local host_name

  require_cmd conda
  require_cmd git
  require_cmd srun
  require_cmd timeout
  normalize_integer_arg "--repro-runs" "${REPRO_RUNS}"
  normalize_integer_arg "--trace-interval" "${TRACE_INTERVAL_SECONDS}"
  normalize_integer_arg "--rss-interval" "${RSS_INTERVAL_SECONDS}"

  [[ "${REPRO_RUNS}" -ge 1 ]] || die "--repro-runs must be at least 1"
  env_exists "${ENV_NAME}" || die "Conda environment not found: ${ENV_NAME}"

  host_name=$(hostname -s 2>/dev/null || hostname)
  case "${host_name}" in
    chrlogin*|nid*|chrysalis*) ;;
    *)
      log "Host ${host_name} does not look like a Chrysalis/LCRC node."
      log "The srun --pty launch shape may not work on this host."
      ;;
  esac
}


write_command_script() {
  local expected_python

  expected_python=$(python_version_for_env "${ENV_NAME}")

  cat > "${COMMAND_SCRIPT}" <<EOF
#!/usr/bin/env bash

set -euo pipefail

readonly ENV_NAME="${ENV_NAME}"
readonly EXPECTED_PYTHON="${expected_python}"
readonly EXPECTED_XARRAY="${XARRAY_VERSION}"
readonly REPO_ROOT="${REPO_ROOT}"
readonly QA_SCRIPT="${QA_SCRIPT}"
readonly RUN_LOG="${RUN_LOG}"
readonly STATUS_FILE="${STATUS_FILE}"
readonly PROBE_DIR="${PROBE_DIR}"
readonly RUN_TIMEOUT="${RUN_TIMEOUT}"
readonly REPRO_RUNS="${REPRO_RUNS}"
readonly TRACE_INTERVAL_SECONDS="${TRACE_INTERVAL_SECONDS}"
readonly RSS_INTERVAL_SECONDS="${RSS_INTERVAL_SECONDS}"

exec > >(tee -a "\${RUN_LOG}") 2>&1

source "${CONDA_BASE}/etc/profile.d/conda.sh"
set +u
conda activate "\${ENV_NAME}"
set -u
cd "\${REPO_ROOT}"

mkdir -p "\${PROBE_DIR}"
cat > "\${PROBE_DIR}/sitecustomize.py" <<'PY'
from __future__ import annotations

import faulthandler
import gc
import os
import signal
import sys
import threading
import time

_TRACE_INTERVAL = int(os.environ.get("E3SM_DIAGS_FAULTHANDLER_INTERVAL", "120"))
_RSS_INTERVAL = int(os.environ.get("E3SM_DIAGS_RSS_INTERVAL", "30"))
_PREFIX = "[py314-gc-trace]"
_ARMED_PIDS: set[int] = set()
_PROCESS_START_BY_PID: dict[int, float] = {}


def _proc_status_values() -> dict[str, str]:
    values: dict[str, str] = {}
    try:
        with open("/proc/self/status", encoding="utf-8") as status_file:
            for line in status_file:
                if line.startswith(("VmRSS:", "VmHWM:", "VmSize:", "Threads:")):
                    key, value = line.split(":", 1)
                    values[key] = " ".join(value.split())
    except OSError as exc:
        values["status_error"] = repr(exc)
    return values


def _compact_gc_stats() -> str:
    parts = []
    for generation, stats in enumerate(gc.get_stats()):
        collected = stats.get("collected", "NA")
        collections = stats.get("collections", "NA")
        uncollectable = stats.get("uncollectable", "NA")
        parts.append(
            f"g{generation}:collections={collections},"
            f"collected={collected},uncollectable={uncollectable}"
        )
    return ";".join(parts)


def _sampler() -> None:
    pid = os.getpid()
    while True:
        status = _proc_status_values()
        elapsed = time.monotonic() - _PROCESS_START_BY_PID.get(pid, time.monotonic())
        print(
            _PREFIX,
            f"elapsed_s={elapsed:.1f}",
            f"pid={pid}",
            f"ppid={os.getppid()}",
            f"VmRSS={status.get('VmRSS', 'NA')}",
            f"VmHWM={status.get('VmHWM', 'NA')}",
            f"VmSize={status.get('VmSize', 'NA')}",
            f"Threads={status.get('Threads', 'NA')}",
            f"gc_count={gc.get_count()}",
            f"gc_threshold={gc.get_threshold()}",
            f"gc_stats={_compact_gc_stats()}",
            flush=True,
            file=sys.stderr,
        )
        time.sleep(_RSS_INTERVAL)


def _enable_for_current_process(reason: str) -> None:
    pid = os.getpid()
    if pid in _ARMED_PIDS:
        return

    _ARMED_PIDS.add(pid)
    _PROCESS_START_BY_PID[pid] = time.monotonic()
    faulthandler.enable(all_threads=True)

    try:
        faulthandler.register(signal.SIGUSR1, file=sys.stderr, all_threads=True, chain=False)
        signal_note = "sigusr1=registered"
    except (AttributeError, OSError, RuntimeError, ValueError) as exc:
        signal_note = f"sigusr1=unavailable:{type(exc).__name__}"

    print(
        _PREFIX,
        f"enabled reason={reason}",
        f"pid={pid}",
        f"ppid={os.getppid()}",
        f"trace_interval_s={_TRACE_INTERVAL}",
        f"rss_interval_s={_RSS_INTERVAL}",
        f"gc_threshold={gc.get_threshold()}",
        signal_note,
        flush=True,
        file=sys.stderr,
    )

    if _TRACE_INTERVAL > 0:
        faulthandler.dump_traceback_later(
            _TRACE_INTERVAL,
            repeat=True,
            file=sys.stderr,
        )

    if _RSS_INTERVAL > 0:
        threading.Thread(
            target=_sampler,
            name=f"e3sm-diags-gc-rss-sampler-{pid}",
            daemon=True,
        ).start()


def _after_fork_in_child() -> None:
    # Keep the fork hook tiny. Starting Python threads or timers here can
    # perturb ProcessPool workers before they start running diagnostics.
    pid = os.getpid()
    message = (
        f"{_PREFIX} forked_child pid={pid} ppid={os.getppid()} "
        "sigusr1=inherited\n"
    )
    try:
        os.write(2, message.encode("utf-8", errors="replace"))
    except OSError:
        pass


_enable_for_current_process("site_import")

if hasattr(os, "register_at_fork"):
    os.register_at_fork(after_in_child=_after_fork_in_child)
PY

export CARTOPY_DATA_DIR="\${CARTOPY_DATA_DIR:-}"
export E3SM_DIAGS_DISABLE_CLIMO_LOCK_WORKAROUND=1
export E3SM_DIAGS_REPRO_RUNS="\${REPRO_RUNS}"
export E3SM_DIAGS_FAULTHANDLER_INTERVAL="\${TRACE_INTERVAL_SECONDS}"
export E3SM_DIAGS_RSS_INTERVAL="\${RSS_INTERVAL_SECONDS}"
export PYTHONPATH="\${PROBE_DIR}:\${PYTHONPATH:-}"

echo "Started: \$(date --iso-8601=seconds)"
echo "Host: \$(hostname)"
echo "SLURM_JOB_ID: \${SLURM_JOB_ID:-unknown}"
echo "Conda env: \${CONDA_DEFAULT_ENV:-unknown}"
echo "QA script: \${QA_SCRIPT}"
echo "Probe dir: \${PROBE_DIR}"
echo "Timeout: \${RUN_TIMEOUT}"
echo "Trace interval: \${TRACE_INTERVAL_SECONDS}s"
echo "RSS interval: \${RSS_INTERVAL_SECONDS}s"
echo "Immediate stack dump: scancel --signal=USR1 \${SLURM_JOB_ID:-<job_id>}"
echo "E3SM_DIAGS_DISABLE_CLIMO_LOCK_WORKAROUND=\${E3SM_DIAGS_DISABLE_CLIMO_LOCK_WORKAROUND}"
echo "E3SM_DIAGS_REPRO_RUNS=\${E3SM_DIAGS_REPRO_RUNS}"

if ! python - <<'PY'
import platform
import sys

import xarray as xr

expected_python = "${expected_python}"
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

  chmod +x "${COMMAND_SCRIPT}"
}


run_confirmation() {
  local srun_args=(--pty --nodes=1 --time="${SRUN_TIME}")

  if [[ -n "${ACCOUNT}" ]]; then
    srun_args+=(--account="${ACCOUNT}")
  fi

  log "Run directory: ${RUN_DIR}"
  log "Starting ${ENV_NAME}"
  srun "${srun_args[@]}" /bin/bash "${COMMAND_SCRIPT}" 2>"${SRUN_STDERR_LOG}" || true
}


write_summary() {
  local job_id="NA"
  local status="failed"
  local exit_code="NA"
  local error_note="NA"
  local key
  local value

  if [[ -f "${STATUS_FILE}" ]]; then
    while IFS='=' read -r key value; do
      case "${key}" in
        job_id) job_id=${value:-NA} ;;
        status) status=${value:-failed} ;;
        exit_code) exit_code=${value:-NA} ;;
        error_note) error_note=${value:-NA} ;;
      esac
    done < "${STATUS_FILE}"
  elif [[ -s "${SRUN_STDERR_LOG}" ]]; then
    status="submit_failed"
    error_note=$(tail -n 5 "${SRUN_STDERR_LOG}" | tr '\n' ' ' | sed 's/[[:space:]]\+/ /g; s/^ //; s/ $//')
  fi

  printf 'env\tpython\txarray\tjob_id\tstatus\texit_code\trun_log\terror_note\n' > "${SUMMARY_FILE}"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${ENV_NAME}" "$(python_version_for_env "${ENV_NAME}")" "${XARRAY_VERSION}" \
    "${job_id}" "${status}" "${exit_code}" "${RUN_LOG}" "${error_note}" \
    >> "${SUMMARY_FILE}"

  log "Summary: ${SUMMARY_FILE}"
  log "Run log: ${RUN_LOG}"
}


main() {
  parse_args "$@"
  init_paths
  preflight
  write_command_script
  run_confirmation
  write_summary
}


main "$@"

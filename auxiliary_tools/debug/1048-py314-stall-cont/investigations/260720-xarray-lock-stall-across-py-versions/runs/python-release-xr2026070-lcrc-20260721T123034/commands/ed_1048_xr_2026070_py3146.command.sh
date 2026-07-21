#!/usr/bin/env bash

set -euo pipefail

readonly ENV_NAME="ed_1048_xr_2026070_py3146"
readonly EXPECTED_PYTHON="3.14.6"
readonly EXPECTED_XARRAY="2026.7.0"
readonly REPO_ROOT="/gpfs/fs1/home/ac.tvo/E3SM-Project/e3sm_diags"
readonly QA_SCRIPT="/gpfs/fs1/home/ac.tvo/E3SM-Project/e3sm_diags/auxiliary_tools/debug/1048-py314-stall-cont/parallel-lcrc/min-scripts/qa.py"
readonly RUN_LOG="/home/ac.tvo/E3SM-Project/e3sm_diags/auxiliary_tools/debug/1048-py314-stall-cont/investigations/260720-xarray-lock-stall-across-py-versions/runs/python-release-xr2026070-lcrc-20260721T123034/logs/ed_1048_xr_2026070_py3146.run.log"
readonly STATUS_FILE="/home/ac.tvo/E3SM-Project/e3sm_diags/auxiliary_tools/debug/1048-py314-stall-cont/investigations/260720-xarray-lock-stall-across-py-versions/runs/python-release-xr2026070-lcrc-20260721T123034/status/ed_1048_xr_2026070_py3146.env_status"
readonly RUN_TIMEOUT="75m"
readonly REPRO_RUNS="5"
readonly CLIMO_LOCK_WORKAROUND="enabled"
readonly CLIMO_LOCK_WORKAROUND_OVERRIDE="0"

exec > >(tee -a "${RUN_LOG}") 2>&1

source "/home/ac.tvo/miniforge3/etc/profile.d/conda.sh"
set +u
conda activate "${ENV_NAME}"
set -u
cd "${REPO_ROOT}"

export CARTOPY_DATA_DIR="${CARTOPY_DATA_DIR:-}"
export E3SM_DIAGS_DISABLE_CLIMO_LOCK_WORKAROUND="${CLIMO_LOCK_WORKAROUND_OVERRIDE}"
export E3SM_DIAGS_REPRO_RUNS="${REPRO_RUNS}"

echo "Started: $(date --iso-8601=seconds)"
echo "Host: $(hostname)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-unknown}"
echo "Conda env: ${CONDA_DEFAULT_ENV:-unknown}"
echo "QA script: ${QA_SCRIPT}"
echo "Timeout: ${RUN_TIMEOUT}"
echo "Climo lock workaround: ${CLIMO_LOCK_WORKAROUND}"
echo "E3SM_DIAGS_DISABLE_CLIMO_LOCK_WORKAROUND=${E3SM_DIAGS_DISABLE_CLIMO_LOCK_WORKAROUND}"
echo "E3SM_DIAGS_REPRO_RUNS=${E3SM_DIAGS_REPRO_RUNS}"

if ! python - <<'PY'
import platform
import sys

import xarray as xr

expected_python = "3.14.6"
expected_xarray = "2026.7.0"
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
  cat > "${STATUS_FILE}" <<STATUS_EOF
env_name=${ENV_NAME}
python_version=${EXPECTED_PYTHON}
job_id=${SLURM_JOB_ID:-}
status=failed
exit_code=2
run_log=${RUN_LOG}
error_note=runtime version validation failed
ended_at=$(date --iso-8601=seconds)
STATUS_EOF
  exit 2
fi

exit_code=0
status="completed"
timeout --signal=TERM --kill-after=2m "${RUN_TIMEOUT}"   python "${QA_SCRIPT}" || exit_code=$?

case "${exit_code}" in
  0) status="completed" ;;
  124|137) status="timeout_stall" ;;
  *) status="failed" ;;
esac

cat > "${STATUS_FILE}" <<STATUS_EOF
env_name=${ENV_NAME}
python_version=${EXPECTED_PYTHON}
job_id=${SLURM_JOB_ID:-}
status=${status}
exit_code=${exit_code}
run_log=${RUN_LOG}
error_note=NA
ended_at=$(date --iso-8601=seconds)
STATUS_EOF

echo "Final status: ${status}"
echo "Exit code: ${exit_code}"
exit "${exit_code}"

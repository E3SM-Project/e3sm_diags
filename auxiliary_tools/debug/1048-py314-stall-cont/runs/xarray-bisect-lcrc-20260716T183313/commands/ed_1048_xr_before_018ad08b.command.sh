#!/usr/bin/env bash

set -euo pipefail
set -o pipefail

readonly ENV_NAME="ed_1048_xr_before_018ad08b"
readonly REPO_ROOT="/gpfs/fs1/home/ac.tvo/E3SM-Project/e3sm_diags"
readonly QA_SCRIPT="/gpfs/fs1/home/ac.tvo/E3SM-Project/e3sm_diags/auxiliary_tools/debug/1048-py314-stall-cont/parallel-lcrc/min-scripts/qa.py"
readonly RUN_LOG="/home/ac.tvo/E3SM-Project/e3sm_diags/auxiliary_tools/debug/1048-py314-stall-cont/runs/xarray-bisect-lcrc-20260716T183313/logs/ed_1048_xr_before_018ad08b.run.log"
readonly STATUS_FILE="/home/ac.tvo/E3SM-Project/e3sm_diags/auxiliary_tools/debug/1048-py314-stall-cont/runs/xarray-bisect-lcrc-20260716T183313/status/ed_1048_xr_before_018ad08b.env_status"
readonly RUN_TIMEOUT="30m"
readonly REPRO_RUNS="10"
readonly XARRAY_COMMIT_SHA="8d271fb393372bcd2ed6ab60c9f469a1625a4aed"
readonly XARRAY_COMMIT_SUBJECT="Make parallel documentation builds threadsafe (#11009)"

mkdir -p "$(dirname "${RUN_LOG}")" "$(dirname "${STATUS_FILE}")"

exec > >(tee -a "${RUN_LOG}") 2>&1

export CARTOPY_DATA_DIR="${CARTOPY_DATA_DIR:-}"
source "${HOME}/miniforge3/etc/profile.d/conda.sh"
set +u
conda activate "ed_1048_xr_before_018ad08b"
set -u
cd "${REPO_ROOT}"

export E3SM_DIAGS_DISABLE_CLIMO_LOCK_WORKAROUND=1
export E3SM_DIAGS_REPRO_RUNS="${E3SM_DIAGS_REPRO_RUNS:-${REPRO_RUNS}}"

echo "Started: $(date --iso-8601=seconds)"
echo "Host: $(hostname)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-unknown}"
echo "Conda env: ${CONDA_DEFAULT_ENV:-unknown}"
echo "QA script: ${QA_SCRIPT}"
echo "Timeout: ${RUN_TIMEOUT}"
echo "E3SM_DIAGS_DISABLE_CLIMO_LOCK_WORKAROUND=${E3SM_DIAGS_DISABLE_CLIMO_LOCK_WORKAROUND}"
echo "E3SM_DIAGS_REPRO_RUNS=${E3SM_DIAGS_REPRO_RUNS}"
echo "Expected xarray commit: ${XARRAY_COMMIT_SHA}"
echo "Expected xarray commit subject: ${XARRAY_COMMIT_SUBJECT}"

python - <<'PY'
import platform
import xarray as xr

print(f"Python version: {platform.python_version()}")
print(f"Python executable: {platform.python_implementation()}")
print(f"xarray version: {xr.__version__}")
print(f"xarray module path: {xr.__file__}")
PY

exit_code=0
status="completed"

timeout --signal=TERM --kill-after=2m "${RUN_TIMEOUT}" python "${QA_SCRIPT}" || exit_code=$?

case "${exit_code}" in
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

cat > "${STATUS_FILE}" <<STATUS_EOF
env_name=${ENV_NAME}
job_id=${SLURM_JOB_ID:-}
status=${status}
exit_code=${exit_code}
run_log=${RUN_LOG}
ended_at=$(date --iso-8601=seconds)
STATUS_EOF

echo "Final status: ${status}"
echo "Exit code: ${exit_code}"
exit "${exit_code}"

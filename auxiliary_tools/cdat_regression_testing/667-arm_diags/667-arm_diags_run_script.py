# python -m auxiliary_tools.cdat_regression_testing.667-arm_diags.667-arm_diags_run_script
from auxiliary_tools.cdat_regression_testing.base_run_script import run_set

SET_NAME = "arm_diags"
SET_DIR = "667-arm_diags-final"
CFG_PATH: str | None = None
# CFG_PATH = (
#     "./auxiliary_tools/cdat_regression_testing/667-arm_diags/arm_diags_model_vs_obs.cfg"
# )
MULTIPROCESSING = True

run_set(SET_NAME, SET_DIR, CFG_PATH, MULTIPROCESSING)

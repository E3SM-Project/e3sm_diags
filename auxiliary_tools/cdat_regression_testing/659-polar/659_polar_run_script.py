# python -m auxiliary_tools.cdat_regression_testing.792-lat-lon-run-script.792_lat_lon_run_script
from auxiliary_tools.cdat_regression_testing.base_run_script import run_set

SET_NAME = "polar"
SET_DIR = "659-polar"
CFG_PATH: str | None = None
# CFG_PATH: str | None = "/global/u2/v/vo13/E3SM-Project/e3sm_diags/auxiliary_tools/cdat_regression_testing/659-polar/659_polar.cfg"
MULTIPROCESSING = False

# %%
run_set(SET_NAME, SET_DIR, CFG_PATH, MULTIPROCESSING)

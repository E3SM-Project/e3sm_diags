# python -m auxiliary_tools.cdat_regression_testing.792-lat-lon-run-script.792_lat_lon_run_script
from auxiliary_tools.cdat_regression_testing.base_run_script import run_set

SET_NAME = "lat_lon"
SET_DIR = "792-lat-lon-debug"
# CFG_PATH: str | None = None
CFG_PATH: str | None = "/global/u2/v/vo13/E3SM-Project/e3sm_diags/auxiliary_tools/cdat_regression_testing/792-lat-lon-run-script/792_lat_lon.cfg"
MULTIPROCESSING = False

# %%
run_set(SET_NAME, SET_DIR, CFG_PATH, MULTIPROCESSING)

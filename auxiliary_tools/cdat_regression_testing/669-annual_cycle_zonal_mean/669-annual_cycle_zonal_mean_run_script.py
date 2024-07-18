# %%
from auxiliary_tools.cdat_regression_testing.base_run_script import run_set

SET_NAME = "annual_cycle_zonal_mean"
SET_DIR = (
    "/global/cfs/cdirs/e3sm/www/cdat-migration-fy24/669-annual_cycle_zonal_mean-debug"
)
CFG_PATH = "auxiliary_tools/cdat_regression_testing/669-annual_cycle_zonal_mean/669-annual_cycle_zonal_mean.cfg"
MULTIPROCESSING = False

run_set(SET_NAME, SET_DIR, CFG_PATH, MULTIPROCESSING)

# %%

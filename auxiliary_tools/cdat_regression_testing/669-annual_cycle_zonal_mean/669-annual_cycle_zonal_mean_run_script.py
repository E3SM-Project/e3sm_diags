# %%
from auxiliary_tools.cdat_regression_testing.base_run_script import run_set

SET_NAME = "annual_cycle_zonal_mean"
SET_DIR = (
    "/global/cfs/cdirs/e3sm/www/cdat-migration-fy24/669-annual_cycle_zonal_mean"
)
# CFG_PATH = "auxiliary_tools/cdat_regression_testing/669-annual_cycle_zonal_mean/669-annual_cycle_zonal_mean.cfg"
CFG_PATH = "e3sm_diags/driver/default_diags/annual_cycle_zonal_mean_model_vs_obs.cfg"
MULTIPROCESSING = True

run_set(SET_NAME, SET_DIR, CFG_PATH, MULTIPROCESSING)

# %%

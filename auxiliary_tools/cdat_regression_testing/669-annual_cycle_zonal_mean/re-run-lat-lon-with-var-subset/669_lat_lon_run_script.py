# python -m auxiliary_tools.cdat_regression_testing.669-annual_cycle_zonal_mean.re-run-lat-lon-with-var-subset.669_lat_lon_run_script
from auxiliary_tools.cdat_regression_testing.base_run_script import run_set

SET_NAME = "lat_lon"
SET_DIR = "669-annual_cycle_zonal_mean_lat_lon_test"
CFG_PATH: str | None = None
# CFG_PATH: str | None = "auxiliary_tools/cdat_regression_testing/669-annual_cycle_zonal_mean/re-run-lat-lon-with-var-subset/669_lat_lon_run_script.py"
MULTIPROCESSING = True

# %%
run_set(SET_NAME, SET_DIR, CFG_PATH, MULTIPROCESSING)

# %%
# python -m auxiliary_tools.cdat_regression_testing.663-enso-diags.662-area-mean-time-series-rerun.run_script
# chmod -R o=rx /global/cfs/cdirs/e3sm/www/cdat-migration-fy24/663-area-mean-time-series-rerun
from auxiliary_tools.cdat_regression_testing.base_run_script import run_set

SET_NAME = "area_mean_time_series"
SET_DIR = "663-area-mean-time-series-rerun"
CFG_PATH: str | None = None
# CFG_PATH = "/global/u2/v/vo13/E3SM-Project/e3sm_diags/auxiliary_tools/cdat_regression_testing/663-enso-diags/662-area-mean-time-series-rerun/run.cfg"
MULTIPROCESSING = True

run_set(SET_NAME, SET_DIR, CFG_PATH, MULTIPROCESSING)

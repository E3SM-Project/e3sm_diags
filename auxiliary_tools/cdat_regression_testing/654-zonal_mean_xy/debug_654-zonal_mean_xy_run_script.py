from auxiliary_tools.cdat_regression_testing.base_run_script import run_set

SET_NAME = "zonal_mean_xy"
SET_DIR = "debug-654-zonal_mean_xy"
# CFG_PATH = "auxiliary_tools/cdat_regression_testing/654-zonal_mean_xy/debug_zonal_mean_xy_model_vs_obs.cfg"

run_set(SET_NAME, SET_DIR)
# run_set(SET_NAME, SET_DIR, CFG_PATH, multiprocessing=False)

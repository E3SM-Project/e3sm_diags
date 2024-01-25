from auxiliary_tools.cdat_regression_testing.base_run_script import run_set

SET_NAME = "cosp_histogram"
SET_DIR = "660-cosp-histogram"
SAVE_NETCDF = True

run_set(SET_NAME, SET_DIR, save_netcdf=SAVE_NETCDF)

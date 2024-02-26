from auxiliary_tools.cdat_regression_testing.base_run_script import run_set

SET_NAME = "lat_lon"
SET_DIR = "792-lat-lon"
CFG_PATH: str | None = None
MULTIPROCESSING = True

run_set(SET_NAME, SET_DIR, CFG_PATH, MULTIPROCESSING)

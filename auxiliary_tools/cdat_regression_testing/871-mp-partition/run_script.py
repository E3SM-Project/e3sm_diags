from auxiliary_tools.cdat_regression_testing.base_run_script import run_set

SET_NAME = "mp_partition"
SET_DIR = "871-mp-partition"
CFG_PATH: str | None = None
MULTIPROCESSING = True

run_set(SET_NAME, SET_DIR, CFG_PATH, MULTIPROCESSING)

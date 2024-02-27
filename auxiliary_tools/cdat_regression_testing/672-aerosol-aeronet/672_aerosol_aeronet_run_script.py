from auxiliary_tools.cdat_regression_testing.base_run_script import run_set

SET_NAME = "aerosol_aeronet"
SET_DIR = "672-aerosol-aeronet"
CFG_PATH: str | None = None
MULTIPROCESSING = True

run_set(SET_NAME, SET_DIR, CFG_PATH, MULTIPROCESSING)

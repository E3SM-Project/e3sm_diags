from auxiliary_tools.cdat_regression_testing.base_run_script import run_set

# TODO: Update SETS_TO_RUN to the single set you are refactoring.
# Example: "lat_lon"
SET_NAME = "lat_lon"
# TODO: Update SET_DIR to <ISSUE-SET_NAME>. This string gets appended
# to the base results_dir, "/global/cfs/cdirs/e3sm/www/cdat-migration-fy24/".
# Example: "671-lat-lon"
SET_DIR = "792-lat-lon"

# TODO: <OPTIONAL> UPDATE CFG_PATH if using a custom cfg file for debugging.
# Example: "auxiliary_tools/cdat_regression_testing/654_zonal_mean_xy.cfg"
CFG_PATH: str | None = None

# TODO: <OPTIONAL> Update MULTIPROCESSING based on whether to run in parallel or
# serial. For debugging purposes, set to False to run serially.
MULTIPROCESSING = True

run_set(SET_NAME, SET_DIR, CFG_PATH, MULTIPROCESSING)

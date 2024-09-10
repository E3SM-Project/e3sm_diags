from auxiliary_tools.cdat_regression_testing.base_run_script import run_set

SET_NAME = [
    "lat_lon",
    "zonal_mean_xy",
    "zonal_mean_2d",
    "zonal_mean_2d_stratosphere",
    "polar",
    "cosp_histogram",
    "meridional_mean_2d",
    "annual_cycle_zonal_mean",
    "enso_diags",
    "qbo",
    "area_mean_time_series",
    "diurnal_cycle",
    "streamflow",
    "arm_diags",
    "tc_analysis",
    "aerosol_aeronet",
    "aerosol_budget",
    "mp_partition",
]

# TODO: Update SET_DIR to <ISSUE-SET_NAME>. This string gets appended
# to the base results_dir, "/global/cfs/cdirs/e3sm/www/cdat-migration-fy24/".
# Example: "671-lat-lon"
SET_DIR = ""

# TODO: <OPTIONAL> UPDATE CFG_PATH if using a custom cfg file for debugging.
# Example: "auxiliary_tools/cdat_regression_testing/654_zonal_mean_xy.cfg"
CFG_PATH: str | None = None

# TODO: <OPTIONAL> Update MULTIPROCESSING based on whether to run in parallel or
# serial. For debugging purposes, set to False to run serially.
MULTIPROCESSING = True

RUN_TYPE = "model_vs_obs"

run_set(SET_NAME, SET_DIR, CFG_PATH, MULTIPROCESSING, RUN_TYPE)  # type: ignore

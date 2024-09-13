# python -m auxiliary_tools.cdat_regression_testing.843-migration-phase3.model_vs_model.run_script
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

SET_DIR = "843-migration-phase3"

CFG_PATH: str | None = None
MULTIPROCESSING = True

RUN_TYPE = "model_vs_model"

run_set(SET_NAME, SET_DIR, CFG_PATH, MULTIPROCESSING, RUN_TYPE)  # type: ignore

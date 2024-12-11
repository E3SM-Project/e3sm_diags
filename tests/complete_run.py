from e3sm_diags.run import runner

from tests.complete_run_params import params

# Run
runner.sets_to_run = [
    "lat_lon",
    "zonal_mean_xy",
    "zonal_mean_2d",
    "polar",
    "cosp_histogram",
    "meridional_mean_2d",
    "enso_diags",
    "qbo",
    "diurnal_cycle",
    "annual_cycle_zonal_mean",
    "streamflow",
    "zonal_mean_2d_stratosphere",
    "arm_diags",
    "tc_analysis",
    "aerosol_aeronet",
    "aerosol_budget",
    "tropical_subseasonal",
]

runner.run_diags(params)

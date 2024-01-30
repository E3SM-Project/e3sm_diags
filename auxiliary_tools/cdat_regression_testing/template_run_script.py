"""
The template run script used for generating results on your development branch.

Steps:
1. Activate your conda dev env for your branch
2. `make install` to install the latest version of your branch code into the env
3. Copy this script into `auxiliary_tools/cdat_regression_testing/<ISSUE>-<SET_NAME>`
4. Update `SET_DIR` string variable
5. Update `SET_NAME` string variable.
   - Options include: "lat_lon", "zonal_mean_xy", "zonal_mean_2d",
     "zonal_mean_2d_stratosphere", "polar", "cosp_histogram",
     "meridional_mean_2d", "annual_cycle_zonal_mean", "enso_diags", "qbo",
     "area_mean_time_series", "diurnal_cycle", "streamflow", "arm_diags",
     "tc_analysis", "aerosol_aeronet", "aerosol_budget", "mp_partition",
6. Run this script
   - Make sure to run this command on NERSC perlmutter cpu:
    `salloc --nodes 1 --qos interactive --time 01:00:00 --constraint cpu --account=e3sm
    conda activate <NAME-OF-DEV-ENV>`
   - python auxiliary_tools/cdat_regression_testing/<ISSUE-<SET_NAME>
7. Make a copy of the CDAT regression testing notebook in the same directory
   as this script and follow the instructions there to start testing.
"""
from auxiliary_tools.cdat_regression_testing.base_run_script import run_set

# TODO: Update SETS_TO_RUN to the single set you are refactoring.
# Example: "lat_lon"
SET_NAME = ""
# TODO: Update SET_DIR to <ISSUE-SET_NAME>. This string gets appended
# to the base results_dir, "/global/cfs/projectdirs/e3sm/e3sm_diags_cdat_test/".
# Example: "671-lat-lon"
SET_DIR = ""

run_set(SET_NAME, SET_DIR)

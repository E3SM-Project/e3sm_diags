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
6. Run this script as a Python module
   - `auxiliary_tools` is not included in `setup.py`, so `-m` is required
     to run the script as a Python module
   - Command: python -m auxiliary_tools.cdat_regression_testing.<ISSUE>-<SET_NAME>.<SCRIPT-NAME>
   - Example: python -m auxiliary_tools.cdat_regression_testing.660_cosp_histogram.run_script
7. Run `chmod -R o=rx <SET_DIR>` to allow public access to viewer outputs on the NERSC webserver
  - Example: `chmod -R o=rx /global/cfs/cdirs/e3sm/www/cdat-migration-fy24/654-zonal_mean_xy`
  - https://portal.nersc.gov/project/e3sm/cdat-migration-fy24/
8. Make a copy of the CDAT regression testing notebook in the same directory
   as this script and follow the instructions there to start testing.
9. <OPTIONAL> Update `CFG_PATH` to a custom cfg file to debug specific variables.
   - It is useful to create a custom cfg based on the default diags to debug
     specific variables that are running into problems.
   - For example, copy `zonal_mean_xy_model_vs_model.cfg` into the same directory
     as the copy of this script, then modify it to specific variables. Afterwards
     update `CFG_PATH` to the path of that .cfg file.
   - Tip: Use VS Code to step through the code with the Python debugger.
"""
from auxiliary_tools.cdat_regression_testing.base_run_script import run_set

# TODO: Update SETS_TO_RUN to the single set you are refactoring.
# Example: "lat_lon"
SET_NAME = "tc_analysis"
# TODO: Update SET_DIR to <ISSUE-SET_NAME>. This string gets appended
# to the base results_dir, "/global/cfs/cdirs/e3sm/www/cdat-migration-fy24/".
# Example: "671-lat-lon"
SET_DIR = "688-tc_analysis"

# TODO: <OPTIONAL> UPDATE CFG_PATH if using a custom cfg file for debugging.
# Example: "auxiliary_tools/cdat_regression_testing/654_zonal_mean_xy.cfg"
CFG_PATH: str | None = None

# TODO: <OPTIONAL> Update MULTIPROCESSING based on whether to run in parallel or
# serial. For debugging purposes, set to False to run serially.
MULTIPROCESSING = True

run_set(SET_NAME, SET_DIR, CFG_PATH, MULTIPROCESSING)

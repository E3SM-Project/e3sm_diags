from auxiliary_tools.cdat_regression_testing.base_run_script import run_set

SET_NAME = "lat_lon"
SET_DIR = "792-lat-lon"
# CFG_PATH: str | None = "/global/u2/v/vo13/E3SM-Project/e3sm_diags/auxiliary_tools/cdat_regression_testing/792-lat-lon-run-script/792_lat_lon.cfg"
CFG_PATH: str | None = None
MULTIPROCESSING = True

# %%
run_set(SET_NAME, SET_DIR, CFG_PATH, MULTIPROCESSING)


# Debugging
# 1. ValueError: Multiple 'Z' axis dims were found in this dataset, ['ilev', 'lev'].
# Please drop the unused dimension(s) before performing grid operations.
#   - Fix: Add `_drop_unused_ilev_axis`
# 2. KeyError: No attribute "units"
#   - Root case: In `cosp_bin_sum()`, calling `.sum()` results in Xarray
#     dropping the original units attribute. The original units are required
#     for the "convert_units()`" function to correctly convert units to "%".
#     Otherwise, units are None and never set to %.
#   - Fix: Add `keep_attrs=True`
# 3. Missing ERA5-TAUXY (ANN/JJA, ocean/land), ERA5-TREFHT (ANN/JJA, lnd),
# MERRA2-TAUXY (ANN/JJA, ocean),MERRA2-TREFHT (ANN/JJA, land), HadISST-SST (ANN, global)

# %%
# python -m auxiliary_tools.cdat_regression_testing.666-diurnal-cycle.run_script
from auxiliary_tools.cdat_regression_testing.base_run_script import run_set

SET_NAME = "diurnal_cycle"
SET_DIR = "666-diurnal-cycle"
CFG_PATH: str | None = None
# CFG_PATH: str | None = "/global/u2/v/vo13/E3SM-Project/e3sm_diags/auxiliary_tools/cdat_regression_testing/666-diurnal-cycle/run.cfg"
MULTIPROCESSING = True

# %%
run_set(SET_NAME, SET_DIR, CFG_PATH, MULTIPROCESSING)


# /global/cfs/cdirs/e3sm/www/cdat-migration-fy24/666-diurnal-cycle/diurnal_cycle/TRMM-3B43v-7_3hr/TRMM-3B43v-7_3hr-PRECT-ANN-20S20N.png
# /global/cfs/cdirs/e3sm/www/cdat-migration-fy24/666-diurnal-cycle/diurnal_cycle/TRMM-3B43v-7_3hr_diff/TRMM-3B43v-7_3hr-PRECT-ANN-20S20N.png

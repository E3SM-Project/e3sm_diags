# python -m auxiliary_tools.cdat_regression_testing.657-meridional-mean-2d.657_meridional_mean_2d_run_script
from auxiliary_tools.cdat_regression_testing.base_run_script import run_set

SET_NAME = "meridional_mean_2d"
SET_DIR = "657-meridional-mean-2d"
CFG_PATH: str | None = None
# CFG_PATH: str | None = "/global/u2/v/vo13/E3SM-Project/e3sm_diags/auxiliary_tools/cdat_regression_testing/657-meridional-mean-2d/657_meridional_mean_2d.cfg"
MULTIPROCESSING = True

# %%
run_set(SET_NAME, SET_DIR, CFG_PATH, MULTIPROCESSING)

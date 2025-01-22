from auxiliary_tools.cdat_regression_testing.base_run_script import run_set

SET_NAME = "polar"
SET_DIR = "909-prov-log"

CFG_PATH: str | None = "/global/u2/v/vo13/E3SM-Project/e3sm_diags/auxiliary_tools/cdat_regression_testing/909-prov-log/run.cfg"
MULTIPROCESSING = True

# %%
run_set(SET_NAME, SET_DIR, CFG_PATH, MULTIPROCESSING)

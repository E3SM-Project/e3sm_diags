# %%
# python -m auxiliary_tools.cdat_regression_testing.673-aerosol-budget.run_script
from auxiliary_tools.cdat_regression_testing.base_run_script import run_set

SET_NAME = "aerosol_budget"
SET_DIR = "673-aerosol-budget"
CFG_PATH: str | None = (
    "auxiliary_tools/cdat_regression_testing/673-aerosol-budget/run.cfg"
)
# CFG_PATH: str | None = None
MULTIPROCESSING = True

run_set(SET_NAME, SET_DIR, CFG_PATH, MULTIPROCESSING)

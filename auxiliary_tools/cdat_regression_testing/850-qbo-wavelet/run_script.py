# %%
from auxiliary_tools.cdat_regression_testing.base_run_script import run_set

SET_NAME = "qbo"
SET_DIR = "850-qbo-wavelet"
CFG_PATH: str | None = None
MULTIPROCESSING = False

# %%
run_set(SET_NAME, SET_DIR, CFG_PATH, MULTIPROCESSING)

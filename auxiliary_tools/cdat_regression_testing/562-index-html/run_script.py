# python -m auxiliary_tools.cdat_regression_testing.792-lat-lon-run-script.792_lat_lon_run_script
from auxiliary_tools.cdat_regression_testing.base_run_script import run_set

SET_NAME = "polar"
SET_DIR = "562-index-html"
# CFG_PATH: str | None = None
CFG_PATH: str | None = "auxiliary_tools/cdat_regression_testing/562-index-html/562_index_html.cfg"
MULTIPROCESSING = False

# %%
run_set(SET_NAME, SET_DIR, CFG_PATH, MULTIPROCESSING)

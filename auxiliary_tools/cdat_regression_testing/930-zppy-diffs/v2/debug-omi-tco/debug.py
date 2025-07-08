#%%
import glob

import numpy as np
import xarray as xr

from e3sm_diags.derivations.derivations import DERIVED_VARIABLES
import pandas as pd

DEV_DIR = "25-02-04-branch-930-zppy-diffs"
DEV_PATH = f"/lcrc/group/e3sm/public_html/cdat-migration-fy24/{DEV_DIR}/"

DEV_GLOB_ALL = sorted(glob.glob(DEV_PATH + "**/**/*.nc"))
DEV_NUM_FILES = len(DEV_GLOB_ALL)

MAIN_DIR = "25-02-04-main-zppy-diffs"
MAIN_PATH = f"/lcrc/group/e3sm/public_html/cdat-migration-fy24/{MAIN_DIR}/"
MAIN_GLOB_ALL = sorted(glob.glob(MAIN_PATH + "**/**/*.nc"))
MAIN_NUM_FILES = len(MAIN_GLOB_ALL)

#%%
KEEP_VARS = ["OMI-MLS"]

DEV_GLOB = [fp for fp in DEV_GLOB_ALL if "diff.nc" not in fp and any(var in fp for var in KEEP_VARS)]
MAIN_GLOB = [fp for fp in MAIN_GLOB_ALL if "diff.nc" not in fp and any(var in fp for var in KEEP_VARS)]

DEV_GLOB_DIFF = [fp for fp in DEV_GLOB_ALL if "diff.nc" in fp and any(var in fp for var in KEEP_VARS)]
MAIN_GLOB_DIFF = [fp for fp in MAIN_GLOB_ALL if "diff.nc" in fp and any(var in fp for var in KEEP_VARS)]

# %%
'/lcrc/group/e3sm/public_html/cdat-migration-fy24/25-02-04-branch-930-zppy-diffs/lat_lon/OMI-MLS/OMI-MLS-TCO-ANN-60S60N_ref.nc',

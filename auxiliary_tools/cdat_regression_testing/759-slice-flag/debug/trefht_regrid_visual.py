"""This script compares the xCDAT and CDAT regridded reference data results.

The `debug_TREFHT_diff.png` shows there there is a difference along the coastlines
of the continents. CDAT doesn't seem to actually regrid the data, which is
noted in `trefht_cdat_only_visual.py`.
"""
# %%
import os

import numpy as np
import xarray as xr

from auxiliary_tools.cdat_regression_testing.utils import get_image_diffs

dir_path = "/global/u2/v/vo13/E3SM-Project/e3sm_diags/auxiliary_tools/cdat_regression_testing/759-slice-flag/debug/"
fp_a = os.path.join(dir_path, "ds_b_regrid.nc")
fp_b = os.path.join(dir_path, "mv2_regrid.nc")


# %%
ds1 = xr.open_dataset(fp_a)
ds2 = xr.open_dataset(fp_b)

var_key = "TREFHT"
ds2["TREFHT"] = ds2["variable_36"].copy()
ds2 = ds2.drop_vars("variable_36")

# %%
try:
    np.testing.assert_allclose(ds1[var_key], ds2[var_key])
except AssertionError as e:
    print(e)

# %%
# Check the sum values -- close
# array(213927.85820162)
np.abs(ds1[var_key]).sum()

# %%
# array(228804.60960445)
np.abs(ds2[var_key]).sum()

# Check the mean values -- close
# array(-7.8665246)
ds1[var_key].mean()

# array(-6.34208681)
ds2[var_key].mean()


# %%
# Check the plots and their diffs
root_dir = "auxiliary_tools/cdat_regression_testing/759-slice-flag"
actual_path = os.path.join(root_dir, "debug_TREFHT_actual.png")
expected_path = os.path.join(root_dir, "debug_TREFHT_expected.png")

ax1 = ds1[var_key].plot()
ax1.figure.savefig(actual_path)

# %%
ax2 = ds2[var_key].plot()
ax2.figure.savefig(expected_path)

# %%
get_image_diffs(actual_path, expected_path)

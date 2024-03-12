# %%
import os

import numpy as np
import xarray as xr

from auxiliary_tools.cdat_regression_testing.utils import get_image_diffs

# %%
ds1 = xr.open_dataset(
    "/global/cfs/cdirs/e3sm/www/cdat-migration-fy24/792-lat-lon/lat_lon/MERRA2/MERRA2-TREFHT-ANN-land_ref.nc"
)
ds2 = xr.open_dataset(
    "/global/cfs/cdirs/e3sm/www/cdat-migration-fy24/main/lat_lon/MERRA2/MERRA2-TREFHT-ANN-land_ref.nc"
)

var_key = "TREFHT"

# %%
try:
    np.testing.assert_allclose(ds1[var_key], ds2[var_key])
except AssertionError as e:
    print(e)

# %%
# Get the nan count -- close
# 706526
nan_count_ds1 = ds1[var_key].isnull().sum()
# 702926
nan_count_ds2 = ds2[var_key].isnull().sum()

# %%
# Check the absolute sum values -- close
# 1470137.
np.abs(ds1[var_key]).sum()
# 1502139.5
np.abs(ds2[var_key]).sum()

# %%
# Check the mean values -- close
# -6.186151
ds1[var_key].mean()

# -6.546063
ds2[var_key].mean()

# %%
# Check the plots and their diffs
root_dir = "auxiliary_tools/cdat_regression_testing/792-lat-lon-run-script"
actual_path = os.path.join(root_dir, "debug_merra2_trefht_actual.png")
expected_path = os.path.join(root_dir, "debug_merra2_trefht_expected.png")

ax1 = ds1[var_key].plot()
ax1.figure.savefig(actual_path)

# %%
ax2 = ds2[var_key].plot()
ax2.figure.savefig(expected_path)

# %%
get_image_diffs(actual_path, expected_path)

# %%

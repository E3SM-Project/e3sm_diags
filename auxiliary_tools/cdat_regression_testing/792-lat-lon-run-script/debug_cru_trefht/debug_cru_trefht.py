# %%
import os

import numpy as np
import xarray as xr

from auxiliary_tools.cdat_regression_testing.utils import get_image_diffs

# %%
ds1 = xr.open_dataset(
    "/global/cfs/cdirs/e3sm/www/cdat-migration-fy24/792-lat-lon-debug/lat_lon/CRU_IPCC/CRU-TREFHT-ANN-land_60S90N_ref.nc"
)
ds2 = xr.open_dataset(
    "/global/cfs/cdirs/e3sm/www/cdat-migration-fy24/main/lat_lon/CRU_IPCC/CRU-TREFHT-ANN-land_60S90N_ref.nc"
)

var_key = "TREFHT"

# %%
np.testing.assert_allclose(ds1[var_key], ds2[var_key])

# %%
# Check the sum values -- close
# 119496.125
np.abs(ds1[var_key]).sum()
# 119496.11
np.abs(ds2[var_key]).sum()

# %%
# Check the mean values -- close
# 8.3522835
ds1[var_key].mean()

# 8.352283
ds2[var_key].mean()

# %%
# Check the plots and their diffs
root_dir = "auxiliary_tools/cdat_regression_testing/792-lat-lon-run-script"
actual_path = os.path.join(root_dir, "debug_trefht_actual.png")
expected_path = os.path.join(root_dir, "debug_trefht_expected.png")

ax1 = ds1[var_key].plot()
ax1.figure.savefig(actual_path)

# %%
ax2 = ds2[var_key].plot()
ax2.figure.savefig(expected_path)

# %%
get_image_diffs(actual_path, expected_path)

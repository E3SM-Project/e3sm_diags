# %%
import os

import numpy as np
import xarray as xr

from auxiliary_tools.cdat_regression_testing.utils import get_image_diffs

ds1 = xr.open_dataset(
    "/global/cfs/cdirs/e3sm/www/cdat-migration-fy24/792-lat-lon/lat_lon/ERA5/ERA5-TAUXY-ANN-ocean_ref.nc"
)
ds2 = xr.open_dataset(
    "/global/cfs/cdirs/e3sm/www/cdat-migration-fy24/main/lat_lon/ERA5/ERA5-TAUXY-ANN-ocean_ref.nc"
)

var_key = "TAUXY"
# %%
try:
    np.testing.assert_allclose(ds1[var_key], ds2[var_key])
except AssertionError as e:
    print(e)

# Check the sum values -- close
# 37431.402
np.abs(ds1[var_key]).sum()
# 37508.312
np.abs(ds2[var_key]).sum()

# Check the mean values -- close
# 0.07276739
ds1[var_key].mean()

# 0.07276866
ds2[var_key].mean()


# %%
# Check the plots and their diffs
root_dir = "auxiliary_tools/cdat_regression_testing/792-lat-lon-run-script"
actual_path = os.path.join(root_dir, "debug_tauxy_actual.png")
expected_path = os.path.join(root_dir, "debug_tauxy_expected.png")

ax1 = ds1[var_key].plot()
ax1.figure.savefig(actual_path)

# %%
ax2 = ds2[var_key].plot()
ax2.figure.savefig(expected_path)

# %%

get_image_diffs(actual_path, expected_path)

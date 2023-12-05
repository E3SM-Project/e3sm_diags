"""This script compares the time series, climatology, and spatial average of the
climatology between the cdat-migration-fy24 and main branches.

CONCLUSION:
 - In the CDAT version of the Dataset class, cdms2.open is being called
   with a slice flag (either "co"/"ccb"). In the case of NETFLUX_SRF, "co" is
   being set because time coordinates start at the beginning of the month.
  - This adds an additional time coordinate point to the end of the time series
    file, which affects the subsequen climatology and spatial averaging
    calculations.

I found omitting the extra time coordinate point before calculating the
climatology and spatial averaging results in an identical result to the xCDAT
version.

Next options:
1. Add a feature to the Dataset class that adds an extra coordinate point using
the slice flag conditional

Related lines of code on `main`:
  - https://github.com/E3SM-Project/e3sm_diags/blob/633b52c314325e605fe7f62687cc4d00e5a0a3d5/e3sm_diags/driver/utils/dataset.py#L665-L672
  - https://github.com/E3SM-Project/e3sm_diags/blob/633b52c314325e605fe7f62687cc4d00e5a0a3d5/e3sm_diags/driver/utils/dataset.py#L699-L700


"""
import cdms2
import numpy as np
import xarray as xr

from e3sm_diags.driver.utils.climo import climo
from e3sm_diags.metrics import mean
from e3sm_diags.metrics.metrics import spatial_avg

# Path to the netcdf files for NET_FLUX_SRF generated on `cdat-migration-fy24`
# and `main` using `ex1.py`.
DIR_PATH = "/global/u2/v/vo13/E3SM-Project/e3sm_diags/auxiliary_tools/cdat_regression_testing/671-lat-lon"


# Compare time series -- identical only if extra time coordinate is removed
# ------------------------------------------------------------------------------
ds_ts1 = xr.open_dataset(f"{DIR_PATH}/671-ts-input.nc")
ds_ts2 = xr.open_dataset(f"{DIR_PATH}/main-ts-input.nc")
ds_ts2 = ds_ts2.rename({"variable_58": "NET_FLUX_SRF"})

# Extra time coordinate becaues of the cdms2 slice flag ("co" for beginning of
# the month/"ccb" for mid-month)
ds_ts2_sub = ds_ts2.isel(time=slice(0, -1))
np.testing.assert_allclose(ds_ts1["NET_FLUX_SRF"], ds_ts2_sub["NET_FLUX_SRF"])

# Result: True

# Compare climatologies -- not identical (due to extra coordinate point)
# ------------------------------------------------------------------------------
ds_climo1 = xr.open_dataset(f"{DIR_PATH}/671-climo.nc")
ds_climo2 = xr.open_dataset(f"{DIR_PATH}/main-climo.nc")
ds_climo2 = ds_climo2.rename({"variable_58": "NET_FLUX_SRF"})

np.testing.assert_allclose(ds_climo1["NET_FLUX_SRF"], ds_climo2["NET_FLUX_SRF"])

# Result: AssertionError:
# Not equal to tolerance rtol=1e-07, atol=0

# Mismatched elements: 33024 / 33024 (100%)
# Max absolute difference: 17.20686057
# Max relative difference: 16932.010712
#  x: array([[-0.439175, -0.439179, -0.439189, ..., -0.439205, -0.439189,
#         -0.439179],
#        [-0.435518, -0.435518, -0.435514, ..., -0.435508, -0.435514,...
#  y: array([[ 0.026507,  0.026507,  0.026507, ...,  0.026508,  0.026507,
#          0.026507],
#        [-0.020339, -0.02031 , -0.020224, ..., -0.02008 , -0.020224,...


# Compare spatial averages -- not identical (due to extra coordinate point)
# ------------------------------------------------------------------------------
ds_avg1 = spatial_avg(ds_climo1, "NET_FLUX_SRF")

ds_avg2 = cdms2.open(f"{DIR_PATH}/main-climo.nc")["variable_58"]
ds_avg2 = mean(ds_avg2)

np.testing.assert_allclose(ds_avg1, ds_avg2.data)
# Mismatched elements: 1 / 1 (100%)
# Max absolute difference: 0.12231419
# Max relative difference: 0.23689147
#  x: array(0.394016)
#  y: array(0.51633)


# Now let's try removing that extra coordinate point from the cdms2 time series
# -----------------------------------------------------------------------------
ds_ts3 = cdms2.open(f"{DIR_PATH}/main-ts-input.nc")["variable_58"]
ds_ts3_sub = ds_ts3(time=("2011-2-1", "2013-12-1"))

# Climatologies are now identical!
ds_climo3 = climo(ds_ts3_sub, "ANN")
np.testing.assert_allclose(ds_climo1["NET_FLUX_SRF"].data, ds_climo3.data)

# Spatial averages are now identical!
ds_avg3 = mean(ds_climo3)
np.testing.assert_allclose(ds_avg1, ds_avg3)

# %%
import xcdat as xc

# Parameters
# Time coordinates at the first day of the month (1850-01-01)
fn1 = "/global/cfs/cdirs/e3sm/e3sm_diags/test_model_data_for_acme_diags/time-series/E3SM_v1/FSNS_185001_201312.nc"
var = "FSNS"
start_time = "2011-01-15"
end_time = "2013-12-15"

slice_flag = "ccb"

# %%
ds_xcdat = xc.open_dataset(fn1, decode_times=True)

ds_xcdat_s = ds_xcdat.sel(time=slice(start_time, end_time))

# %%
# 35 2013-12-01 00:00:00
print(len(ds_xcdat_s.time), ds_xcdat_s.time.values[-1])

# [cftime.DatetimeNoLeap(2013, 11, 1, 0, 0, 0, 0, has_year_zero=True)
# cftime.DatetimeNoLeap(2013, 12, 1, 0, 0, 0, 0, has_year_zero=True)]
print(ds_xcdat_s.time_bnds.values[-1])


"""Possible solutions:
1. Use last time coordinate (2014, 1, 1) with subsetting
2. Use start value of the last bounds
  - (2013, 12, 1), (2014, 1, 1)
"""

# 1. Get the last time coordinate as ds1
# 2. Combine with ds2 (subsetted dataset)

ds1 = ds_xcdat.isel(time=-1)

# %%
import xarray as xr

ds_final = xr.concat([ds_xcdat_s, ds1], dim="time")

# %%
import cdms2
import xcdat as xc

# Parameters
# Time coordinates at the first day of the month (1850-01-01)
fn1 = "/global/cfs/cdirs/e3sm/e3sm_diags/test_model_data_for_acme_diags/time-series/E3SM_v1/FSNS_185001_201312.nc"
var = "FSNS"
start_time = "2011-01-15"
end_time = "2013-12-15"

slice_flag = "ccb"


ds_cdat2 = cdms2.open(fn1)
var_s = ds_cdat2(var, time=(start_time, end_time, slice_flag))(squeeze=1)

# %%
import numpy as np

np.testing.assert_allclose(ds_final["FSNS"].values, var_s.data, rtol=0, atol=0)
# %%

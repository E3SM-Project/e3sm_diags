# %%
from glob import glob

import cftime
import xcdat as xc

args = {
    "add_bounds": ["X", "Y"],
    "coords": "minimal",
    "compat": "override",
    "chunks": "auto",
}

filepath = "/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/climatology/MERRA2_Aerosols/MERRA2_Aerosols_[0-1][0-9]_*climo.nc"
paths = sorted(glob(filepath))

# Check time coordinates are in ascending order with wildcard and sorted glob
ds_mf = xc.open_mfdataset(filepath, decode_times=True, **args)
ds_mf_raw = xc.open_mfdataset(paths, **args, decode_times=True)

# %%

# filepath 1: '/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/climatology/MERRA2_Aerosols/MERRA2_Aerosols_01_198001_202101_climo.nc'
ds_raw_time = xc.open_mfdataset(paths[0], **args, decode_times=False)

# 10782720
time_int = ds_raw_time.time.values.item()
# 'minutes since 1980-01-01 00:30:00'
units = ds_raw_time.time.units
# None, so "standard"
calendar = ds_raw_time.time.attrs.get("calendar", "standard")

# cftime.DatetimeGregorian(2000, 7, 2, 0, 30, 0, 0, has_year_zero=False)
cftime.num2date(time_int, units, calendar=calendar)

# %%
import datetime

first_step = datetime.datetime(1980, 1, 1, hour=0, minute=30)
time_delta = datetime.timedelta(minutes=10782720)

# datetime.datetime(2000, 7, 2, 0, 30)
print(first_step + time_delta)

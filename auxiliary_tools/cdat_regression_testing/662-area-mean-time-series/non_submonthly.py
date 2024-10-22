# %%
from datetime import datetime

import pandas as pd
import xarray as xr
import xcdat as xc

# %%
fn1 = "/global/cfs/cdirs/e3sm/e3sm_diags/test_model_data_for_acme_diags/time-series/E3SM_v1/FSNS_185001_201312.nc"
var = "FSNS"
slice_flag = "ccb"

start_time = "2011-01-15"
end_time = "2013-12-15"
time_slice = slice(start_time, end_time, None)

# %%
ds = xr.open_dataset(fn1)


# %%
def get_new_end_time(ds: xr.Dataset, end_time: str) -> str:
    # Get the dimension coordinates and bounds key.
    time_coords = xc.get_dim_coords(ds, axis="T")
    time_dim = time_coords.name
    time_bnds_key = time_coords.attrs["bounds"]

    # Extract the sub-dataset for all data at the last time coordinate.
    ds_last_time = ds.isel({time_dim: -1})

    # Get the delta between time coordinates using the difference between
    # bounds values.
    time_bnds = ds_last_time[time_bnds_key]
    time_delta = time_bnds[-1] - time_bnds[0]
    time_delta_py = pd.to_timedelta(time_delta.values).to_pytimedelta()

    old_end_time = datetime.strptime(end_time, "%Y-%m-%d")
    new_end_time = old_end_time + time_delta_py
    new_end_time_str = new_end_time.strftime("%Y-%m-%d")

    return new_end_time_str


new_end_time = get_new_end_time(ds, end_time)

ds.sel(time=slice(start_time, new_end_time, None)).squeeze()

# %%
import timeit

import xarray as xr

filepaths = [
    "/lcrc/group/e3sm/diagnostics/observations/Atm/time-series/ERA5/ua_197901_201912.nc"
]

ds = xr.open_mfdataset(filepaths)

ds_sub = ds.sel(time=slice("1996-01-15", "1997-01-15", None))

# %%
start_time = timeit.default_timer()
ds_sub.ua.load()
elapsed = timeit.default_timer() - start_time
print(f"Time taken to load ds_xc_sub: {elapsed} seconds")

# %%

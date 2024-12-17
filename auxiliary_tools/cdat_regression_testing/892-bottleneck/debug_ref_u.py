"""
This script is used to debug the bottleneck issue in the reference u variable.
"""

# %%
import timeit

import xcdat as xc

filepaths = [
    "/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/time-series/ERA5/ua_197901_201912.nc"
]
time_slice = slice("1996-01-15", "1997-01-15", None)

# %%
# Test case 1 - OPEN_MFDATASET() + "ua" dataset (76 GB) + subsetting + `.load()`
# Result: .load() hangs when using `open_mfdataset`
# ------------------------------------------------------------------------------
ds_ua_omfd = xc.open_mfdataset(
    filepaths[0],
    add_bounds=["X", "Y", "T"],
    decode_times=True,
    use_cftime=True,
    coords="minimal",
    compat="override",
)
ds_ua_omfd_sub = ds_ua_omfd.sel(time=time_slice)

# %%
start_time = timeit.default_timer()
ds_ua_omfd_sub.load()
elapsed = timeit.default_timer() - start_time
print(f"Time taken to load ds_xc_sub: {elapsed} seconds")

# %%
# Test case 2 - OPEN_DATASET() + "ua" dataset (76 GB) + subsetting + `.load()`
# Result: load() works fine when using `open_dataset`
# ------------------------------------------------------------------------------
ds_ua_od = xc.open_dataset(
    filepaths[0],
    add_bounds=["X", "Y", "T"],
    decode_times=True,
    use_cftime=True,
    # coords="minimal",
    # compat="override",
)
ds_ua_od_sub = ds_ua_od.sel(time=time_slice)

# %%
start_time = timeit.default_timer()
ds_ua_od_sub.load()
elapsed = timeit.default_timer() - start_time
print(f"Time taken to load ds_xc_sub: {elapsed} seconds")

# %%
# Test case 3 - OPEN_MFDATASET() + "pr" dataset (2 GB) + subsetting + `.load()`
# Result: ds.load() works fine with pr variable, but not with ua variable
# Notes: pr is 3D variable (time, lat, lon), ua is a 4D variable (time, lat, lon, plev).
# ------------------------------------------------------------------------------
filepaths_pr = [
    "/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/time-series/ERA5/pr_197901_201912.nc"
]
ds_pr = xc.open_mfdataset(
    filepaths_pr,
    add_bounds=["X", "Y", "T"],
    decode_times=True,
    use_cftime=True,
    coords="minimal",
    compat="override",
)

# %%
# pr dataset is ~2 GB without subsetting. There is no need to subset.
start_time = timeit.default_timer()
ds_pr.load()
elapsed = timeit.default_timer() - start_time
print(f"Time taken to load ds_xc_sub_0: {elapsed} seconds")
# %%

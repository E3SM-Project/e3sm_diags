"""
This dataset has a single time coordinate (2005, 1, 1) with units month. Attempting to add
bounds results and loading the dataset into memory results in:
`ValueError: Non-integer years and months are ambiguous and not currently supported.`.

Instead, we need to squeeze the time dimension to drop coordinates bounds before
loading the dataset into memory.
"""
# %%
import xcdat as xc

from e3sm_diags.driver.utils.dataset_xr import squeeze_time_dim

filepaths = [
    "/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/climatology/OMI-MLS/OMI-MLS_ANN_200501_201712_climo.nc"
]

# %%
ds = xc.open_mfdataset(filepaths, decode_times=True)

# ds = squeeze_time_dim(ds)
# %%

# %%
ds.load()

# %%

# %%
import numpy as np
import pandas as pd
import xarray as xr
import timeit

import dask.array as da

# %%
# Define the dimensions
time = 12
plev = 37
lat = 721
lon = 1440

# Create the data arrays using dask.
data = da.random.random(size=(time, plev, lat, lon), chunks=(12, 37, 721, 1440)).astype(
    np.float32
)

# Create the coordinates.
times = pd.date_range("2000-01-01", periods=time)
plevs = np.linspace(100000, 10, plev)
lats = np.linspace(-90, 90, lat)
lons = np.linspace(0, 360, lon, endpoint=False)

# Create the dataset and write out to a file.
ds = xr.Dataset(
    {"data": (["time", "plev", "lat", "lon"], data)},
    coords={"time": times, "plev": plevs, "lat": lats, "lon": lons},
)
# %%
ds.to_netcdf("dask_bottleneck.nc")

# %%
# Open the dataset.
ds_open = xr.open_mfdataset("dask_bottleneck.nc")

# %%
# Load the dataset into memory
start_time = timeit.default_timer()
ds.load()
end_time = timeit.default_timer()

print(f"Time taken to load the dataset: {end_time - start_time} seconds")


# %%

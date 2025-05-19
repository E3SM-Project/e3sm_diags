#%%
import numpy as np
import xarray as xr
import xcdat as xc

ds_b = xr.open_dataset("/lcrc/group/e3sm/diagnostics/observations/Atm/climatology/GPCP_v3.2/GPCP_v3.2_ANN_198301_202112_climo.nc")

# Check for existing lat bounds -- False
print("lat_bnds" in ds_b)  # False

# Add lat bounds and check the orientation
ds_b = ds_b.bounds.add_bounds(axis="Y")

#%%
# Print the first 5 latitude coordinates of ds_b
latitudes = ds_b["lat"].values
print("\nLatitude Coordinates (First 5):")
print(latitudes[:5])

# Check if latitude coordinates are ascending or descending
lat_order = "ascending" if np.all(latitudes[:-1] <= latitudes[1:]) else "descending"
print(f"Latitude coordinates are {lat_order}.")

# Print the first 5 latitude bounds of ds_b
lat_bnds = ds_b["lat_bnds"].values
print("\nLatitude Bounds (First 5):")
for i, bounds in enumerate(lat_bnds[:5]):
    print(f"Bounds {i + 1}: {bounds}")

# Check if latitude bounds are ascending or descending
bnds_order = "ascending" if np.all(lat_bnds[:, 0] <= lat_bnds[:, 1]) else "descending"
print(f"Latitude bounds are {bnds_order}.")

# %%

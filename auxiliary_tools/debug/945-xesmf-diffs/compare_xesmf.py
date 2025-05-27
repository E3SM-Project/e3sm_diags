
import cdms2
import numpy as np
import pandas as pd
import xarray as xr
import xcdat as xc
import xesmf as xe

from regrid2 import Regridder


def print_stats(*arrays, labels=None):
    """Prints statistical comparison of multiple arrays using a pandas DataFrame."""
    if labels is None:
        labels = [f"Array {i + 1}" for i in range(len(arrays))]
    elif len(labels) != len(arrays):
        raise ValueError("Number of labels must match the number of arrays.")

    stats = {
        "Min": [np.min(arr) for arr in arrays],
        "Max": [np.max(arr) for arr in arrays],
        "Mean": [np.mean(arr) for arr in arrays],
        "Sum": [np.sum(arr) for arr in arrays],
        "Std": [np.std(arr) for arr in arrays],
    }

    # Create a DataFrame from the stats dictionary
    df = pd.DataFrame(stats, index=labels)

    # Print the DataFrame
    print("\nStatistical Comparison:")
    print(df)

def regrid_with_xesmf(ds, output_grid, method="conservative_normed"):
    """Regrid a dataset using xESMF."""
    print("Regridding with xESMF...")
    print("--" * 20)
    # Print lat and lat_bnds ascending or descending with first 5 values
    lat_ascending = np.all(np.diff(ds["lat"].values) > 0)
    lat_bnds_ascending = np.all(np.diff(ds["lat_bnds"].values, axis=-1) > 0)

    print("Latitude is ascending:", lat_ascending)
    print("First 5 lat values:", ds["lat"].values[:5])
    print("Latitude bounds are ascending:", lat_bnds_ascending)
    print("First 5 lat_bnds values:", ds["lat_bnds"].values[:5])

    regridder = xe.Regridder(ds, output_grid, method=method)
    return regridder(ds["PRECT"])


# %%
# Load datasets
ds_a = xr.open_dataset(
    "/lcrc/group/e3sm/public_html/cdat-migration-fy24/25-02-14-branch-930-polar-after/polar/GPCP_v3.2/test.nc"
)
ds_b = xr.open_dataset(
    "/lcrc/group/e3sm/public_html/cdat-migration-fy24/25-02-14-branch-930-polar-after/polar/GPCP_v3.2/ref.nc"
)

# %%
# 1. xCDAT + xESMF (ascending lat, descending lat_bnds -- default)
# ----------------------------------------------------
output_grid_xesmf = ds_a.regridder.grid
ds_b1 = ds_b.copy(deep=True)
unaligned1 = regrid_with_xesmf(ds_b1, output_grid_xesmf)

# 2. xCDAT + xESMF (descending lat, ascending lat_bnds)
# ----------------------------------------------------
ds_b2 = ds_b.copy(deep=True).sortby("lat", ascending=False)
ds_b2["lat_bnds"].values = np.sort(ds_b2["lat_bnds"], axis=-1)
unaligned2 = regrid_with_xesmf(ds_b2, output_grid_xesmf)

# 3. xCDAT + xESMF (ascending lat, ascending lat_bnds)
# ----------------------------------------------------
ds_b3 = ds_b.copy(deep=True)
ds_b3["lat_bnds"].values = np.sort(ds_b3["lat_bnds"], axis=-1)
aligned1 = regrid_with_xesmf(ds_b3, output_grid_xesmf)

# 4. xCDAT + xESMF (descending lat, descending lat_bnds)
# ----------------------------------------------------
ds_b4 = ds_b.copy(deep=True)
ds_b4 = ds_b4.sortby("lat", ascending=False)
ds_b4["lat_bnds"].values = np.sort(ds_b4["lat_bnds"], axis=-1)[:, ::-1]
aligned2 = regrid_with_xesmf(ds_b4, output_grid_xesmf)

# 5. Desc lat, no lat_bnds
# ----------------------------------------------------
ds_b5 = ds_b.copy(deep=True).sortby("lat", ascending=False)
ds_b5 = ds_b5.drop_vars("lat_bnds")
ds_b5 = ds_b5.bounds.add_bounds(axis="Y")
nobounds1 = regrid_with_xesmf(ds_b5, output_grid_xesmf)

# 6. Asc lat, no lat_bnds
# ----------------------------------------------------
ds_b6 = ds_b.copy(deep=True).sortby("lat", ascending=True)
ds_b6 = ds_b6.drop_vars("lat_bnds")
ds_b6 = ds_b6.bounds.add_bounds(axis="Y")
nobounds2 = regrid_with_xesmf(ds_b6, output_grid_xesmf)

# %%
# Compare statistics
# ----------------------------------------------------
print_stats(
    unaligned1.values,
    unaligned2.values,
    aligned1.values,
    aligned2.values,
    nobounds1.values,
    nobounds2.values,
    labels=[
        "asc lat, desc lat_bnds",
        "desc lat, asc lat_bnds",
        "asc lat, asc lat_bnds",
        "desc lat, desc lat_bnds",
        "desc lat, no lat_bnds",
        "asc lat, no lat_bnds",

    ],
)

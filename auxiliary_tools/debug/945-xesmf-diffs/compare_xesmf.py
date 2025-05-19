
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
unaligned1 = regrid_with_xesmf(ds_b, output_grid_xesmf)

# 2. xCDAT + xESMF (descending lat, ascending lat_bnds)
# ----------------------------------------------------
ds_b2 = ds_b.copy(deep=True)
ds_b2 = ds_b2.sortby("lat", ascending=False)
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


# %%
# Compare statistics
# ----------------------------------------------------
print_stats(
    unaligned1.values,
    unaligned2.values,
    aligned1.values,
    aligned2.values,
    labels=[
        "asc lat, desc lat_bnds",
        "desc lat, asc lat_bnds",
        "asc lat, asc lat_bnds",
        "desc lat, desc lat_bnds",

    ],
)
"""
Statistical Comparison:
                              Min        Max      Mean           Sum       Std
Original                 0.162348  11.557144  1.398501  80553.664062  1.045462
PRECT (xesmf, unsorted)  0.165452   9.341431  1.423206  20494.160156  1.052469
PRECT(xesmf, sorted)     0.164364  10.046874  1.398573  20139.449219  1.039437
PRECT (cdat esmf,)       0.164364  10.046873  1.398573  20139.449219  1.039437
PRECT (regrid2)          0.164363  10.047869  1.398592  20139.730469  1.039466
"""

# %%

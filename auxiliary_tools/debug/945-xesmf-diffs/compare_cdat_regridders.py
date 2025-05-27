# %%
"""
This script compares regridding results between two methods:
1. `xesmf` (xarray-based regridding library).
2. `regrid2` (cdms2-based regridding library).

Key Steps:
1. Load datasets.
2. Perform regridding using `xesmf` with unsorted and sorted latitude bounds.
3. Perform regridding using `regrid2`.
4. Compare statistical differences in results.

Findings:
- Regridding results differ between `xesmf` and `regrid2` due to algorithmic differences.
- `xesmf` depends on having coordinates and coordinate bounds aligned.
- Statistical differences (e.g., min, max, mean, sum, std) highlight sensitivity to grid preparation and implementation.

conda create -n xcdat_cdat latest python xcdat=0.8.0 cdms2=3.1.5 ipykernel
conda activate xcdat_cdat
"""

# %%
import cdms2
import numpy as np
import pandas as pd
from regrid2 import Regridder
from regrid2.horizontal import extractBounds


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

#%%
def make_lat_descending(var):
    lat = var.getLatitude()
    lat_index = next(i for i, ax in enumerate(var.getAxisList()) if ax.id == lat.id)

    # Reverse latitude values
    lat_vals = lat[:][::-1]
    lat_reversed = cdms2.createAxis(lat_vals)
    lat_reversed.id = lat.id
    lat_reversed.units = lat.units
    lat_reversed.designateLatitude()

    # Reverse data along latitude axis
    slicer = [slice(None)] * var.ndim
    slicer[lat_index] = slice(None, None, -1)
    data_reversed = var[tuple(slicer)]

    # Replace the latitude axis in the axis list
    new_axes = list(var.getAxisList())
    new_axes[lat_index] = lat_reversed

    # Create new variable with updated latitude axis
    var_reversed = cdms2.createVariable(data_reversed, axes=new_axes, id=var.id)

    return var_reversed

def drop_bounds(var, axis_ids=("latitude",)):
    """
    Returns a copy of `var` with bounds removed from specified axes.
    """
    axes = []
    for ax in var.getAxisList():
        ax_copy = cdms2.createAxis(ax[:])
        ax_copy.id = ax.id
        ax_copy.units = getattr(ax, "units", "")
        if ax.id.lower() in axis_ids or ax.isLatitude() or ax.isLongitude():
            ax_copy.setBounds(None)
        axes.append(ax_copy)

    new_var = cdms2.createVariable(var[:], axes=axes, id=var.id)
    return new_var

# %%
# 1. CDAT + Regrid2 (ascending latitude, descending latitude bounds) -- -- default values, automatically sorted
# --------------------------------------------------------------------
# Convert xarray datasets to cdms2 variables
with (
    cdms2.open(
        "/lcrc/group/e3sm/public_html/cdat-migration-fy24/25-02-14-branch-930-polar-after/polar/GPCP_v3.2/test.nc"
    ) as f_a,
    cdms2.open(
        "/lcrc/group/e3sm/public_html/cdat-migration-fy24/25-02-14-branch-930-polar-after/polar/GPCP_v3.2/ref.nc"
    ) as f_b,
):
    var_a1 = f_a("PRECT")
    var_b1 = f_b("PRECT")

# Create regridder using regrid2
misaligned1 = Regridder(var_b1.getGrid(), var_a1.getGrid())(var_b1)

#%%
# 2. CDAT + Regrid2 (descending latitude, ascending latitude bounds)
# --------------------------------------------------------------------
# Convert xarray datasets to cdms2 variables
with (
    cdms2.open(
        "/lcrc/group/e3sm/public_html/cdat-migration-fy24/25-02-14-branch-930-polar-after/polar/GPCP_v3.2/test.nc"
    ) as f_a,
    cdms2.open(
        "/lcrc/group/e3sm/public_html/cdat-migration-fy24/25-02-14-branch-930-polar-after/polar/GPCP_v3.2/ref.nc"
    ) as f_b,
):
    var_a2 = f_a("PRECT")
    var_b2 = f_b("PRECT")

    var_a2 = make_lat_descending(var_a2)
    var_b2 = make_lat_descending(var_b2)


# Create regridder using regrid2
aligned = Regridder(var_b2.getGrid(), var_a2.getGrid())(var_b2)


# %%
# 3. CDAT + Regrid2 (ascending latitude, no latitude bounds)
# --------------------------------------------------------------------
# Convert xarray datasets to cdms2 variables
with (
    cdms2.open(
        "/lcrc/group/e3sm/public_html/cdat-migration-fy24/25-02-14-branch-930-polar-after/polar/GPCP_v3.2/test.nc"
    ) as f_a,
    cdms2.open(
        "/lcrc/group/e3sm/public_html/cdat-migration-fy24/25-02-14-branch-930-polar-after/polar/GPCP_v3.2/ref.nc"
    ) as f_b,
):
    var_a3 = f_a("PRECT")
    var_b3 = f_b("PRECT")

    var_a3 = drop_bounds(var_a3)
    var_b3 = drop_bounds(var_a3)


# Create regridder using regrid2
no_bnds1 = Regridder(var_b3.getGrid(), var_a3.getGrid())(var_b3)


# %%
# 4. CDAT + Regrid2 (ascending latitude, no latitude bounds)
# --------------------------------------------------------------------
# Convert xarray datasets to cdms2 variables
with (
    cdms2.open(
        "/lcrc/group/e3sm/public_html/cdat-migration-fy24/25-02-14-branch-930-polar-after/polar/GPCP_v3.2/test.nc"
    ) as f_a,
    cdms2.open(
        "/lcrc/group/e3sm/public_html/cdat-migration-fy24/25-02-14-branch-930-polar-after/polar/GPCP_v3.2/ref.nc"
    ) as f_b,
):
    var_a4 = f_a("PRECT")
    var_b4 = f_b("PRECT")

    var_a4 = make_lat_descending(var_a4)
    var_b4 = make_lat_descending(var_a4)

    var_a4 = drop_bounds(var_a4)
    var_b4 = drop_bounds(var_a4)


# Create regridder using regrid2
no_bnds2 = Regridder(var_b4.getGrid(), var_a4.getGrid())(var_b4)


# %%
# Compare statistics
# ----------------------------------------------------
print_stats(
    misaligned1,
    aligned,
    no_bnds1,
    no_bnds2,
    labels=[
        "asc lat, desc lat_bnds",
        "desc lat, desc lat_bnds",
        "asc lat, no lat_bnds",
        "desc lat, no lat_bnds",
    ],
)

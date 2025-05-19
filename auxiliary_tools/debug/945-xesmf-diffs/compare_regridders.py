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
- Sorting latitude bounds before regridding with `xesmf` impacts results.
- Statistical differences (e.g., min, max, mean, sum, std) highlight sensitivity to grid preparation and implementation.

conda create -n xcdat_cdat latest python xcdat=0.8.0 cdms2=3.1.5 ipykernel
conda activate xcdat_cdat
"""

# %%
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

# %%
# Load datasets
ds_a = xr.open_dataset(
    "/lcrc/group/e3sm/public_html/cdat-migration-fy24/25-02-14-branch-930-polar-after/polar/GPCP_v3.2/test.nc"
)
ds_b = xr.open_dataset(
    "/lcrc/group/e3sm/public_html/cdat-migration-fy24/25-02-14-branch-930-polar-after/polar/GPCP_v3.2/ref.nc"
)

# ds_b = xr.open_dataset("/lcrc/group/e3sm/diagnostics/observations/Atm/climatology/GPCP_v3.2/GPCP_v3.2_ANN_198301_202112_climo.nc")

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
# 1. xCDAT + xESMF (unsorted latitude bounds)
# ----------------------------------------------------
output_grid_xesmf = ds_a.regridder.grid
regridder_xesmf = xe.Regridder(ds_b, output_grid_xesmf, method="conservative_normed")
prect_xesmf_unsorted = regridder_xesmf(ds_b["PRECT"])

# 2. xCDAT + xESMF (sorted latitude bounds)
# ----------------------------------------------------
# Sort lat bounds in ascending order for ds_b
ds_b_sorted = ds_b.copy(deep=True)
ds_b_sorted["lat_bnds"].values = np.sort(ds_b_sorted["lat_bnds"], axis=-1)

output_grid_xesmf = ds_a.regridder.grid
regridder_xesmf = xe.Regridder(
    ds_b_sorted, output_grid_xesmf, method="conservative_normed"
)
prect_xesmf_sorted = regridder_xesmf(ds_b_sorted["PRECT"])


# %%
# 3. CDAT + ESMF (unsorted latitude bounds)
# ----------------------------------------------------
# Convert xarray datasets to cdms2 variables
with (
    cdms2.open(
        "/lcrc/group/e3sm/public_html/cdat-migration-fy24/25-02-14-branch-930-polar-after/polar/GPCP_v3.2/test.nc"
    ) as f_a,
    cdms2.open(
        "/lcrc/group/e3sm/public_html/cdat-migration-fy24/25-02-14-branch-930-polar-after/polar/GPCP_v3.2/ref.nc"
    ) as f_b,
):
    var_a = f_a("PRECT")
    var_b = f_b("PRECT")

# Create regridder using regrid2
prect_cdat_esmf = var_b.regrid(var_a.getGrid(), regridTool="esmf", regridMethod="conservative")

# %%
# 4. CDAT + Regrid2 (unsorted latitude bounds) -- automatically sorted
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
    var_a = f_a("PRECT")
    var_b = f_b("PRECT")

# Create regridder using regrid2
regrid2 = Regridder(var_b.getGrid(), var_a.getGrid())
prect_regrid2 = regrid2(var_b)


# %%
# Compare statistics
print_stats(
    ds_b["PRECT"].values,
    prect_xesmf_unsorted.values,
    prect_xesmf_sorted.values,
    prect_cdat_esmf.data,
    prect_regrid2,
    labels=["Original", "PRECT (xesmf, unsorted)", "PRECT(xesmf, sorted)", "PRECT (cdat esmf,)", "PRECT (regrid2)"],
)

#%%

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
# Compare floating point differences
def compare_arrays(array1, array2, rtol, atol, err_msg):
    """Compares two arrays using np.testing.assert_allclose and prints the result."""
    try:
        np.testing.assert_allclose(
            array1, array2, rtol=rtol, atol=atol, err_msg=err_msg
        )
        print(f"{err_msg}: Arrays are close within tolerance.")
    except AssertionError as e:
        print(e)


# %%
compare_arrays(
    prect_xesmf_sorted.values,
    prect_regrid2.data,
    rtol=1e-5,
    atol=0,
    err_msg="Differences found between xesmf sorted and regrid2 results",
)

"""
Not equal to tolerance rtol=1e-05, atol=0
Differences found between xesmf sorted and regrid2 results
Mismatched elements: 9004 / 14400 (62.5%)
Max absolute difference among violations: 0.00323868
Max relative difference among violations: 0.00062978
 ACTUAL: array([[2.356618, 2.652923, 2.578529, ..., 2.912792, 2.402574, 2.316819],
       [2.049907, 2.211725, 2.076397, ..., 2.646537, 2.309698, 2.113865],
       [1.979461, 2.010734, 1.869982, ..., 2.347431, 2.147887, 1.941646],...
 DESIRED: array([[2.35669 , 2.653266, 2.578602, ..., 2.913047, 2.402567, 2.316867],
       [2.049983, 2.21183 , 2.076715, ..., 2.646721, 2.309823, 2.113984],
       [1.97932 , 2.010875, 1.869887, ..., 2.347089, 2.14766 , 1.941615],...
"""

# %%
# 4. Regrid using xcdat's regrid2 implementation
# ----------------------------------------------------
# Perform regridding using xcdat's regridder
prect_xcdat_regrid2 = ds_b.regridder.horizontal(
    "PRECT",
    ds_a.regridder.grid,
    tool="regrid2",
)["PRECT"]

# %%
# Compare statistics including xcdat's regrid2
print_stats(
    prect_xesmf_unsorted.values,
    prect_xesmf_sorted.values,
    prect_regrid2,
    prect_xcdat_regrid2.values,
    labels=[
        "PRECT (xesmf, unsorted)",
        "PRECT (xesmf, sorted)",
        "PRECT (regrid2)",
        "PRECT (xcdat regrid2)",
    ],
)

"""
Stat       PRECT (xesmf, unsorted) PRECT (xesmf, sorted) PRECT (regrid2) PRECT (xcdat regrid2)
--------------------------------------------------------------------------------------------
Min        ...                     ...                   ...             ...
Max        ...                     ...                   ...             ...
Mean       ...                     ...                   ...             ...
Sum        ...                     ...                   ...             ...
Std        ...                     ...                   ...             ...
"""

# %%
# Compare floating point differences with xcdat's regrid2
compare_arrays(
    prect_xcdat_regrid2.values,
    prect_regrid2.data,
    rtol=1e-5,
    atol=0,
    err_msg="Differences found between xcdat regrid2 and regrid2 results",
)

compare_arrays(
    prect_xcdat_regrid2.values,
    prect_xesmf_sorted.values,
    rtol=1e-5,
    atol=0,
    err_msg="Differences found between xcdat regrid2 and xesmf sorted results",
)

# %%

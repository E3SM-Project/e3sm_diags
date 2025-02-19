"""
This script demonstrates the difference between regridding with and without
pre-existing lat/lon bounds.

It uses two datasets, one with lat/lon bounds and one without. The script
finds that the regridded arrays are not exactly the same, even though the
input arrays are the same.

The name of the bounds dimensions is different in the two datasets, which
causes the regridding to be different.
"""
#%%
import xcdat as xc
import numpy as np

ds_a = xc.open_dataset("/lcrc/group/e3sm/public_html/cdat-migration-fy24/25-02-14-branch-930-polar-after/polar/GPCP_v3.2/test.nc")
ds_b = xc.open_dataset("/lcrc/group/e3sm/public_html/cdat-migration-fy24/25-02-14-branch-930-polar-after/polar/GPCP_v3.2/ref.nc")

ds_a = ds_a.rename({'nbnd': 'bnds'})

# 1. Check the bounds dimension names
print("Bounds in ds_a:", ds_a["lon_bnds"].dims, ds_a["lat_bnds"].dims)
print("Bounds in ds_b:", ds_b["lon_bnds"].dims, ds_b["lat_bnds"].dims)

"""
Bounds in ds_a: ('lon', 'nbnd') ('lat', 'nbnd')
Bounds in ds_b: ('lon', 'bnds') ('lat', 'bnds')

Notice how the bounds dimensions are different ('nbnd' vs. 'bnd'.). This
seems to affect results when regridding.
"""
# 2. Check if bounds values are equal
np.testing.assert_allclose(ds_a["lon_bnds"].values, ds_b["lon_bnds"].values)
np.testing.assert_allclose(ds_a["lat_bnds"].values, ds_b["lat_bnds"].values)

# 2. Regrid with pre-existing lat/lon bounds. -- This one produces larger diffs compared to CDAT
output_grid = ds_a.regridder.grid
prect_a = ds_b.regridder.horizontal("PRECT", output_grid, tool="xesmf", method="conservative_normed")["PRECT"]

# 3. Regrid without pre-existng bounds. -- This one is close to CDAT.
ds_a_no_bnds = ds_a.drop_vars(["lon_bnds", "lat_bnds"])
ds_b_no_bnds = ds_b.drop_vars(["lon_bnds", "lat_bnds"])
output_grid_no_bnds = ds_a_no_bnds.regridder.grid
prect_b = ds_b_no_bnds.regridder.horizontal("PRECT", output_grid_no_bnds, tool="xesmf", method="conservative_normed")["PRECT"]

# 4. Print the stats -- they are "relatively" close but still different (which we shouldn't expect).
def print_stats(arr1, arr2, label1="Array 1", label2="Array 2"):
    stats = {
        "Min": (np.min(arr1), np.min(arr2)),
        "Max": (np.max(arr1), np.max(arr2)),
        "Mean": (np.mean(arr1), np.mean(arr2)),
        "Sum": (np.sum(arr1), np.sum(arr2)),
        "Std": (np.std(arr1), np.std(arr2)),
    }

    print(f"\n{'Stat':<10} {label1:<15} {label2:<15}")
    print("-" * 40)
    for stat, values in stats.items():
        print(f"{stat:<10} {values[0]:<15.6f} {values[1]:<15.6f}")

print_stats(prect_a.values, prect_b.values, "PRECT (with bounds)", "PRECT (no bounds)")

"""
Stat       PRECT (with bounds) PRECT (no bounds)
----------------------------------------
Min        0.165452        0.164364
Max        9.341431        10.046874
Mean       1.423206        1.398573
Sum        20494.160156    20139.449219
Std        1.052469        1.039437
"""

# 5. Compare arrays for floating point differences -- notice how they don't align
try:
    np.testing.assert_allclose(prect_a.values, prect_b.values, rtol=1e-5, atol=0)
except AssertionError as e:
    print("Arrays are not within rtol 1e-5")
    print(e)
else:
    print("Arrays are within rtol 1e-5")

"""
Arrays are not within rtol 1e-5

Not equal to tolerance rtol=1e-05, atol=0

Mismatched elements: 14392 / 14400 (99.9%)
Max absolute difference among violations: 3.5014582
Max relative difference among violations: 1.2256833
 ACTUAL: array([[2.257148, 2.883649, 2.568689, ..., 3.145134, 2.34432 , 2.263298],
       [2.288256, 2.331007, 2.508708, ..., 2.673229, 2.40844 , 2.271105],
       [1.978262, 2.111453, 1.774192, ..., 2.471611, 2.190455, 2.001195],...
 DESIRED: array([[2.356618, 2.652923, 2.578529, ..., 2.912792, 2.402574, 2.316819],
       [2.049907, 2.211725, 2.076397, ..., 2.646537, 2.309698, 2.113865],
       [1.979461, 2.010734, 1.869982, ..., 2.347431, 2.147887, 1.941646],...

"""
# %%

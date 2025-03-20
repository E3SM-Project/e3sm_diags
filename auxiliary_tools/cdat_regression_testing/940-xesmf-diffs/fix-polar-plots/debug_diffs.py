"""
conda create -n xcdat_cdat -c conda-forge -c xcdat cdms2 ipykernel
"""

#%%
import xarray as xr
import xcdat as xc # noqa: F401

# The variable in the dataset to compare.
VAR_KEY = "PRECT"

# Regrid using xCDAT + xESMF
# --------------------------------
# 1. Open the datasets with Xarray.
ds_test = xr.open_dataset("/lcrc/group/e3sm/public_html/cdat-migration-fy24/25-02-14-branch-930-polar/polar/GPCP_v3.2/GPCP_v3.2-PRECT-ANN-polar_N_test.nc")
ds_ref = xr.open_dataset("/lcrc/group/e3sm/public_html/cdat-migration-fy24/25-02-14-branch-930-polar/polar/GPCP_v3.2/GPCP_v3.2-PRECT-ANN-polar_N_ref.nc")

# 2. Regrid the reference variable to the test variable using xESMF and "conservative".
test_grid_xc = ds_test.regridder.grid
ds_ref_regrid = ds_ref.regridder.horizontal(
    VAR_KEY, test_grid_xc, tool="xesmf", method="conservative"
)

# 3. Get the difference between the test and regridded reference datasets.
ds_diff = ds_test.copy()
ds_diff[VAR_KEY] = ds_test[VAR_KEY] - ds_ref_regrid[VAR_KEY]

# 4. Extract the variables for easier comparison.
test_var_xc = ds_test[VAR_KEY].values
ref_var_xc_reg = ds_ref_regrid[VAR_KEY].values
diff_var_xc = ds_diff[VAR_KEY].values

#%%
##% Regrid using CDAT + ESMF
# --------------------------------
import cdms2 as cd

# 1. Open the datasets with cdms2.
test_var_cd = cd.open("/lcrc/group/e3sm/public_html/cdat-migration-fy24/25-02-14-branch-930-polar-main/polar/GPCP_v3.2/GPCP_v3.2-PRECT-ANN-polar_N_test.nc")(VAR_KEY)
ref_var_cd = cd.open("/lcrc/group/e3sm/public_html/cdat-migration-fy24/25-02-14-branch-930-polar-main/polar/GPCP_v3.2/GPCP_v3.2-PRECT-ANN-polar_N_ref.nc")(VAR_KEY)

# 2. Regrid the reference variable to the test variable using ESMF and "conservative".
test_grid_cd = test_var_cd.getGrid()
ref_var_cd_reg = ref_var_cd.regrid(
    test_grid_cd, regridTool="esmf", regridMethod="conservative"
)

diff_var_cd = test_var_cd - ref_var_cd_reg

#%%
def print_stats(arr1, arr2, label1="Array 1", label2="Array 2"):
    stats = {
        "Min": (np.min(arr1), np.min(arr2)),
        "Max": (np.max(arr1), np.max(arr2)),
        "Mean": (np.mean(arr1), np.mean(arr2)),
        "Std": (np.std(arr1), np.std(arr2)),
    }

    print(f"{'Stat':<10} {label1:<15} {label2:<15}")
    print("-" * 40)
    for stat, values in stats.items():
        print(f"{stat:<10} {values[0]:<15.6f} {values[1]:<15.6f}")

# %%
### Compare the results
import numpy as np
# 1. Check the grid to use for regridding are the same.
# Result: Both are the same
np.testing.assert_allclose(test_grid_xc.lat.values, test_grid_cd.getLatitude(), rtol=1e-5, atol=0)
np.testing.assert_allclose(test_grid_xc.lon.values, test_grid_cd.getLongitude(), rtol=1e-5, atol=0)

print("xCDAT Grid Latitude Values:", test_grid_xc.lat.values)
print("CDAT Grid Latitude Values:", test_grid_cd.getLatitude()[:])

print("xCDAT Grid Longitude Values:", test_grid_xc.lon.values)
print("CDAT Grid Longitude Values:", test_grid_cd.getLongitude()[:])

#%%
# 2. Check the test variables are the same
try:
    np.testing.assert_allclose(test_var_xc, test_var_cd.data, rtol=1e-5, atol=0)
except AssertionError as e:
    print("Arrays are not within relative tolerance (1e-5).")
    print(e)
else:
    print("Arrays are the within relative tolerance (1e-5).")


# Example usage
print_stats(test_var_xc, test_var_cd.data, label1="xCDAT Test", label2="CDAT Test")

#%%
# 3. Check the regridded ref variables are the same
try:
    np.testing.assert_allclose(ref_var_xc_reg, ref_var_cd_reg.data, rtol=1e-5, atol=0)
except AssertionError as e:
    print("Arrays are not within relative tolerance (1e-5).")
    print(e)
else:
    print("Arrays are the within relative tolerance (1e-5).")

"""
Mismatched elements: 14392 / 14400 (99.9%)
Max absolute difference among violations: 3.6998303
Max relative difference among violations: 1.2256833
 ACTUAL: array([[1.123808, 1.435736, 1.278921, ..., 1.565927, 1.16721 , 1.126871],
       [2.288256, 2.331007, 2.508708, ..., 2.673229, 2.40844 , 2.271105],
       [1.978262, 2.111453, 1.774192, ..., 2.471611, 2.190455, 2.001195],...
 DESIRED: array([[2.356618, 2.652923, 2.578529, ..., 2.912792, 2.402574, 2.316819],
       [2.049907, 2.211725, 2.076397, ..., 2.646537, 2.309698, 2.113865],
       [1.979461, 2.010734, 1.869982, ..., 2.347431, 2.147887, 1.941646],...
"""
print_stats(ref_var_xc_reg, ref_var_cd_reg.data, label1="xCDAT Ref", label2="CDAT Ref")

# %%
# 4. Check the regridded ref variables are the same
try:
    np.testing.assert_allclose(diff_var_xc, diff_var_cd.data, rtol=1e-5, atol=0)
except AssertionError as e:
    print("Arrays are not within relative tolerance (1e-5).")
    print(e)
else:
    print("Arrays are the within relative tolerance (1e-5).")

print_stats(diff_var_xc, diff_var_cd.data, label1="xCDAT Diff", label2="CDAT Diff")






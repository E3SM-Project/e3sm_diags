#%%
"""
Debug differences between Xarray and CDAT RMSE results for the SST variable.

Environment setup:
    conda create -n 940-xesmf-diffs-sst-rmse xarray xcdat cdms2 cdutil genutil xskillscore ipykernel
"""
import cdms2 as cd
import genutil
import numpy as np
import xarray as xr
import xcdat as xc
import xskillscore as xs

BASE_PATH = "/lcrc/group/e3sm/public_html/cdat-migration-fy24/25-02-19-branch-940-xesmf-diffs-mask-fix-sst-rmse/prov"
VAR_KEY ="SST"

ds1 = xr.open_dataset(f"{BASE_PATH}/sst_test_xcdat_conservative_normed.nc")
ds2 = xr.open_dataset(f"{BASE_PATH}/sst_ref_xcdat_conservative_normed.nc")

sst1 = cd.open(f"{BASE_PATH}/sst_test_cdat_bilinear.nc")["variable_16"]
sst2 = cd.open(f"{BASE_PATH}/sst_ref_cdat_bilinear.nc")["sst_CdmsRegrid"]

#%%
# Compare arrays for closeness, treating np.nan values as equal
try:
    np.testing.assert_allclose(ds1[VAR_KEY].values, sst1[:].filled(np.nan))
except AssertionError as e:
    print(e)
    # Get nan mismatch location counts between ds1 and sst1
    nan_mismatch_ds1_sst1 = np.sum(np.isnan(ds1[VAR_KEY].values) != np.isnan(sst1[:].filled(np.nan)))
    print(f"NaN mismatch count between ds1 and sst1: {nan_mismatch_ds1_sst1}")


try:
    np.testing.assert_allclose(ds2[VAR_KEY].values, sst2[:].filled(np.nan))
except AssertionError as e:
    print(e)

    # Get nan mismatch location counts between ds2 and sst2
    nan_mismatch_ds2_sst2 = np.sum(np.isnan(ds2[VAR_KEY].values) != np.isnan(sst2[:].filled(np.nan)))
    total_elements_ds2_sst2 = np.prod(ds2[VAR_KEY].values.shape)
    print(f"NaN mismatch count between ds2 and sst2: {nan_mismatch_ds2_sst2} out of {total_elements_ds2_sst2} elements")

"""
Not equal to tolerance rtol=1e-07, atol=0

nan location mismatch:
 ACTUAL: array([[nan, nan, nan, ..., nan, nan, nan],
       [nan, nan, nan, ..., nan, nan, nan],
       [nan, nan, nan, ..., nan, nan, nan],...
 DESIRED: array([[nan, nan, nan, ..., nan, nan, nan],
       [nan, nan, nan, ..., nan, nan, nan],
       [nan, nan, nan, ..., nan, nan, nan],...
NaN mismatch count between ds2 and sst2: 2522 out of 64800 elements
"""

#%%
# Compare statistics between ds1 and sst1
print("Comparing ds1 and sst1:")
print(f"Min: {np.nanmin(ds1[VAR_KEY].values):.6e} vs {np.nanmin(sst1[:].filled(np.nan)):.6e}")
print(f"Max: {np.nanmax(ds1[VAR_KEY].values):.6e} vs {np.nanmax(sst1[:].filled(np.nan)):.6e}")
print(f"Mean: {np.nanmean(ds1[VAR_KEY].values):.6e} vs {np.nanmean(sst1[:].filled(np.nan)):.6e}")
print(f"Sum: {np.nansum(ds1[VAR_KEY].values):.6e} vs {np.nansum(sst1[:].filled(np.nan)):.6e}")
print(f"Std: {np.nanstd(ds1[VAR_KEY].values):.6e} vs {np.nanstd(sst1[:].filled(np.nan)):.6e}")

# Compare statistics between ds2 and sst2
print("\nComparing ds2 and sst2:")
print(f"Min: {np.nanmin(ds2[VAR_KEY].values):.6e} vs {np.nanmin(sst2[:].filled(np.nan)):.6e}")
print(f"Max: {np.nanmax(ds2[VAR_KEY].values):.6e} vs {np.nanmax(sst2[:].filled(np.nan)):.6e}")
print(f"Mean: {np.nanmean(ds2[VAR_KEY].values):.6e} vs {np.nanmean(sst2[:].filled(np.nan)):.6e}")
print(f"Sum: {np.nansum(ds2[VAR_KEY].values):.6e} vs {np.nansum(sst2[:].filled(np.nan)):.6e}")
print(f"Std: {np.nanstd(ds2[VAR_KEY].values):.6e} vs {np.nanstd(sst2[:].filled(np.nan)):.6e}")

"""
Comparing ds1 and sst1:
Min: 1.035706e-01 vs 1.035706e-01
Max: 2.969695e+01 vs 2.969695e+01
Mean: 1.824388e+01 vs 1.824388e+01
Sum: 5.513300e+05 vs 5.513300e+05
Std: 8.615571e+00 vs 8.615571e+00

Comparing ds2 and sst2:
Min: -1.714309e+00 vs -1.714309e+00
Max: 2.963584e+01 vs 2.963584e+01
Mean: 1.535821e+01 vs 1.545739e+01
Sum: 6.200109e+05 vs 6.160696e+05
Std: 1.057837e+01 vs 1.054536e+01
"""

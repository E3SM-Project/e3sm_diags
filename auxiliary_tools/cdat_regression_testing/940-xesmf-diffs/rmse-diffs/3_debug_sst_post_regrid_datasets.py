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

"""
nan location mismatch:
 ACTUAL: array([[nan, nan, nan, ..., nan, nan, nan],
       [nan, nan, nan, ..., nan, nan, nan],
       [nan, nan, nan, ..., nan, nan, nan],...
 DESIRED: array([[nan, nan, nan, ..., nan, nan, nan],
       [nan, nan, nan, ..., nan, nan, nan],
       [nan, nan, nan, ..., nan, nan, nan],...
"""
try:
    np.testing.assert_allclose(ds2[VAR_KEY].values, sst2[:].filled(np.nan))
except AssertionError as e:
    print(e)

# Get nan mismatch location counts between ds1 and sst1
nan_mismatch_ds1_sst1 = np.sum(np.isnan(ds1[VAR_KEY].values) != np.isnan(sst1[:].filled(np.nan)))
print(f"NaN mismatch count between ds1 and sst1: {nan_mismatch_ds1_sst1}")

# Get nan mismatch location counts between ds2 and sst2
nan_mismatch_ds2_sst2 = np.sum(np.isnan(ds2[VAR_KEY].values) != np.isnan(sst2[:].filled(np.nan)))
total_elements_ds2_sst2 = np.prod(ds2[VAR_KEY].values.shape)
print(f"NaN mismatch count between ds2 and sst2: {nan_mismatch_ds2_sst2} out of {total_elements_ds2_sst2} elements")

#%%
# Compare statistics between ds1 and sst1
print("Comparing ds1 and sst1:")
print(f"Min: {np.nanmin(ds1[VAR_KEY].values)} vs {np.nanmin(sst1[:].filled(np.nan))}")
print(f"Max: {np.nanmax(ds1[VAR_KEY].values)} vs {np.nanmax(sst1[:].filled(np.nan))}")
print(f"Mean: {np.nanmean(ds1[VAR_KEY].values)} vs {np.nanmean(sst1[:].filled(np.nan))}")
print(f"Sum: {np.nansum(ds1[VAR_KEY].values)} vs {np.nansum(sst1[:].filled(np.nan))}")
print(f"Std: {np.nanstd(ds1[VAR_KEY].values)} vs {np.nanstd(sst1[:].filled(np.nan))}")

# Compare statistics between ds2 and sst2
print("\nComparing ds2 and sst2:")
print(f"Min: {np.nanmin(ds2[VAR_KEY].values)} vs {np.nanmin(sst2[:].filled(np.nan))}")
print(f"Max: {np.nanmax(ds2[VAR_KEY].values)} vs {np.nanmax(sst2[:].filled(np.nan))}")
print(f"Mean: {np.nanmean(ds2[VAR_KEY].values)} vs {np.nanmean(sst2[:].filled(np.nan))}")
print(f"Sum: {np.nansum(ds2[VAR_KEY].values)} vs {np.nansum(sst2[:].filled(np.nan))}")
print(f"Std: {np.nanstd(ds2[VAR_KEY].values)} vs {np.nanstd(sst2[:].filled(np.nan))}")

"""
Comparing ds1 and sst1:
Min: 0.10357055813074112 vs 0.10357055664064774
Max: 29.69695472717285 vs 29.696954345703148
Mean: 18.243877410888672 vs 18.243878141513044
Sum: 551330.0 vs 551329.9974365241
Std: 8.615571022033691 vs 8.615570625389605

Comparing ds2 and sst2:
Min: -1.7143093347549438 vs -1.7143093347549438
Max: 29.63583755493164 vs 29.63583755493164
Mean: 15.358208656311035 vs 15.358208656311035
Sum: 620010.875 vs 620010.875
Std: 10.578368186950684 vs 10.578368186950684
"""
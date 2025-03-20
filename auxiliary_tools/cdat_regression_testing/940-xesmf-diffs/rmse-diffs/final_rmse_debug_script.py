#%%
"""
This script compares the xskillscore and genutil RMSE functions to determine if
they produce similar results. Specifically, it fixes the weight generation
in xCDAT to produce the same weights as CDAT (normalizing and masking).

Results: After normalizing and masking the weight, the RMSE values are close.
    - Xarray RMSE (weighted): 1.7713504288524176
    - CDAT RMSE (weighted): 1.7713504288524178

Environment Setup:
  - conda create -n 940-xesmf-diffs-sst-rmse xarray xcdat cdms2 cdutil genutil xskillscore ipykernel
"""
import cdms2 as cd
import genutil
import numpy as np
import xarray as xr
import xcdat as xc
import xskillscore as xs
from genutil.averager import area_weights

BASE_PATH = "/lcrc/group/e3sm/public_html/cdat-migration-fy24/25-02-19-branch-940-xesmf-diffs-mask-fix-sst-rmse/prov"

VAR_KEY = "SST"
TEST_VAR_KEY ="variable_16"
REF_VAR_KEY ="sst_CdmsRegrid"

#%%
# CDAT setup.
# ----------
sst1_cd = cd.open(f"{BASE_PATH}/sst_test_cdat_conservative.nc")[TEST_VAR_KEY][:]
sst2_cd = cd.open(f"{BASE_PATH}/sst_ref_cdat_conservative.nc")[REF_VAR_KEY][:]

sst3_cd = cd.open(f"{BASE_PATH}/sst_test_cdat_bilinear.nc")[TEST_VAR_KEY][:]
sst4_cd = cd.open(f"{BASE_PATH}/sst_ref_cdat_bilinear.nc")[REF_VAR_KEY][:]

#%%
# Xarray setup.
# ------------
ds1_xr = xr.open_dataset(f"{BASE_PATH}/sst_test_xcdat_conservative_normed.nc")
ds2_xr = xr.open_dataset(f"{BASE_PATH}/sst_ref_xcdat_conservative_normed.nc")
ds3_xr = xr.open_dataset(f"{BASE_PATH}/sst_test_xcdat_bilinear.nc")
ds4_xr = xr.open_dataset(f"{BASE_PATH}/sst_ref_xcdat_bilinear.nc")
ds5_xr = xr.open_dataset(f"{BASE_PATH}/sst_test_xcdat_conservative_normed_bnds_fix.nc")
ds6_xr = xr.open_dataset(f"{BASE_PATH}/sst_ref_xcdat_conservative_normed_bnds_fix.nc")

ds1_xr = ds1_xr.bounds.add_missing_bounds()
ds2_xr = ds2_xr.bounds.add_missing_bounds()
ds3_xr = ds3_xr.bounds.add_missing_bounds()
ds4_xr = ds4_xr.bounds.add_missing_bounds()
ds5_xr = ds5_xr.bounds.add_missing_bounds()
ds6_xr = ds6_xr.bounds.add_missing_bounds()

# Extract the variables for re-use.
sst1_xr = ds1_xr[VAR_KEY]
sst2_xr = ds2_xr[VAR_KEY]
sst3_xr = ds3_xr[VAR_KEY]
sst4_xr = ds4_xr[VAR_KEY]
sst5_xr = ds5_xr[VAR_KEY]
sst6_xr = ds6_xr[VAR_KEY]

# Get CDAT weights.
weights_cd_bilinear = area_weights(sst3_cd, "xy")
weights_cd_cons = area_weights(sst1_cd, "xy")

# Get xCDAT weights for each regrid method.
weights_bilinear = ds3_xr.spatial.get_weights(["X", "Y"], data_var=VAR_KEY)
weights_cons = ds1_xr.spatial.get_weights(["X", "Y"], data_var=VAR_KEY)
weights_norm = (weights_cons * 0.5)/360
weights_norm = (weights_norm.where(ds1_xr[VAR_KEY].notnull(), np.nan))
weights_cons_bnds_fix = ds5_xr.spatial.get_weights(["X", "Y"], data_var=VAR_KEY)
#%%
# Results
# -------
result_cd_bilinear = genutil.statistics.rms(sst3_cd, sst4_cd, axis="xy", weights=weights_cd_bilinear)
result_cd_cons = genutil.statistics.rms(sst1_cd, sst2_cd, axis="xy", weights=weights_cd_cons)

result_xr_bilinear = xs.rmse(sst3_xr, sst4_xr, dim=["lat", "lon"], weights=weights_bilinear, skipna=True)
result_xr_cons_norm = xs.rmse(sst1_xr, sst2_xr, dim=["lat", "lon"], weights=weights_cons, skipna=True)
result_xr_cons_norm_bnds_fix = xs.rmse(sst5_xr, sst6_xr, dim=["lat", "lon"], weights=weights_cons_bnds_fix, skipna=True)

# Calculate unweighted RMSE
result_cd_bilinear_unweighted = genutil.statistics.rms(sst3_cd, sst4_cd, axis="xy")
result_cd_cons_unweighted = genutil.statistics.rms(sst1_cd, sst2_cd, axis="xy")

result_xr_bilinear_unweighted = xs.rmse(sst3_xr, sst4_xr, dim=["lat", "lon"], skipna=True)
result_xr_cons_norm_unweighted = xs.rmse(sst1_xr, sst2_xr, dim=["lat", "lon"], skipna=True)
result_xr_cons_norm_bnds_fix_unweighted = xs.rmse(sst5_xr, sst6_xr, dim=["lat", "lon"], skipna=True)

print("v2.12.1 CDAT RMSE (unweighted, bilinear):", result_cd_bilinear_unweighted)
print("v2.12.1 CDAT RMSE (unweighted, conservative_normed):", result_cd_cons_unweighted)
print("v3.0.0 Xarray RMSE (unweighted, bilinear):", result_xr_bilinear_unweighted.values)
print("v3.0.0 Xarray RMSE (unweighted, conservative_normed):", result_xr_cons_norm_unweighted.values)
print("v3.0.0 Xarray RMSE (unweighted, conservative_normed, bnds fix):", result_xr_cons_norm_bnds_fix_unweighted.values)

"""
v2.12.1 CDAT RMSE (unweighted, bilinear): 1.601377718177987
v2.12.1 CDAT RMSE (unweighted, conservative_normed): 1.6013035930408892
v3.0.0 Xarray RMSE (unweighted, bilinear): 1.6013004
v3.0.0 Xarray RMSE (unweighted, conservative_normed): 1.9159542
v3.0.0 Xarray RMSE (unweighted, conservative_normed, bnds fix): 1.6013036
"""

print("\nv2.12.1 CDAT RMSE (weighted, bilinear):", result_cd_bilinear)
print("v2.12.1 CDAT RMSE (weighted, conservative):", result_cd_cons)
print("v3.0.0 Xarray RMSE (weighted, bilinear):", result_xr_bilinear.values)
print("v3.0.0 Xarray RMSE (weighted normalized, conservative_normed):", result_xr_cons_norm.values)
print("v3.0.0 Xarray RMSE (weighted normalized, conservative_normed, bnds fix):", result_xr_cons_norm_bnds_fix.values)
"""
v2.12.1 CDAT RMSE (weighted, bilinear): 1.476555924046852
v2.12.1 CDAT RMSE (weighted, conservative): 1.4764821102067773
v3.0.0 Xarray RMSE (weighted, bilinear): 1.4763235405423747
v3.0.0 Xarray RMSE (weighted normalized, conservative_normed): 1.7713504288524176
v3.0.0 Xarray RMSE (weighted normalized, conservative_normed, bnds fix): 1.4764820110242953
"""

# %%

# Add mask variable to Xarray datasets.
ds1_xr["mask"] = xr.where(~np.isnan(ds1_xr[VAR_KEY]), 1, 0)
ds2_xr["mask"] = xr.where(~np.isnan(ds2_xr[VAR_KEY]), 1, 0)
ds3_xr["mask"] = xr.where(~np.isnan(ds3_xr[VAR_KEY]), 1, 0)
ds4_xr["mask"] = xr.where(~np.isnan(ds4_xr[VAR_KEY]), 1, 0)
ds5_xr["mask"] = xr.where(~np.isnan(ds5_xr[VAR_KEY]), 1, 0)
ds6_xr["mask"] = xr.where(~np.isnan(ds6_xr[VAR_KEY]), 1, 0)

# Get xCDAT weights for each regrid method with mask applied.
weights_bilinear_masked = ds3_xr.spatial.get_weights(["X", "Y"], data_var=VAR_KEY) * ds3_xr["mask"]
weights_cons_masked = ds1_xr.spatial.get_weights(["X", "Y"], data_var=VAR_KEY) * ds1_xr["mask"]
weights_cons_bnds_fix_masked = ds5_xr.spatial.get_weights(["X", "Y"], data_var=VAR_KEY) * ds5_xr["mask"]

# Normalize the weights.
weights_bilinear_masked = weights_bilinear_masked / weights_bilinear_masked.sum()
weights_cons_masked = weights_cons_masked / weights_cons_masked.sum()
weights_cons_bnds_fix_masked = weights_cons_bnds_fix_masked / weights_cons_bnds_fix_masked.sum()

# Calculate RMSE with masked weights.
result_xr_bilinear_masked = xs.rmse(sst3_xr, sst4_xr, dim=["lat", "lon"], weights=weights_bilinear_masked, skipna=True)
result_xr_cons_masked = xs.rmse(sst1_xr, sst2_xr, dim=["lat", "lon"], weights=weights_cons_masked, skipna=True)
result_xr_cons_bnds_fix_masked = xs.rmse(sst5_xr, sst6_xr, dim=["lat", "lon"], weights=weights_cons_bnds_fix_masked, skipna=True)

# Calculate unweighted RMSE with masked weights.
result_xr_bilinear_masked_unweighted = xs.rmse(sst3_xr, sst4_xr, dim=["lat", "lon"], skipna=True)
result_xr_cons_masked_unweighted = xs.rmse(sst1_xr, sst2_xr, dim=["lat", "lon"], skipna=True)
result_xr_cons_bnds_fix_masked_unweighted = xs.rmse(sst5_xr, sst6_xr, dim=["lat", "lon"], skipna=True)

print("\nv3.0.0 Xarray RMSE (weighted, bilinear, masked):", result_xr_bilinear_masked.values)
print("v3.0.0 Xarray RMSE (weighted, conservative_normed, masked):", result_xr_cons_masked.values)
print("v3.0.0 Xarray RMSE (weighted normalized, conservative_normed, bnds fix, masked):", result_xr_cons_bnds_fix_masked.values)
"""
v3.0.0 Xarray RMSE (weighted, bilinear, masked): 1.4763235405423747
v3.0.0 Xarray RMSE (weighted, conservative_normed, masked): 1.7713504288524176
v3.0.0 Xarray RMSE (weighted normalized, conservative_normed, bnds fix, masked): 1.4764820110242953
"""

print("\nv3.0.0 Xarray RMSE (unweighted, bilinear, masked):", result_xr_bilinear_masked_unweighted.values)
print("v3.0.0 Xarray RMSE (unweighted, conservative_normed, masked):", result_xr_cons_masked_unweighted.values)
print("v3.0.0 Xarray RMSE (unweighted, conservative_normed, bnds fix, masked):", result_xr_cons_bnds_fix_masked_unweighted.values)
"""
v3.0.0 Xarray RMSE (unweighted, bilinear, masked): 1.6013004
v3.0.0 Xarray RMSE (unweighted, conservative_normed, masked): 1.9159542
v3.0.0 Xarray RMSE (unweighted, conservative_normed, bnds fix, masked): 1.6013036
"""

# %%

"""
This script compares the xskillscore and genutil RMSE functions to determine if
they produce similar results.

Results:
  - RMSE functions generate similar results for both Xarray and CDAT, weighted and unweighted.
  - Xarray RMSE (unweighted): 1.601377718177987
    CDAT RMSE (unweighted): 1.601377718177987
  - Xarray RMSE (weighted): 1.4765559240468518
    CDAT RMSE (weighted): 1.476555924046852

Environment Setup:
  - conda create -n 940-xesmf-diffs-sst-rmse xarray xcdat cdms2 cdutil genutil xskillscore ipykernel
"""
#%%
import cdms2 as cd
import genutil
import numpy as np
import xarray as xr
import xcdat as xc
import xskillscore as xs
from genutil.averager import area_weights

BASE_PATH = "/lcrc/group/e3sm/public_html/cdat-migration-fy24/25-02-19-branch-940-xesmf-diffs-mask-fix-sst-rmse/prov"
VAR_KEY ="SST"

# Xarray setup.
ds1_xr = xr.open_dataset(f"{BASE_PATH}/sst_test_cdat_bilinear.nc")
ds2_xr = xr.open_dataset(f"{BASE_PATH}/sst_ref_cdat_bilinear.nc")
sst1_xr = ds1_xr["variable_16"]
sst2_xr = ds2_xr["sst_CdmsRegrid"]

#%%
# CDAT setup.
sst1_cd = cd.open(f"{BASE_PATH}/sst_test_cdat_bilinear.nc")["variable_16"][:]
sst2_cd = cd.open(f"{BASE_PATH}/sst_ref_cdat_bilinear.nc")["sst_CdmsRegrid"][:]

#%%
# Unweighted
# --------
result_xr = xs.rmse(sst1_xr, sst2_xr, dim=["lat", "lon"], weights=None, skipna=True)
result_cd = genutil.statistics.rms(sst1_cd, sst2_cd, axis="xy", weights=None)

np.testing.assert_allclose(result_xr, result_cd, rtol=1e-5, atol=0) # True
print("Xarray RMSE (unweighted):", result_xr.values)
print("CDAT RMSE (unweighted):", result_cd)
"""
Xarray RMSE (unweighted): 1.601377718177987
CDAT RMSE (unweighted): 1.601377718177987
"""

#%%
# Weighted
# -----------
ds1_xr = ds1_xr.bounds.add_missing_bounds()
ds2_xr = ds2_xr.bounds.add_missing_bounds()
weights = ds1_xr.spatial.get_weights(["X","Y"], data_var="variable_16")
weights = weights.fillna(0)

result_xr_w = xs.rmse(sst1_xr, sst2_xr, dim=["lat", "lon"], weights=weights, skipna=True)
result_cd_w = genutil.statistics.rms(sst1_cd, sst2_cd, axis="xy", weights="generate")

np.testing.assert_allclose(result_xr, result_cd, rtol=1e-5, atol=0) # True
print("Xarray RMSE (weighted):", result_xr_w.values)
print("CDAT RMSE (weighted):", result_cd_w)
"""
Xarray RMSE (weighted): 1.4765559240468518
CDAT RMSE (weighted): 1.476555924046852
"""
# %%

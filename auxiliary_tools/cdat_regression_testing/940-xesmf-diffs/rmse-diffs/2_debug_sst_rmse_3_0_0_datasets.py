"""
This script compares the xskillscore and genutil RMSE functions to identify
any discrepancies in their results.

Purpose:
  - Determine if the RMSE function is causing differences in results.
  - If results are similar, investigate earlier operations for discrepancies.

Results:
  - Unweighted RMSE values are similar for both Xarray and CDAT:
    - Xarray: 1.9159542
    - CDAT: 1.9159544
  - Weighted RMSE values differ:
    - Xarray: 1.7713504288524176
    - CDAT: 2.267541010250707
  - The difference in weighted RMSE values may be due to different `np.nan`
    positioning in the post-regridded datasets, affecting the weights.

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
# ------------
ds1_xr = xr.open_dataset(f"{BASE_PATH}/sst_test_xcdat_conservative_normed.nc")
ds2_xr = xr.open_dataset(f"{BASE_PATH}/sst_ref_xcdat_conservative_normed.nc")
ds1_xr = ds1_xr.bounds.add_missing_bounds()
ds2_xr = ds2_xr.bounds.add_missing_bounds()

# Extract the variables for re-use.
sst1_xr = ds1_xr[VAR_KEY]
sst2_xr = ds2_xr[VAR_KEY]

# CDAT setup.
# ----------
sst1_cd = cd.open(f"{BASE_PATH}/sst_test_xcdat_conservative_normed.nc")[VAR_KEY][:]
sst2_cd = cd.open(f"{BASE_PATH}/sst_ref_xcdat_conservative_normed.nc")[VAR_KEY][:]

# Update the mask and set the fill_value to 1e+20
sst1_cd = np.ma.masked_invalid(sst1_cd)
sst2_cd = np.ma.masked_invalid(sst2_cd)
sst1_cd.set_fill_value(1e+20)
sst2_cd.set_fill_value(1e+20)

# Update the name of the x and y axis for the TransientVariable sst1 and sst2 to lat and lon
sst1_cd.getAxis(0).id = 'lat'
sst1_cd.getAxis(1).id = 'lon'
sst2_cd.getAxis(0).id = 'lat'
sst2_cd.getAxis(1).id = 'lon'

#%%
# Unweighted RMSE
# ---------------
result_xr = xs.rmse(sst1_xr, sst2_xr, dim=["lat", "lon"], weights=None, skipna=True)
result_cd = genutil.statistics.rms(sst1_cd, sst2_cd, axis="xy", weights=None)

np.testing.assert_allclose(result_xr, result_cd, rtol=1e-5, atol=0) # True
print("Xarray RMSE (unweighted):", result_xr.values)
print("CDAT RMSE (unweighted):", result_cd)
"""
Xarray RMSE (unweighted): 1.9159542
CDAT RMSE (unweighted): 1.9159544
"""

#%%
# Weighted RMSE
# --------------
weights = ds1_xr.spatial.get_weights(["X","Y"], data_var=VAR_KEY)
weights = weights.fillna(0)

result_xr_w = xs.rmse(sst1_xr, sst2_xr, dim=["lat", "lon"], weights=weights, skipna=True)
result_cd_w = genutil.statistics.rms(sst1_cd, sst2_cd, axis="xy", weights="generate")

np.testing.assert_allclose(result_xr, result_cd, rtol=1e-5, atol=0) # True
print("Xarray RMSE (weighted):", result_xr_w.values)
print("CDAT RMSE (weighted):", result_cd_w)
"""
Xarray RMSE (weighted): 1.7713504288524176
CDAT RMSE (weighted): 2.267541010250707
"""


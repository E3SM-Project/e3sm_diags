#%%
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
import cdms2 as cd
import genutil
import numpy as np
import xarray as xr
import xcdat as xc
import xskillscore as xs
from genutil.averager import area_weights

BASE_PATH = "/lcrc/group/e3sm/public_html/cdat-migration-fy24/25-02-19-branch-940-xesmf-diffs-mask-fix-sst-rmse/prov"

TEST_VAR_KEY ="variable_16"
REF_VAR_KEY ="sst_CdmsRegrid"

# Xarray setup.
# ------------
ds1_xr = xr.open_dataset(f"{BASE_PATH}/sst_test_cdat_bilinear.nc")
ds2_xr = xr.open_dataset(f"{BASE_PATH}/sst_ref_cdat_bilinear.nc")
ds1_xr = ds1_xr.bounds.add_missing_bounds()
ds2_xr = ds2_xr.bounds.add_missing_bounds()

# Extract the variables for re-use.
sst1_xr = ds1_xr[TEST_VAR_KEY]
sst2_xr = ds2_xr[REF_VAR_KEY]

# CDAT setup.
# ----------
sst1_cd = cd.open(f"{BASE_PATH}/sst_test_cdat_bilinear.nc")[TEST_VAR_KEY][:]
sst2_cd = cd.open(f"{BASE_PATH}/sst_ref_cdat_bilinear.nc")[REF_VAR_KEY][:]

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

print("Xarray RMSE (unweighted):", result_xr.values)
print("CDAT RMSE (unweighted):", result_cd)
"""
Xarray RMSE (unweighted): 1.601377718177987
CDAT RMSE (unweighted): 1.601377718177987
"""

#%%
# Weighted RMSE (different weights)
# ------------------------------
weights = ds1_xr.spatial.get_weights(["X","Y"], data_var=TEST_VAR_KEY)
weights = weights.fillna(0)

result_xr_w = xs.rmse(sst1_xr, sst2_xr, dim=["lat", "lon"], weights=weights, skipna=True)
result_cd_w = genutil.statistics.rms(sst1_cd, sst2_cd, axis="xy", weights="generate")

print("Xarray RMSE (weighted):", result_xr_w.values)
print("CDAT RMSE (weighted):", result_cd_w)
"""
Xarray RMSE (weighted): 1.4765559240468518
CDAT RMSE (weighted): 1.8924341633197894
"""

#%%
# Weighted RMSE (use same weights from xCDAT)
# -------------------------------------------
weights = ds1_xr.spatial.get_weights(["X","Y"], data_var=TEST_VAR_KEY)
weights = weights.fillna(0)

result_xr_w = xs.rmse(sst1_xr, sst2_xr, dim=["lat", "lon"], weights=weights, skipna=True)
weights_mv2 = cd.MV2.array(weights.values.T)
weights_mv2.getAxis(0).id = 'lat'
weights_mv2.getAxis(1).id = 'lon'
weights_mv2.getAxis(0).id = 'lat'
weights_mv2.getAxis(1).id = 'lon'

result_cd_w = genutil.statistics.rms(sst1_cd, sst2_cd, axis="xy", weights=weights_mv2)

print("Xarray RMSE (weighted):", result_xr_w.values)
print("CDAT RMSE (weighted):", result_cd_w)
"""
Xarray RMSE (weighted): 1.4765559240468518
CDAT RMSE (weighted): 1.4765559240468518
"""

#%%
# Weighted RMSE (use same weights from CDAT)
# -------------------------------------------
weights_cd = area_weights(sst1_cd, "xy")
weights_cd_xr = xr.DataArray(weights_cd, dims=["lat", "lon"], coords={"lat": sst1_xr.lat, "lon": sst1_xr.lon})

result_xr_w = xs.rmse(sst1_xr, sst2_xr, dim=["lat", "lon"], weights=weights_cd_xr, skipna=True)
result_cd_w = genutil.statistics.rms(sst1_cd, sst2_cd, axis="xy", weights=weights_cd)

print("Xarray RMSE (weighted):", result_xr_w.values)
print("CDAT RMSE (weighted):", result_cd_w)
"""
Xarray RMSE (weighted): 1.8924341633197892
CDAT RMSE (weighted): 1.8924341633197894
"""


# %%

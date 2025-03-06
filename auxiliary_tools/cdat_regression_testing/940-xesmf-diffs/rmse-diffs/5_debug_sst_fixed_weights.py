#%%
"""
This script compares the xskillscore and genutil RMSE functions to determine if
they produce similar results. Specifically, it fixes the weight generation
in xCDAT to produce the same weights as CDAT (normalizing and masking).

Results: After normalizing and masking the weight, the RMSE values are close.
    - Xarray RMSE (weighted): 1.4765559240468518
    - CDAT RMSE (weighted): 1.476555924046852

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

# Weighted RMSE (normalize the xCDAT weights with masking)
# --------------------------------------------------------
# Get area weights.
weights_fixed = ds1_xr.spatial.get_weights(["X", "Y"], data_var=TEST_VAR_KEY)
# Normalize the weights.
weights_norm = (weights_fixed * 0.5)/360
# Mask the weights.
weights_norm_masked = (weights_norm.where(ds1_xr[TEST_VAR_KEY].notnull(), np.nan)).T

result_xr_w = xs.rmse(sst1_xr, sst2_xr, dim=["lat", "lon"], weights=weights_norm_masked.T, skipna=True)

weights_cd = area_weights(sst1_cd, "yx")
result_cd_w = genutil.statistics.rms(sst1_cd, sst2_cd, axis="xy", weights=weights_cd)

print("Xarray RMSE (weighted):", result_xr_w.values)
print("CDAT RMSE (weighted):", result_cd_w)
"""
Xarray RMSE (weighted): 1.4765559240468518
CDAT RMSE (weighted): 1.476555924046852
"""

# %%

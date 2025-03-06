"""
This script compares the weights generated by xCDAT and CDAT.

xCDAT: `spatial.get_weights()`
    - computes the relative weight of each rectilinear grid cell, but doesn't do any
      normalizing or masking.
CDAT: `genutil.averager.area_weights()`
    - computes the relative weight of each rectilinear grid cell with normalizing
      and masking.
    - https://github.com/CDAT/cdms/blob/3f8c7baa359f428628a666652ecf361764dc7b7a/Lib/grid.py#L445-L463
"""
#%%
import cdms2 as cd
import numpy as np
import xarray as xr
import xcdat as xc # noqa: F401
from genutil.averager import area_weights

BASE_PATH = "/lcrc/group/e3sm/public_html/cdat-migration-fy24/25-02-19-branch-940-xesmf-diffs-mask-fix-sst-rmse/prov"
# VAR_KEY = "SST"
VAR_KEY = "variable_16"

# CDAT Weights via `area_weights()`
# --------------------------------------
sst1_cd = cd.open(f"{BASE_PATH}/sst_test_cdat_bilinear.nc")[VAR_KEY][:]
weights_cd = area_weights(sst1_cd, "xy")


#%% Breakdown of `area_weights` function
# --------------------------------------
dsgr = sst1_cd.getGrid()

# Return normalized area weights
# https://github.com/CDAT/cdms/blob/3f8c7baa359f428628a666652ecf361764dc7b7a/Lib/grid.py#L445-L463
latwts, lonwts = dsgr.getWeights()

# np.outer computes the outer product of two vectors. The outer product of vectors a and b results in a matrix where each element (i, j) is the product of a[i] and b[j].
wt = np.outer(np.array(latwts), np.array(lonwts))
wt_var = cd.createVariable(
    np.ma.masked_array(wt,np.ma.getmask(sst1_cd)),
    axes=sst1_cd.getAxisList()
)

initial_order = "yx"
result = wt_var(order=initial_order)

#%%
# xCDAT Weights to align with `area_weights()`
# --------------------------------------------
ds1_xr = xr.open_dataset(f"{BASE_PATH}/sst_test_cdat_bilinear.nc")
ds1_xr = ds1_xr.bounds.add_missing_bounds()

# Get area weights.
weights = ds1_xr.spatial.get_weights(["X", "Y"], data_var=VAR_KEY)

# Normalize the weights.
weights_norm = (weights * 0.5)/360

# Mask the weights.
weights_norm_masked = weights_norm.where(ds1_xr[VAR_KEY].notnull(), np.nan)

# %%
import copy

import cdms2
import numpy as np
import xarray as xr
import xcdat as xc

from e3sm_diags.derivations.acme import mask_by  # noqa: F401
from e3sm_diags.driver.utils.regrid import _drop_unused_ilev_axis

VAR_KEY = "ts"
MASK_VAR_KEY = "LANDFRAC"
LOWER_MASK_LIMIT = 0.65


# %%
# 1. Get the mask file being used on both branches
mask_path = "/global/cfs/cdirs/e3sm/e3sm_diags/postprocessed_e3sm_v2_data_for_e3sm_diags/20210528.v2rc3e.piControl.ne30pg2_EC30to60E2r2.chrysalis/climatology/rgr/20210528.v2rc3e.piControl.ne30pg2_EC30to60E2r2.chrysalis_ANN_005101_006012_climo.nc"
data_path = "/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/climatology/ERA5/ERA5_ANN_197901_201912_climo.nc"

# Data
ds = xr.open_dataset(data_path)
ds = _drop_unused_ilev_axis(ds)
output_grid = ds.regridder.grid

# Mask
ds_mask = xr.open_dataset(mask_path)
ds_mask_new = _drop_unused_ilev_axis(ds_mask)

ds_mask_regrid = ds_mask_new.regridder.horizontal(
    MASK_VAR_KEY,
    output_grid,
    tool="xesmf",
    method="bilinear",
)

# Get the land masking condition and apply it
land_sea_mask = ds_mask_regrid[MASK_VAR_KEY]
cond = land_sea_mask > LOWER_MASK_LIMIT

masked_var = ds[VAR_KEY].where(cond=cond, drop=False)

# %%
var = cdms2.open(data_path)(VAR_KEY)
var_mask = cdms2.open(mask_path)(MASK_VAR_KEY)
var_mask_rg = var_mask.regrid(
    var.getGrid(),
    regridTool="esmf",
    regridMethod="bilinear",
)

# Apply the land masking condition
var_new = copy.deepcopy(var)
var_new.mask = LOWER_MASK_LIMIT > var_mask_rg

# %%

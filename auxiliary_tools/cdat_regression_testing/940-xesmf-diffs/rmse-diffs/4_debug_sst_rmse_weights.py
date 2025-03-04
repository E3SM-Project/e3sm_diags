#%%
import cdms2 as cd
import numpy as np
import xarray as xr
import xcdat as xc # noqa: F401
from genutil.averager import area_weights

BASE_PATH = "/lcrc/group/e3sm/public_html/cdat-migration-fy24/25-02-19-branch-940-xesmf-diffs-mask-fix-sst-rmse/prov"
# VAR_KEY = "SST"
VAR_KEY = "variable_16"

# CDAT Weights
sst1_cd = cd.open(f"{BASE_PATH}/sst_test_cdat_bilinear.nc")[VAR_KEY][:]
weights_cd = area_weights(sst1_cd, "xy")
regular_array = weights_cd.filled(np.nan)
regular_array = regular_array[~np.isnan(regular_array)]


# xCDAT Weights
ds1_xr = xr.open_dataset(f"{BASE_PATH}/sst_test_cdat_bilinear.nc")
ds1_xr = ds1_xr.bounds.add_missing_bounds()
weights_xc = ds1_xr.spatial.get_weights(["X", "Y"], data_var=VAR_KEY)
weights_xc = weights_xc.where(ds1_xr[VAR_KEY].notnull())
# weights_xc = weights_xc.fillna(0)

weights_xc.values[~np.isnan(weights_xc.values)]

# %%

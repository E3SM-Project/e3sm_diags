# %%
from glob import glob

import xarray as xr
import xcdat as xc

args = {
    "decode_times": False,
    "add_bounds": ["X", "Y"],
    "coords": "minimal",
    "compat": "override",
    "chunks": "auto",
}

# %%
filepath = "/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/climatology/AOD_550/AOD_550_[0-1][0-9]_*climo.nc"
filepaths = glob(filepath)
ds = xc.open_mfdataset(filepaths, **args)

# array([ 15.5101,  46.03  ,  76.5498, 107.57  , 137.5895, 168.6097,
#       198.6292, 229.6494, 259.6689, 290.6891, 321.7142, 351.2334])
ds.time.values

# %%
ds1 = xc.open_mfdataset(filepaths[0], **args)
ds2 = xc.open_mfdataset(filepaths[1], **args)

# 'days since 2000-03-01'
ds1.time.units
# 'days since 2000-03-01'
ds2.time.units

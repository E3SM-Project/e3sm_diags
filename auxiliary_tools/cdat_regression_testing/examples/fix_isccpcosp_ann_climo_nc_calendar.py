# %%


import xarray as xr
import xcdat as xc


filepath = (
    "/global/cfs/cdirs/e3sm/e3sm_diags/obs_for_e3sm_diags/climatology/ISCCPCOSP/"
    "ISCCPCOSP_ANN_climo.nc"
)

ds = xc.open_dataset(filepath, decode_times=True)

# %%
ds.time.values[:] = 150

# %%
ds.to_netcdf("ISCCPCOSP_ANN_climo.nc")

# %%
ds = xc.open_dataset("ISCCPCOSP_ANN_climo.nc")
# %%
from datetime import datetime
from dateutil import parser
from dateutil import relativedelta as rd

import numpy as np

flat_offsets = np.array([150.5])
ref_date = "1983-06"
units_type = "months"

ref_datetime: datetime = parser.parse(ref_date, default=datetime(2000, 1, 1))


# ValueError: Non-integer years and months are ambiguous and not currently supported.
times = np.array(
    [
        ref_datetime + rd.relativedelta(**{units_type: offset})
        for offset in flat_offsets
    ],
    dtype="object",
)

# %%

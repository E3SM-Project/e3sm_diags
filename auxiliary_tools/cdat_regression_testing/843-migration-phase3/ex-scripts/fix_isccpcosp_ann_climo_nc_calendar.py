# %%
from datetime import datetime

import numpy as np
import xarray as xr
import xcdat as xc
from dateutil import parser
from dateutil import relativedelta as rd

# %%
# Reproduce error: ValueError: Non-integer years and months are ambiguous and not currently supported.
# -------------------------------------------------------------------------------------------

filepath = "/global/cfs/cdirs/e3sm/e3sm_diags/obs_for_e3sm_diags/climatology/ISCCPCOSP/ISCCPCOSP_ANN_climo.nc"
ds = xc.open_dataset(filepath, decode_times=True)


# %%
# Debug -- issue is that the time value is a float representing middle of the month (150.5).
# relativedelta expects a time value at the beginning of the month (e.g,. 150).
# -----

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
# Fix error -- replace 150.5 (middle month float) with 150 (integer month)
# -----------------------------------------------------------------------------
ds = xc.open_dataset(filepath, decode_times=False)
ds.time.values[:] = 150
ds.to_netcdf(filepath)

# %%

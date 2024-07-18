# %%
from glob import glob

import xarray as xr
import xcdat as xc

args = {
    "decode_times": True,
    "add_bounds": ["X", "Y"],
    "coords": "minimal",
    "compat": "override",
    "chunks": "auto",
}

filepath = "/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/climatology/MERRA2_Aerosols/MERRA2_Aerosols_[0-1][0-9]_*climo.nc"
paths = sorted(glob(filepath))

# filepath 1: '/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/climatology/MERRA2_Aerosols/MERRA2_Aerosols_01_198001_202101_climo.nc'
ds_fp1 = xc.open_mfdataset(paths[0], **args)
# filepath 2: '/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/climatology/MERRA2_Aerosols/MERRA2_Aerosols_02_198002_202102_climo.nc'
ds_fp2 = xc.open_mfdataset(paths[1], **args)

# Check the units -- different
# ----------------------------
# 'minutes since 1980-01-01 00:30:00'
ds_fp1.time.attrs["units"]
# 'minutes since 1980-02-01 00:30:00'
ds_fp2.time.attrs["units"]

# Check the time values -- same
# ----------------------------
# 10782720
ds_fp1.time.values[0]
# 10782720
ds_fp2.time.values[0]


# %%
import cftime


def _get_encoded_time(var: xr.DataArray) -> xr.DataArray:
    """Convert `cftime` datetime objects to encoded float objects.

    This function is useful for plotters that plot the time axis, which
    requires time coordinates to be encoded as floats if arithmetic is
    performed.

    Parameters
    ----------
    var : xr.DataArray
        The variable.

    Returns
    -------
    xr.DataArray
        The encoded time coordinates.
    """
    time_coords = xc.get_dim_coords(var, axis="T")

    for index, time_val in enumerate(time_coords):
        new_time_val = cftime.date2num(
            time_val.item(),
            time_val.encoding["units"],
            calendar=time_val.encoding["calendar"],
        )
        time_coords.values[index] = new_time_val
        time_coords[time_coords.name].values[index] = new_time_val

    return time_coords


# %%
new_time = _get_encoded_time(ds_fp1)

# %%

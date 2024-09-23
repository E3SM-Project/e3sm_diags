import xarray as xr
import xcdat as xc

filepath = "/global/cfs/cdirs/e3sm/e3sm_diags/obs_for_e3sm_diags/climatology/ISCCPCOSP/ISCCPCOSP_ANN_climo.nc"

# %%
ds_xc = xc.open_dataset(filepath)
# ValueError: Non-integer years and months are ambiguous and not currently supported.

ds_xr = xr.open_dataset(filepath)
# ValueError: Failed to decode variable 'time': unable to decode time units 'months since 1983-06' # with 'the default calendar'. Try opening your dataset with decode_times=False or installing
# cftime if it is not installed.

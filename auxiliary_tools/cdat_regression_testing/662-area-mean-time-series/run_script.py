# python -m auxiliary_tools.cdat_regression_testing.792-lat-lon-run-script.792_lat_lon_run_script
from auxiliary_tools.cdat_regression_testing.base_run_script import run_set

SET_NAME = "area_mean_time_series"
SET_DIR = "662-area-mean-time-series"
# CFG_PATH: str | None = None
CFG_PATH: str | None = "/global/u2/v/vo13/E3SM-Project/e3sm_diags/auxiliary_tools/cdat_regression_testing/662-area-mean-time-series/run.cfg"
MULTIPROCESSING = False

# %%
run_set(SET_NAME, SET_DIR, CFG_PATH, MULTIPROCESSING)


# # %%
# import xarray as xr
# import xcdat as xc

# filepath = "/global/cfs/cdirs/e3sm/e3sm_diags/postprocessed_e3sm_v2_data_for_e3sm_diags/20210528.v2rc3e.piControl.ne30pg2_EC30to60E2r2.chrysalis/time-series/rgr/FLUT_005101_006012.nc"
# time_slice = slice("51-01-15", "60-12-15", None)

# # %%
# ds = xr.open_dataset(filepath, decode_times=True, use_cftime=True)
# time_dim = xc.get_dim_keys(ds, axis="T")

# # %%
# ds_subset = ds.sel({time_dim: time_slice}).squeeze()

# # %%
# from xarray.coding.cftimeindex import parse_iso8601_like

# # %%
# parse_iso8601_like("51-01-01")

# # %%
# from datetime import datetime

# # %%
# year = "51"

# year.zfill(1)

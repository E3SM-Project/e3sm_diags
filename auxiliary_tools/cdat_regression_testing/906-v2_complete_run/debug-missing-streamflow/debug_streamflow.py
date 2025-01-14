# %%
import xarray as xr

fp = "/global/cfs/cdirs/e3sm/e3sm_diags/postprocessed_e3sm_v2_data_for_e3sm_diags/20210528.v2rc3e.piControl.ne30pg2_EC30to60E2r2.chrysalis/time-series/rgr/RIVER_DISCHARGE_OVER_LAND_LIQ_005101_006012.nc"

ds1 = xr.open_mfdataset(fp)

# %%

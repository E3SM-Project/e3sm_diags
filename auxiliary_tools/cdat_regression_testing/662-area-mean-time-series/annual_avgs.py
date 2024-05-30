# %%
import cdms2
import cdutil
import xcdat as xc

# %%
# CDAT
# Time coordinates at the first day of the month (1850-01-01)
fn1 = "/global/cfs/cdirs/e3sm/e3sm_diags/test_model_data_for_acme_diags/time-series/E3SM_v1/FSNS_185001_201312.nc"
var_key = "FSNS"

ds_cdat2 = cdms2.open(fn1)
var = ds_cdat2(var_key)

var_avg = cdutil.averager(var, axis="xy")
var_avg_year = cdutil.YEAR(var_avg)

print(var_avg_year.getTime().asComponentTime()[-3:-1])
# [2011-7-2 12:0:0.0, 2012-7-2 12:0:0.0]

# %%
# xCDAT
ds_xcdat = xc.open_dataset(fn1, decode_times=True)
ds_xcdat_avg = ds_xcdat.spatial.average(var_key, axis=["X", "Y"])
ds_xcdat_avg_year = ds_xcdat_avg.temporal.group_average(var_key, freq="year")

var_avg_year_xc = ds_xcdat_avg_year[var_key]

print(var_avg_year_xc.time.values[-3:-1])
# [cftime.DatetimeNoLeap(2012, 1, 1, 0, 0, 0, 0, has_year_zero=True)
#  cftime.DatetimeNoLeap(2013, 1, 1, 0, 0, 0, 0, has_year_zero=True)]

# %%
import cdms2
import xcdat as xc

# Parameters
# Time coordinates at the first day of the month (1850-01-01)
fn1 = "/global/cfs/cdirs/e3sm/e3sm_diags/test_model_data_for_acme_diags/time-series/E3SM_v1/FSNS_185001_201312.nc"
var = "FSNS"
start_time = "2011-01-01"
end_time = "2014-01-01"

# "co" - YYYY-01-01 +1 for the end year
slice_flag = "co"


ds_cdat = cdms2.open(fn1)

# With slice -- nothing changes because "co" allows the right side to be open,
# The actual end time (2013-12-01) is within the slice end time (2014-01-01).
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
var_s = ds_cdat(var, time=(start_time, end_time, slice_flag))(squeeze=1)

time_s = var_s.getTime()
time_s_dt = time_s.asComponentTime()

"""
[2011-1-1 0:0:0.0,
 2011-2-1 0:0:0.0,
 2011-3-1 0:0:0.0,
 2011-4-1 0:0:0.0,
 2011-5-1 0:0:0.0,
 2011-6-1 0:0:0.0,
 2011-7-1 0:0:0.0,
 2011-8-1 0:0:0.0,
 2011-9-1 0:0:0.0,
 2011-10-1 0:0:0.0,
 2011-11-1 0:0:0.0,
 2011-12-1 0:0:0.0,
 2012-1-1 0:0:0.0,
 2012-2-1 0:0:0.0,
 2012-3-1 0:0:0.0,
 2012-4-1 0:0:0.0,
 2012-5-1 0:0:0.0,
 2012-6-1 0:0:0.0,
 2012-7-1 0:0:0.0,
 2012-8-1 0:0:0.0,
 2012-9-1 0:0:0.0,
 2012-10-1 0:0:0.0,
 2012-11-1 0:0:0.0,
 2012-12-1 0:0:0.0,
 2013-1-1 0:0:0.0,
 2013-2-1 0:0:0.0,
 2013-3-1 0:0:0.0,
 2013-4-1 0:0:0.0,
 2013-5-1 0:0:0.0,
 2013-6-1 0:0:0.0,
 2013-7-1 0:0:0.0,
 2013-8-1 0:0:0.0,
 2013-9-1 0:0:0.0,
 2013-10-1 0:0:0.0,
 2013-11-1 0:0:0.0,
 2013-12-1 0:0:0.0]
"""


# 36 2013-12-1 0:0:0.0
print(len(time_s), time_s_dt[-1])
# [59799. 59829.]
print(time_s.getBounds()[-1])


# ~~~~~~~

# xCDAT
ds_xcdat = xc.open_dataset(fn1, decode_times=True)

ds_xcdat_s = ds_xcdat.sel(time=slice(start_time, end_time))

# 35 2013-12-01 00:00:00
print(len(ds_xcdat_s.time), ds_xcdat_s.time.values[-1])

# [cftime.DatetimeNoLeap(2013, 11, 1, 0, 0, 0, 0, has_year_zero=True)
# cftime.DatetimeNoLeap(2013, 12, 1, 0, 0, 0, 0, has_year_zero=True)]
print(ds_xcdat_s.time_bnds.values[-1])


# No slice
# RESULT: Adds a coordinate to the end ("ccn" is default)
# ~~~~~~~
var_ns = ds_cdat(var, time=(start_time, end_time))(squeeze=1)
time_ns = var_ns.getTime()
time_ns_dt = time_ns.asComponentTime()

"""
[2011-1-1 0:0:0.0,
 2011-2-1 0:0:0.0,
 2011-3-1 0:0:0.0,
 2011-4-1 0:0:0.0,
 2011-5-1 0:0:0.0,
 2011-6-1 0:0:0.0,
 2011-7-1 0:0:0.0,
 2011-8-1 0:0:0.0,
 2011-9-1 0:0:0.0,
 2011-10-1 0:0:0.0,
 2011-11-1 0:0:0.0,
 2011-12-1 0:0:0.0,
 2012-1-1 0:0:0.0,
 2012-2-1 0:0:0.0,
 2012-3-1 0:0:0.0,
 2012-4-1 0:0:0.0,
 2012-5-1 0:0:0.0,
 2012-6-1 0:0:0.0,
 2012-7-1 0:0:0.0,
 2012-8-1 0:0:0.0,
 2012-9-1 0:0:0.0,
 2012-10-1 0:0:0.0,
 2012-11-1 0:0:0.0,
 2012-12-1 0:0:0.0,
 2013-1-1 0:0:0.0,
 2013-2-1 0:0:0.0,
 2013-3-1 0:0:0.0,
 2013-4-1 0:0:0.0,
 2013-5-1 0:0:0.0,
 2013-6-1 0:0:0.0,
 2013-7-1 0:0:0.0,
 2013-8-1 0:0:0.0,
 2013-9-1 0:0:0.0,
 2013-10-1 0:0:0.0,
 2013-11-1 0:0:0.0,
 2013-12-1 0:0:0.0,
 2014-1-1 0:0:0.0]
"""

# 37 2014-1-1 0:0:0.0
print(len(time_ns), time_ns_dt[-1])
# [59829. 59860.]
print(time_ns.getBounds()[-1])

import os
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.parameter.diurnal_cycle_parameter import DiurnalCycleParameter
from e3sm_diags.parameter.tropical_subseasonal_parameter import (
    TropicalSubseasonalParameter,
)
from e3sm_diags.run import runner

param = CoreParameter()

param.reference_data_path = (
    "/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/climatology"
)

test_base_path = "/pscratch/sd/c/chengzhu/ne256pg2_ne256pg2.F20TR-SCREAMv1.rainfrac1.spanc1000.auto2700.acc150.n0128"
ref_climo = "/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/climatology/"
ref_ts = "/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/time-series"
param.test_data_path = f"{test_base_path}/rgr/climo"
param.test_name = "1ma_ne30pg2.AVERAGE.nmonths_x1"
param.seasons = ["ANN", "DJF", "MAM", "JJA", "SON"]
# param.save_netcdf = True

prefix = "/global/cfs/cdirs/e3sm/www/zhang40/tests/eamxx"
param.results_dir = os.path.join(prefix, "eamxx_ne256_0520_trop")
params = [param]

trop_param = TropicalSubseasonalParameter()
trop_param.test_data_path = f"{test_base_path}/rgr/ts_daily"
trop_param.short_test_name = "Daily_3hi_ne30pg2"
trop_param.test_start_yr = 1996
trop_param.test_end_yr = 2001

# Obs
trop_param.reference_data_path = ref_ts
trop_param.ref_start_yr = 2001
trop_param.ref_end_yr = 2010

params.append(trop_param)

dc_param = DiurnalCycleParameter()
dc_param.test_data_path = f"{test_base_path}/rgr/climo_diurnal_3hrly"
dc_param.test_name = "3hi_ne30pg2.INSTANT.nhours_x3"
dc_param.short_test_name = "3hi_ne30pg2"
# Plotting diurnal cycle amplitude on different scales. Default is True
dc_param.normalize_test_amp = False

# Obs
dc_param.reference_data_path = ref_climo

params.append(dc_param)

runner.sets_to_run = ["tropical_subseasonal"]

runner.run_diags(params)

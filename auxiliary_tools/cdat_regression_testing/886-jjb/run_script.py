import os
from e3sm_diags.parameter.tropical_subseasonal_parameter import (
    TropicalSubseasonalParameter,
)
from e3sm_diags.run import runner

param = TropicalSubseasonalParameter()

param.reference_data_path = (
    "/global/cfs/cdirs/e3sm/e3sm_diags/obs_for_e3sm_diags/time-series"
)
# param.reference_data_path = '/global/cfs/cdirs/e3sm/chengzhu/e3sm_diags_zppy_test_complete_run_output/v2.LR.historical_0101_20240130/post/atm/180x360_aave/ts/daily/15yr'
param.test_data_path = "/global/cfs/cdirs/e3sm/chengzhu/e3sm_diags_zppy_test_complete_run_output/v2.LR.historical_0101_20240130/post/atm/180x360_aave/ts/daily/15yr"
param.test_name = "E3SMv2"
param.results_dir = "/global/cfs/cdirs/e3sm/www/cdat-migration-fy24/886-jjb"
# param.run_type = "model_vs_model"
# param.ref_name = 'E3SMv2'
param.test_start_yr = "2000"
param.test_end_yr = "2000"
param.ref_start_yr = "2001"
param.ref_end_yr = "2001"
param.save_netcdf = True

runner.sets_to_run = ["tropical_subseasonal"]
runner.run_diags([param])

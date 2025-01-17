import os

from e3sm_diags.parameter.tropical_subseasonal_parameter import (
    TropicalSubseasonalParameter,
)
from e3sm_diags.run import runner

param = TropicalSubseasonalParameter()

param.reference_data_path = '/global/cfs/cdirs/e3sm/e3sm_diags/obs_for_e3sm_diags/time-series'
#param.reference_data_path = '/global/cfs/cdirs/e3sm/chengzhu/e3sm_diags_zppy_test_complete_run_output/v2.LR.historical_0101_20240130/post/atm/180x360_aave/ts/daily/15yr'
param.test_data_path = '/global/cfs/cdirs/e3sm/chengzhu/e3sm_diags_zppy_test_complete_run_output/v2.LR.historical_0101_20240130/post/atm/180x360_aave/ts/daily/15yr'
param.test_name = 'E3SMv2'
prefix = '/Users/zhang40/Documents/repos/e3sm_diags/auxiliary_tools/tropical_subseasonal_diags/data1'
prefix = '/global/cfs/cdirs/e3sm/www/chengzhu/tests/tropical_diags_subsetting'

param.results_dir = os.path.join(prefix, 'tropical_variability_model_obs_refine')
#param.run_type = "model_vs_model"
#param.ref_name = 'E3SMv2'
param.test_start_yr = '2000'
param.test_end_yr = '2000'
param.ref_start_yr = '2001'
param.ref_end_yr = '2001'
param.save_netcdf = True

runner.sets_to_run = ['tropical_subseasonal']
runner.run_diags([param])



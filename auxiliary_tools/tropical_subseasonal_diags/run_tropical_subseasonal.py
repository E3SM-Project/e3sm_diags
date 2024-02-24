import os
from e3sm_diags.parameter.tropical_subseasonal_parameter import TropicalSubseasonalParameter
from e3sm_diags.run import runner

param = TropicalSubseasonalParameter()

param.reference_data_path = '/Users/zhang40/Documents/e3sm_diags_data/e3sm_diags_test_data/E3SM_v2_daily'
param.test_data_path = '/Users/zhang40/Documents/e3sm_diags_data/e3sm_diags_test_data/E3SM_v2_daily'
param.test_name = 'E3SMv2'
prefix = '/Users/zhang40/Documents/repos/e3sm_diags/auxiliary_tools/tropical_subseasonal_diags/data1'
param.results_dir = os.path.join(prefix, 'tropical_variability')
param.test_start_yr = '2010'
param.test_end_yr = '2014'
param.ref_start_yr = '2010'
param.ref_end_yr = '2014'

runner.sets_to_run = ['tropical_subseasonal']
runner.run_diags([param])

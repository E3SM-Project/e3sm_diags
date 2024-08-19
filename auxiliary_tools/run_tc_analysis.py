import os
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.parameter.tc_analysis_parameter import TCAnalysisParameter
from e3sm_diags.run import runner

param = CoreParameter()

#param.test_data_path = '/Users/zhang40/Documents/ACME_simulations/E3SM_v1_H1_6hourly_TC/test'
param.test_data_path = '/global/cfs/cdirs/e3sm/e3sm_diags/postprocessed_e3sm_v2_data_for_e3sm_diags/20210528.v2rc3e.piControl.ne30pg2_EC30to60E2r2.chrysalis/1_2/1_2'
param.test_name = '20210528.v2rc3e.piControl.ne30pg2_EC30to60E2r2.chrysalis' 
param.test_data_path = '/global/cfs/cdirs/e3sm/e3sm_diags/postprocessed_e3sm_v2_data_for_e3sm_diags/extendedOutput.v3.LR.historical_0101/tc-analysis'
param.test_name = 'extendedOutput.v3.LR.historical_0101' 
param.short_test_name = 'v3.LR.historical_threshold_1.0'
param.reference_data_path = '/global/cfs/cdirs/e3sm/e3sm_diags/postprocessed_e3sm_v2_data_for_e3sm_diags/20210528.v2rc3e.piControl.ne30pg2_EC30to60E2r2.chrysalis/0_3/0_3'
param.ref_name = '20210528.v2rc3e.piControl.ne30pg2_EC30to60E2r2.chrysalis'
param.short_ref_name = 'threshold_0.3'
param.run_type = 'model_vs_model'
prefix = '/global/cfs/cdirs/e3sm/www/zhang40/tests'
#param.multiprocessing = True
#param.num_workers = 4
param.results_dir = os.path.join(prefix, 'tc_analysis_model_model')
#param.test_timeseries_input = True
#param.ref_timeseries_input = True
param.test_start_yr = '2000'
param.test_end_yr = '2014'
param.ref_start_yr = '0051'
param.ref_end_yr = '0060'

runner.sets_to_run = ['tc_analysis']
runner.run_diags([param])


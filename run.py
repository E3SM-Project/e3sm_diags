# TODO: Remove before merging
import os
from acme_diags.parameter.core_parameter import CoreParameter
from acme_diags.run import runner

# `pip install . && python run.py`

param = CoreParameter()
# # Acme1
# param.test_data_path = '/Users/zhang40/Documents/ACME_simulations/E3SM_v1'
# param.reference_data_path = '/Users/zhang40/Documents/ACME_simulations/obs_for_e3sm_diags/time-series'
# # Cori
# param.test_data_path = '/global/project/projectdirs/acme/acme_diags/test_model_data_for_acme_diags/time-series/E3SM_v1'
# param.reference_data_path = '/global/project/projectdirs/acme/acme_diags/obs_for_e3sm_diags/time-series'
# Mac
# p = '/Users/forsyth2/projectAIMS/data/'
# p = '/Users/forsyth2/projectAIMS/data/test_model_data_for_acme_diags'
# p = '/Users/forsyth2/projectAIMS/data/test_model_data_for_acme_diags/climatology'
# p = '/Users/forsyth2/projectAIMS/data/test_model_data_for_acme_diags/time-series'
# p = '/Users/forsyth2/projectAIMS/data/test_model_data_for_acme_diags/time-series/E3SM_v1'
p = '/Users/forsyth2/projectAIMS/data/test_model_data_for_acme_diags/time-series/E3SM_v1_amip'
# p = '/Users/forsyth2/projectAIMS/data/obs_for_e3sm_diags'
# p = '/Users/forsyth2/projectAIMS/data/obs_for_e3sm_diags/climatology'
# p = '/Users/forsyth2/projectAIMS/data/obs_for_e3sm_diags/climatology/ERA-Interim'
# p = '/Users/forsyth2/projectAIMS/data/obs_for_e3sm_diags/time-series'
# p = '/Users/forsyth2/projectAIMS/data/obs_for_e3sm_diags/time-series/ERA-Interim'
p = '/Users/forsyth2/projectAIMS/data/wind_vector_data'
p = '/Users/forsyth2/projectAIMS/data/test_model_data_for_acme_diags/climatology/'
# param.test_timeseries_input = True
# param.test_start_yr = 2008
# param.test_end_yr = 2013
# param.ref_timeseries_input = True
# param.ref_start_yr = 2008
# param.ref_end_yr = 2013
param.test_data_path = p
param.reference_data_path = p
param.test_data_name = '20161118.beta0.F1850COSP.ne30_ne30.edison'
param.ref_data_name = '20161118.beta0.F1850COSP.ne30_ne30.edison'
param.test_name = '20161118.beta0.F1850COSP.ne30_ne30.edison'
param.ref_name = '20161118.beta0.F1850COSP.ne30_ne30.edison'

#param.test_name = 'e3sm_v1'
#param.test_name = 'E3SM_v1'
#param.ref_name = 'E3SM_v1 (ref)'
prefix = '/Users/forsyth2/projectAIMS/E3SM-Project/e3sm_diags'
param.results_dir = os.path.join(prefix, 'wind_vectors_results')
param.output_file = 'wind_vectors'
#param.output_format_subplot = ['png']
param.save_netcdf = True

runner.sets_to_run = ['lat_lon_vector']
runner.run_diags([param])

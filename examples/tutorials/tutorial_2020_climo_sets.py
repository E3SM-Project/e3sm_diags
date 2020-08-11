import os
from acme_diags.run import runner
from acme_diags.parameter.core_parameter import CoreParameter

param = CoreParameter()

param.reference_data_path = '/global/cfs/cdirs/e3sm/acme_diags/obs_for_e3sm_diags/climatology'
param.test_data_path = '/global/cfs/cdirs/e3sm/acme_diags/test_model_data_for_acme_diags/climatology/'
param.test_name = '20161118.beta0.FC5COSP.ne30_ne30.edison'
param.seasons = ["ANN"]
param.multiprocessing = True
param.num_workers = 24

prefix = '/global/cfs/cdirs/e3sm/www/zhang40/tutorial2020'
param.results_dir = os.path.join(prefix, 'climo_sets_haswell_unified')

runner.sets_to_run = ['lat_lon','zonal_mean_xy', 'zonal_mean_2d', 'polar', 'cosp_histogram', 'meridional_mean_2d']
runner.run_diags([param])


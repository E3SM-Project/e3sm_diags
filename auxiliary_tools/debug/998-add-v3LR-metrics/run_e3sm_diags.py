import os
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.run import runner

param = CoreParameter()

param.reference_data_path = '/global/cfs/cdirs/e3sm/e3sm_diags/obs_for_e3sm_diags/climatology/'
param.test_data_path = '/global/cfs/cdirs/e3sm/e3sm_diags/test_model_data_for_acme_diags/climatology/'

param.reference_data_path = '/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/climatology/'
param.test_data_path = '/global/cfs/cdirs/e3sm/chengzhu/tests/zppy_example_v3/v3.LR.amip_0101/post/atm/180x360_aave/clim/10yr/'
param.test_name = 'v3.LR.amip_0101'
param.seasons = ["ANN","DJF", "MAM", "JJA", "SON"] #will run,if comment out"

prefix = '/global/cfs/cdirs/e3sm/www/chengzhu/tests/test_zppy/'
param.results_dir = os.path.join(prefix, 'lat_lon_v3_metrics')
# Use the following if running in parallel:
param.multiprocessing = True
param.num_workers = 24

# Use below to run all core sets of diags:
#runner.sets_to_run = ['lat_lon','zonal_mean_xy', 'zonal_mean_2d', 'polar', 'cosp_histogram', 'meridional_mean_2d']
# Use below to run lat_lon map only:
runner.sets_to_run = ['lat_lon']
runner.run_diags([param])

import os
from e3sm_diags.parameter.lat_lon_land_parameter import LatLonLandParameter # Note the change.
from e3sm_diags.run import runner

param = LatLonLandParameter()
param.test_data_path = '/lcrc/group/e3sm2/ac.zhang40/E3SMv3/v3.LR.piControl_land_ilamb/post/lnd/native/clim/50yr'
param.test_name = 'v3.LR.piControl' 
#param.short_test_name = 'alpha20_rrm_test'
param.reference_data_path = '/lcrc/group/e3sm2/ac.zhang40/E3SMv3/v3.LR.piControl_land_ilamb/post/lnd/native/clim/50yr'
param.ref_name = 'v3.LR.piControl' 

param.run_type = 'model_vs_model'
prefix = '/lcrc/group/e3sm/public_html/diagnostic_output/ac.zhang40/tests/999-all-land-var'
param.seasons = ["ANN"]
param.multiprocessing = True
#param.num_workers = 16
param.results_dir = os.path.join(prefix, 'lat_lon_test_land_model_vs_model')
runner.sets_to_run = ['lat_lon_land']   # Note the change
runner.run_diags([param])

import os

import numpy

from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.parameter.lat_lon_land_parameter import LatLonLandParameter
from e3sm_diags.run import runner

short_name = "v3.LR.historical_0051"

param = CoreParameter()

# Model
param.test_name = "v3.LR.historical_0051"


# Ref

# Output dir
param.results_dir = "/lcrc/group/e3sm/public_html/ac.tvo/zppy_weekly_comprehensive_v3_output/test_zppy_pr651_20250115/v3.LR.historical_0051/post/"

# Additional settings
param.run_type = "model_vs_model"
param.diff_title = "Difference"
param.multiprocessing = True
param.num_workers = 8
param.seasons = ["ANN"]
# param.fail_on_incomplete = True
params = [param]

# Model land
land_param = LatLonLandParameter()
land_param.test_data_path = "/lcrc/group/e3sm/ac.zhang40/zppy_weekly_comprehensive_v3_output/test_zppy_pr651_20250115/v3.LR.historical_0051/post/lnd/180x360_aave/clim/2yr"

# Reference
land_param.reference_data_path = "/lcrc/group/e3sm/ac.zhang40/zppy_weekly_comprehensive_v3_output/test_zppy_pr651_20250115/v3.LR.historical_0051/post/lnd/180x360_aave/clim/2yr"
land_param.ref_name = "v3.LR.historical_0051"
land_param.short_ref_name = "same simulation"
# Optionally, swap test and reference model
if False:
    land_param.test_data_path, param.reference_data_path = (
        param.reference_data_path,
        param.test_data_path,
    )
    land_param.test_name, param.ref_name = param.ref_name, param.test_name
    land_param.short_test_name, param.short_ref_name = (
        param.short_ref_name,
        param.short_test_name,
    )
params.append(land_param)
# Run
runner.sets_to_run = ["lat_lon_land"]
runner.run_diags(params)

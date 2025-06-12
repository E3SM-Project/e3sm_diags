import os
import sys

from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.run import runner

# Once the run is done You can navigate to https://portal.nersc.gov/cfs/e3sm/<username> for results viewer.

param = CoreParameter()

param.reference_data_path = "/lcrc/group/e3sm/diagnostics/observations/Atm/climatology/"
param.test_data_path = "/lcrc/group/e3sm/ac.zhang40/example_v3/v3.LR.historical_0051/post/atm/180x360_aave/clim/30yr"
param.test_name = "v3.LR.historical_0051"
param.short_test_name = "v3.LR.historical_0051"
seasons = ["ANN", "DJF", "MAM", "JJA", "SON"]
param.multiprocessing = True
param.num_workers = 8

# Additional parameters:
# param.short_test_name = 'beta0.FC5COSP.ne30'
# param.run_type = 'model_vs_model'
# param.diff_title = 'Difference'
# param.output_format = ['png']
# param.output_format_subplot = ['pdf']
# param.save_netcdf = True

param.results_dir = '/lcrc/group/e3sm/public_html/cdat-migration-fy24/985-perf-degrade-compute-node-esmf-nompi'

# cfg_path = "auxiliary_tools/debug/985-perf-degrade/2-reproduce/core_set.cfg"
# sys.argv.extend(["--diags", cfg_path])

runner.sets_to_run = [
    "lat_lon",
    "zonal_mean_xy",
    "zonal_mean_2d",
    "polar",
    "cosp_histogram",
    "meridional_mean_2d",
]
runner.run_diags([param])

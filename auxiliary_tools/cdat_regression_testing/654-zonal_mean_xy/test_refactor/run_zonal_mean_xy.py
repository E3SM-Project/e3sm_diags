import os
import sys

from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.run import runner

param = CoreParameter()

# Location of the data.
param.reference_data_path = (
    "/global/cfs/cdirs/e3sm/e3sm_diags/obs_for_e3sm_diags/climatology/"
)
param.test_data_path = (
    "/global/cfs/cdirs/e3sm/e3sm_diags/postprocessed_e3sm_v2_data_for_e3sm_diags/20210528.v2rc3e.piControl.ne30pg2_EC30to60E2r2.chrysalis/climatology/rgr"
)
# Name of the test model data, used to find the climo files.
param.test_name = "20210528.v2rc3e.piControl.ne30pg2_EC30to60E2r2.chrysalis"
# An optional, shorter name to be used instead of the test_name.
param.short_test_name = "v2rc3e"

# What plotsets to run the diags on.
#param.sets = ["lat_lon"]
# Name of the folder where the results are stored.
# Change `prefix` to use your directory.
prefix = "/global/cfs/cdirs/e3sm/www/chengzhu/test_e3sm_refactor"
param.results_dir = os.path.join(prefix, "ex5_model_to_obs")

# Below are more optional arguments.

# Title of the difference plots.
param.diff_title = "Model - Obs."
# Save the netcdf files for each of the ref, test, and diff plot.
param.save_netcdf = True
# For running with multiprocessing.
# param.multiprocessing = True
# param.num_workers = 32
# Use the specified `.cfg` file for debugging
CFG_PATH = "examples/test_refactor/diags.cfg"
sys.argv.extend(["-d", CFG_PATH])

runner.sets_to_run = ["zonal_mean_xy"]
runner.run_diags([param])

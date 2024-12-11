"""
QA diffs

"""
# %%
import os
import sys

from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.run import runner

param = CoreParameter()

# Location of the data.
param.test_data_path = "/global/cfs/cdirs/e3sm/chengzhu/tutorial2024/v3.LR.historical_0101/post/atm/180x360_aave/clim/15yr"
param.reference_data_path = '/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/climatology/'


# Set this parameter to True.
# By default, e3sm_diags expects the test data to be climo data.

# Set this parameter to True.
# By default, e3sm_diags expects the ref data to be climo data.

# When running with time-series data, you don't need to specify the name of the data.
# But you should, otherwise nothing is displayed when the test/ref name is needed.

# This parameter modifies the software to accommodate model vs model runs.
# The default setting for run_type is 'model_vs_obs'.
param.run_type = "model_vs_obs"
# Name of the folder where the results are stored.
# Change `prefix` to use your directory.
prefix = "/global/cfs/cdirs/e3sm/www/chengzhu/e3sm_diags_tests/'"
param.results_dir = os.path.join(prefix, "12-11-v3-main")

# Below are more optional arguments.

# What plotsets to run the diags on.
# If not defined, then all available sets are used.
param.sets = ["lat_lon"]
# What seasons to run the diags on.
# If not defined, diags are run on ['ANN', 'DJF', 'MAM', 'JJA', 'SON'].
param.seasons = ["ANN"]

# For running with multiprocessing.
param.multiprocessing = False
# param.num_workers = 24

# %%
DIR_PATH = "/global/u2/c/chengzhu/e3sm_diags/auxiliary_tools/cdat_regression_testing/906-v3_complete_run"
CFG_PATH = os.path.join(DIR_PATH, "906-diags.cfg")
sys.argv.extend(["-d", CFG_PATH])
runner.run_diags([param])
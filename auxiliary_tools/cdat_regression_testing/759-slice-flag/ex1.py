# %%
import os
import sys

from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.run import runner

param = CoreParameter()

# Location of the data.
param.test_data_path = "/global/cfs/cdirs/e3sm/e3sm_diags/test_model_data_for_acme_diags/time-series/E3SM_v1"
param.reference_data_path = "/global/cfs/cdirs/e3sm/e3sm_diags/test_model_data_for_acme_diags/time-series/E3SM_v1"

# Variables
# param.variables = ["PRECT"]

# Set this parameter to True.
# By default, e3sm_diags expects the test data to be climo data.
param.test_timeseries_input = True
# Years to slice the test data, base this off the years in the filenames.
param.test_start_yr = "2011"
param.test_end_yr = "2013"

# Set this parameter to True.
# By default, e3sm_diags expects the ref data to be climo data.
param.ref_timeseries_input = True
# Years to slice the ref data, base this off the years in the filenames
param.ref_start_yr = "1850"
param.ref_end_yr = "1852"

# When running with time-series data, you don't need to specify the name of the data.
# But you should, otherwise nothing is displayed when the test/ref name is needed.
param.short_test_name = "historical_H1"
param.short_ref_name = "historical_H1"

# This parameter modifies the software to accommodate model vs model runs.
# The default setting for run_type is 'model_vs_obs'.
param.run_type = "model_vs_model"
# Name of the folder where the results are stored.
# Change `prefix` to use your directory.
prefix = "/global/cfs/cdirs/e3sm/www/cdat-migration-fy24"
param.results_dir = os.path.join(prefix, "759-slice-flag")

# Below are more optional arguments.

# What plotsets to run the diags on.
# If not defined, then all available sets are used.
param.sets = ["lat_lon"]
# What seasons to run the diags on.
# If not defined, diags are run on ['ANN', 'DJF', 'MAM', 'JJA', 'SON'].
param.seasons = ["ANN"]
# Title of the difference plots.
param.diff_title = "Model (2011-2013) - Model (1850-1852)"

# For running with multiprocessing.
param.multiprocessing = True
# param.num_workers = 24

# %%
cfg_path = "/global/u2/v/vo13/E3SM-Project/e3sm_diags/auxiliary_tools/cdat_regression_testing/759-slice-flag/671-diags.cfg"
sys.argv.extend(["--diags", cfg_path])
runner.run_diags([param])

# %%

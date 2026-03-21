import os

from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.run import runner

param = CoreParameter()
# Location of the data.
param.test_data_path = "/pscratch/sd/c/chengzhu/WCYCL20TR_ne30pg2_r05_IcoswISC30E3r5_JRA55_FOSIRL_1980050100.EN00/ts/rgr"
param.reference_data_path = "/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/climatology/"

# Set this parameter to True.
# By default, e3sm_diags expects the test data to be climo data.
param.test_timeseries_input = True
# Years to slice the test data, base this off the years in the filenames.
# With start_month=5, each cycle year runs May-April. Two cycle years:
# cycle year 1980: May 1980 - April 1981
# cycle year 1981: May 1981 - April 1982
param.test_start_yr = 1980
param.test_end_yr = 1981


# The starting month of the annual cycle for computing climatologies.
param.start_month = 5

# When running with time-series data, you don't need to specify the name of the data.
# But you should, otherwise nothing is displayed when the test/ref name is needed.
param.short_test_name = "FOSIRL_1980050100.EN00"

# Name of the folder where the results are stored.
# Change `prefix` to use your directory.
prefix = "/global/cfs/cdirs/e3sm/www/chengzhu/examples"
param.results_dir = os.path.join(prefix, "e3sm_diags_for_s2d_w_3dvar")

# What plotsets to run the diags on.
# If not defined, then all available sets are used.
param.sets = ["lat_lon"]
# What seasons to run the diags on.
# If not defined, diags are run on ['ANN', 'DJF', 'MAM', 'JJA', 'SON'].
#param.seasons = ["ANN"]
# For running with multiprocessing.
param.multiprocessing = True
param.num_workers = 32

runner.sets_to_run = ["lat_lon"]
runner.run_diags([param])

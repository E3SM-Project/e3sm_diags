import os

from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.run import runner

param = CoreParameter()

# Test data: time-series input from the original ncclimo run.
param.test_data_path = "/pscratch/sd/c/chengzhu/WCYCL20TR_ne30pg2_r05_IcoswISC30E3r5_JRA55_FOSIRL_1980050100.EN00/ts/rgr"
param.test_timeseries_input = True
param.test_start_yr = 1980
param.test_end_yr = 1981

# Reference data: climo from flexible-months ncclimo.
param.reference_data_path = "/pscratch/sd/c/chengzhu/WCYCL20TR_ne30pg2_r05_IcoswISC30E3r5_JRA55_FOSIRL_1980050100.EN00/ts/rgr_ncclimo_flex_month"

# The starting month of the annual cycle for computing climatologies.
param.start_month = 5

param.ref_name = "WCYCL20TR_ne30pg2_r05_IcoswISC30E3r5_JRA55_FOSIRL_1980050100.EN00"
param.short_test_name = "FOSIRL_1980050100.EN00 (ts)"
param.short_ref_name = "FOSIRL_1980050100.EN00 (ncclimo climo)"

# Name of the folder where the results are stored.
prefix = "/global/cfs/cdirs/e3sm/www/chengzhu/examples"
param.results_dir = os.path.join(prefix, "e3sm_diags_for_s2d_model_vs_model")

param.multiprocessing = True
param.num_workers = 32

runner.sets_to_run = ["lat_lon"]
runner.run_diags([param])

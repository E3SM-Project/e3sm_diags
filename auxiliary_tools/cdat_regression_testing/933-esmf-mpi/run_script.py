
import os
import sys

import numpy

from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.run import runner

import debugpy
debugpy.listen(("0.0.0.0", 5678))
print("Waiting for debugger attach...")
debugpy.wait_for_client()

short_name = 'v2.LR.historical_0201'
test_ts = '/lcrc/group/e3sm/ac.forsyth2/zppy_weekly_comprehensive_v2_output/test_pr651_both_commits_20250117/v2.LR.historical_0201/post/atm/180x360_aave/ts/monthly/2yr'

param = CoreParameter()

# Model
param.test_data_path = '/lcrc/group/e3sm/ac.forsyth2/zppy_weekly_comprehensive_v2_output/test_pr651_both_commits_20250117/v2.LR.historical_0201/post/atm/180x360_aave/clim/2yr'
param.test_name = 'v2.LR.historical_0201'
param.short_test_name = short_name

# Ref

# Obs
param.reference_data_path = '/lcrc/group/e3sm/diagnostics/observations/Atm/climatology/'


# Output dir
param.results_dir = 'model_vs_obs_1982-1983'

# Additional settings
param.run_type = 'model_vs_obs'
param.diff_title = 'Model - Observations'
param.multiprocessing = True
param.num_workers = 8
#param.fail_on_incomplete = True
params = [param]

# Run
cfg_path = "auxiliary_tools/cdat_regression_testing/933-esmf-mpi/v2_run.cfg"
sys.argv.extend(["--diags", cfg_path])

runner.sets_to_run = ['lat_lon', ]
runner.run_diags(params)

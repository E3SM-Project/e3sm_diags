import os
from acme_diags.parameter.core_parameter import CoreParameter
from acme_diags.parameter.enso_diags_parameter import EnsoDiagsParameter
from acme_diags.run import runner

# Generic parameters. See enso.cfg for specific parameters.
# Run with `python run_enso_diags.py -d enso.cfg`.
param = CoreParameter()
machine_path_prefix = '/global/project/projectdirs/acme/acme_diags'
# /global/project/projectdirs/acme/acme_diags/obs_for_e3sm_diags/time-series contains:
# GPCP_v2.2/PRECT_197901_201412.nc
# HadISST/sst_187001_201712.nc
param.reference_data_path = os.path.join(machine_path_prefix, 'obs_for_e3sm_diags/time-series')
# /global/project/projectdirs/acme/acme_diags/test_model_data_for_acme_diags/time-series/E3SM_v1 contains:
# PRECC_185001_201312.nc
# PRECL_185001_201312.nc
# TS_185001_201312.nc
# OCNFRAC_185001_201312.nc
# The first two are used to derive PRECT
# The last two are used to derive sst
param.test_data_path = os.path.join(machine_path_prefix, 'test_model_data_for_acme_diags/time-series/E3SM_v1')
# html_prefix = '/global/project/projectdirs/acme/www/<username>'
html_prefix = '/global/project/projectdirs/acme/www/forsyth'
param.results_dir = os.path.join(html_prefix, 'enso_diags_1_8')
param.test_name = 'e3sm_v1'

# We're passing in this new object as well, in
# addtion to the CoreParameter object.
ts_param = EnsoDiagsParameter()
ts_param.test_name = 'e3sm_v1'
ts_param.start_yr = '2000'
ts_param.end_yr = '2004'
#ts_param.end_yr = '2010'

runner.sets_to_run = ['enso_diags']
# Call acme_diags.run.runner, which calls
# creates an acme_diags.run.Run() object and runs
# acme_diags.run.Run.run_diags, which calls
# acme_diags.acme_diags_driver.main, which
# 1) uses cdp.cdp_run.{multiprocess, distribute, serial} to call
# acme_diags.acme_diags.run_diag, which calls
# acme_diags.driver.{set_name}_driver.run_diag
# e.g. acme_diags.driver.enso_diags_driver.run_diag is called
# 2) calls acme_diags.viewer.main.create_viewer, which calls
# acme_diags.viewer.main.SET_TO_VIEWER[set_name], which is
# acme_diags.viewer.{set_name}_viewer.create_viewer
# e.g. acme_diags.viewer.enso_diags_viewer.create_viewer is called
runner.run_diags([param, ts_param])

"""Visual test for the new enso_diags index_timeseries plot_type.

Run with the companion cfg (limits the run to a single index_timeseries block):

    python auxiliary_tools/debug/enso_index_timeseries/run_index_timeseries.py \
        -d auxiliary_tools/debug/enso_index_timeseries/run_index_timeseries.cfg

Uses a real E3SM v2 piControl TS time series (years 0051-0060) with the built-in
HadISST nino index as the reference.
"""

from e3sm_diags.parameter.enso_diags_parameter import EnsoDiagsParameter
from e3sm_diags.run import runner

base = "/global/cfs/cdirs/e3sm/e3sm_diags/postprocessed_e3sm_v2_data_for_e3sm_diags"
test_ts = f"{base}/20210528.v2rc3e.piControl.ne30pg2_EC30to60E2r2.chrysalis/time-series/rgr"
obs_ts = "/global/cfs/cdirs/e3sm/e3sm_diags/obs_for_e3sm_diags/time-series"

param = EnsoDiagsParameter()
param.reference_data_path = obs_ts
param.test_data_path = test_ts
param.test_name = "v2rc3e.piControl"
param.test_start_yr = "0051"
param.test_end_yr = "0060"
# Obs nino index built-in record spans 1870-2018.
param.ref_start_yr = "2001"
param.ref_end_yr = "2010"
param.save_netcdf = True
param.results_dir = "/global/cfs/cdirs/e3sm/www/chengzhu/tests/enso_index_timeseries_test"

runner.sets_to_run = ["enso_diags"]
runner.run_diags([param])

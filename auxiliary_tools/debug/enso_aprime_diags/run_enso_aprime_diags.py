"""Visual test for the ported a-prime enso_diags plot types.

Covers all five ported plot types -- nino_index_timeseries, seasonality,
interannual_variability, equatorial_soi, and lead_lag -- in a single run. Run
with the companion cfg (the ``-d`` cfg limits the run to the ported blocks):

    python auxiliary_tools/debug/enso_aprime_diags/run_enso_aprime_diags.py \
        -d auxiliary_tools/debug/enso_aprime_diags/run_enso_aprime_diags.cfg

Uses a real E3SM v2 piControl time series (years 0051-0060). Observations come
from the standard e3sm_diags obs directory
(/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/time-series). References:

* nino_index_timeseries / seasonality: built-in HadISST nino index.
* interannual_variability: gridded HadISST SST observations.
* equatorial_soi: observed PSL (ERA5 time series) for the EQSOI index and the
  built-in HadISST record for the Nino3.4 index.
* lead_lag: observed fields for each variable (ERA5, and GPCP_v3.2 for PRECT)
  and the built-in HadISST record for the Nino3.4 index.

The lead_lag set is the slowest -- a global regression at every grid point for
five lags and both cases, over the same variables as the regression map
(TS, PRECT, TAUX, TAUY, LHFLX, SHFLX, NET_FLUX_SRF); a full run takes a while.
"""

from e3sm_diags.parameter.enso_diags_parameter import EnsoDiagsParameter
from e3sm_diags.run import runner

base = "/global/cfs/cdirs/e3sm/e3sm_diags/postprocessed_e3sm_v2_data_for_e3sm_diags"
test_ts = f"{base}/20210528.v2rc3e.piControl.ne30pg2_EC30to60E2r2.chrysalis/time-series/rgr"
obs_ts = "/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/time-series"

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
param.results_dir = "/global/cfs/cdirs/e3sm/www/chengzhu/tests/enso_aprime_diags_test"

runner.sets_to_run = ["enso_diags"]
runner.run_diags([param])

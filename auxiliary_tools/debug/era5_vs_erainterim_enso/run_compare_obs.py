"""Run the enso_diags response (map) and feedback (scatter) plots twice -- once
with the original observations (ERA-Interim / GPCP_v2.3) and once with the newer
ones (ERA5 / GPCP_v3.2) -- so the two references can be compared in the viewer.

Anchored on the v3.LR.amip_0101 run over its real years (1995-2004), so the test
column is identical between the paired blocks and only the reference changes.

    PYTHONPATH=/global/u2/c/chengzhu/e3sm_diags:$PYTHONPATH \
      python auxiliary_tools/debug/era5_vs_erainterim_enso/run_compare_obs.py \
        -d auxiliary_tools/debug/era5_vs_erainterim_enso/compare_obs.cfg
"""

from e3sm_diags.parameter.enso_diags_parameter import EnsoDiagsParameter
from e3sm_diags.run import runner

test_ts = (
    "/global/cfs/cdirs/e3sm/chengzhu/tests/zppy_example_v3/"
    "v3.LR.amip_0101/post/atm/180x360_aave/ts/monthly/10yr"
)
obs_ts = "/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/time-series"

param = EnsoDiagsParameter()
param.reference_data_path = obs_ts
param.test_data_path = test_ts
param.test_name = "v3.LR.amip_0101"
param.test_start_yr = "1995"
param.test_end_yr = "2004"
param.ref_start_yr = "1995"
param.ref_end_yr = "2004"
param.save_netcdf = True
param.results_dir = "/global/cfs/cdirs/e3sm/www/chengzhu/tests/enso_obs_compare"

runner.sets_to_run = ["enso_diags"]
runner.run_diags([param])

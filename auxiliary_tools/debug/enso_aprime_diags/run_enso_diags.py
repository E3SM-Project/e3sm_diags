"""Run the enso_diags set against the v3.LR.amip_0101 simulation (1995-2004).

Run without ``-d`` for the full default enso_diags set; pass
``-d run_enso_aprime_diags.cfg`` to limit the run to the ported plot types.
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
param.results_dir = "/global/cfs/cdirs/e3sm/www/chengzhu/tests/enso_aprime_diags_test"

# Run cfg blocks in parallel (one process per block). Use a compute node, not
# the login node (30 GB/user cgroup cap); 8 workers peak ~20-28 GB.
param.multiprocessing = True
param.num_workers = 8

runner.sets_to_run = ["enso_diags"]
runner.run_diags([param])

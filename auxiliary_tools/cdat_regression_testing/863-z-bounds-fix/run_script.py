import os
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.run import runner

param = CoreParameter()

# param.reference_data_path = '/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/climatology'
# param.test_data_path = '/global/cfs/cdirs/e3sm/zhang40/e3sm_diags_for_EAMxx/data/Cess'
param.reference_data_path = (
    "/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/climatology"
)
param.test_data_path = "/global/cfs/cdirs/e3sm/chengzhu/eamxx/post/data/rgr"
param.test_name = "eamxx_decadal"
param.seasons = ["ANN"]
# param.save_netcdf = True

prefix = "/global/cfs/cdirs/e3sm/www/cdat-migration-fy24"
param.results_dir = os.path.join(prefix, "889-z-axis-bnds")

runner.sets_to_run = ["zonal_mean_2d"]

runner.run_diags([param])

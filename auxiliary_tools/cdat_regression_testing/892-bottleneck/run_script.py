import os
import sys

from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.run import runner

param = CoreParameter()


param.reference_data_path = (
    "/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/time-series"
)
param.test_data_path = "/global/cfs/cdirs/e3sm/chengzhu/eamxx/post/data/rgr"
param.test_name = "eamxx_decadal"
param.seasons = ["ANN"]
# param.save_netcdf = True

param.ref_timeseries_input = True
# Years to slice the ref data, base this off the years in the filenames.
param.ref_start_yr = "1996"
param.ref_end_yr = "1996"

prefix = "/global/cfs/cdirs/e3sm/www/cdat-migration-fy24/892-bottleneck"
param.results_dir = os.path.join(prefix, "eamxx_decadal_1996_1107_edv3")

cfg_path = "auxiliary_tools/cdat_regression_testing/892-bottleneck/run_script.cfg"
sys.argv.extend(["--diags", cfg_path])

runner.sets_to_run = [
    "lat_lon",
    "zonal_mean_xy",
    "zonal_mean_2d",
    "zonal_mean_2d_stratosphere",
    "polar",
    "cosp_histogram",
    "meridional_mean_2d",
    "annual_cycle_zonal_mean",
]

runner.run_diags([param])

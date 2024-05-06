import os

from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.run import runner

data_prefix = "/global/cfs/cdirs/e3sm/e3sm_diags"
html_prefix = "/global/cfs/cdirs/e3sm/www/<username>"  # Change <username>
html_prefix = "/global/cfs/cdirs/e3sm/www/chengzhu"  # Change <username>
# Once the run is done You can navigate to https://portal.nersc.gov/cfs/e3sm/<username> for results viewer.

param = CoreParameter()

param.reference_data_path = "/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/climatology/"
param.test_data_path = "/global/cfs/cdirs/e3sm/chengzhu/e3sm_diags_zppy_test_complete_run_output/v3.LR.historical_0101_20240423/post/atm/180x360_aave/clim/15yr"
param.test_name = "extendedOutput.v3.LR.historical_0101"
param.short_test_name = "v3.LR.historical_0101"
param.seasons = [
    "JJA",
]  # Default setting: seasons = ["ANN", "DJF", "MAM", "JJA", "SON"]
param.multiprocessing = True
param.num_workers = 24

# Additional parameters:
# param.short_test_name = 'beta0.FC5COSP.ne30'
# param.run_type = 'model_vs_model'
# param.diff_title = 'Difference'
# param.output_format = ['png']
# param.output_format_subplot = ['pdf']
# param.save_netcdf = True

param.results_dir = os.path.join(html_prefix, "tutorial_2024/e3sm_diags_core_sets")

runner.sets_to_run = [
    "lat_lon",
    "zonal_mean_xy",
    "zonal_mean_2d",
    "polar",
    "cosp_histogram",
    "meridional_mean_2d",
]
runner.run_diags([param])

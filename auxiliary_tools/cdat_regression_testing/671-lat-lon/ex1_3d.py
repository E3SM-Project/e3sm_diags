# %%
import os

from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.run import runner

param = CoreParameter()

# %%
param.sets = ["lat_lon"]
param.case_id = "ERA-Interim"
param.variables = ["T"]
param.seasons = ["ANN"]
param.plevs = [850.0]
param.contour_levels = [240, 245, 250, 255, 260, 265, 270, 275, 280, 285, 290, 295]
param.diff_levels = [-10, -7.5, -5, -4, -3, -2, -1, -0.5, 0.5, 1, 2, 3, 4, 5, 7.5, 10]

param.test_name = "system tests"
param.short_test_name = "short_system tests"
param.ref_name = "ERA-Interim"
param.reference_name = "ERA-Interim Reanalysis 1979-2015"
param.reference_data_path = (
    "/global/u2/v/vo13/E3SM-Project/e3sm_diags/tests/integration/integration_test_data"
)
param.ref_file = "ta_ERA-Interim_ANN_198001_201401_climo.nc"
param.test_data_path = (
    "/global/u2/v/vo13/E3SM-Project/e3sm_diags/tests/integration/integration_test_data"
)
param.test_file = "T_20161118.beta0.FC5COSP.ne30_ne30.edison_ANN_climo.nc"

param.backend = "mpl"
prefix = "/global/cfs/cdirs/e3sm/www/vo13/examples"
param.results_dir = os.path.join(prefix, "lat_lon_3d_var_test")
param.debug = True
param.multiprocessing = False


# %%
runner.run_diags([param])

# %%

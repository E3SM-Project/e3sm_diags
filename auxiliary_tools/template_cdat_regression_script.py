# %%
import os

from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.run import runner

param = CoreParameter()

# %%
param.sets = ["cosp_histogram"]
param.case_id = "MISR-COSP"
param.variables = ["COSP_HISTOGRAM_MISR"]
param.seasons = ["ANN"]
param.contour_levels = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5]
param.diff_levels = [
    -3.0,
    -2.5,
    -2.0,
    -1.5,
    -1.0,
    -0.5,
    0,
    0.5,
    1.0,
    1.5,
    2.0,
    2.5,
    3.0,
]

param.test_name = "system tests"
param.short_test_name = "short_system tests"
param.ref_name = "MISRCOSP"
param.reference_name = "MISR COSP (2000-2009)"
param.reference_data_path = "tests/integration"
param.ref_file = "CLDMISR_ERA-Interim_ANN_198001_201401_climo.nc"
param.test_data_path = "tests/integration"
param.test_file = "CLD_MISR_20161118.beta0.FC5COSP.ne30_ne30.edison_ANN_climo.nc"

param.backend = "mpl"
prefix = "/global/cfs/cdirs/e3sm/www/vo13/examples"
param.results_dir = os.path.join(prefix, "cdat_regression_tests/", param.sets[0])
param.multiprocessing = False

# %%
runner.run_diags([param])

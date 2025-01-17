# %%
import os

import numpy

from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.parameter.enso_diags_parameter import EnsoDiagsParameter
from e3sm_diags.parameter.streamflow_parameter import StreamflowParameter
from e3sm_diags.parameter.tc_analysis_parameter import TCAnalysisParameter
from e3sm_diags.run import runner

short_name = "v3.LR.historical_0051"
# test_ts = "ts"
test_ts = "/lcrc/group/e3sm/ac.forsyth2/zppy_min_case_e3sm_diags_cdat_migrated_output/test-diags-no-cdat-20240917/v3.LR.historical_0051/post/atm/180x360_aave/ts/monthly/2yr"
start_yr = int("1987")
end_yr = int("1988")
num_years = end_yr - start_yr + 1
ref_start_yr = 1985

param = CoreParameter()

# Model
param.test_data_path = "climo"
param.test_name = "v3.LR.historical_0051"
param.short_test_name = short_name

# Ref

# Obs
param.reference_data_path = "/lcrc/group/e3sm/diagnostics/observations/Atm/climatology/"

# Output dir
param.results_dir = "model_vs_obs_1987-1988"

# Additional settings
param.run_type = "model_vs_obs"
param.diff_title = "Model - Observations"
param.output_format = ["png"]
param.output_format_subplot = []
param.multiprocessing = True
param.num_workers = 8
# param.fail_on_incomplete = True
params = [param]

# Model land
enso_param = EnsoDiagsParameter()
enso_param.test_data_path = test_ts
enso_param.test_name = short_name
enso_param.test_start_yr = start_yr
enso_param.test_end_yr = end_yr

# Obs
enso_param.reference_data_path = (
    "/lcrc/group/e3sm/diagnostics/observations/Atm/time-series/"
)
enso_param.ref_start_yr = ref_start_yr
enso_param.ref_end_yr = ref_start_yr + 10

params.append(enso_param)
streamflow_param = StreamflowParameter()
streamflow_param.reference_data_path = (
    "/lcrc/group/e3sm/diagnostics/observations/Atm/time-series/"
)
streamflow_param.test_data_path = "/lcrc/group/e3sm/ac.forsyth2/zppy_min_case_e3sm_diags_cdat_migrated_output/test-diags-no-cdat-20240917/v3.LR.historical_0051/post/rof/native/ts/monthly/2yr"
streamflow_param.test_name = short_name
streamflow_param.test_start_yr = start_yr
streamflow_param.test_end_yr = end_yr

# Obs
streamflow_param.reference_data_path = (
    "/lcrc/group/e3sm/diagnostics/observations/Atm/time-series/"
)
streamflow_param.ref_start_yr = (
    "1986"  # Streamflow gauge station data range from year 1986 to 1995
)
streamflow_param.ref_end_yr = "1995"

params.append(streamflow_param)
tc_param = TCAnalysisParameter()
tc_param.test_data_path = "/lcrc/group/e3sm/ac.forsyth2/zppy_min_case_e3sm_diags_cdat_migrated_output/test-diags-no-cdat-20240917/v3.LR.historical_0051/post/atm/tc-analysis_1987_1988"
tc_param.short_test_name = short_name
tc_param.test_start_yr = "1987"
tc_param.test_end_yr = "1988"

# Obs
tc_param.reference_data_path = (
    "/lcrc/group/e3sm/diagnostics/observations/Atm/tc-analysis/"
)
# For model vs obs, the ref start and end year can be any four digit strings
# For now, use all available years from obs by default
tc_param.ref_start_yr = "1979"
tc_param.ref_end_yr = "2018"

params.append(tc_param)

# Run
runner.sets_to_run = ["lat_lon", "tc_analysis", "enso_diags", "streamflow"]
# runner.sets_to_run = ["tc_analysis", "enso_diags", "streamflow"]
# runner.sets_to_run = ["enso_diags", "streamflow"]
# runner.sets_to_run = ["enso_diags"]
runner.run_diags(params)

# %%

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
streamflow_param.test_data_path = "/lcrc/group/e3sm/ac.forsyth2/zppy_min_case_e3sm_diags_cdat_migrated_output/test-diags-no-cdat-20240917/v3.LR.historical_0051/post/rof"
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
# runner.sets_to_run = ["lat_lon", "tc_analysis", "enso_diags", "streamflow"]
# runner.sets_to_run = ["tc_analysis", "enso_diags", "streamflow"]
runner.sets_to_run = ["tc_analysis"]
runner.run_diags(params)

"""
Generating TC Metrics from TE Stitch Files
2024-10-09 16:20:01,321 [INFO]: tc_analysis_driver.py(generate_tc_metrics_from_te_stitch_file:174) >> ============================================
2024-10-09 16:20:01,330 [INFO]: tc_analysis_driver.py(_calc_num_storms_and_max_len:235) >> Number of storms: 0
2024-10-09 16:20:01,331 [INFO]: tc_analysis_driver.py(_calc_num_storms_and_max_len:236) >> Max length of storms: 0
2024-10-09 16:20:01,332 [ERROR]: core_parameter.py(_run_diag:343) >> Error in e3sm_diags.driver.tc_analysis_driver
Traceback (most recent call last):
  File "/gpfs/fs1/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/parameter/core_parameter.py", line 340, in _run_diag
    single_result = module.run_diag(self)
  File "/gpfs/fs1/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/driver/tc_analysis_driver.py", line 91, in run_diag
    test_data["metrics"] = generate_tc_metrics_from_te_stitch_file(test_te_file)
  File "/gpfs/fs1/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/driver/tc_analysis_driver.py", line 181, in generate_tc_metrics_from_te_stitch_file
    te_stitch_vars = _get_vars_from_te_stitch(lines, max_len, num_storms)
  File "/gpfs/fs1/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/driver/tc_analysis_driver.py", line 258, in _get_vars_from_te_stitch
    year_start = int(lines[0].split("\t")[2])
IndexError: list index out of range
2024-10-09 16:20:01,360 [WARNING]: e3sm_diags_driver.py(main:378) >> There was not a single valid diagnostics run, no viewer created.
2024-10-09 16:20:01,361 [ERROR]: run.py(run_diags:91) >> Error traceback:
Traceback (most recent call last):
  File "/gpfs/fs1/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/run.py", line 89, in run_diags
    params_results = main(params)
  File "/gpfs/fs1/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/e3sm_diags_driver.py", line 397, in main
    if parameters_results[0].fail_on_incomplete and (
IndexError: list index out of range
2024-10-09 16:20:01,368 [INFO]: logger.py(move_log_to_prov_dir:106) >> Log file saved in model_vs_obs_1987-1988/prov/e3sm_diags_run.log
"""

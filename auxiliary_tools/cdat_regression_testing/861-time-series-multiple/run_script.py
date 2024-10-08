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
streamflow_param.test_data_path = "rof"
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
runner.sets_to_run = ["enso_diags"]
runner.run_diags(params)

# https://github.com/E3SM-Project/zppy/pull/598#issuecomment-2384025994
"""
OSError: No files found for target variable TAUX or derived variables ([('tauu',), ('surf_mom_flux_\
U',)]) in ts.
OSError: No files found for target variable LHFLX or derived variables ([('hfls',), ('QFLX',), ('su\
rface_upward_latent_heat_flux',)]) in ts.
OSError: No files found for target variable SHFLX or derived variables ([('hfss',), ('surf_sens_flu\
x',)]) in ts
OSError: No files found for target variable LHFLX or derived variables ([('hfls',), ('QFLX',), ('su\
rface_upward_latent_heat_flux',)]) in ts.
OSError: No files found for target variable SHFLX or derived variables ([('hfss',), ('surf_sens_flu\
x',)]) in ts.
OSError: No files found for target variable NET_FLUX_SRF or derived variables ([('FSNS', 'FLNS', 'Q\
FLX', 'PRECC', 'PRECL', 'PRECSC', 'PRECSL', 'SHFLX'), ('FSNS', 'FLNS', 'LHFLX', 'SHFLX'), ('FSNS', \
'FLNS', 'QFLX', 'SHFLX'), ('rsds', 'rsus', 'rlds', 'rlus', 'hfls', 'hfss')]) in ts.
OSError: No files found for target variable PRECT or derived variables ([('PRECT',), ('pr',), ('PRE\
CC', 'PRECL'), ('sat_gauge_precip',), ('PrecipLiqSurfMassFlux', 'PrecipIceSurfMassFlux')]) in ts.
OSError: No files found for target variable TAUX or derived variables ([('tauu',), ('surf_mom_flux_\
U',)]) in ts.
OSError: No files found for target variable TAUY or derived variables ([('tauv',), ('surf_mom_flux_\
V',)]) in ts.
OSError: No files found for target variable FSNS or derived variables ([('sfc_net_sw_all_mon',), ('\
rsds', 'rsus')]) in ts.
OSError: No files found for target variable NET_FLUX_SRF or derived variables ([('FSNS', 'FLNS', 'Q\
FLX', 'PRECC', 'PRECL', 'PRECSC', 'PRECSL', 'SHFLX'), ('FSNS', 'FLNS', 'LHFLX', 'SHFLX'), ('FSNS', \
'FLNS', 'QFLX', 'SHFLX'), ('rsds', 'rsus', 'rlds', 'rlus', 'hfls', 'hfss')]) in ts.
"""


# %%
# import xcdat as xc

# filepath = "/lcrc/group/e3sm/ac.forsyth2/zppy_min_case_e3sm_diags_cdat_migrated_output/test-diags-no-cdat-20240917/v3.LR.historical_0051/post/atm/180x360_aave/ts/monthly/2yr/TS*.nc"

# ds = xc.open_mfdataset(
#     filepath, add_bounds=["X", "Y", "T"], decode_times=True, use_cftime=True
# )

# # %%

# 2024-10-08 14:54:33,975 [ERROR]: core_parameter.py(_run_diag:343) >> Error in e3sm_diags.driver.enso_diags_driver
# Traceback (most recent call last):
#   File "/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/parameter/core_parameter.py", line 340, in _run_diag
#     single_result = module.run_diag(self)
#   File "/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/driver/enso_diags_driver.py", line 57, in run_diag
#     return run_diag_map(parameter)
#   File "/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/driver/enso_diags_driver.py", line 101, in run_diag_map
#     ds_test_reg_coe, da_test_conf_lvls = calc_linear_regression(
#   File "/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/driver/enso_diags_driver.py", line 449, in calc_linear_regression
#     reg_coe = xs.linslope(independent_var, anomaly_var, keep_attrs=True)
#   File "/gpfs/fs1/home/ac.tvo/mambaforge/envs/e3sm_diags_dev_673/lib/python3.10/site-packages/xskillscore/core/deterministic.py", line 143, in linslope
#     return xr.apply_ufunc(
#   File "/gpfs/fs1/home/ac.tvo/mambaforge/envs/e3sm_diags_dev_673/lib/python3.10/site-packages/xarray/core/computation.py", line 1270, in apply_ufunc
#     return apply_dataarray_vfunc(
#   File "/gpfs/fs1/home/ac.tvo/mambaforge/envs/e3sm_diags_dev_673/lib/python3.10/site-packages/xarray/core/computation.py", line 316, in apply_dataarray_vfunc
#     result_var = func(*data_vars)
#   File "/gpfs/fs1/home/ac.tvo/mambaforge/envs/e3sm_diags_dev_673/lib/python3.10/site-packages/xarray/core/computation.py", line 771, in apply_variable_ufunc
#     raise ValueError(
# ValueError: dimension time on 0th function argument to apply_ufunc with dask='parallelized' consists of multiple chunks, but is also a core dimension. To fix, either rechunk into a single array chunk along this dimension, i.e., ``.chunk(dict(time=-1))``, or pass ``allow_rechunk=True`` in ``dask_gufunc_kwargs`` but beware that this may significantly increase memory usage.

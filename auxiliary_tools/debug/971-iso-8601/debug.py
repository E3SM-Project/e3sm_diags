"""
This script is designed to debug a specific issue in the E3SM Diagnostics tool
(e3sm_diags) related to the tropical_subseasonal set. The bug occurs when
processing year values less than 1000, resulting in the error:
"ValueError: no ISO-8601 or cftime-string-like match for string: 1-01-01".
This issue was reported in PR #971.

The script replicates the behavior of the following command-line invocation
of e3sm_diags:

e3sm_diags tropical_subseasonal --no_viewer --reference_data_path '/lcrc/soft/climate/e3sm_diags_data/obs_for_e3sm_diags/time-series' --test_data_path '/lcrc/group/e3sm2/ac.xzheng/E3SMv3_dev/20250404.wcycl1850.ne120pg2_r025_RRSwISC6to18E3r5.test4.chrysalis//post/atm/180x360_traave/ts/daily/1yr' --results_dir '/lcrc/group/e3sm/public_html/diagnostic_output/ac.zhang40/tests/tropical_subseasonal_time_fix' --case_id 'wavenumber-frequency' --ref_timeseries_input --test_timeseries_input --run_type 'model_vs_obs' --sets 'tropical_subseasonal' --variables 'PRECT' --seasons 'ANN' 'DJF' 'MAM' 'JJA' 'SON' --regions '15S15N' --regrid_tool 'xesmf' --regrid_method 'conservative_normed' --multiprocessing --num_workers '32' --backend 'cartopy' --output_format 'png' --output_format_subplot 'pdf' --canvas_size_w '1212' --canvas_size_h '1628' --figsize '8.5' '11.0' --dpi '150' --arrows --test_name 'v3.HR_test4' --short_test_name 'v3.HR_test4' --test_colormap 'cet_rainbow.rgb' --ref_name 'IMERG_Daily' --reference_name 'IMERG Daily' --reference_colormap 'cet_rainbow.rgb' --diff_title 'percent difference' --diff_colormap 'diverging_bwr.rgb' --granulate 'variables' 'plevs' 'regions' --selectors 'sets' 'seasons' --test_start_yr 2 --test_end_yr 18 --ref_start_yr 2001 --ref_end_yr 2010

The script uses the e3sm_diags Python API to configure and run the diagnostics
with the same parameters as the command-line invocation. It is intended to
help identify and resolve the issue with year values less than 1000.
"""

from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.run import runner
from e3sm_diags.parameter.tropical_subseasonal_parameter import TropicalSubseasonalParameter

# Set up parameters
param = CoreParameter()
param.no_viewer = True
param.reference_data_path = '/lcrc/soft/climate/e3sm_diags_data/obs_for_e3sm_diags/time-series'
param.test_data_path = '/lcrc/group/e3sm2/ac.xzheng/E3SMv3_dev/20250404.wcycl1850.ne120pg2_r025_RRSwISC6to18E3r5.test4.chrysalis//post/atm/180x360_traave/ts/daily/1yr'
param.results_dir = '/lcrc/group/e3sm/public_html/diagnostic_output/ac.tvo/tests/tropical_subseasonal_time_fix'
param.case_id = 'wavenumber-frequency'
param.ref_timeseries_input = True
param.test_timeseries_input = True
param.run_type = 'model_vs_obs'
param.sets = ['tropical_subseasonal']
param.variables = ['PRECT']
param.seasons = ['ANN', 'DJF', 'MAM', 'JJA', 'SON']
param.regions = ['15S15N']
param.regrid_tool = 'xesmf'
param.regrid_method = 'conservative_normed'
param.multiprocessing = True
param.num_workers = 32
param.output_format_subplot = ['pdf']
param.canvas_size_w = 1212
param.canvas_size_h = 1628
param.figsize = [8.5, 11.0]
param.dpi = 150
param.arrows = True
param.short_test_name = 'v3.HR_test4'
param.test_colormap = 'cet_rainbow.rgb'
param.ref_name = 'IMERG_Daily'
param.reference_name = 'IMERG Daily'
param.reference_colormap = 'cet_rainbow.rgb'
param.diff_title = 'percent difference'
param.diff_colormap = 'diverging_bwr.rgb'
param.granulate = ['variables', 'plevs', 'regions']
param.selectors = ['sets', 'seasons']


trop_param = TropicalSubseasonalParameter()
trop_param.test_start_yr = 2
trop_param.test_name = 'v3.HR_test4'

trop_param.test_end_yr = 18
trop_param.ref_start_yr = 2001
trop_param.ref_end_yr = 2010

# Run the diagnostics
runner.sets_to_run = ['tropical_subseasonal']
runner.run_diags([param, trop_param])
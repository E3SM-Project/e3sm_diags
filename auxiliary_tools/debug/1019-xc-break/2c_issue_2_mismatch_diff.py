"""
- 2c - Test Set 3
- Issue 2 — Diff Plots Appear Inverted for lat_lon_land

Debug mismatching third diff plots between expected and actual output
specifically for lat_lon_land.

Source: https://web.lcrc.anl.gov/public/e3sm/diagnostic_output/ac.forsyth2/zppy_weekly_comprehensive_v3_www/test-zppy-diags-1019-xc-break/v3.LR.historical_0051/e3sm_diags/lnd_monthly_mvm_lnd/model_vs_model_1987-1988/prov/e3sm.py
Log file: https://web.lcrc.anl.gov/public/e3sm/diagnostic_output/ac.forsyth2/zppy_weekly_comprehensive_v3_www/test-zppy-diags-1019-xc-break/v3.LR.historical_0051/e3sm_diags/lnd_monthly_mvm_lnd/model_vs_model_1987-1988/prov/e3sm_diags_run.log
zppy prov: https://web.lcrc.anl.gov/public/e3sm/diagnostic_output/ac.forsyth2/zppy_weekly_comprehensive_v3_www/test-zppy-diags-1019-xc-break/v3.LR.historical_0051/provenance.20251205_185051_134335.cfg
Viewer: https://web.lcrc.anl.gov/public/e3sm/diagnostic_output/ac.forsyth2/zppy_weekly_comprehensive_v3_www/test-zppy-diags-1019-xc-break/v3.LR.historical_0051/e3sm_diags/lnd_monthly_mvm_lnd/model_vs_model_1987-1988/viewer/index.html
"""
# source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_chrysalis.sh
# Reproduce e3sm_diags dependencies in E3SM Unified 1.12.0
# conda create -n e3sm_diags_unified_1120 -c conda-forge e3sm-diags=3.12.0 xcdat=0.10.1 xarray=2025.09.0 matplotlib=3.10.6 xesmf=0.8.8 esmpy=8.9.0 esmf=8.9.0 cartopy=0.24.0 cartopy_offlinedata=0.24.0
#  python auxiliary_tools/debug/1019-xc-break/2c_issue_2_mismatch_diff.py

import os
import sys
import numpy
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.parameter.lat_lon_land_parameter import LatLonLandParameter


from e3sm_diags.run import runner

short_name = 'v3.LR.historical_0051'
test_ts = "/lcrc/group/e3sm/ac.forsyth2/zppy_weekly_comprehensive_v3_output/test-zppy-diags-1019-xc-break/v3.LR.historical_0051/post/lnd/180x360_aave/ts/monthly/2yr"
start_yr = int('1987')
end_yr = int('1988')
num_years = end_yr - start_yr + 1

param = CoreParameter()

# Model
param.test_name = 'v3.LR.historical_0051'
param.short_test_name = short_name

# Ref

# Output dir
param.results_dir = '/lcrc/group/e3sm/public_html/diagnostic_output/ac.tvo/tests/annual_cycle_xr_2025.11/2c_issue_2_mismatch_diff_cartopy_025_mvce/model_vs_model_1987-1988'

# Additional settings
param.run_type = 'model_vs_model'
param.diff_title = 'Difference'
param.output_format = ['png']
param.output_format_subplot = []
param.multiprocessing = True
param.num_workers = 8
#param.fail_on_incomplete = True
params = [param]

# Model land
land_param = LatLonLandParameter()
land_param.test_data_path = "climo_test_land"


# Reference
land_param.reference_data_path = "/lcrc/group/e3sm/ac.forsyth2/zppy_weekly_comprehensive_v3_output/zppy-diags-1019-xc-break-test4/v3.LR.historical_0051/post/lnd/180x360_aave/clim/2yr"
land_param.ref_name = 'v3.LR.historical_0051'
land_param.short_ref_name = 'same simulation'

# Optionally, swap test and reference model
if False:
   land_param.test_data_path, param.reference_data_path = param.reference_data_path, param.test_data_path
   land_param.test_name, param.ref_name = param.ref_name, param.test_name
   land_param.short_test_name, param.short_ref_name = param.short_ref_name, param.short_test_name

params.append(land_param)

cfg_path = "auxiliary_tools/debug/1019-xc-break/2c_issue_2_mismatch_diff.cfg"
sys.argv.extend(["--diags", cfg_path])

# Run
runner.sets_to_run = ['lat_lon_land']
runner.run_diags(params)


#!/bin/bash

# Running on pm-cpu

#SBATCH  --job-name=e3sm_diags_extended
#SBATCH  --account=e3sm
#SBATCH  --nodes=1
#SBATCH  --exclusive
#SBATCH  --time=1:00:00
#SBATCH  --output=/global/cfs/cdirs/e3sm/www/chengzhu/tutorial2024/logs/e3sm_diags_extended.o%j


#SBATCH  --qos=regular
#SBATCH  --constraint=cpu


source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh

# Run E3SM Diags
echo
echo ===== RUN E3SM DIAGS =====
echo

# Prepare configuration file
cat > e3sm.py << EOF
import os
import numpy
from e3sm_diags.parameter.arm_diags_parameter import ARMDiagsParameter
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.parameter.diurnal_cycle_parameter import DiurnalCycleParameter
from e3sm_diags.parameter.enso_diags_parameter import EnsoDiagsParameter
from e3sm_diags.parameter.qbo_parameter import QboParameter
from e3sm_diags.parameter.streamflow_parameter import StreamflowParameter
from e3sm_diags.parameter.tc_analysis_parameter import TCAnalysisParameter
from e3sm_diags.parameter.tropical_subseasonal_parameter import TropicalSubseasonalParameter


from e3sm_diags.run import runner

case = 'extendedOutput.v3.LR.historical_0101'
short_name = 'v3.LR.historical_0101'
results_dir = '/global/cfs/cdirs/e3sm/www/chengzhu/tutorial2024/e3sm_diags_extended'

test_climo = '/global/cfs/cdirs/e3sm/chengzhu/tutorial2024/v3.LR.historical_0101/post/atm/180x360_aave/clim/15yr'
test_ts = '/global/cfs/cdirs/e3sm/chengzhu/tutorial2024/v3.LR.historical_0101/post/atm/180x360_aave/ts/monthly/15yr'
test_ts_daily_dir = '/global/cfs/cdirs/e3sm/chengzhu/tutorial2024/v3.LR.historical_0101/post/atm/180x360_aave/ts/daily/15yr'

ref_climo = '/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/climatology/'
ref_ts = '/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/time-series'

start_yr = 2000
end_yr = 2014

param = CoreParameter()

# Model
param.test_data_path = test_climo
param.test_name = case
param.short_test_name = short_name

# Ref/Obs
param.reference_data_path = ref_climo

# Output dir
param.results_dir = results_dir

# Additional settings
param.run_type = 'model_vs_obs'
param.diff_title = 'Model - Observations'
param.output_format = ['png']
param.output_format_subplot = []
param.multiprocessing = True
param.num_workers = 24
param.seasons = ["ANN"]
params = [param]

# Model 
enso_param = EnsoDiagsParameter()
enso_param.test_data_path = test_ts
#enso_param.test_name = short_name
enso_param.test_start_yr = start_yr
enso_param.test_end_yr = end_yr

# Obs
enso_param.reference_data_path = ref_ts
enso_param.ref_start_yr = start_yr
enso_param.ref_end_yr = end_yr

params.append(enso_param)

trop_param = TropicalSubseasonalParameter()
trop_param.test_data_path = test_ts_daily_dir
#trop_param.test_name = short_name
trop_param.test_start_yr = start_yr
trop_param.test_end_yr = end_yr

# Obs
trop_param.reference_data_path = ref_ts
trop_param.ref_start_yr = 2001
trop_param.ref_end_yr = 2010

params.append(trop_param)
qbo_param = QboParameter()
qbo_param.test_data_path = test_ts
#qbo_param.test_name = short_name
qbo_param.test_start_yr = start_yr
qbo_param.test_end_yr = end_yr
qbo_param.ref_start_yr = start_yr
qbo_param.ref_end_yr = end_yr

# Obs
qbo_param.reference_data_path = ref_ts 

params.append(qbo_param)
dc_param = DiurnalCycleParameter()
dc_param.test_data_path = '/global/cfs/cdirs/e3sm/chengzhu/tutorial2024/v3.LR.historical_0101/post/atm/180x360_aave/clim_diurnal_8xdaily/15yr/'
#dc_param.short_test_name = short_name
# Plotting diurnal cycle amplitude on different scales. Default is True
dc_param.normalize_test_amp = False

# Obs
dc_param.reference_data_path = ref_climo

params.append(dc_param)
streamflow_param = StreamflowParameter()
streamflow_param.reference_data_path = ref_ts 
streamflow_param.test_data_path = '/global/cfs/cdirs/e3sm/chengzhu/tutorial2024/v3.LR.historical_0101/post/rof/native/ts/monthly/15yr/'
#streamflow_param.test_name = short_name
streamflow_param.test_start_yr = start_yr
streamflow_param.test_end_yr = end_yr

# Obs
streamflow_param.reference_data_path = ref_ts 
streamflow_param.ref_start_yr = "1986" # Streamflow gauge station data range from year 1986 to 1995
streamflow_param.ref_end_yr = "1995"

params.append(streamflow_param)
tc_param = TCAnalysisParameter()
tc_param.test_data_path = "/global/cfs/cdirs/e3sm/chengzhu/tutorial2024/v3.LR.historical_0101/post/atm/tc-analysis_2000_2014"
#tc_param.short_test_name = short_name
tc_param.test_start_yr = start_yr
tc_param.test_end_yr = end_yr

# Obs
tc_param.reference_data_path = '/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/tc-analysis/'
# For model vs obs, the ref start and end year can be any four digit strings
# For now, use all available years from obs by default
tc_param.ref_start_yr = 1979
tc_param.ref_end_yr = 2018

params.append(tc_param)

arm_param = ARMDiagsParameter()
arm_param.reference_data_path = '/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/arm-diags-data'
arm_param.ref_name = 'armdiags'
arm_param.test_data_path = '/global/cfs/cdirs/e3sm/chengzhu/tutorial2024/v3.LR.historical_0101/post/atm/site'
arm_param.test_name = short_name
arm_param.test_start_yr = start_yr
arm_param.test_end_yr = end_yr
# For model vs obs, the ref start and end year can be any four digit strings.
# For now, will use all available years form obs
arm_param.ref_start_yr = "0001"
arm_param.ref_end_yr = "0001"

params.append(arm_param)

# Run
runner.sets_to_run = ['lat_lon', 'zonal_mean_xy', 'zonal_mean_2d', 'polar', 'cosp_histogram', 'meridional_mean_2d', 'enso_diags', 'qbo', 'diurnal_cycle', 'streamflow', 'arm_diags', 'aerosol_aeronet', 'aerosol_budget', 'tropical_subseasonal']
runner.run_diags(params)

EOF

command="srun -n 1 python -u e3sm.py"


# Run diagnostics
time ${command}
exit 0

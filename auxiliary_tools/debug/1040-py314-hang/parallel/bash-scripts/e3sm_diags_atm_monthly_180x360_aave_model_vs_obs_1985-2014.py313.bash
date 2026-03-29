#!/bin/bash

# Running on chrysalis

#SBATCH  --job-name=e3sm_diags_atm_monthly_180x360_aave_model_vs_obs_1985-2014
#SBATCH  --account=e3sm
#SBATCH  --nodes=1
#SBATCH  --output=/lcrc/group/e3sm/ac.zhang40/zppy_example_v3.2.0/v3.LR.historical_0051/post/scripts/e3sm_diags_atm_monthly_180x360_aave_model_vs_obs_1985-2014.o%j
#SBATCH  --exclusive
#SBATCH  --time=4:00:00


#SBATCH  --partition=compute


# Turn on debug output if needed
debug=False
if [[ "${debug,,}" == "true" ]]; then
  set -x
fi

# Script dir
cd /lcrc/group/e3sm/ac.zhang40/zppy_example_v3.2.0/v3.LR.historical_0051/post/scripts

# Get jobid
id=${SLURM_JOBID}

# Update status file
STARTTIME=$(date +%s)
echo "RUNNING ${id}" > e3sm_diags_atm_monthly_180x360_aave_model_vs_obs_1985-2014.status
set -e
source "${HOME}/miniforge3/etc/profile.d/conda.sh"; conda activate ed_1040_py313
set +e

# Make sure UVCDAT doesn't prompt us about anonymous logging
export UVCDAT_ANONYMOUS_LOG=False

# Basic definitions
case="v3.LR.historical_0051"
short="v3.LR.historical_0051"
www="/lcrc/group/e3sm/public_html/diagnostic_output/zppy_example/v3.2.0"
y1=1985
y2=2014
Y1="1985"
Y2="2014"

run_type="model_vs_obs"
tag="model_vs_obs"

results_dir=/lcrc/group/e3sm/public_html/ac.tvo/1040-py314-hang-tom-py313/${tag}_${Y1}-${Y2}

# Create temporary workdir
hash=`mktemp --dry-run -d XXXX`
workdir=tmp.e3sm_diags_atm_monthly_180x360_aave_model_vs_obs_1985-2014.${id}.${hash}
mkdir ${workdir}
cd ${workdir}

create_links_climo()
{
  climo_dir_source=$1
  climo_dir_destination=$2
  nc_prefix=$3
  begin_year=$4
  end_year=$5
  error_num=$6
  mkdir -p ${climo_dir_destination}
  cd ${climo_dir_destination}
  cp -s ${climo_dir_source}/${nc_prefix}_*_${begin_year}??_${end_year}??_climo.nc .
  if [ $? != 0 ]; then
    cd /lcrc/group/e3sm/ac.zhang40/zppy_example_v3.2.0/v3.LR.historical_0051/post/scripts
    echo "ERROR (${error_num})" > e3sm_diags_atm_monthly_180x360_aave_model_vs_obs_1985-2014.status
    exit ${error_num}
  fi
  cd ..
}

create_links_climo_diurnal()
{
  climo_diurnal_dir_source=$1
  climo_diurnal_dir_destination=$2
  nc_prefix=$3
  begin_year=$4
  end_year=$5
  error_num=$6
  mkdir -p ${climo_diurnal_dir_destination}
  cd ${climo_diurnal_dir_destination}
  cp -s ${climo_diurnal_dir_source}/${nc_prefix}.*_*_${begin_year}??_${end_year}??_climo.nc .
  if [ $? != 0 ]; then
    cd /lcrc/group/e3sm/ac.zhang40/zppy_example_v3.2.0/v3.LR.historical_0051/post/scripts
    echo "ERROR (${error_num})" > e3sm_diags_atm_monthly_180x360_aave_model_vs_obs_1985-2014.status
    exit ${error_num}
  fi
  cd ..
}

climo_dir_primary=climo

# Create local links to input climo files
climo_dir_source=/lcrc/group/e3sm/ac.zhang40/zppy_example_v3.2.0/v3.LR.historical_0051/post/atm/180x360_aave/clim/30yr
create_links_climo ${climo_dir_source} ${climo_dir_primary} ${case} ${Y1} ${Y2} 1


climo_diurnal_dir_primary=climo_diurnal_8xdaily

# Create local links to input diurnal cycle climo files
climo_diurnal_dir_source=/lcrc/group/e3sm/ac.zhang40/zppy_example_v3.2.0/v3.LR.historical_0051/post/atm/180x360_aave/clim_diurnal_8xdaily/30yr
create_links_climo_diurnal ${climo_diurnal_dir_source} ${climo_diurnal_dir_primary} ${case} ${Y1} ${Y2} 3

# Create xml files for time series variables
ts_dir_primary=/lcrc/group/e3sm/ac.zhang40/zppy_example_v3.2.0/v3.LR.historical_0051/post/atm/180x360_aave/ts/monthly/5yr

ts_daily_dir=/lcrc/group/e3sm/ac.zhang40/zppy_example_v3.2.0/v3.LR.historical_0051/post/atm/180x360_aave/ts/daily/5yr

ts_rof_dir_primary="/lcrc/group/e3sm/ac.zhang40/zppy_example_v3.2.0/v3.LR.historical_0051/post/rof/native/ts/monthly/5yr"





# Run E3SM Diags
echo
echo ===== RUN E3SM DIAGS =====
echo

# Prepare configuration file
cat > e3sm.py << EOF
import os
import numpy
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.parameter.diurnal_cycle_parameter import DiurnalCycleParameter
from e3sm_diags.parameter.enso_diags_parameter import EnsoDiagsParameter
from e3sm_diags.parameter.qbo_parameter import QboParameter
from e3sm_diags.parameter.streamflow_parameter import StreamflowParameter
from e3sm_diags.parameter.mp_partition_parameter import MPpartitionParameter
from e3sm_diags.parameter.precip_pdf_parameter import PrecipPDFParameter
from e3sm_diags.parameter.tc_analysis_parameter import TCAnalysisParameter
from e3sm_diags.parameter.tropical_subseasonal_parameter import TropicalSubseasonalParameter


from e3sm_diags.run import runner

short_name = '${short}'
test_ts = '${ts_dir_primary}'
start_yr = int('${Y1}')
end_yr = int('${Y2}')
num_years = end_yr - start_yr + 1
ref_start_yr = 1985

param = CoreParameter()

# Model
param.test_data_path = '${climo_dir_primary}'
param.test_name = '${case}'
param.short_test_name = short_name

# Ref

# Obs
param.reference_data_path = '/lcrc/group/e3sm/diagnostics/observations/Atm/climatology/'


# Output dir
param.results_dir = '${results_dir}_units'

# Additional settings
param.run_type = 'model_vs_obs'
param.diff_title = 'Model - Observations'
param.output_format = ['png']
param.output_format_subplot = []
param.multiprocessing = True
param.num_workers = 8
#param.fail_on_incomplete = True
params = [param]

# Model land
enso_param = EnsoDiagsParameter()
enso_param.test_data_path = test_ts
enso_param.test_name = short_name
enso_param.test_start_yr = start_yr
enso_param.test_end_yr = end_yr

# Obs
enso_param.reference_data_path = '/lcrc/group/e3sm/diagnostics/observations/Atm/time-series/'
enso_param.ref_start_yr = ref_start_yr
enso_param.ref_end_yr = ref_start_yr + 10

params.append(enso_param)
trop_param = TropicalSubseasonalParameter()
trop_param.test_data_path = '${ts_daily_dir}'
trop_param.test_name = short_name
trop_param.test_start_yr = start_yr
trop_param.test_end_yr = end_yr

# Obs
trop_param.reference_data_path = '/lcrc/group/e3sm/diagnostics/observations/Atm/time-series/'
trop_param.ref_start_yr = 2001
trop_param.ref_end_yr = 2010

params.append(trop_param)
qbo_param = QboParameter()
qbo_param.test_data_path = test_ts
qbo_param.test_name = short_name
qbo_param.test_start_yr = start_yr
qbo_param.test_end_yr = end_yr
qbo_param.ref_start_yr = ref_start_yr
ref_end_yr = ref_start_yr + num_years - 1
if (ref_end_yr <= 2014):
  qbo_param.ref_end_yr = ref_end_yr
else:
  qbo_param.ref_end_yr = 2014

# Obs
qbo_param.reference_data_path = '/lcrc/group/e3sm/diagnostics/observations/Atm/time-series/'

params.append(qbo_param)
mp_param = MPpartitionParameter()
mp_param.test_data_path = test_ts
mp_param.test_name = short_name
mp_param.short_test_name = short_name
mp_param.test_start_yr = start_yr
mp_param.test_end_yr = end_yr

params.append(mp_param)
precip_pdf_param = PrecipPDFParameter()
precip_pdf_param.test_data_path = '${ts_daily_dir}'
precip_pdf_param.test_name = short_name
precip_pdf_param.short_test_name = short_name
precip_pdf_param.test_start_yr = start_yr
precip_pdf_param.test_end_yr = end_yr

# Obs
precip_pdf_param.reference_data_path = '/lcrc/group/e3sm/diagnostics/observations/Atm/time-series/'
precip_pdf_param.ref_start_yr = 2001
precip_pdf_param.ref_end_yr = 2010

params.append(precip_pdf_param)
dc_param = DiurnalCycleParameter()
dc_param.test_data_path = '${climo_diurnal_dir_primary}'
dc_param.short_test_name = short_name
# Plotting diurnal cycle amplitude on different scales. Default is True
dc_param.normalize_test_amp = False

# Obs
dc_param.reference_data_path = '/lcrc/group/e3sm/diagnostics/observations/Atm/climatology/'

params.append(dc_param)
streamflow_param = StreamflowParameter()
streamflow_param.test_data_path = '${ts_rof_dir_primary}'
streamflow_param.test_name = short_name
streamflow_param.test_start_yr = start_yr
streamflow_param.test_end_yr = end_yr

# Obs
streamflow_param.reference_data_path = '/lcrc/group/e3sm/diagnostics/observations/Atm/time-series/'
streamflow_param.ref_start_yr = "1986" # Streamflow gauge station data range from year 1986 to 1995
streamflow_param.ref_end_yr = "1995"

params.append(streamflow_param)
tc_param = TCAnalysisParameter()
tc_param.test_data_path = "/lcrc/group/e3sm/ac.zhang40/zppy_example_v3.2.0/v3.LR.historical_0051/post/atm/tc-analysis_${Y1}_${Y2}"
tc_param.short_test_name = short_name
tc_param.test_start_yr = "${Y1}"
tc_param.test_end_yr = "${Y2}"

# Obs
tc_param.reference_data_path = '/lcrc/group/e3sm/diagnostics/observations/Atm/tc-analysis/'
# For model vs obs, the ref start and end year can be any four digit strings
# For now, use all available years from obs by default
tc_param.ref_start_yr = "1979"
tc_param.ref_end_yr = "2018"

params.append(tc_param)

# Run
runner.sets_to_run = ['lat_lon', 'zonal_mean_xy', 'zonal_mean_2d', 'polar', 'cosp_histogram', 'meridional_mean_2d', 'annual_cycle_zonal_mean', 'enso_diags', 'qbo', 'diurnal_cycle', 'zonal_mean_2d_stratosphere', 'aerosol_aeronet', 'mp_partition', 'tropical_subseasonal', 'precip_pdf', 'tc_analysis', 'streamflow']
runner.run_diags(params)

EOF

# Handle cases when cfg file is explicitly provided

command="srun -n 1 python -u e3sm.py"


# Run diagnostics
time ${command}
if [ $? != 0 ]; then
  cd /lcrc/group/e3sm/ac.zhang40/zppy_example_v3.2.0/v3.LR.historical_0051/post/scripts
  echo 'ERROR (9)' > e3sm_diags_atm_monthly_180x360_aave_model_vs_obs_1985-2014.status
  exit 9
fi

# Copy output to web server
echo
echo ===== COPY FILES TO WEB SERVER =====
echo

# Create top-level directory
web_dir=${www}/${case}/e3sm_diags/atm_monthly_180x360_aave
mkdir -p ${web_dir}
if [ $? != 0 ]; then
  cd /lcrc/group/e3sm/ac.zhang40/zppy_example_v3.2.0/v3.LR.historical_0051/post/scripts
  echo 'ERROR (10)' > e3sm_diags_atm_monthly_180x360_aave_model_vs_obs_1985-2014.status
  exit 10
fi



# Copy files
rsync -a --delete ${results_dir} ${web_dir}/
if [ $? != 0 ]; then
  cd /lcrc/group/e3sm/ac.zhang40/zppy_example_v3.2.0/v3.LR.historical_0051/post/scripts
  echo 'ERROR (11)' > e3sm_diags_atm_monthly_180x360_aave_model_vs_obs_1985-2014.status
  exit 11
fi




# For LCRC, change permissions of new files
pushd ${web_dir}/
chmod -R go+rX,go-w ${results_dir}
popd


# Delete temporary workdir
cd ..
if [[ "${debug,,}" != "true" ]]; then
  rm -rf ${workdir}
fi

# Update status file and exit

ENDTIME=$(date +%s)
ELAPSEDTIME=$(($ENDTIME - $STARTTIME))

echo ==============================================
echo "Elapsed time: $ELAPSEDTIME seconds"
echo ==============================================
rm -f e3sm_diags_atm_monthly_180x360_aave_model_vs_obs_1985-2014.status
echo 'OK' > e3sm_diags_atm_monthly_180x360_aave_model_vs_obs_1985-2014.status
exit 0

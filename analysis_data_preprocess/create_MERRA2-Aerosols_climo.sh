#!/bin/bash
# Script Generated by Jill Zhang (zhang40@llnl.gov)
# Data downloaded from https://disc.gsfc.nasa.gov/ 
#MERRA-2 tavgM_2d_aer_Nx: 2d,Monthly mean,Time-averaged,Single-Level,Assimilation,Aerosol Diagnostics V5.12.4 (M2TMNXAER 5.12.4)
#File name: MERRA2_400.tavgM_2d_aer_Nx.202112.nc4 and are being renamed to MERRA2_202112.nc4 to facilitate post-processing
path='/p/user_pub/e3sm/zhang40/analysis_data_e3sm_diags/MERRA2_Aerosols/'

original_data_path=$path'original_data/'
time_series_output_path=$path'time_series/'
climo_output_path=$path'climatology/'
tmp=$path'tmp/'


mkdir $time_series_output_path
mkdir $climo_output_path

cd ${original_data_path}
ncrcat ${original_data_path}MERRA2_*nc4 ${time_series_output_path}MERRA2_Aerosols_198001_202112.nc
cd $climo_output_path
ncclimo -a sdd -c MERRA2_198001.nc4 -s 1980 -e 2021
for i in *.nc; do mv "$i" "${i/MERRA2_/MERRA2_Aerosols_}" ; done

exit


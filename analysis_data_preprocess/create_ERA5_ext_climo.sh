#!/bin/bash

# This script is to generate climatology and time-series based on ERA5 variables not included in obs4mip archive. pr and et are duplicated for cross-validation.
path='/p/user_pub/e3sm/zhang40/analysis_data_e3sm_diags/ERA5/'

ext='ext_1'
# variables include:
# si10: 10 metre wind speed
# d2m: 2 metre dewpoint temperature
# sp: Surface pressure

filename='adaptor.mars.internal-1683916975.6853056-27160-9-445faae7-c61c-4ee6-8cb2-32d1c02b0839.nc'
#ext='ext'
# variables include:
# t2m: 2 meter temp"
# cp: Convective precipitation
# e: Evaporation
# lsp: Large-scale precipitation
# ro: Runoff
# tp: Total precipitation
# vimd: Vertically integrated moisture divergence
#filename='adaptor.mars.internal-1649447170.9796565-18358-9-69c2693a-cb8f-49a8-8f09-36e7e03c7239.nc'
original_data_path=$path'original_'${ext}'/'
time_series_output_path=$path'time_series_'${ext}'/'
climo_output_path=$path'climatology_'${ext}'/'
tmp=$path'tmp_'${ext}'/'

mkdir $time_series_output_path
mkdir $climo_output_path
mkdir $tmp


start_yr=1979
end_yr=2019

## reduce expver dimension (only needed for recent data), reference: https://confluence.ecmwf.int/display/CUSF/ERA5+CDS+requests+which+return+a+mixture+of+ERA5+and+ERA5T+data 
#cdo --reduce_dim -copy ${original_data_path}${filename} ${tmp}reduce_expver.nc

#switch latitude to N-to-S
ncpdq -a time,-lat,lon ${original_data_path}${filename} ${tmp}N-to-S.nc
#uncompress variable
ncpdq -U ${tmp}N-to-S.nc ${tmp}N-to-S_uncompress.nc

cdo splityear ${tmp}N-to-S_uncompress.nc ${tmp}ERA5_ext

for yr in $(eval echo "{$start_yr..$end_yr}"); do
    yyyy=`printf "%04d" $yr`
    echo $yyyy

    for mth in {1..12}; do
        mm=`printf "%02d" $mth`
        ncks -O -F -d time,${mth} ${tmp}ERA5_ext${yyyy}.nc ${tmp}ERA5_ext_${yyyy}${mm}.nc
        done
done

cd ${tmp}

ncclimo -a sdd -c ${tmp}ERA5_ext_${start_yr}01.nc -s $start_yr -e $end_yr
mv *climo.nc $climo_output_path

ncrcat ${tmp}ERA5_ext_*nc ${time_series_output_path}ERA5_ext_${start_yr}01_${end_yr}12.nc

# time series of variables are splitted into one variable each file, ex:
#ncks -v sp ERA5_ext_197901_201912.nc sp_197901_201912.nc

# climatology are appended
#declare -a sn=("ANN" "DJF" "MAM" "JJA" "SON" "01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12")
#for j in "${sn[@]}"
#do
#    ncks -A /p/user_pub/e3sm/zhang40/analysis_data_e3sm_diags/ERA5/climatology_ext_1/ERA5_ext_${j}_*nc /p/user_pub/e3sm/zhang40/analysis_data_e3sm_diags/ERA5/climatology_ext/ERA5_ext_${j}_*nc
#done
#

exit

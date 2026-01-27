#!/bin/bash
# Script generated to process GPCPDAY v3.3 0.5Degree(360x720) into time_series files as input of e3sm_diags by Jill Zhang (zhang40@llnl.gov)
# GPCPDAY v3.3 related URL:https://disc.gsfc.nasa.gov/datasets/GPCPDAY_3.3/summary 

path='/global/cfs/cdirs/e3sm/zhang40/analysis_data_e3sm_diags/GPCPDAY_0.5D_3.3/'

time_series_output_path=$path'time_series/'
original_path=$path'original/'
tmp=$path'tmp/'

mkdir $time_series_output_path
mkdir $tmp

cd $original_path
echo $path
start_yr=1998
end_yr=2023

##
for yr in $(eval echo "{$start_yr..$end_yr}"); do
    yyyy=`printf "%04d" $yr`
    echo $yyyy
    for mth in {1..12}; do
        mm=`printf "%02d" $mth`
	# number of days in month (handles leap years)
        ndays=$(cal "$mth" "$yr" | awk 'NF {d=$NF} END {print d}')
    
        for day in $(seq 1 "$ndays"); do
          dd=$(printf "%02d" "$day")
    
          in_file="${original_path}GPCPDAY_L3_${yyyy}${mm}${dd}_V3.3.nc4"
          out_file="${tmp}time_rec_${yyyy}${mm}${dd}.nc"
    
          # assign time as record dimension
          ncks --mk_rec_dmn time "$in_file" "$out_file"
    
          # If you later want to rotate lon, you can adapt the commented block per-file here.
          # Example (uncomment/adapt if needed):
          # ncks -O --msa_usr_rdr -d lon,0.25,179.75 -d lon,-179.75,-0.25 "$out_file" "${tmp}GPCP_v3.3_${yyyy}${mm}${dd}.nc"
          # ncap2 -O -s 'where(lon < 0) lon=lon+360' "${tmp}GPCP_v3.3_${yyyy}${mm}${dd}.nc" "${tmp}GPCP_v3.3_${yyyy}${mm}${dd}.nc"
    
        done
    done
done

cd ${tmp};eval ls time_rec_{${start_yr}..${end_yr}}*.nc | ncrcat ${tmp}time_rec_*nc ${time_series_output_path}precip_${start_yr}01_${end_yr}12.nc


# Regrid from 0.5 degree to 1 degree
# ncremap -m /global/cfs/cdirs/e3sm/diagnostics/maps/map_r05_to_cmip6_180x360_esmfaave.20231110.nc /global/cfs/cdirs/e3sm/diagnostics/observations/Atm/time-series/GPCPDAY_0.5D_3.3/precip_199801_202312.nc /global/cfs/cdirs/e3sm/diagnostics/observations/Atm/time-series/GPCPDAY_3.3/precip_199801_202312.nc

exit



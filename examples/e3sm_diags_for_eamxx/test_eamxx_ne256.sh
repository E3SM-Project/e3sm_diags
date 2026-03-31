#!/bin/bash                                 

source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh
#
# Test case 1
#drc_in=/global/cfs/cdirs/e3sm/beydoun/ne256pg2_ne256pg2.F20TR-SCREAMv1.rainfrac1.spanc1000.auto2700.acc150.n0128/run
#drc_out=/pscratch/sd/c/chengzhu/ne256pg2_ne256pg2.F20TR-SCREAMv1.rainfrac1.spanc1000.auto2700.acc150.n0128

# Test case 2
drc_in=/pscratch/sd/z/zhan391/EAMxx/ne256pg2_ne256pg2.F20TR-SCREAMv1.July-1.spanc800.2xauto.acc150.n0032.test2.1/run
drc_out=/pscratch/sd/c/chengzhu/ne256pg2_ne256pg2.F20TR-SCREAMv1.July-1.spanc800.2xauto.acc150.n0032.test2.1
map_file=/global/cfs/projectdirs/e3sm/zender/maps/map_ne30pg2_to_cmip6_180x360_traave.20231201.nc
#Monthly averaged: 1ma_ne30pg2.AVERAGE.nmonths_x1.
#3hourly instantaneous: 3hi_ne30pg2.INSTANT.nhours_x3.
#1hourly instantaneous: 1hi_ne30pg2.INSTANT.nhours_x1.
start=1996
end=2005

echo "Generating Climatology files"
mkdir -p $drc_out
drc_rgr=${drc_out}/rgr/climo
drc=${drc_out}/native/climo
echo $drc_out
#cd ${drc_in};eval ls 1ma_ne30pg2.AVERAGE.nmonths_x1.*{${start}..${end}}*.nc | ncclimo -P eamxx -p serial --fml_nm=1ma_ne30pg2.AVERAGE.nmonths_x1 --yr_srt=${start} --yr_end=${end} --drc_out=$drc -O $drc_rgr --map=${map_file}

# Diurnal Cycle Climo: commented out pending ncclimo update to correctly handle
# EAMxx instantaneous 3-hourly output without time bounds.
# https://github.com/E3SM-Project/e3sm_diags/pull/1046/changes#r2991031741
#drc_rgr=${drc_out}/rgr/climo_diurnal_3hrly
#drc=${drc_out}/native/climo_diurnal_3hrly
#cd ${drc_in};eval ls 3hi_ne30pg2.INSTANT.nhours_x3.*{${start}..${end}}*-10800.nc | ncclimo -P eamxx --clm_md=hfc --caseid=3hi_ne30pg2.INSTANT.nhours_x3 -v precip_liq_surf_mass_flux,precip_ice_surf_mass_flux --yr_srt=${start} --yr_end=${end} --drc_out=${drc} -O $drc_rgr --map=${map_file} #--tpd=8
##
echo "Generating per-variable time-series from 3hourly instantaneous output."
drc_rgr=${drc_out}/rgr/ts_3hrly
drc=${drc_out}/native/ts_3hrly

## Make sure to include --tpd=8 explicitly when 3hi_ne30pg2.INSTANT.nhours_x3 stream does not include time bounds
cd ${drc_in};eval ls 3hi_ne30pg2.INSTANT.nhours_x3.*{${start}..${end}}*-10800.nc | ncclimo -P eamxx --clm_md=hfs --caseid=3hi_ne30pg2.INSTANT.nhours_x3 --var=precip_liq_surf_mass_flux,precip_ice_surf_mass_flux,U_at_850hPa,LW_flux_up_at_model_top --yr_srt=$start --yr_end=$end --drc_out=${drc} -O $drc_rgr --map=${map_file} --split --tpd=8 -p serial #--thr_nbr=1

echo "Average 3hrly to daily"
drc_rgr=${drc_out}/rgr/ts_daily
mkdir -p $drc_rgr
cd ${drc_out}/rgr/ts_3hrly
for i in *.nc; do ncra --mro -O -d time,0,,8,8 "$i" "${drc_rgr}/$i";done

exit



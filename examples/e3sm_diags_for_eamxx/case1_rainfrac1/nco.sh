#!/bin/bash                                 

source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh

drc_in=/global/cfs/cdirs/e3sm/chengzhu/eamxx/run
drc_out=/global/cfs/cdirs/e3sm/chengzhu/eamxx/post/data
caseid=output.scream.decadal.monthlyAVG_ne30pg2.AVERAGE.nmonths_x1

# spoofed climatology files with data from 1995-09 to 1996-08

# create climatology files
cd ${drc_in};ls ${caseid}*1996-0[1-8]*.nc ${caseid}*1995-09*.nc ${caseid}*1995-1[0-2]*.nc | ncclimo -P eamxx --fml_nm=eamxx_decadal --yr_srt=1996 --yr_end=1996 --drc_out=$drc_out


map=/global/cfs/projectdirs/e3sm/zender/maps/map_ne30pg2_to_cmip6_180x360_traave.20231201.nc
# remaping climo files to regular lat-lon
cd $drc_out;ls *.nc | ncremap -P eamxx --prm_opt=time,lwband,swband,ilev,lev,plev,cosp_tau,cosp_cth,cosp_prs,dim2,ncol --map=${map} --drc_out=${drc_out}/rgr

exit


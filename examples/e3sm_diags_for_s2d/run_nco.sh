#!/bin/bash
#
# Generate per-variable monthly time series from E3SM S2D output.
# This script uses ncclimo to split, remap, and create time series
# files suitable for e3sm_diags with a non-January start month.
#
# Usage: bash run_nco.sh

source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh

# --- Configuration ---
caseid="WCYCL20TR_ne30pg2_r05_IcoswISC30E3r5_JRA55_FOSIRL_1980050100.EN00"
atm_name="eam"       # Use "cam" for v1 or "eam" for v2 production simulations

start='1980'
end='1982'
# Define the annual cycle: May (month 5) through April (month 4).
# This aligns with the start_month=5 setting in run_e3sm_diags_s2d.py.
month_start='5'
month_end='4'

input_path=/pscratch/sd/z/zhan391/e3sm_project/E3SMv3_S2D/WCYCL20TR_ne30pg2_r05_IcoswISC30E3r5_JRA55_FOSIRL_1980050100/EN00/archive/atm/hist
result_dir=/pscratch/sd/c/chengzhu/${caseid}
map_file=/global/cfs/projectdirs/e3sm/zender/maps/map_ne30pg2_to_cmip6_180x360_traave.20231201.nc

# --- Generate time series ---
drc_rgr=${result_dir}/ts/rgr
drc_out=${result_dir}/ts/native

echo "Generating per-variable monthly time series."
cd ${input_path}
eval ls ${caseid}.${atm_name}.h0.*{${start}..${end}}*.nc \
    | ncclimo -P ${atm_name} \
        --caseid=${caseid} \
        --mth_srt=${month_start} \
        --mth_end=${month_end} \
        --var=U,CLDICE,CLDLIQ,CLDHGH,CLDLOW,CLDMED,CLDTOT,FLNS,FLUT,FSNS,FSNT,FSNTOA,LANDFRAC,LHFLX,LWCF,OCNFRAC,PRECC,PRECL,PSL,QFLX,SHFLX,SWCF,T,TAUX,TAUY,TREFHT,TS \
        --yr_srt=$start \
        --yr_end=$end \
        --drc_out=${drc_out} \
        -O $drc_rgr \
        --map=${map_file} \
        --split

# --- Alternative: Generate climatologies with flexible-months ncclimo ---
# The latest NCO snapshot supports generating climatologies starting from
# any month (not just January). Use --mth_srt and --mth_end to define the
# annual cycle. For climos, mth_end should precede mth_srt by one month so
# that an integral number of years are supplied.
# On Chrysalis, use ~ac.zender/bin/ncclimo instead of ~zender/bin/ncclimo.
drc_rgr_ncclimo_flex_month=${result_dir}/ts/rgr_ncclimo_flex_month
drc_out_ncclimo_flex_month=${result_dir}/ts/native_ncclimo_flex_month

echo "Generating climatologies with flexible-months ncclimo."
cd ${input_path}
eval ls ${caseid}.${atm_name}.h0.*{${start}..${end}}*.nc | /global/cfs/cdirs/e3sm/zender/bin/ncclimo --npo -P ${atm_name} \
    --caseid=${caseid} \
    --mth_srt=${month_start} \
    --mth_end=${month_end} \
    --yr_srt=$start \
    --yr_end=$end \
    -o $drc_out_ncclimo_flex_month \
    -O $drc_rgr_ncclimo_flex_month \
    --map=${map_file}


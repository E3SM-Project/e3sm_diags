#!/bin/bash
mdl="-P eam " # Model-specific options, with out -v to include all variables
clm="-c extendedOutput.v3.LR.historical_0101 -s 2000 -e 2014" # Climo. params  
hrz="--map=/global/homes/z/zender/data/maps/map_ne30pg2_to_cmip6_180x360_traave.20231201.nc"  
drc_in=/global/cfs/cdirs/e3sm/www/Tutorials/2024/simulations/extendedOutput.v3.LR.historical_0101/archive/atm/hist  
in="-i ${drc_in}" # Input directory  
out="-o ${SCRATCH}/clm -O ${SCRATCH}/rgr" # Output directories

ncclimo ${mdl} ${clm} ${vrt} ${hrz} ${in} ${out}  
exit 0

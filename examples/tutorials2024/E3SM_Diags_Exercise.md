# E3SM Diags Exercise

This exercise will help you learn two ways for running E3SM Diags.
Try to solve issues by looking through the [documentation](https://docs.e3sm.org/e3sm_diags).

## Run E3SM Diags standalone with a run script

**Step 0**: Obtain an interactive compute node and then activate `E3SM Unified` on Perlmutter:
```
salloc --nodes 1 --qos interactive --time 1:00:00 --constraint cpu --reservation=e3sm_day2 -A ntrain6

source /global/common/software/e3sm/anaconda_envs
/load_latest_e3sm_unified_pm-cpu.sh
```
**Step 1**: Generate climatology with ncclimo, an example bash script (Based on Charlie's Day 1 lecture) as following:
```
#!/bin/bash
source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh
mdl="-P eam " # Model-specific options, with out -v to include all variables
clm="-c extendedOutput.v3.LR.historical_0101 -s 2000 -e 2014" # Climo. params  
hrz="--map=/global/homes/z/zender/data/maps/map_ne30pg2_to_cmip6_180x360_traave.20231201.nc"  
drc_in=/global/cfs/cdirs/e3sm/www/Tutorials/2024/simulations/extendedOutput.v3.LR.historical_0101/archive/atm/hist  
in="-i ${drc_in}" # Input directory  
out="-o ${SCRATCH}/clm -O ${SCRATCH}/rgr" # Output directories

ncclimo ${mdl} ${clm} ${vrt} ${hrz} ${in} ${out}  
exit 0
```
The ncclimo run takes about 7 mins.

**Step 2**: Generate a Python script to config an E3SM Diags run.  Get an example run script using: `wget https://github.com/E3SM-Project/e3sm_diags/examples/tutorials2024/examples/tutorials2024/e3sm_diags_core_sets.py`

Update the `html_prefix` with your user name and `test_data_path`  with the data path generated form step 1 (${SCRATCH}/rgr). 

You can run it directly with `srun -n 1 python e3sm_diags_core_sets.py`
or wrap it in a batch script ([example](https://github.com/E3SM-Project/e3sm_diags/examples/tutorials2024/examples/tutorials2024/e3sm_diags_core_sets.bash)) and submit with `sbatch`.

The E3SM Diags run takes about 18 mins.

**Step 3**: Go over results at https://portal.nersc.gov/cfs/ntrain6/your_user_name/tutorial2024/e3sm_diags_core_sets/viewer/
Other than the diagnostics plots, notice some useful metadata to be saved:

-  In the main viewer look for a provenance folder where you can locate the the run script and log file.
- Under each set of plots,  look for `Hide Output Metadata`, click to show the command line to produce the plots.

**Step 4(Optional, advanced)**: If you have another longer than one year simulation available, repeat step 1 with modification to generate climatology files for your run. Generate a new Python script based on `e3sm_diags_core_sets.py` (looking for changes needed for model vs model). Run the new script to generate a model vs model comparison.    

## Run E3SM Diags with zppy ( Homework)
If you are running with extended diagnostics sets and with long simulation, in which for E3SM Diags to run every 10 20, and 50 years simulation interval. The recommended way is to set up a zppy configuration file, which takes care of pre-processing, and can launch many runs automatically at once.

**Step 0**: Activate  `E3SM Unified` on Perlmutter to get zppy:
```
source /global/common/software/e3sm/anaconda_envs
/load_latest_e3sm_unified_pm-cpu.sh
```
**Step 1**: Generate a `zppy` configuration file following the example, 
wget https://github.com/E3SM-Project/e3sm_diags/examples/tutorials2024/examples/tutorials2024/zppy_cfg_e3sm_diags_v3.cfg

Update `output` and `www` paths with your username.

**Step2**: Run `zppy` with `zppy -c zppy_cfg_e3sm_diags_v3.cfg`

**Step3**: View results at https://portal.nersc.gov/cfs/ntrain6/your_user_name/tutorial2024/v3.LR.historical_0101

More details about zppy will be covered in zppy post-processing practicum in Day 3.

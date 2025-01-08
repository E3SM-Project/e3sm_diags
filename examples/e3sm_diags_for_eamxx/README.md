# Initial Instruction to Run E3SM Diags on EAMxx output (e.g. monthly ne30pg2 output)

0. Secure an interactive compute node and to activate the E3SM-Unified enviroment:

salloc --nodes 1 --qos interactive --time 02:00:00 --constraint cpu --account e3sm

source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh

(The version of E3SM Diags (v3) that has EAMxx variable support is available in E3SM-Unified v1.11 (Mid Feb 2025 release. ) 

1. To remap  monthly ne30pg2 data to regular lat-lon data to prepare for E3SM Diags run. An example usage based on a EAMxx decadal run is provided in following script ``nco.sh``. To run the script:

bash nco.sh

2. Generate a python script for running E3SM Diags. Two example is provided here:

python run_e3sm_diags_1996.py: to compare 1996 climatology from EAMxx to available 1990 obs climatology

python run_e3sm_diags_climo.py: to compare 1996 climatology from EAMxx to pre-calculated obs climatology 





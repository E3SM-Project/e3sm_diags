#!/bin/bash
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --constraint=cpu
#SBATCH --time=00:30:00
#SBATCH --reservation=e3sm_dryrun
#SBATCH --account=ntrain6

source /global/common/software/e3sm/anaconda_envs/test_e3sm_unified_1.10.0rc2_pm-cpu.sh

cd /global/homes/c/chengzhu/e3sm_diags/examples/tutorials2024
time python e3sm_diags_core_sets.py




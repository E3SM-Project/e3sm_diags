#!/bin/bash
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --constraint=cpu
#SBATCH --time=00:30:00
#SBATCH --reservation=e3sm_day2
#SBATCH --account=ntrain6

source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh

cd /global/homes/c/chengzhu/e3sm_diags/examples/tutorials2024
srun -n 1 python e3sm_diags_core_sets.py




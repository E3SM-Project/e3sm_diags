#!/usr/bin/env python3
"""
Simple run script to test the updated mp_partition driver.
Based on auxiliary_tools/cdat_regression_testing/871-mp-partition/run_script.py
"""

import os

from e3sm_diags.parameter.mp_partition_parameter import MPpartitionParameter
from e3sm_diags.run import runner

# Test data path
TEST_DATA_PATH = "/lcrc/group/e3sm2/ac.zhang40/E3SMv3/v3.LR.historical_0101/post/atm/180x360_aave/ts/monthly/5yr"

# Create the mp_partition parameter object
mp_param = MPpartitionParameter()

# Set up basic parameters
mp_param.test_data_path = TEST_DATA_PATH
#mp_param.test_name = "20210528.v2rc3e.piControl.ne30pg2_EC30to60E2r2.chrysalis"
#mp_param.short_test_name = "e3sm_v2"
#mp_param.test_start_yr = "0051"
#mp_param.test_end_yr = "0060"
mp_param.test_name = "v3.LR.historical_0101"
mp_param.short_test_name = "e3sm_v3"
mp_param.test_start_yr = "2000"
mp_param.test_end_yr = "2014"

# Set results directory
mp_param.results_dir = "/lcrc/group/e3sm/public_html/diagnostic_output/ac.zhang40/tests/1028-LCF_mixed-phase_partition_update"

# Set run type to model-vs-obs (uses benchmark data)
mp_param.run_type = "model_vs_obs"

# Enable multiprocessing for faster execution
mp_param.multiprocessing = True
mp_param.num_workers = 2

# Save netCDF files for debugging
mp_param.save_netcdf = True

# Optional: Enable debugging output
mp_param.debug = False

print("=" * 60)
print("MP_PARTITION TEST CONFIGURATION")
print("=" * 60)
print(f"Test data path: {mp_param.test_data_path}")
print(f"Test name: {mp_param.test_name}")
print(f"Test years: {mp_param.test_start_yr}-{mp_param.test_end_yr}")
print(f"Results directory: {mp_param.results_dir}")
print(f"Run type: {mp_param.run_type}")
print(f"Multiprocessing: {mp_param.multiprocessing}")
print("=" * 60)

# Verify test data path exists
if not os.path.exists(mp_param.test_data_path):
    print(f"ERROR: Test data path does not exist: {mp_param.test_data_path}")
    print("Please update TEST_DATA_PATH in this script to point to your E3SM data.")
    exit(1)

# Create results directory if it doesn't exist
os.makedirs(mp_param.results_dir, exist_ok=True)

# Set which diagnostic set to run
runner.sets_to_run = ["mp_partition"]

# Run the diagnostics
try:
    print("Starting mp_partition diagnostic run...")
    runner.run_diags([mp_param])
    print(f"Results saved to: {mp_param.results_dir}")
except Exception as e:
    print("=" * 60)
    print(f"ERROR: MP_PARTITION diagnostic failed: {e}")
    print("=" * 60)
    raise

#!/usr/bin/env python3
"""
This script runs e3sm_diags with the lat_lon_native set to visualize native grid data.
"""

import os
import sys

from e3sm_diags.parameter.lat_lon_native_parameter import LatLonNativeParameter
from e3sm_diags.run import runner

# Auto-detect username
username = os.environ.get('USER', 'unknown_user')

# Create parameter objects for 3 different runs
params = []

## (1) First test configuration
#param1 = LatLonNativeParameter()
#param1.results_dir = f"/lcrc/group/e3sm/public_html/diagnostic_output/{username}/tests/lat_lon_native_test_1"
#param1.test_data_path = "/lcrc/group/e3sm/public_html/e3sm_diags_test_data/native_grid"
#param1.test_name = "v3.LR.amip_0101"
#param1.short_test_name = "v3.LR.amip_0101"
#param1.reference_data_path = "/lcrc/group/e3sm/public_html/e3sm_diags_test_data/native_grid"
#param1.ref_name = "v3.HR.test4"
#param1.short_ref_name = "v3.HR.test4"
#param1.seasons = ["DJF"]
#param1.test_grid_file = "/lcrc/group/e3sm/diagnostics/grids/ne30pg2.nc"
#param1.ref_grid_file = "/lcrc/group/e3sm/diagnostics/grids/ne120pg2.nc"
#param1.case_id = "model_vs_model"
#param1.run_type = "model_vs_model"
#params.append(param1)
#
## (2) Second test configuration
#param2 = LatLonNativeParameter()
#param2.results_dir = f"/lcrc/group/e3sm/public_html/diagnostic_output/{username}/tests/lat_lon_native_test_2"
#param2.test_data_path = "/lcrc/group/e3sm/public_html/e3sm_diags_test_data/native_grid"
#param2.test_file = "v3.LR.amip_0101_DJF_climo.nc"
#param2.short_test_name = "v3.LR.amip_0101"
#param2.reference_data_path = "/lcrc/group/e3sm/public_html/e3sm_diags_test_data/native_grid"
#param2.ref_file = "v3.HR.test4_DJF_climo.nc"
#param2.short_ref_name = "v3.HR.test4"
#param2.seasons = ["DJF"]
#param2.test_grid_file = "/lcrc/group/e3sm/diagnostics/grids/ne30pg2.nc"
#param2.ref_grid_file = "/lcrc/group/e3sm/diagnostics/grids/ne120pg2.nc"
#param2.case_id = "model_vs_model"
#param2.run_type = "model_vs_model"
#params.append(param2)

# (3) Third test configuration
param3 = LatLonNativeParameter()
param3.results_dir = f"/lcrc/group/e3sm/public_html/diagnostic_output/{username}/tests/lat_lon_native_test_3"
param3.test_data_path = "/lcrc/group/e3sm/public_html/e3sm_diags_test_data/native_grid"
param3.test_file = "v3.LR.amip_0101.eam.h0.1989-12.nc"
param3.reference_data_path = "/lcrc/group/e3sm/public_html/e3sm_diags_test_data/native_grid"
param3.ref_file = "v3.LR.amip_0101.eam.h0.1989-12.nc"
param3.short_ref_name = "v3.HR.test4"
param3.time_slices = ["0"]
param3.test_grid_file = "/lcrc/group/e3sm/diagnostics/grids/ne30pg2.nc"
param3.ref_grid_file = "/lcrc/group/e3sm/diagnostics/grids/ne30pg2.nc"
param3.case_id = "model_vs_model"
param3.run_type = "model_vs_model"
params.append(param3)

# Run the single diagnostic, comment out for complete diagnostics.
cfg_path = "auxiliary_tools/debug/968-native-grid-vis/TGCLDLWP.cfg"
sys.argv.extend(["--diags", cfg_path])

runner.sets_to_run = ["lat_lon_native"]

# Run each test sequentially
for i, param in enumerate(params, 1):
    print(f"\n{'='*60}")
    print(f"Running Test {i}: {param.results_dir}")
    print(f"{'='*60}")

    # Create results directory
    if not os.path.exists(param.results_dir):
        os.makedirs(param.results_dir)

    # Run the diagnostic
    runner.run_diags([param])
    print(f"Test {i} completed!")

print(f"\n{'='*60}")
print("All tests completed!")
print(f"{'='*60}")



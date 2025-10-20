#!/usr/bin/env python3
"""
Example 8: Native Grid Visualization

This example demonstrates how to visualize model data on its native grid
(e.g., cubed-sphere, unstructured grids) using UXarray, without regridding
to a regular lat-lon grid.

This preserves native grid features and is particularly useful for:
- High-resolution models with complex grid structures
- Comparing models with different native grids
- Analyzing grid-dependent features

This feature was introduced in E3SM Diags v3.1.0.
"""

import os

from e3sm_diags.parameter.lat_lon_native_parameter import LatLonNativeParameter
from e3sm_diags.run import runner

# Auto-detect username
username = os.environ.get('USER', 'unknown_user')

# Create parameter object
param = LatLonNativeParameter()

# Location of the data
param.test_data_path = "/lcrc/group/e3sm/public_html/e3sm_diags_test_data/native_grid"
param.test_file = "v3.LR.amip_0101.eam.h0.1989-12.nc"

param.reference_data_path = "/lcrc/group/e3sm/public_html/e3sm_diags_test_data/native_grid"
param.ref_file = "v3.LR.amip_0101.eam.h0.1989-12.nc"

# Short names for display
param.short_test_name = "v3.LR.amip_0101"
param.short_ref_name = "v3.HR.test4"

# Native grid files
# These specify the grid structure for native visualization
param.test_grid_file = "/lcrc/group/e3sm/diagnostics/grids/ne30pg2.nc"
param.ref_grid_file = "/lcrc/group/e3sm/diagnostics/grids/ne30pg2.nc"

# Time selection: use snapshot instead of climatology
# Use time_slices to analyze a specific time index
param.time_slices = ["0"]  # First time step

# Comparison settings
param.case_id = "model_vs_model"
param.run_type = "model_vs_model"

# Antialiasing setting
param.antialiased = False

# Name of the folder where the results are stored.
prefix = f"/lcrc/group/e3sm/public_html/diagnostic_output/{username}/e3sm_diags_examples"
param.results_dir = os.path.join(prefix, "ex8_native_grid")

# Below are more optional arguments.

# For running with multiprocessing.
# param.multiprocessing = True
# param.num_workers = 32

# Run the diagnostic
runner.sets_to_run = ["lat_lon_native"]
runner.run_diags([param])

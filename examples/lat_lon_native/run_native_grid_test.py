#!/usr/bin/env python3
"""
This script runs e3sm_diags with the lat_lon_native set to visualize native grid data.
"""

import os

from e3sm_diags.parameter.lat_lon_native_parameter import LatLonNativeParameter
from e3sm_diags.run import runner

# Create parameter object
param = LatLonNativeParameter()

# Auto-detect username
username = os.environ.get('USER', 'unknown_user')

# Basic parameters
param.results_dir = f"/lcrc/group/e3sm/public_html/diagnostic_output/{username}/tests/lat_lon_native_test_TGCLDLWP"

# Create results directory if it doesn't exist
if not os.path.exists(param.results_dir):
    os.makedirs(param.results_dir)

# Model data
param.test_data_path = "/lcrc/group/e3sm/public_html/e3sm_diags_test_data/native_grid"
param.test_file = "v3.LR.amip_0101.eam.h0.1989-12.nc"

param.reference_data_path = "/lcrc/group/e3sm/public_html/e3sm_diags_test_data/native_grid"
param.ref_file = "v3.LR.amip_0101.eam.h0.1989-12.nc"
param.short_ref_name = "v3.HR.test4"

param.case_id = "model_vs_model"

param.variables = ["TGCLDLWP"]
param.regions = ["global"]
## Variables to plot
# param.variables = ["PRECC"]
# param.seasons = ["ANN"]
# param.seasons = ["DJF"]
# param.regions = ["global"]
# param.regions = ["60S60N"]

# Time slices for snapshot-based analysis
param.time_slices = ["0"]
param.seasons = ["ANN"]

# Native grid settings
param.test_grid_file = "/lcrc/group/e3sm/diagnostics/grids/ne30pg2.nc"
param.ref_grid_file = "/lcrc/group/e3sm/diagnostics/grids/ne30pg2.nc"

param.antialiased = False

## No reference data for this test - model only
# param.model_only = True
param.run_type = "model_vs_model"

# Plot settings
# param.contour_levels = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15, 20]
# param.test_colormap = 'WhiteBlueGreenYellowRed.rgb'
# param.test_colormap = "plasma"

# Run the diagnostic
runner.sets_to_run = ["lat_lon_native"]
runner.run_diags([param])

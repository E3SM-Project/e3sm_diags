#!/usr/bin/env python3
"""
This script runs e3sm_diags with the lat_lon_native set to visualize native grid data.
"""

import os

from e3sm_diags.parameter.lat_lon_native_parameter import LatLonNativeParameter
from e3sm_diags.run import runner

# Create parameter object
param = LatLonNativeParameter()

# Basic parameters
param.results_dir = os.path.expanduser(
    "~/Documents/repos/e3sm_diags/results_native_grid"
)

# Create results directory if it doesn't exist
if not os.path.exists(param.results_dir):
    os.makedirs(param.results_dir)


# Model data
param.test_data_path = os.path.expanduser(
    "~/Documents/ACME_simulations/E3SM_v2/native_grid_data/"
)
param.test_name = "v3.LR.historical_0051"

param.reference_data_path = os.path.expanduser(
    "~/Documents/ACME_simulations/E3SM_v2/native_grid_data/"
)
param.ref_name = "v3.HR.test4"

param.case_id = "model_vs_model"

## Variables to plot
# param.variables = ["PRECC"]
# param.seasons = ["ANN"]
# param.regions = ["global"]
# param.regions = ["60S60N"]
# param.regions = ["global", "60S60N"]

# Native grid settings
param.test_grid_file = os.path.expanduser(
    "~/Documents/ACME_simulations/E3SM_v2/native_grid_data/ne30pg2.nc"
)

param.ref_grid_file = os.path.expanduser(
    "~/Documents/ACME_simulations/E3SM_v2/native_grid_data/ne120pg2.nc"
)

param.split_periodic_elements = True
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

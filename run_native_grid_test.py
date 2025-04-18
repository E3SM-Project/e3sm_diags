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
param.results_dir = os.path.expanduser('~/Documents/repos/e3sm_diags/results_native_grid')

# Create results directory if it doesn't exist
if not os.path.exists(param.results_dir):
    os.makedirs(param.results_dir)


# Model data
param.test_data_path = os.path.expanduser('~/Documents/ACME_simulations/E3SM_v2/native_grid_data/')
param.test_name = 'v3.LR.historical_0051'
param.case_id = "model_vs_model"

# Variables to plot
param.variables = ['PRECC']
param.seasons = ['ANN']
param.regions = ['global']

# Native grid settings
param.test_grid_file = os.path.expanduser('~/Documents/ACME_simulations/E3SM_v2/native_grid_data/ne30pg2.nc')
param.split_periodic_elements = True
param.antialiased = False

# No reference data for this test - model only
param.model_only = True

# Plot settings
param.contour_levels = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
param.test_colormap = 'WhiteBlueGreenYellowRed.rgb'

# Run the diagnostic
runner.sets_to_run = ["lat_lon_native"]
runner.run_diags([param])

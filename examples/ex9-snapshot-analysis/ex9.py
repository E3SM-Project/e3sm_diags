#!/usr/bin/env python3
"""
Example 9: Snapshot Analysis for Core Sets

This example demonstrates time slice analysis on core diagnostic sets.
Instead of computing climatological seasonal means, this analyzes individual
time steps from model output using index-based time selection.

This is useful for:
- Analyzing specific events or time periods
- Comparing model states at specific time points
- Understanding temporal evolution without time averaging
- Event-based or process-oriented diagnostics

This feature was introduced in E3SM Diags v3.1.0.

Supported diagnostic sets:
- lat_lon
- lat_lon_native
- polar
- zonal_mean_2d
- meridional_mean_2d
- zonal_mean_2d_stratosphere
"""

import os

from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.run import runner

# Create parameter object
param = CoreParameter()

# Location of the data
param.test_data_path = "/lcrc/group/e3sm/public_html/e3sm_diags_test_data/postprocessed_e3sm_v2_data_for_e3sm_diags/20210528.v2rc3e.piControl.ne30pg2_EC30to60E2r2.chrysalis/time-series/rgr"
param.test_file = "T_005101_006012.nc"

param.reference_data_path = "/lcrc/group/e3sm/public_html/e3sm_diags_test_data/postprocessed_e3sm_v2_data_for_e3sm_diags/20210528.v2rc3e.piControl.ne30pg2_EC30to60E2r2.chrysalis/time-series/rgr"
param.ref_file = "T_005101_006012.nc"

# Short names for display
param.short_test_name = "v2 test"
param.short_ref_name = "v2 test"

# Key difference: Use time_slices instead of seasons
# Time slices are zero-based indices into the time dimension
# Examples:
#   ["0"] = first time step
#   ["5"] = 6th time step
#   ["0", "1", "2"] = first 3 time steps (analyzed separately)
param.time_slices = ["0", "1"]

# IMPORTANT: time_slices and seasons are mutually exclusive
# When using time_slices, do NOT set param.seasons

# Comparison settings
param.case_id = "model_vs_model"
param.run_type = "model_vs_model"

# Name of the folder where the results are stored.
# Change `prefix` to use your directory.
prefix = "/global/cfs/cdirs/e3sm/www/<your directory>/examples"
param.results_dir = os.path.join(prefix, "ex9_snapshot_analysis")

# Below are more optional arguments.

# For running with multiprocessing.
# param.multiprocessing = True
# param.num_workers = 32

# Run the diagnostics on multiple sets
runner.sets_to_run = [
    "lat_lon",
    "zonal_mean_xy",
    "zonal_mean_2d",
    "zonal_mean_2d_stratosphere",
    "polar",
    "meridional_mean_2d",
]

runner.run_diags([param])

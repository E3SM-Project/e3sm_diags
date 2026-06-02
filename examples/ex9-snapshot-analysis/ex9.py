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

# Auto-detect username
username = os.environ.get('USER', 'unknown_user')

# Create parameter object
param = CoreParameter()

# Location of the data
param.test_data_path = "/lcrc/group/e3sm/public_html/e3sm_diags_test_data/postprocessed_e3sm_v2_data_for_e3sm_diags/20210528.v2rc3e.piControl.ne30pg2_EC30to60E2r2.chrysalis/time-series/rgr"
param.test_file = "T_005101_006012.nc"

param.reference_data_path = "/lcrc/group/e3sm/public_html/e3sm_diags_test_data/postprocessed_e3sm_v2_data_for_e3sm_diags/20210528.v2rc3e.piControl.ne30pg2_EC30to60E2r2.chrysalis/time-series/rgr"
param.ref_file = "T_005101_006012.nc"

# To compare the model against observations at the same calendar time, point the
# reference at an observational time-series file (with a time dimension and the
# same variable name) and use date-based time_slices. Because each date selects
# the nearest step independently for the test and reference, the comparison is
# aligned even if the model and obs time axes differ. For example:
#
#   param.run_type = "model_vs_obs"
#   param.reference_data_path = "/path/to/obs/time-series"
#   param.ref_file = "T_obs_time_series.nc"
#   param.reference_name = "ERA5"
#   param.time_slices = ["2010-01", "2010-07"]

# Short names for display
param.short_test_name = "v2 test"
param.short_ref_name = "v2 test"

# Key difference: Use time_slices instead of seasons.
#
# time_slices entries are either:
#   - Positional indices (zero-based) into the time dimension:
#       ["0"]            = first time step
#       ["0", "1", "2"]  = first 3 time steps (analyzed separately)
#     The same index is applied to both the test and reference datasets.
#   - Dates in "YYYY-MM" or "YYYY-MM-DD" format:
#       ["0051-01"]      = the step nearest to Jan 0051
#       ["0051-01-15"]   = same step as "0051-01" for monthly data
#     Each date selects the time step *nearest* to that date, applied
#     independently to the test and reference datasets. This keeps the test and
#     reference aligned on the same calendar time even when their time axes
#     differ in start date, cadence, or length -- the recommended mode for
#     model_vs_obs snapshot comparisons (see the README).
#
# This file (T_005101_006012.nc) holds monthly data for years 0051-0060, so we
# select two months by date here.
param.time_slices = ["0051-01", "0051-07"]

# IMPORTANT: time_slices and seasons are mutually exclusive
# When using time_slices, do NOT set param.seasons

# Comparison settings
param.case_id = "model_vs_model"
param.run_type = "model_vs_model"

# Name of the folder where the results are stored.
# Change `prefix` to use your directory.
prefix = f"/lcrc/group/e3sm/public_html/diagnostic_output/{username}/e3sm_diags_examples"
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

#!/usr/bin/env python3
"""
Example 9 (model_vs_obs): Snapshot analysis against observations.

Companion to ``ex9.py`` (which runs model_vs_model). This example compares a
model against an observational time series at the same calendar time using
date-based ``time_slices`` instead of climatological seasons.

Why date-based time_slices for model_vs_obs:
    Each date selects the time step *nearest* to that date, applied
    independently to the test and reference datasets. This keeps the model and
    observations aligned on the same calendar time even when their time axes
    differ in start date, cadence, or length -- the recommended mode for
    snapshot comparisons.

This example also exercises two convenience features:
    - Directory mode: ``test_file`` / ``ref_file`` are NOT set, so the per-
      variable time series files are discovered under ``test_data_path`` /
      ``reference_data_path``.
    - Multi-file derivation: derived variables (e.g. ``PRECT``) are derived from
      their separate source-variable files found in the directory.

Supported diagnostic sets:
    lat_lon, lat_lon_native, zonal_mean_xy, zonal_mean_2d,
    zonal_mean_2d_stratosphere, polar, meridional_mean_2d
"""

import os

from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.run import runner

# Auto-detect username for the output directory.
username = os.environ.get("USER", "unknown_user")

param = CoreParameter()

# Directory mode: point at the directories that hold the per-variable monthly
# time series. Do NOT set test_file / ref_file -- the files are discovered by
# variable name (and, for derived variables, by their source-variable files).
param.test_data_path = (
    "/pscratch/sd/c/chengzhu/tests/"
    "ne256pg2_ne256pg2.F20TR-SCREAMv1.July-1.spanc800.2xauto.acc150.n0032.test2.1_05182026"
    "/post/atm/180x360_aave/ts/monthly/5yr"
)
param.reference_data_path = (
    "/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/time-series/"
)

param.short_test_name = "screamv1 test"

param.run_type = "model_vs_obs"

# Date-based snapshots (nearest monthly step). Each entry is a "YYYY-MM" or
# "YYYY-MM-DD" date; clear the default seasons so the run does not expand over
# seasons x time_slices.
#
# NOTE: time_slices and seasons are mutually exclusive. When time_slices is set
# it takes precedence, but clearing seasons keeps the intent explicit.
param.seasons = []
param.time_slices = ["2000-01"]

# Name of the folder where the results are stored.
# Change `prefix` to use your directory.
prefix = f"/global/cfs/cdirs/e3sm/www/{username}/examples"
param.results_dir = os.path.join(prefix, "ex9_model_obs_snapshot")

# All sets that support time_slices and run on this regridded lat-lon data.
# (lat_lon_native also supports time_slices but requires native-grid input, so
# it is not included here.)
runner.sets_to_run = [
    "lat_lon",
    "zonal_mean_xy",
    "zonal_mean_2d",
    "zonal_mean_2d_stratosphere",
    "polar",
    "meridional_mean_2d",
]
runner.run_diags([param])

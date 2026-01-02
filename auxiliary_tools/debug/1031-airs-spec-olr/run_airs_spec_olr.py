#!/usr/bin/env python
"""
Standalone run script for testing AIRS spectral OLR diagnostics.
Usage: python run_airs_spec_olr.py
"""

import os
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.run import runner

# =============================================================================
# USER CONFIGURATION - Update these paths
# =============================================================================

# Path to test model data (climatology files)
test_data_path = "/global/cfs/cdirs/e3sm/zhang40/analysis_data_e3sm_diags/AIRS_specOLR/processing/climo_2003-2021"

# Path to reference AIRS spectral OLR data
reference_data_path = "/global/cfs/cdirs/e3sm/zhang40/analysis_data_e3sm_diags/AIRS_specOLR/climatology"

# Results output directory
# Change prefix to use your directory
prefix = "/global/cfs/cdirs/e3sm/www/chengzhu/tests"
results_dir = os.path.join(prefix, "1031-airs-spec-olr")

test_name = "v3.LR.F20TR_spectral.pm-cpu"

# Short names for plot titles
short_test_name = "v3.LR.F20TR_spectral.pm-cpu"

# =============================================================================
# PARAMETER SETUP
# =============================================================================

param = CoreParameter()

# Data paths
param.test_data_path = test_data_path
param.reference_data_path = reference_data_path

param.short_test_name = short_test_name

# Results directory
param.results_dir = results_dir

# Run configuration
param.run_type = "model_vs_obs"
param.sets = ["lat_lon"]

# Use the custom cfg file for AIRS spectral OLR variables
param.cfg_path = os.path.join(
    os.path.dirname(__file__), "1031-airs-spec-olr.cfg"
)

# Optional: Enable multiprocessing for faster execution
# param.multiprocessing = True
# param.num_workers = 8

# Optional: Customize seasons (default is all seasons if not specified)
param.seasons = ["ANN"]

# =============================================================================
# RUN DIAGNOSTICS
# =============================================================================

runner.sets_to_run = ["lat_lon"]
runner.run_diags([param])

print(f"\nDiagnostics complete! Results saved to: {results_dir}")

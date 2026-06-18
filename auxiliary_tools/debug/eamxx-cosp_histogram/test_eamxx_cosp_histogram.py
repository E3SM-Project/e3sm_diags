#!/usr/bin/env python
"""
Test script for EAMxx COSP histogram support.

This script tests the new EAMxx COSP histogram variables:
- isccp_ctptau (replaces FISCCP1_COSP)
- modis_ctptau (replaces CLMODIS)
- misr_cthtau (replaces CLD_MISR)

It runs both the `cosp_histogram` set and the `lat_lon` set (cloud-fraction
maps derived from the COSP histograms via cosp_bin_sum), driven by the
`cosp_model_vs_obs.cfg` file in this directory.
"""

import os
import sys
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.run import runner

# Load the custom COSP diagnostics defined in cosp_model_vs_obs.cfg. CALIPSO
# variables are intentionally excluded because EAMxx output does not provide
# them.
cfg_path = os.path.join(os.path.dirname(__file__), "cosp_model_vs_obs.cfg")
sys.argv.extend(["--diags", cfg_path])

# Set up parameters
param = CoreParameter()

# Test case name
param.case_id = "EAMxx_cosp_histogram_test"
param.short_name = "EAMxx_COSP"

# Test data path - updated EAMxx output that already ships COSP coordinate
# values (cosp_tau, cosp_prs, cosp_cth), so no preprocessing step is needed.
# Original preprocessed-data path (preprocessing no longer required):
# test_data_path = "/pscratch/sd/c/chengzhu/EAMxx/ne256pg2_ne256pg2.F20TR-SCREAMv1.July-1.spanc800.2xauto.acc150.n0032.test2.1/rgr/climo"
test_data_path = "/pscratch/sd/y/yuying/e3sm_scratch/pm-gpu/ne30pg2_ne30pg2.F20TR-SCREAMv1.260525.branch25/run/climo"
param.test_data_path = test_data_path
# param.test_name = "1ma_ne30pg2.AVERAGE"
param.test_name = "eamxx_ne30"

# Reference data (use obs data for COSP)
param.reference_data_path = "/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/climatology"
param.results_dir = os.path.join(
    '/global/cfs/cdirs/e3sm/www/chengzhu/tests',
    'eamxx_cosp_histogram'
)


# Run lat_lon (COSP cloud-fraction maps) and cosp_histogram sets
param.sets = ["lat_lon", "cosp_histogram"]
#param.seasons = ["ANN"]

# Run the diagnostics
runner.sets_to_run = param.sets
runner.run_diags([param])


#!/usr/bin/env python
"""
Test script for EAMxx COSP histogram support.

This script tests the new EAMxx COSP histogram variables:
- isccp_ctptau (replaces FISCCP1_COSP)
- modis_ctptau (replaces CLMODIS)
- misr_cthtau (replaces CLD_MISR)
"""

import os
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.run import runner

# Set up parameters
param = CoreParameter()

# Test case name
param.case_id = "EAMxx_cosp_histogram_test"
param.short_name = "EAMxx_COSP"

# Test data path
test_data_path = "/pscratch/sd/t/terai/EAMxx/ne256pg2_ne256pg2.F20TR-SCREAMv1.July-1.spanc800.2xauto.acc150.n0032.test2.1/rgr/climo"
test_data_path = "/global/cfs/cdirs/e3sm/chengzhu/tests/zppy_example_v3/v3.LR.amip_0101/post/atm/180x360_aave/clim/10yr"
param.test_data_path = test_data_path
param.test_name = "1ma_ne30pg2.AVERAGE"
param.test_name = "v3.LR.amip_0101"

# Reference data (use obs data for COSP)
param.reference_data_path = "/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/climatology"
param.results_dir = os.path.join(
    '/global/cfs/cdirs/e3sm/www/chengzhu/tests',
    'eamxx_cosp_histogram'
)


# Run cosp_histogram set
param.sets = ["cosp_histogram"]
param.seasons = ["ANN"]

# Run the diagnostics
runner.sets_to_run = param.sets
runner.run_diags([param])


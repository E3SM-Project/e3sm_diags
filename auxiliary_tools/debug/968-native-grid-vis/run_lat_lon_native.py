#!/usr/bin/env python3
"""
This script runs e3sm_diags with the lat_lon_native set to visualize native grid data.
"""

import os
import sys

from e3sm_diags.parameter.lat_lon_native_parameter import LatLonNativeParameter
from e3sm_diags.run import runner

# Create parameter object
param = LatLonNativeParameter()

# Basic parameters
param.results_dir = "/lcrc/group/e3sm/public_html/diagnostic_output/ac.tvo/tests/lat_lon_native_file"
# Create results directory if it doesn't exist
if not os.path.exists(param.results_dir):
    os.makedirs(param.results_dir)


# Model data

##(1)
#param.test_data_path = "/home/ac.zhang40/test"
#param.test_name = "v3.LR.amip_0101"
#param.short_test_name = "v3.LR.amip_0101"
#param.reference_data_path = "/home/ac.zhang40/test"
#param.ref_name = "v3.HR.test4"
#param.short_ref_name = "v3.HR.test4"
#param.seasons = ["DJF"]


##(2)
#param.test_data_path = "/home/ac.zhang40/test"
#param.test_file = "v3.LR.amip_0101_DJF_climo.nc"
#param.short_test_name = "v3.LR.amip_0101"
#param.reference_data_path = "/home/ac.zhang40/test"
#param.ref_file = "v3.HR.test4_DJF_climo.nc"
#param.short_ref_name = "v3.HR.test4"
#param.seasons = ["DJF"]
##(3)
param.test_data_path = "/lcrc/group/e3sm2/ac.wlin/E3SMv3/v3.LR.historical_0051/archive/atm/hist"
param.test_file = "v3.LR.historical_0051.eam.h0.1989-12.nc"
#param.short_test_name = "v3.LR.amip_0101"
param.reference_data_path = "/lcrc/group/e3sm2/ac.wlin/E3SMv3/v3.LR.historical_0051/archive/atm/hist"
param.ref_file = "v3.LR.historical_0051.eam.h0.1989-12.nc"
#param.reference_data_path = "/lcrc/group/e3sm2/ac.jwolfe/E3SMv3_dev/20250404.wcycl1850.ne120pg2_r025_RRSwISC6to18E3r5.test4.chrysalis/archive/atm/hist/"
#param.ref_file = "20250404.wcycl1850.ne120pg2_r025_RRSwISC6to18E3r5.test4.chrysalis.eam.h0.0018-12.nc"
param.short_ref_name = "v3.HR.test4"
#param.seasons = ["DJF"]
param.time_slices=["0"]


param.case_id = "model_vs_model"


# Native grid settings
param.test_grid_file = "/lcrc/group/e3sm/diagnostics/grids/ne30pg2.nc"
param.ref_grid_file = "/lcrc/group/e3sm/diagnostics/grids/ne30pg2.nc"
#param.ref_grid_file = "/lcrc/group/e3sm/diagnostics/grids/ne120pg2.nc"

param.split_periodic_elements = True
param.antialiased = False

# param.model_only = True
param.run_type = "model_vs_model"

# Run the single diagnostic, comment out for complete diagnostics.
cfg_path = "auxiliary_tools/debug/968-native-grid-vis/run_lat_lon_native.cfg"

sys.argv.extend(["--diags", cfg_path])

runner.sets_to_run = ["lat_lon_native"]
runner.run_diags([param])



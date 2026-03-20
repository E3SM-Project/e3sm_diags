#!/usr/bin/env python
"""Test e3sm_diags GPCP PDF calculation"""
import sys
sys.path.insert(0, '/global/u2/c/chengzhu/e3sm_diags')

import numpy as np
from e3sm_diags.driver.precip_pdf_driver import calculate_precip_pdf, extract_regional_pdf
from e3sm_diags.parameter.precip_pdf_parameter import PrecipPDFParameter
from e3sm_diags.driver.utils.dataset_xr import Dataset

# Create parameter for GPCP reference data
param = PrecipPDFParameter()
param.reference_data_path = '/global/cfs/cdirs/e3sm/zhang40/analysis_data_e3sm_diags/GPCP_Daily_1DD/time_series/'
param.ref_name = 'GPCP_Daily_1DD'
param.ref_start_yr = '2001'
param.ref_end_yr = '2010'
param.ref_timeseries_input = True
param.run_type = 'model_vs_obs'

# Create dataset object for reference
ref_data = Dataset(param, data_type="ref")

# Get the dataset
print("Loading GPCP data from e3sm_diags reference path...")
print(f"Path: {param.reference_data_path}")
print(f"Years: {param.ref_start_yr}-{param.ref_end_yr}")

try:
    ref_ds = ref_data.get_time_series_dataset('PRECT', single_point=True)

    # Calculate PDF
    print("Calculating GPCP PDF...")
    ref_pdf, ref_start, ref_end = calculate_precip_pdf(ref_ds, 'PRECT')

    # Extract global mean
    print("Extracting global mean...")
    global_pdf = extract_regional_pdf(ref_pdf, 'global')

    # Get values
    bc = global_pdf['bin_centers'].values
    freq = global_pdf['FREQPDF'].values
    amnt = global_pdf['AMNTPDF'].values

    print('\n=== e3sm_diags GPCP Calculation ===')
    print(f'Time range: {ref_start}-{ref_end}')
    print(f'GPCP Freq max: {np.max(freq):.6f} at {bc[np.argmax(freq)]:.2f} mm/day')
    print(f'GPCP Amnt max: {np.max(amnt):.6f} at {bc[np.argmax(amnt)]:.2f} mm/day')

    print('\n=== Comparison with Terai GPCP ===')
    print('Terai GPCP Freq max: 0.254749')
    print('Terai GPCP Amnt max: 2.962387')
    print(f'Freq difference: {np.max(freq) - 0.254749:.6f}')
    print(f'Amnt difference: {np.max(amnt) - 2.962387:.6f}')

except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()

#!/usr/bin/env python
"""Test e3sm_diags PDF calculation against Terai's method"""
import sys
sys.path.insert(0, '/global/u2/c/chengzhu/e3sm_diags')

import numpy as np
import xarray as xr
from e3sm_diags.driver.precip_pdf_driver import calculate_precip_pdf, extract_regional_pdf
from e3sm_diags.parameter.precip_pdf_parameter import PrecipPDFParameter
from e3sm_diags.driver.utils.dataset_xr import Dataset

# Create parameter
param = PrecipPDFParameter()
param.test_data_path = '/global/cfs/cdirs/e3sm/chengzhu/tests/zppy_example_v3/v3.LR.amip_0101/post/atm/180x360_aave/ts/daily/10yr'
param.test_start_yr = '1995'
param.test_end_yr = '2004'
param.test_timeseries_input = True

# Create dataset object
test_data = Dataset(param, data_type="test")

# Get the dataset
print("Loading data...")
test_ds = test_data.get_time_series_dataset('PRECT', single_point=True)

# Calculate PDF
print("Calculating PDF...")
test_pdf, test_start, test_end = calculate_precip_pdf(test_ds, 'PRECT')

# Extract global mean
print("Extracting global mean...")
global_pdf = extract_regional_pdf(test_pdf, 'global')

# Get values
bc = global_pdf['bin_centers'].values
freq = global_pdf['FREQPDF'].values
amnt = global_pdf['AMNTPDF'].values

print('\n=== e3sm_diags Method ===')
print(f'Number of bins: {len(bc)}')
print(f'Bin range: {bc[0]:.3f} to {bc[-1]:.2f} mm/day')
print(f'Freq PDF max: {np.max(freq):.6f} at {bc[np.argmax(freq)]:.2f} mm/day')
print(f'Amnt PDF max: {np.max(amnt):.6f} at {bc[np.argmax(amnt)]:.2f} mm/day')
print(f'Freq PDF at 1 mm/day (bin ~20): {freq[20]:.6f}')
print(f'Amnt PDF at 10 mm/day (bin ~50): {amnt[50]:.6f}')

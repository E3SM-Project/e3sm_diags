"""
Run script for precipitation PDF diagnostics (model vs obs).

This script demonstrates how to run the precip_pdf diagnostic set
comparing model output against observational datasets (IMERG, GPCP).

With season_subset=True, this will generate PDFs for:
- ANN (all months)
- DJF (December, January, February)
- MAM (March, April, May)
- JJA (June, July, August)
- SON (September, October, November)
"""
import os

from e3sm_diags.parameter.precip_pdf_parameter import PrecipPDFParameter
from e3sm_diags.run import runner

# Create parameter object
param = PrecipPDFParameter()

# Set data paths
# Reference data: observational daily precipitation (IMERG, GPCP, etc.)
param.reference_data_path = '/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/time-series/'

# Test data: model daily time series output
param.test_data_path = '/global/cfs/cdirs/e3sm/chengzhu/tests/zppy_example_v3/v3.LR.amip_0101/post/atm/180x360_aave/ts/daily/10yr'

# Set test model name
param.test_name = 'E3SMv3.LR.amip_0101'

# Set results output directory
param.results_dir = os.path.join(
    '/global/cfs/cdirs/e3sm/www/chengzhu/tests',
    'precip_pdf_test_multi_ref'
)

# Set time ranges
param.test_start_yr = '1995'
param.test_end_yr = '2004'
param.ref_start_yr = '2001'
param.ref_end_yr = '2010'

# Save global PDF netCDF files for offline use
param.save_netcdf = True

# Enable seasonal subsetting (default is False)
# When True: generates PDFs for ANN, DJF, MAM, JJA, SON
# When False: generates PDF for ANN only
param.season_subset = False

# Run the diagnostics
runner.sets_to_run = ['precip_pdf']
runner.run_diags([param])

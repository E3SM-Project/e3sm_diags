# Example 10: Precipitation PDF Diagnostics

This example demonstrates how to run precipitation PDF (Probability Density Function) diagnostics comparing model output against observational datasets.

## What it does

The `precip_pdf` diagnostic set:
- Calculates frequency and amount PDFs of daily precipitation
- Compares model output against observational datasets (GPCP, IMERG)
- Analyzes multiple regions: TROPICS, CONUS, and global
- Supports seasonal subsetting (ANN, DJF, MAM, JJA, SON)
- Generates interactive HTML viewer with detailed plots

## Files

- `ex10.py`: Main run script with parameter configuration
- `diags.cfg`: Configuration file specifying regions and reference datasets

## Key Parameters

### Data Paths
- `reference_data_path`: Location of observational daily precipitation data
- `test_data_path`: Location of model daily time series output

### Time Ranges
- `test_start_yr`, `test_end_yr`: Years to analyze from model data
- `ref_start_yr`, `ref_end_yr`: Years to analyze from observational data

### Analysis Options
- `season_subset`: When `True`, analyzes all seasons (DJF, MAM, JJA, SON) plus annual
- `save_netcdf`: When `True`, saves calculated PDFs to NetCDF for later use
- `regions`: List of regions to analyze ("global", "TROPICS", "CONUS")
- `ref_name`: Reference dataset(s) - can be single string or list

## Running the Example

### Method 1: Using Python Script

```bash
# Allocate a node (or use batch job)
salloc --nodes 1 --qos interactive --time 01:00:00 --constraint cpu --account=e3sm

# Load E3SM Unified environment
source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh

# Run the diagnostics
cd examples/ex10-precip-pdf
python ex10.py -d diags.cfg --multiprocessing --num_workers=32

# Adjust permissions to view results
chmod -R 755 <your web directory>
```

### Method 2: Using Command Line

```bash
e3sm_diags precip_pdf --no_viewer \
  --reference_data_path '/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/time-series/' \
  --test_data_path '/global/cfs/cdirs/e3sm/chengzhu/tests/zppy_example_v3/v3.LR.amip_0101/post/atm/180x360_aave/ts/daily/10yr' \
  --results_dir '/global/cfs/cdirs/e3sm/www/chengzhu/tests/precip_pdf_test_seasons_multi_ref' \
  --case_id 'precipitation-pdf' \
  --sets 'precip_pdf' \
  --variables 'PRECT' \
  --regions 'TROPICS' \
  --save_netcdf \
  --test_name 'E3SMv3.LR.amip_0101' \
  --ref_name 'GPCP_1DD_Daily' 'IMERG_Daily' \
  --season_subset \
  --test_start_yr 1995 \
  --test_end_yr 2004 \
  --ref_start_yr 2001 \
  --ref_end_yr 2010
```

## Output

The diagnostics generate:
- **Plots**: PNG files showing frequency and amount PDFs for each region/season
- **NetCDF files**: Cached global gridded PDFs (if `save_netcdf=True`)
- **HTML Viewer**: Interactive viewer at `<results_dir>/viewer/precip_pdf/index.html`

### Viewer Structure
- Single "Precipitation PDF" group containing all results
- Each row represents a unique region + reference dataset combination
- Detail pages for each season showing:
  - Frequency PDF (df/dlog(P))
  - Amount PDF (dA/dlog(P))
  - Links to NetCDF data files
  - Metadata and command to recreate

## Supported Reference Datasets

- **GPCP_1DD_Daily**: GPCP 1-Degree Daily v1.3 (1996-2017)
- **IMERG_Daily**: GPM IMERG Daily (2001-2021)

## Supported Regions

- **global**: Global mean
- **TROPICS**: 30°S - 30°N
- **CONUS**: Continental United States

## Notes

- Daily precipitation data is required (not monthly climatologies)
- Model variable should be `PRECT` (total precipitation rate)
- Reference data is automatically subset to match test data time range
- PDFs are calculated on exponentially-spaced bins (0.1 to ~900 mm/day)
- The viewer fix (v3.1.0+) ensures that all regions generate distinct detail HTML pages

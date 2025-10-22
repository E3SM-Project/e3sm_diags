# Example 9: Snapshot Analysis for Core Sets

This example demonstrates the **snapshot analysis** feature introduced in **E3SM Diags v3.1.0**.

## What This Example Does

- Analyzes individual time steps instead of seasonal climatological means
- Uses index-based time selection to examine specific time points
- Demonstrates time_slices parameter on multiple core diagnostic sets
- Compares model states at specific indices without temporal averaging

## Key Features

The snapshot analysis capability:
- Enables event-based or process-oriented diagnostics
- Analyzes specific time points without climatological averaging
- Supports multiple time indices analyzed separately
- Works across multiple core diagnostic sets

## Key Parameters

- `time_slices` - List of time indices to analyze (e.g., ["0", "1", "2"])
  - Time slices are zero-based indices into the time dimension
  - ["0"] = first time step
  - ["5"] = 6th time step
  - ["0", "1", "2"] = first 3 time steps (each analyzed separately)
- **IMPORTANT**: `time_slices` and `seasons` are mutually exclusive
  - When using `time_slices`, do NOT set `seasons`

## Supported Diagnostic Sets

The following core diagnostic sets support snapshot analysis:
- `lat_lon` - Latitude-Longitude contour maps
- `lat_lon_native` - Native grid visualization
- `polar` - Polar contour maps
- `zonal_mean_2d` - Pressure-Latitude zonal mean contour plots
- `meridional_mean_2d` - Pressure-Longitude meridional mean contour plots
- `zonal_mean_2d_stratosphere` - Stratospheric zonal mean plots

## Data Requirements

This example uses test data located at LCRC:
- Data path: `/lcrc/group/e3sm/public_html/e3sm_diags_test_data/postprocessed_e3sm_v2_data_for_e3sm_diags/`

For your own data, ensure:
1. Model output files contain time dimension
2. Files have sufficient time steps for requested indices
3. Data files are accessible from the diagnostic runs

## Running This Example

### Using the Python Script

```bash
# Edit ex9.py to set your output directory
# Update the `prefix` variable to point to your web directory

# Run with default settings
python ex9.py

# Run with custom configuration file
python ex9.py -d diags.cfg 
```

### Using Command-Line Interface (example)

```bash
e3sm_diags zonal_mean_2d \
  --no_viewer \
  --reference_data_path '/lcrc/group/e3sm/public_html/e3sm_diags_test_data/postprocessed_e3sm_v2_data_for_e3sm_diags/20210528.v2rc3e.piControl.ne30pg2_EC30to60E2r2.chrysalis/time-series/rgr' \
  --test_data_path '/lcrc/group/e3sm/public_html/e3sm_diags_test_data/postprocessed_e3sm_v2_data_for_e3sm_diags/20210528.v2rc3e.piControl.ne30pg2_EC30to60E2r2.chrysalis/time-series/rgr' \
  --results_dir '/lcrc/group/e3sm/public_html/diagnostic_output/$USER/e3sm_diags_examples/ex9_snapshot_analysis' \
  --case_id 'model_vs_model' \
  --run_type 'model_vs_model' \
  --sets 'zonal_mean_2d' \
  --variables 'T' \
  --time_slices 0 \
  --multiprocessing \
  --main_title 'T 0 global' \
  --contour_levels '180' '185' '190' '200' '210' '220' '230' '240' '250' '260' '270' '280' '290' '295' '300' \
  --short_test_name 'v2 test' \
  --ref_file 'T_005101_006012.nc' \
  --diff_levels '-3.0' '-2.5' '-2' '-1.5' '-1' '-0.5' '-0.25' '0.25' '0.5' '1' '1.5' '2' '2.5' '3.0' \
  --test_file 'T_005101_006012.nc'
```

**Note:** Use `--no_viewer` for command-line usage to avoid directory creation issues. For HTML viewer output, use the Python script approach instead.

**Important**: Do not use both `--time_slices` and `--seasons` in the same command!

## Configuration File

The `diags.cfg` file allows you to customize settings for each diagnostic set:
- Variables to plot (e.g., T for temperature)
- Pressure levels for 3D variables (e.g., 850.0 mb)
- Regions of interest (e.g., polar_S, polar_N)
- Colormap settings
- Contour levels for test, reference, and difference plots
- Regridding method (for lat_lon set)

## Expected Output

The diagnostic will generate:
- Plots for each time slice specified
- Test model plot at each time index
- Reference model plot at each time index
- Difference plots (Test - Reference) at each time index
- HTML viewer with columns for each time slice

Results will be saved in: `<your_directory>/ex9_snapshot_analysis/viewer/`

In the viewer:
- Rows represent different variables/regions/levels
- Columns represent different time slices (0, 1, 2, etc.)
- Click on any cell to see detailed plots

## Notes

- Time slices are **zero-based indices** (0 = first time step, 1 = second, etc.)
- Each time slice is analyzed **separately** (not averaged together)
- The viewer displays time slices as column headers instead of seasons
- Time slices are sorted **numerically** (0, 1, 2, ..., not alphabetically)
- Make sure your data files have enough time steps for the requested indices

## Differences from Seasonal Climatology

Unlike seasonal climatology analysis which uses `seasons = ["ANN", "DJF", "JJA", "SON"]`:
- Snapshot analysis uses `time_slices = ["0", "1", "2", ...]`
- No temporal averaging - analyzes exact time points
- Useful for event-based studies and temporal evolution
- Can analyze any arbitrary time index in your dataset

## Use Cases

Snapshot analysis is particularly useful for:
1. **Event Studies**: Analyzing specific weather events or phenomena
2. **Model Spin-up**: Examining early time steps in model initialization
3. **Temporal Evolution**: Tracking how fields change over successive time steps
4. **Intercomparison**: Comparing models at synchronized time points
5. **Debugging**: Investigating specific time steps with unusual behavior

## Combining with Native Grid Visualization

You can combine snapshot analysis with native grid visualization:

```python
from e3sm_diags.parameter.lat_lon_native_parameter import LatLonNativeParameter

param = LatLonNativeParameter()
param.time_slices = ["0", "5", "10"]  # Snapshot analysis
param.test_grid_file = "/path/to/grid.nc"  # Native grid
# ... other parameters ...

runner.sets_to_run = ["lat_lon_native"]
runner.run_diags([param])
```

This combines both v3.1.0 features for snapshot analysis on native grids!

## More Information

For more details, see:
- [E3SM Diags Documentation](https://e3sm-project.github.io/e3sm_diags)
- [E3SM Diags README - v3.1.0 Features](https://github.com/E3SM-Project/e3sm_diags#new-features-in-v310)

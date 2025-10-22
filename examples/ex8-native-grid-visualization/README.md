# Example 8: Native Grid Visualization

This example demonstrates the **native grid visualization** feature introduced in **E3SM Diags v3.1.0**.

## What This Example Does

- Visualizes model data on native grids (e.g., cubed-sphere, unstructured grids)
- Uses UXarray for grid-aware operations
- Compares two models without regridding to a regular lat-lon grid
- Preserves native grid features and characteristics

## Key Features

The native grid visualization capability:

- Supports various native grid formats (cubed-sphere, unstructured, etc.)
- Eliminates artifacts introduced by regridding
- Enables comparison of models with different native grids
- Particularly useful for high-resolution models

## Key Parameters

- `LatLonNativeParameter` - Required parameter class for native grid visualization
- `test_grid_file` - Path to test model's grid file (UGRID format)
- `ref_grid_file` - Path to reference model's grid file (optional for model-only runs)
- `time_slices` - Use snapshot analysis instead of climatology (e.g., ["0"])
- `antialiased` - Whether to apply antialiasing to the plot

## Data Requirements

This example uses test data located at LCRC:

- Data path: `/lcrc/group/e3sm/public_html/e3sm_diags_test_data/native_grid`
- Grid files: `/lcrc/group/e3sm/diagnostics/grids/`

For your own data, ensure you have:

1. Model output files on native grid
2. Corresponding grid files in UGRID format

## Running This Example

### Using the Python Script

```bash
# Run with default settings (automatically uses your username for output directory)
python ex8.py

# Run with custom configuration file
python ex8.py -d diags.cfg 
```

### Using Command-Line Interface

```bash
e3sm_diags lat_lon_native \
  --no_viewer \
  --reference_data_path '/lcrc/group/e3sm/public_html/e3sm_diags_test_data/native_grid' \
  --test_data_path '/lcrc/group/e3sm/public_html/e3sm_diags_test_data/native_grid' \
  --results_dir '/lcrc/group/e3sm/public_html/diagnostic_output/$USER/e3sm_diags_examples/ex8_native_grid' \
  --case_id 'model_vs_model' \
  --run_type 'model_vs_model' \
  --sets 'lat_lon_native' \
  --variables 'TGCLDLWP' \
  --time_slices 0 \
  --main_title 'TGCLDLWP 0 global' \
  --contour_levels '10' '25' '50' '75' '100' '125' '150' '175' '200' '225' '250' \
  --short_test_name 'v3.LR.amip_0101' \
  --ref_file 'v3.LR.amip_0101.eam.h0.1989-12.nc' \
  --diff_colormap 'RdBu' \
  --diff_levels '-35' '-30' '-25' '-20' '-15' '-10' '-5' '5' '10' '15' '20' '25' '30' '35' \
  --test_grid_file '/lcrc/group/e3sm/diagnostics/grids/ne30pg2.nc' \
  --ref_grid_file '/lcrc/group/e3sm/diagnostics/grids/ne30pg2.nc' \
  --test_file 'v3.LR.amip_0101.eam.h0.1989-12.nc'
```

**Note:** Use `--no_viewer` for command-line usage to avoid directory creation issues. For HTML viewer output, use the Python script approach instead.

## Configuration File

The `diags.cfg` file allows you to customize:

- Variables to plot (e.g., TGCLDLWP)
- Regions of interest
- Colormap settings
- Contour levels for test, reference, and difference plots

## Expected Output

The diagnostic will generate:

- Native grid visualizations for specified variables
- Test model plot
- Reference model plot
- Difference plot (Test - Reference)
- HTML viewer for browsing results

Results will be saved in: `/lcrc/group/e3sm/public_html/diagnostic_output/$USER/e3sm_diags_examples/ex8_native_grid/viewer/`

## Notes

- Native grid visualization requires the UXarray library, which is included as a dependency of E3SM diagnostics and the E3SM Unified environment.
- Grid files must be in UGRID format
- This example uses `time_slices` for snapshot analysis; you can also use `seasons` for climatology
- For model-only runs (no reference data), set `model_only = True` and omit ref_grid_file

## Differences from Regular lat_lon Set

Unlike the standard `lat_lon` set which regrids data to a regular lat-lon grid:

- `lat_lon_native` preserves the original grid structure
- No interpolation artifacts
- Better representation of native grid features
- Requires grid files in UGRID format

## More Information

For more details, see:

- [E3SM Diags Documentation](https://e3sm-project.github.io/e3sm_diags)
- [UXarray Documentation](https://uxarray.readthedocs.io/)

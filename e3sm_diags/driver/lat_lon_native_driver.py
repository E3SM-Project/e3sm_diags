from __future__ import annotations

import os
from typing import TYPE_CHECKING, List, Tuple

import matplotlib.colors as mcolors
import numpy as np
import uxarray as ux
import xarray as xr

from e3sm_diags.driver.utils.climo_xr import ClimoFreq
from e3sm_diags.driver.utils.dataset_xr import Dataset
from e3sm_diags.driver.utils.regrid import (
    _apply_land_sea_mask,
    _subset_on_region,
    get_z_axis,
    has_z_axis,
    regrid_z_axis_to_plevs,
)
from e3sm_diags.driver.utils.type_annotations import MetricsDict
from e3sm_diags.logger import _setup_child_logger
from e3sm_diags.metrics.metrics import spatial_avg, std

logger = _setup_child_logger(__name__)

if TYPE_CHECKING:
    from e3sm_diags.parameter.lat_lon_native_parameter import LatLonNativeParameter

# The default value for metrics if it is not calculated. This value was
# preserved from the legacy CDAT codebase because the viewer expects this
# value for metrics that aren't calculated.
# TODO: Update `lat_lon_viewer.py` to handle missing metrics with None value.
METRICS_DEFAULT_VALUE = 999.999


def run_diag(parameter: LatLonNativeParameter) -> LatLonNativeParameter:
    """Get metrics for the lat_lon_native diagnostic set.

    This function loops over each variable, season, pressure level (if 3-D),
    and region.

    Parameters
    ----------
    parameter : LatLonNativeParameter
        The parameter for the diagnostic.

    Returns
    -------
    LatLonNativeParameter
        The parameter for the diagnostic with the result (completed or failed).

    Raises
    ------
    RuntimeError
        If the dimensions of the test and reference datasets are not aligned
        (e.g., one is 2-D and the other is 3-D).
    """
    variables = parameter.variables
    seasons = parameter.seasons
    ref_name = getattr(parameter, "ref_name", "")
    regions = parameter.regions

    # Variables storing xarray `Dataset` objects start with `ds_` and
    # variables storing e3sm_diags `Dataset` objects end with `_ds`. This
    # is to help distinguish both objects from each other.
    test_ds = Dataset(parameter, data_type="test")
    ref_ds = Dataset(parameter, data_type="ref")

    for var_key in variables:
        logger.info("Variable: {}".format(var_key))
        parameter.var_id = var_key

        for season in seasons:
            parameter._set_name_yrs_attrs(test_ds, ref_ds, season)

            ds_test = test_ds.get_climo_dataset(var_key, season)
            ds_ref = _get_ref_climo_dataset(ref_ds, var_key, season)
            ds_land_sea_mask: xr.Dataset = test_ds._get_land_sea_mask(season)

            # Debug information about the test dataset
            logger.info("Test dataset info:")
            logger.info(f"Variables: {list(ds_test.variables)}")
            if hasattr(ds_test, 'file_path'):
                logger.info(f"Dataset file path: {ds_test.file_path}")
            if hasattr(ds_test, 'filepath'):
                logger.info(f"Dataset filepath: {ds_test.filepath}")
            if hasattr(ds_test, '_file_obj') and hasattr(ds_test._file_obj, 'name'):
                logger.info(f"Dataset file object name: {ds_test._file_obj.name}")
            try:
                for name, var in ds_test.variables.items():
                    if hasattr(var, 'file'):
                        logger.info(f"Variable {name} file: {var.file}")
                        break
            except:
                pass

            # Helper function for variable derivation and transformation
            from e3sm_diags.derivations.derivations import (
                DERIVED_VARIABLES,
                FUNC_NEEDS_TARGET_VAR,
            )

            def process_variable_derivations(dataset, variable_key, dataset_name=''):
                """
                Apply variable derivations and transformations to a dataset.

                Parameters
                ----------
                dataset : ux.dataset.UxDataset
                    The dataset to process
                variable_key : str
                    The target variable key
                dataset_name : str, optional
                    Name of the dataset for logging purposes (e.g., 'test', 'reference')

                Returns
                -------
                bool
                    True if the variable exists in the dataset after processing
                """
                name_suffix = f" in {dataset_name} dataset" if dataset_name else ""

                # Process the variable whether it exists or not
                if variable_key in dataset:
                    logger.info(f"Found variable {variable_key}{name_suffix}")
                    # For existing variables, check if we need to apply unit conversion or other transformations
                    if variable_key in DERIVED_VARIABLES:
                        # Look for a direct transformation of this variable
                        for vars_tuple, func in DERIVED_VARIABLES[variable_key].items():
                            # If there's a direct conversion for this variable (e.g., unit conversion)
                            if len(vars_tuple) == 1 and vars_tuple[0] == variable_key:
                                try:
                                    # Apply the conversion function
                                    original_units = dataset[variable_key].attrs.get('units', 'unknown')
                                    result = func(dataset[variable_key])
                                    if isinstance(result, xr.DataArray):
                                        dataset[variable_key] = result
                                        new_units = dataset[variable_key].attrs.get('units', 'unknown')
                                        logger.info(f"Applied conversion to {variable_key}{name_suffix}: {original_units} -> {new_units}")
                                        break
                                except Exception as e:
                                    logger.warning(f"Failed to apply conversion to {variable_key}{name_suffix}: {e}")
                else:
                    # If not found directly, attempt to derive it
                    logger.info(f"Variable {variable_key} not found directly{name_suffix}, attempting derivation")
                    if variable_key in DERIVED_VARIABLES:
                        # Try each derivation method for this variable
                        derived = False
                        for vars_tuple, func in DERIVED_VARIABLES[variable_key].items():
                            # Skip direct conversions as the variable doesn't exist
                            if len(vars_tuple) == 1 and vars_tuple[0] == variable_key:
                                continue

                            # Check if all source variables are available in the dataset
                            if all(v in dataset for v in vars_tuple):
                                try:
                                    # Apply derivation function
                                    if func in FUNC_NEEDS_TARGET_VAR:
                                        result = func(*[dataset[v] for v in vars_tuple], variable_key)
                                    else:
                                        result = func(*[dataset[v] for v in vars_tuple])

                                    # Add the derived variable to the dataset
                                    if isinstance(result, xr.DataArray):
                                        dataset[variable_key] = result
                                        logger.info(f"Successfully derived {variable_key} from {vars_tuple}{name_suffix}")
                                        derived = True
                                        break
                                except Exception as e:
                                    logger.warning(f"Failed to derive {variable_key} from {vars_tuple}{name_suffix}: {e}")

                        if not derived:
                            logger.warning(f"Could not derive {variable_key}{name_suffix} - required source variables not available")
                    else:
                        logger.warning(f"{variable_key} is not a recognized derivable variable")

                # Verify if variable exists and log possible matches if not
                if variable_key not in dataset:
                    logger.warning(f"Variable {variable_key} not found{name_suffix}!")
                    logger.info(f"Available variables{name_suffix}: {list(dataset.data_vars)}")

                    # Try to find the variable with a different name or naming convention
                    possible_matches = []
                    for data_var in dataset.data_vars:
                        if variable_key.lower() in data_var.lower() or data_var.lower() in variable_key.lower():
                            possible_matches.append(data_var)

                    if possible_matches:
                        logger.info(f"Possible variable matches{name_suffix}: {possible_matches}")

                    return False
                return True

            # Load the native grid information for test data
            uxds_test = None
            if parameter.test_grid_file:
                try:
                    logger.info(f"Loading test native grid from: {parameter.test_grid_file}")
                    # When loading the dataset, include the test data to map it onto the grid
                    uxds_test = ux.open_dataset(parameter.test_grid_file, parameter.test_data_file_path)
                    logger.info("Successfully loaded test native grid data with uxarray")

                    # Process variable derivations for test dataset
                    process_variable_derivations(uxds_test, var_key, 'test')

                except Exception as e:
                    logger.error(f"Failed to load test native grid: {e}")
                    import traceback
                    logger.error(traceback.format_exc())
                    uxds_test = None

            # Load the native grid information for reference data
            uxds_ref = None
            if ds_ref is not None:
                logger.info("=============== REFERENCE DATA LOADING DEBUG INFO ===============")
                logger.info(f"ds_ref exists: {ds_ref is not None}")
                logger.info(f"ds_ref var_key exists: {var_key in ds_ref}")
                if var_key in ds_ref:
                    logger.info(f"ds_ref[{var_key}] shape: {ds_ref[var_key].shape}")
                logger.info(f"ref_grid_file attribute exists: {hasattr(parameter, 'ref_grid_file')}")
                if hasattr(parameter, 'ref_grid_file'):
                    logger.info(f"ref_grid_file value: {parameter.ref_grid_file}")
                    logger.info(f"ref_grid_file exists: {os.path.exists(parameter.ref_grid_file) if parameter.ref_grid_file else False}")

                if not hasattr(parameter, 'ref_grid_file') or parameter.ref_grid_file is None:
                    logger.warning("No ref_grid_file specified in parameter. This is required for model_vs_model native grid visualization.")
                    logger.warning("Make sure your parameter configuration includes 'ref_grid_file'")
                    # Try to use test_grid_file as a fallback for models with the same grid
                    if hasattr(parameter, 'test_grid_file') and parameter.test_grid_file:
                        logger.warning(f"Attempting to use test_grid_file as fallback for reference: {parameter.test_grid_file}")
                        try:
                            # Attempt to use ref_data_file_path if available
                            if hasattr(parameter, 'ref_data_file_path') and parameter.ref_data_file_path:
                                logger.info(f"Using ref_data_file_path with test grid: {parameter.ref_data_file_path}")
                                uxds_ref = ux.open_dataset(parameter.test_grid_file, parameter.ref_data_file_path)
                            else:
                                # Otherwise fall back to ds_ref
                                uxds_ref = ux.open_dataset(parameter.test_grid_file, ds_ref)

                            logger.info("Successfully loaded reference data using test grid file as fallback")
                            logger.info(f"uxds_ref is now {type(uxds_ref)}")
                            logger.info(f"uxds_ref variables: {list(uxds_ref.data_vars)}")
                            process_variable_derivations(uxds_ref, var_key, 'reference')
                        except Exception as e:
                            logger.error(f"Failed to load reference using test grid as fallback: {e}")
                            import traceback
                            logger.error(traceback.format_exc())
                            uxds_ref = None
                elif parameter.ref_grid_file:
                    try:
                        logger.info(f"Loading reference native grid from: {parameter.ref_grid_file}")

                        # Try to use ref_data_file_path if available
                        if hasattr(parameter, 'ref_data_file_path') and parameter.ref_data_file_path:
                            logger.info(f"Using ref_data_file_path: {parameter.ref_data_file_path}")
                            try:
                                uxds_ref = ux.open_dataset(parameter.ref_grid_file, parameter.ref_data_file_path)
                                logger.info("Successfully loaded reference grid with ref_data_file_path")
                            except Exception as e:
                                logger.warning(f"Failed to load reference with ref_data_file_path: {e}")
                                # Fall back to using ds_ref
                                uxds_ref = ux.open_dataset(parameter.ref_grid_file, ds_ref)
                                logger.info("Fallback: Successfully loaded reference using ds_ref")
                        else:
                            # When loading the dataset, include the reference data to map it onto the grid
                            uxds_ref = ux.open_dataset(parameter.ref_grid_file, ds_ref)
                            logger.info("Successfully loaded reference native grid data with uxarray")
                        logger.info(f"uxds_ref is now {type(uxds_ref)}")

                        # Log the variables in the loaded reference dataset
                        logger.info(f"Reference dataset variables: {list(uxds_ref.data_vars)}")
                        if var_key in uxds_ref:
                            logger.info(f"uxds_ref[{var_key}] already exists")
                        else:
                            logger.info(f"uxds_ref does not contain {var_key} yet, will try derivation")

                        # Process variable derivations for reference dataset
                        if process_variable_derivations(uxds_ref, var_key, 'reference'):
                            logger.info(f"Successfully derived/found {var_key} in reference dataset")
                            logger.info(f"uxds_ref[{var_key}] shape: {uxds_ref[var_key].shape}")
                            logger.info(f"uxds_ref[{var_key}] data type: {uxds_ref[var_key].dtype}")
                            if hasattr(uxds_ref[var_key], 'units'):
                                logger.info(f"uxds_ref[{var_key}] units: {uxds_ref[var_key].units}")
                        else:
                            logger.warning(f"Unable to find or derive {var_key} in reference dataset")

                    except Exception as e:
                        logger.error(f"Failed to load reference native grid: {e}")
                        import traceback
                        logger.error(traceback.format_exc())
                        uxds_ref = None

                logger.info("=============== END REFERENCE DATA DEBUG INFO ===============")

            if ds_ref is None:
                is_vars_3d = has_z_axis(ds_test[var_key])

                # Only proceed with native grid diagnostics if uxds_test is available
                if not is_vars_3d and uxds_test is not None:
                    _run_diags_2d_model_only(
                        parameter,
                        season,
                        regions,
                        var_key,
                        ref_name,
                        uxds_test,
                    )
                elif not is_vars_3d:
                    logger.warning("Skipping native grid diagnostics: uxds_test is None")
                else:
                    _run_diags_3d_model_only(
                        parameter,
                        ds_test,
                        ds_land_sea_mask,
                        season,
                        regions,
                        var_key,
                        ref_name,
                        uxds_test,  # Use the consistent naming
                    )
            else:
                is_vars_3d, is_dims_diff = _check_var_dims(ds_test, ds_ref, var_key)

                if is_dims_diff:
                    raise RuntimeError(
                        "Dimensions of the two variables are different. Aborting."
                    )
                elif not is_vars_3d:
                    _run_diags_2d(
                        parameter,
                        ds_test,
                        ds_ref,
                        ds_land_sea_mask,
                        season,
                        regions,
                        var_key,
                        ref_name,
                        uxds_test,  # Use the consistent naming
                        uxds_ref,  # Use the consistent naming
                    )
                elif is_vars_3d:
                    _run_diags_3d(
                        parameter,
                        ds_test,
                        ds_ref,
                        ds_land_sea_mask,
                        season,
                        regions,
                        var_key,
                        ref_name,
                        uxds_test,  # Use the consistent naming
                        uxds_ref,  # Use the consistent naming
                    )

    return parameter


def _save_native_plot(
    parameter: LatLonNativeParameter,
    var_key: str,
    region: str,
    uxds_test: ux.dataset.UxDataset = None,
    uxds_ref: ux.dataset.UxDataset = None,
    uxds_diff: ux.dataset.UxDataset = None,
):
    """Save plots using native grid datasets directly.

    This function creates visualizations of data on native (unstructured) grids using uxarray,
    without regridding to a regular lat-lon grid. The layout matches the standard lat_lon_plot
    with 3 panels (test, ref, diff). This function explicitly uses native grid dataset from uxarray.

    Parameters
    ----------
    parameter : LatLonNativeParameter
        The parameter object.
    var_key : str
        The variable key.
    region : str
        The region name, used to determine map extents.
    uxds_test : ux.dataset.UxDataset, optional
        The test native grid dataset.
    uxds_ref : ux.dataset.UxDataset, optional
        The reference native grid dataset.
    uxds_diff : ux.dataset.UxDataset, optional
        The difference native grid dataset.
    """
    # Debug information for investigating reference data issues
    logger.info("=== SAVE_NATIVE_PLOT DEBUG INFO ===")
    logger.info(f"Function call with var_key={var_key}, region={region}")

    # Check test dataset
    logger.info(f"uxds_test is None: {uxds_test is None}")
    if uxds_test is not None:
        logger.info(f"uxds_test type: {type(uxds_test)}")
        logger.info(f"var_key in uxds_test: {var_key in uxds_test}")

    # Check reference dataset - the key issue we're trying to debug
    logger.info(f"uxds_ref is None: {uxds_ref is None}")
    if uxds_ref is not None:
        logger.info(f"uxds_ref type: {type(uxds_ref)}")
        logger.info(f"var_key in uxds_ref: {var_key in uxds_ref}")
        if var_key in uxds_ref:
            logger.info(f"uxds_ref[{var_key}] shape: {uxds_ref[var_key].shape}")
    logger.info("=== END SAVE_NATIVE_PLOT DEBUG INFO ===")


    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import matplotlib.pyplot as plt
    import numpy as np

    # Import required modules
    from e3sm_diags.derivations.default_regions_xr import REGION_SPECS
    from e3sm_diags.driver.utils.io import _get_output_dir
    from e3sm_diags.plot.utils import DEFAULT_PANEL_CFG, _get_colormap

    # Check if uxarray dataset and variable are available for plotting
    if uxds_test is None or var_key not in uxds_test:
        logger.error(f"Cannot plot native grid data. Either uxds_test is None or {var_key} not in dataset")
        if uxds_test is not None:
            logger.error(f"Available variables: {list(uxds_test.data_vars)}")
        return

    # Check if reference data is available
    has_reference = uxds_ref is not None and var_key in uxds_ref
    if has_reference:
        logger.info(f"Reference data available for {var_key}, implementing model vs model visualization")

    # Set the viewer description to the "long_name" attr of the variable
    if 'long_name' in uxds_test[var_key].attrs:
        parameter.viewer_descr[var_key] = uxds_test[var_key].attrs['long_name']
    else:
        parameter.viewer_descr[var_key] = var_key

    # Create figure with standard layout
    fig = plt.figure(figsize=parameter.figsize, dpi=parameter.dpi)
    # Use the same title formatting as in lat_lon_plot (fontsize=18, y=0.96)
    fig.suptitle(parameter.main_title, x=0.5, y=0.96, fontsize=18)

    # Get region information for setting map extents
    region_specs = REGION_SPECS.get(region, None)

    # Set map bounds based on region
    if region_specs:
        lat_bounds = region_specs.get('lat', (-90, 90))
        lon_bounds = region_specs.get('lon', (0, 360))
        is_global_domain = lat_bounds == (-90, 90) and lon_bounds == (0, 360)
    else:
        lat_bounds = (-90, 90)
        lon_bounds = (0, 360)
        is_global_domain = True

    # Determine projection and extents based on region
    projection = ccrs.PlateCarree()
    if is_global_domain:
        projection = ccrs.PlateCarree(central_longitude=180)

    logger.info(f"Region: {region}, lat_bounds: {lat_bounds}, lon_bounds: {lon_bounds}")

    # Determine X and Y ticks
    lat_south, lat_north = lat_bounds
    lon_west, lon_east = lon_bounds

    # Extract metrics directly from the uxarray dataset
    if uxds_test is not None and var_key in uxds_test:
        logger.info(f"Extracting metrics from uxds_test for variable {var_key}")
        test_min = uxds_test[var_key].min().item()
        test_max = uxds_test[var_key].max().item()
        test_mean = uxds_test[var_key].mean().item()  # For native grid, use simple mean as approximation
        units = uxds_test[var_key].attrs.get('units', '')
        logger.info(f"Test data metrics - min: {test_min}, mean: {test_mean}, max: {test_max}, units: {units}")

        # Check for data issues that might affect visualization
        n_unique_values = len(np.unique(uxds_test[var_key].values))
        if n_unique_values <= 1:
            logger.warning(f"⚠️ Variable {var_key} has only {n_unique_values} unique value(s)!")
    else:
        # This should not happen since we check earlier, but just in case
        logger.error(f"Missing test data for variable {var_key} in native grid dataset")
        return

    # Extract metrics for reference data if available
    ref_min = ref_max = ref_mean = diff_min = diff_max = diff_mean = None
    if has_reference:
        ref_min = uxds_ref[var_key].min().item()
        ref_max = uxds_ref[var_key].max().item()
        ref_mean = uxds_ref[var_key].mean().item()
        ref_units = uxds_ref[var_key].attrs.get('units', '')

        # Check if units match between test and reference
        if units != ref_units:
            logger.warning(f"Units mismatch between test ({units}) and reference ({ref_units})")

        logger.info(f"Reference data metrics - min: {ref_min}, mean: {ref_mean}, max: {ref_max}, units: {ref_units}")

        # Calculate metrics for difference
        diff_mean = test_mean - ref_mean
        # For min/max, we use the range of possible differences
        diff_min = test_min - ref_max  # This may be a simplification
        diff_max = test_max - ref_min  # This may be a simplification
        logger.info(f"Difference metrics - min: {diff_min}, mean: {diff_mean}, max: {diff_max}")

    # Create panels following the lat_lon_plot layout
    # Panel 1: Test data (always created)
    ax1 = fig.add_axes(DEFAULT_PANEL_CFG[0], projection=projection)
    # Use the standard title configuration from utils._configure_titles
    # Format: (years_text on left, main title in center, units on right)
    ax1.set_title(parameter.test_name_yrs, loc="left", fontdict={"fontsize": 9.5})
    ax1.set_title(parameter.test_title, fontdict={"fontsize": 11.5})
    ax1.set_title(units, loc="right", fontdict={"fontsize": 9.5})

    # Initialize ax2 and ax3 as None - they'll only be created when reference data exists
    ax2 = None
    ax3 = None

    # Only create panels 2 and 3 when reference data is available
    if has_reference:
        # Panel 2: Reference data
        ax2 = fig.add_axes(DEFAULT_PANEL_CFG[1], projection=projection)
        ax2.set_title(parameter.ref_name_yrs, loc="left", fontdict={"fontsize": 9.5})
        ax2.set_title(parameter.reference_title, fontdict={"fontsize": 11.5})
        ax2.set_title(units, loc="right", fontdict={"fontsize": 9.5})

        # Panel 3: Difference plot
        ax3 = fig.add_axes(DEFAULT_PANEL_CFG[2], projection=projection)
        ax3.set_title("", loc="left", fontdict={"fontsize": 9.5})
        ax3.set_title(parameter.diff_title, fontdict={"fontsize": 11.5})
        ax3.set_title(units, loc="right", fontdict={"fontsize": 9.5})

    # Configure map settings for all created panels
    panels = [ax for ax in [ax1, ax2, ax3] if ax is not None]
    for ax in panels:
        # Debug extent issues
        logger.info(f"Setting extent: lon=[{lon_west}, {lon_east}], lat=[{lat_south}, {lat_north}]")

        # Handle global domain specially - don't use set_extent for global domain
        if is_global_domain:
            logger.info("Using global view")
            ax.set_global()
        else:
            try:
                # More robust longitude handling for map extents
                lon_west_orig, lon_east_orig = lon_west, lon_east

                # Normalize to [-180, 180] range for consistency with cartopy
                if lon_west > 180:
                    lon_west -= 360
                if lon_east > 180:
                    lon_east -= 360

                # Handle cases where the region crosses the dateline
                if lon_east < lon_west:
                    # This is a dateline-crossing region (e.g., Pacific)
                    logger.info(f"Detected dateline crossing region: lon=[{lon_west}, {lon_east}]")

                    # For dateline-crossing regions, adjust the central longitude of projection
                    center_lon = (lon_west + lon_east + 360) / 2.0
                    if center_lon > 180:
                        center_lon -= 360

                    logger.info(f"Using central longitude: {center_lon}")

                    # Create a new projection with the adjusted central longitude
                    ax.projection = ccrs.PlateCarree(central_longitude=center_lon)

                    # When using a central_longitude, we need to transform our coordinates
                    # Adjust longitudes for the new central longitude
                    if lon_west < 0:
                        lon_west += 360
                    if lon_east < 0:
                        lon_east += 360

                    logger.info(f"Transformed coordinates for central_longitude={center_lon}: lon=[{lon_west}, {lon_east}]")

                # Make sure longitudes are properly ordered
                if lon_east < lon_west:
                    logger.warning("East longitude still less than west after transforms - swapping values")
                    lon_west, lon_east = lon_east, lon_west

                logger.info(f"Final map extent: lon=[{lon_west}, {lon_east}], lat=[{lat_south}, {lat_north}]")

                # Set the extent using the adjusted longitude values
                ax.set_extent([lon_west, lon_east, lat_south, lat_north], crs=ccrs.PlateCarree())

                # Log actual extent for debugging
                actual_extent = ax.get_extent(crs=ccrs.PlateCarree())
                logger.info(f"Actual extent after set_extent: {actual_extent}")

            except Exception as e:
                # Comprehensive error handling
                logger.error(f"Error setting map extent: {e}")
                logger.error(f"Original lon bounds: [{lon_west_orig}, {lon_east_orig}]")
                logger.error(f"Transformed lon bounds: [{lon_west}, {lon_east}]")
                import traceback
                logger.error(traceback.format_exc())

                # Fallback to global view
                logger.info("Falling back to global view due to extent error")
                ax.set_global()

        # Add map features
        ax.coastlines(linewidth=0.5)
        ax.add_feature(cfeature.BORDERS, linewidth=0.3)

        # Configure gridlines and labels
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                         linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False

        # Set aspect ratio only for non-global views
        if not is_global_domain:
            # Ensure we have a valid aspect ratio by clamping to reasonable values
            width = lon_east - lon_west
            height = lat_north - lat_south
            if width <= 0:
                width += 360  # Handle wraparound cases
            aspect_ratio = width / (2 * max(height, 1))  # Avoid division by zero or negative values
            logger.info(f"Setting aspect ratio: {aspect_ratio}")
            ax.set_aspect(aspect_ratio)

    # Debug information
    logger.info(f"Plotting variable {var_key} from uxarray dataset")
    logger.info(f"Variable shape: {uxds_test[var_key].shape}")
    logger.info(f"Variable min/max: {test_min}/{test_max}")
    if 'units' in uxds_test[var_key].attrs:
        logger.info(f"Variable units: {units}")

    # Create the PolyCollection for test data
    logger.info("Creating PolyCollection")

    # Check for valid face coordinates in uxarray data
    if hasattr(uxds_test, 'face') and len(uxds_test.face) > 0:
        logger.info(f"Found valid face mesh with {len(uxds_test.face)} elements")
    else:
        logger.warning("No valid face mesh found in uxarray dataset!")

    # Handle periodic elements setting
    periodic_elements = None
    if hasattr(parameter, 'split_periodic_elements'):
        periodic_elements = "split" if parameter.split_periodic_elements else None
        logger.info(f"Using periodic_elements={periodic_elements}")

    # Squeeze the data array to remove singleton dimensions
    test_var = uxds_test[var_key].squeeze()
    logger.info(f"Variable dimensions after squeeze: {test_var.dims}")

    # Get colormap from parameter
    cmap = _get_colormap(parameter.test_colormap)
    logger.info(f"Using colormap: {parameter.test_colormap}")

    # Create normalization based on contour levels or data range
    if hasattr(parameter, 'contour_levels') and parameter.contour_levels and len(parameter.contour_levels) > 0:
        # Use the contour levels to create a BoundaryNorm
        levels = parameter.contour_levels
        logger.info(f"Using contour levels for normalization: {levels}")
        # Create boundary norm with extended boundaries for values beyond the levels
        boundaries = [-1.0e8] + levels + [1.0e8]
        norm = mcolors.BoundaryNorm(boundaries=boundaries, ncolors=256)
    else:
        # Use a simple min/max normalization based on data range
        vmin, vmax = test_min, test_max
        # Add a small buffer to prevent issues with constant values
        if vmin == vmax:
            buffer = max(0.1, abs(vmin * 0.1))
            vmin -= buffer
            vmax += buffer
            logger.warning(f"Data has constant value {test_min}, adding buffer for visualization: [{vmin}, {vmax}]")
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
        logger.info(f"Using data range for normalization: {vmin} to {vmax}")

    # Create the PolyCollection
    pc = test_var.to_polycollection(periodic_elements=periodic_elements)

    # Set visualization properties
    pc.set_cmap(cmap)
    pc.set_norm(norm)  # Apply the normalization
    pc.set_antialiased(parameter.antialiased)

    # Add to panel 1 (test data panel)
    ax1.add_collection(pc)

    # Add reference data visualization if available
    if has_reference and ax2 is not None:
        # Reference data processing (panel 2)
        logger.info("Creating PolyCollection for reference data")

        # Squeeze the data array to remove singleton dimensions
        ref_var = uxds_ref[var_key].squeeze()
        logger.info(f"Reference variable dimensions: {ref_var.dims}")

        # Get colormap for reference data
        ref_cmap = _get_colormap(parameter.reference_colormap)
        logger.info(f"Using colormap for reference: {parameter.reference_colormap}")

        # Create normalization for reference data
        if hasattr(parameter, 'contour_levels') and parameter.contour_levels and len(parameter.contour_levels) > 0:
            # Use the contour levels to create a BoundaryNorm
            ref_levels = parameter.contour_levels
            logger.info(f"Using contour levels for reference normalization: {ref_levels}")
            # Create boundary norm with extended boundaries for values beyond the levels
            ref_boundaries = [-1.0e8] + ref_levels + [1.0e8]
            ref_norm = mcolors.BoundaryNorm(boundaries=ref_boundaries, ncolors=256)
        else:
            # Use a simple min/max normalization based on data range
            ref_vmin, ref_vmax = ref_min, ref_max
            # Add a small buffer to prevent issues with constant values
            if ref_vmin == ref_vmax:
                buffer = max(0.1, abs(ref_vmin * 0.1))
                ref_vmin -= buffer
                ref_vmax += buffer
                logger.warning(f"Reference data has constant value {ref_min}, adding buffer: [{ref_vmin}, {ref_vmax}]")
            ref_norm = mcolors.Normalize(vmin=ref_vmin, vmax=ref_vmax)
            logger.info(f"Using data range for reference normalization: {ref_vmin} to {ref_vmax}")

        # Create the PolyCollection for reference data
        pc_ref = ref_var.to_polycollection(periodic_elements=periodic_elements)

        # Set visualization properties
        pc_ref.set_cmap(ref_cmap)
        pc_ref.set_norm(ref_norm)
        pc_ref.set_antialiased(parameter.antialiased)

        # Add to panel 2 (reference data panel)
        ax2.add_collection(pc_ref)

        # Add colorbar for reference data panel
        cbax_rect_ref = (
            DEFAULT_PANEL_CFG[1][0] + 0.6635,  # Position relative to the panel
            DEFAULT_PANEL_CFG[1][1] + 0.0215,
            0.0326,   # Width
            0.1792,   # Height
        )
        cbax_ref = fig.add_axes(cbax_rect_ref)
        cbar_ref = fig.colorbar(pc_ref, cax=cbax_ref, extend='both')

        # Configure reference colorbar ticks based on contour levels
        if hasattr(parameter, 'contour_levels') and parameter.contour_levels and len(parameter.contour_levels) > 0:
            cbar_ref.set_ticks(parameter.contour_levels)

            # Format tick labels
            maxval = np.amax(np.absolute(parameter.contour_levels))
            if maxval < 0.01:
                fmt, pad = "%.1e", 35
            elif maxval < 0.2:
                fmt, pad = "%5.3f", 28
            elif maxval < 10.0:
                fmt, pad = "%5.2f", 25
            elif maxval < 100.0:
                fmt, pad = "%5.1f", 25
            elif maxval > 9999.0:
                fmt, pad = "%.0f", 40
            else:
                fmt, pad = "%6.1f", 30

            labels = [fmt % level for level in parameter.contour_levels]
            cbar_ref.ax.set_yticklabels(labels, ha="right")
            cbar_ref.ax.tick_params(labelsize=9.0, pad=pad, length=0)
        else:
            cbar_ref.ax.tick_params(labelsize=9.0, length=0)

        # Add units label
        cbar_ref.set_label(units, fontsize=9.5)

        # Add metrics text for reference panel
        metrics_ref = (ref_max, ref_mean, ref_min)

        # Format metrics based on size
        fmt_m_ref = []
        for i in range(3):
            if metrics_ref[i] < 0.01 and metrics_ref[i] > 0:
                fs = "e"
            elif metrics_ref[i] > 100000.0:
                fs = "e"
            else:
                fs = "2f"
            fmt_m_ref.append(fs)

        fmt_metrics_ref = f"%.{fmt_m_ref[0]}\n%.{fmt_m_ref[1]}\n%.{fmt_m_ref[2]}"

        # Add "Max, Mean, Min" labels for reference panel
        fig.text(
            DEFAULT_PANEL_CFG[1][0] + 0.6635,
            DEFAULT_PANEL_CFG[1][1] + 0.2107,
            "Max\nMean\nMin",
            ha="left",
            fontdict={"fontsize": 9.5},
        )

        # Add the actual metric values for reference panel
        fig.text(
            DEFAULT_PANEL_CFG[1][0] + 0.7635,
            DEFAULT_PANEL_CFG[1][1] + 0.2107,
            fmt_metrics_ref % metrics_ref,
            ha="right",
            fontdict={"fontsize": 9.5},
        )

        # For panel 3 (difference), add a placeholder message until differential comparison is implemented
        if ax3 is not None:
            ax3.text(0.5, 0.5, "Difference calculation not yet implemented for native grid data",
                    transform=ax3.transAxes, ha='center', va='center', fontsize=11)

    # Add colorbar for the test data panel
    # Use the standard layout position from lat_lon_plot
    cbax_rect = (
        DEFAULT_PANEL_CFG[0][0] + 0.6635,  # Position relative to the panel
        DEFAULT_PANEL_CFG[0][1] + 0.0215,
        0.0326,   # Width
        0.1792,   # Height
    )
    # Create the colorbar axes
    cbax = fig.add_axes(cbax_rect)
    # Add the colorbar using the PolyCollection with arrows at both ends
    # The extend='both' parameter adds arrows at both ends of the colorbar
    # This matches the behavior in the standard lat_lon_plot
    cbar = fig.colorbar(pc, cax=cbax, extend='both')

    # Configure colorbar ticks based on contour levels if provided
    if hasattr(parameter, 'contour_levels') and parameter.contour_levels and len(parameter.contour_levels) > 0:
        # Use the specified contour levels for ticks
        levels = parameter.contour_levels
        logger.info(f"Using contour levels for colorbar: {levels}")

        # Add extended boundaries for BoundaryNorm to match _get_c_levels_and_norm function
        c_levels = [-1.0e8] + levels + [1.0e8]
        cbar.set_ticks(levels)  # Only show the actual levels, not the extended boundaries

        # Format tick labels based on the range of values - directly matching utils._get_contour_label_format_and_pad
        maxval = np.amax(np.absolute(levels))
        if maxval < 0.01:
            fmt, pad = "%.1e", 35  # Scientific notation for very small values
        elif maxval < 0.2:
            fmt, pad = "%5.3f", 28  # 3 decimal places for small values
        elif maxval < 10.0:
            fmt, pad = "%5.2f", 25  # 2 decimal places for moderate values
        elif maxval < 100.0:
            fmt, pad = "%5.1f", 25  # 1 decimal place for larger values
        elif maxval > 9999.0:
            fmt, pad = "%.0f", 40   # No decimals for very large values
        else:
            fmt, pad = "%6.1f", 30  # 1 decimal place with padding for other values

        # Format the labels and apply them
        labels = [fmt % level for level in levels]
        cbar.ax.set_yticklabels(labels, ha="right")
        cbar.ax.tick_params(labelsize=9.0, pad=pad, length=0)
    else:
        # Default formatting if no contour levels specified
        cbar.ax.tick_params(labelsize=9.0, length=0)

    # Add the units label with matching style from lat_lon_plot
    cbar.set_label(units, fontsize=9.5)

    # Add the metrics text (min, mean, max) below the panel
    # Same position used in lat_lon_plot
    # Add "Max, Mean, Min" labels
    fig.text(
        DEFAULT_PANEL_CFG[0][0] + 0.6635,  # Left position
        DEFAULT_PANEL_CFG[0][1] + 0.2107,  # Bottom position
        "Max\nMean\nMin",
        ha="left",
        fontdict={"fontsize": 9.5},
    )

    # Format the metric values based on their size, matching utils._get_float_format
    # Print in scientific notation if values are very small or very large
    metrics = (test_max, test_mean, test_min)
    fmt_m = []

    # Determine format for each value (max, mean, min)
    for i in range(3):
        if metrics[i] < 0.01 and metrics[i] > 0:
            # For very small positive values, use scientific notation
            fs = "e"
        elif metrics[i] > 100000.0:
            # For very large values, use scientific notation
            fs = "e"
        else:
            # For regular values, use fixed-point notation
            fs = "2f"
        fmt_m.append(fs)

    # Build format string with appropriate format for each value
    fmt_metrics = f"%.{fmt_m[0]}\n%.{fmt_m[1]}\n%.{fmt_m[2]}"

    # Add the actual metric values
    fig.text(
        DEFAULT_PANEL_CFG[0][0] + 0.7635,  # Right position
        DEFAULT_PANEL_CFG[0][1] + 0.2107,  # Bottom position
        fmt_metrics % metrics,  # Using the metrics tuple from above
        ha="right",
        fontdict={"fontsize": 9.5},
    )
    # No need to add text to panels that don't exist
    # In lat_lon_plot style, model-only runs just have a single panel

    # Save the plot - use the standard output directory structure to ensure viewer compatibility
    output_dir = _get_output_dir(parameter)

    # Log the directory structure for debugging
    logger.info(f"Saving plot to directory: {output_dir}")
    logger.info(f"Output filename base: {parameter.output_file}")

    for fmt in parameter.output_format:
        filename = os.path.join(output_dir, f"{parameter.output_file}.{fmt}")
        fig.savefig(filename)
        logger.info(f"Plot saved to: {filename}")

    # Save individual subplots if requested
    if hasattr(parameter, 'output_format_subplot') and parameter.output_format_subplot:
        # Following the pattern from _save_plot in utils.py
        # For model-only runs, only save the first panel
        if not has_reference:
            panels_to_save = [0]  # Just the test panel (index 0)
        else:
            panels_to_save = range(len(panels))  # All panels

        for i in panels_to_save:
            for fmt in parameter.output_format_subplot:
                subplot_filename = os.path.join(output_dir, f"{parameter.output_file}.{i}.{fmt}")

                # Get the position of the current panel and expand it slightly
                panel_pos = fig.axes[i].get_position()
                extent = panel_pos.expanded(1.1, 1.2)

                # Save the individual panel
                fig.savefig(subplot_filename, bbox_inches=extent)
                logger.info(f"Panel {i} saved to: {subplot_filename}")

    plt.close(fig)

def _run_diags_2d_model_only(
    parameter: LatLonNativeParameter,
    season: str,
    regions: List[str],
    var_key: str,
    ref_name: str,
    uxds_test: ux.dataset.UxDataset,
):
    """Run a model-only diagnostics on a 2D variable using native grid.

    This function plots the native grid data directly using uxarray dataset.

    Parameters
    ----------
    parameter : LatLonNativeParameter
        The parameter object.
    season : str
        The season.
    regions : List[str]
        The list of regions.
    var_key : str
        The key of the variable.
    ref_name : str
        The reference name.
    uxds_test : ux.dataset.UxDataset
        The uxarray dataset containing the test native grid information.
    """
    # Debug the uxarray dataset
    logger.info("===== TEST UXDS DEBUG INFO =====")
    logger.info(f"Type: {type(uxds_test)}")
    logger.info(f"Variables: {list(uxds_test.data_vars)}")
    logger.info(f"Has face attribute: {hasattr(uxds_test, 'face')}")

    # Check if the variable exists in the uxarray dataset
    if var_key not in uxds_test:
        logger.warning(f"Variable {var_key} not found in uxarray dataset!")
        logger.info(f"Available variables: {list(uxds_test.data_vars)}")
        return

    # Check data values for potential issues
    logger.info(f"Variable {var_key} shape: {uxds_test[var_key].shape}")

    # Get min/max values
    var_min = uxds_test[var_key].min().item()
    var_max = uxds_test[var_key].max().item()
    logger.info(f"Variable {var_key} min/max: {var_min}/{var_max}")

    # Check if all values are identical - this would result in a single-color plot
    var_data = uxds_test[var_key].values
    n_unique = len(np.unique(var_data))
    if n_unique <= 1:
        logger.warning(f"⚠️ WARNING: Variable {var_key} contains only {n_unique} unique value!")
        logger.warning("Data will appear as a solid color - possible issue with data or derivation")

    # Log units and other attributes to confirm derivation worked correctly
    if 'units' in uxds_test[var_key].attrs:
        logger.info(f"Variable {var_key} units: {uxds_test[var_key].attrs['units']}")
    else:
        logger.warning(f"Variable {var_key} has no units attribute")

    # Log other important attributes
    for attr_name in ['long_name', 'standard_name', 'cell_methods', 'description']:
        if attr_name in uxds_test[var_key].attrs:
            logger.info(f"Variable {var_key} {attr_name}: {uxds_test[var_key].attrs[attr_name]}")

    logger.info("================================")

    for region in regions:
        parameter._set_param_output_attrs(var_key, season, region, ref_name, ilev=None)
        print(f"Processing {var_key} for region {region}")

        # Use the native grid plotting function with explicit uxarray dataset
        _save_native_plot(
            parameter,
            var_key,
            region,
            uxds_test=uxds_test,
            uxds_ref=None,
            uxds_diff=None,
        )


def _run_diags_3d_model_only(
    parameter: LatLonNativeParameter,
    ds_test: xr.Dataset,
    ds_land_sea_mask: xr.Dataset,
    season: str,
    regions: List[str],
    var_key: str,
    ref_name: str,
    uxds_test: ux.dataset.UxDataset = None,
):
    """Run a model-only diagnostics on a 3D variable using native grid.

    This function gets the variable's metrics by region, then saves the
    metrics, metric plots, and data (optional, `CoreParameter.save_netcdf`).

    Parameters
    ----------
    parameter : LatLonNativeParameter
        The parameter object.
    ds_test : xr.Dataset
        The dataset containing the test variable.
    ds_land_sea_mask : xr.Dataset
        The land sea mask dataset, which is only used for masking if the region
        is "land" or "ocean".
    season : str
        The season.
    regions : List[str]
        The list of regions.
    var_key : str
        The key of the variable.
    ref_name : str
        The reference name.
    uxds_test : ux.dataset.UxDataset, optional
        The uxarray dataset containing the test native grid information.
    """
    plev = parameter.plevs
    logger.info("Selected pressure level(s): {}".format(plev))

    ds_test_rg = regrid_z_axis_to_plevs(ds_test, var_key, parameter.plevs)

    for ilev in plev:
        z_axis_key = get_z_axis(ds_test_rg[var_key]).name
        ds_test_ilev = ds_test_rg.sel({z_axis_key: ilev})

        for region in regions:
            ds_test_region = _process_test_dataset(
                parameter, ds_test_ilev, ds_land_sea_mask, var_key, region
            )
            parameter._set_param_output_attrs(var_key, season, region, ref_name, ilev)
            print(f"Processing {var_key} for region {region} at level {ilev}")

            # Use the native grid plotting function
            _save_native_plot(
                parameter,
                var_key,
                region,
                uxds_test=uxds_test,
                uxds_ref=None,
            )


def _check_var_dims(
    ds_test: xr.Dataset, ds_ref: xr.Dataset, var_key: str
) -> Tuple[bool, bool]:
    """Check if the variables have 3D dimensions and if their dimensions are different.

    Parameters
    ----------
    ds_test : xr.Dataset
        The test dataset.
    ds_ref : xr.Dataset
        The reference dataset.
    var_key : str
        The key of the variable.

    Returns
    -------
    Tuple[bool, bool]
        A tuple containing two boolean values:
        - is_vars_3d: True if both variables have 3D dimensions, False otherwise.
        - is_dims_diff: True if the dimensions of the two variables are different, False otherwise.
    """
    dv_test = ds_test[var_key]
    dv_ref = ds_ref[var_key]

    is_vars_3d = has_z_axis(dv_test) and has_z_axis(dv_ref)
    is_dims_diff = has_z_axis(dv_test) != has_z_axis(dv_ref)

    return is_vars_3d, is_dims_diff


def _run_diags_2d(
    parameter: LatLonNativeParameter,
    ds_test: xr.Dataset,
    ds_ref: xr.Dataset,
    ds_land_sea_mask: xr.Dataset,
    season: str,
    regions: List[str],
    var_key: str,
    ref_name: str,
    uxds_test: ux.dataset.UxDataset = None,
    uxds_ref: ux.dataset.UxDataset = None,
):
    """Run diagnostics on a 2D variable using native grid.

    This function gets the variable's metrics by region, then saves the
    metrics, metric plots, and data (optional, `CoreParameter.save_netcdf`).

    Parameters
    ----------
    parameter : LatLonNativeParameter
        The parameter object.
    ds_test : xr.Dataset
        The dataset containing the test variable.
    ds_ref : xr.Dataset
        The dataset containing the ref variable.
    ds_land_sea_mask : xr.Dataset
        The land sea mask dataset, which is only used for masking if the region
        is "land" or "ocean".
    season : str
        The season.
    regions : List[str]
        The list of regions.
    var_key : str
        The key of the variable.
    ref_name : str
        The reference name.
    uxds_test : ux.dataset.UxDataset, optional
        The uxarray dataset containing the test native grid information.
    uxds_ref : ux.dataset.UxDataset, optional
        The uxarray dataset containing the reference native grid information.
    """
    # Log whether we have reference data in uxarray format
    if uxds_ref is None:
        logger.warning("No uxarray reference dataset (uxds_ref is None) - this will cause model_vs_model to fall back to model-only")
    elif var_key not in uxds_ref:
        logger.warning(f"Variable {var_key} not found in reference uxarray dataset")
        logger.info(f"Available variables in reference dataset: {list(uxds_ref.data_vars)}")
    else:
        logger.info(f"Reference data for {var_key} is available in uxarray format")

    # Check test dataset too
    if uxds_test is None:
        logger.error("No test uxarray dataset (uxds_test is None) - cannot proceed with native grid visualization")
        return
    elif var_key not in uxds_test:
        logger.error(f"Variable {var_key} not found in test uxarray dataset")
        logger.info(f"Available variables in test dataset: {list(uxds_test.data_vars)}")
        return
    # Check if we have valid reference data
    has_valid_ref = uxds_ref is not None and var_key in uxds_ref

    for region in regions:
        parameter._set_param_output_attrs(var_key, season, region, ref_name, ilev=None)

        if has_valid_ref:
            print(f"Processing {var_key} for region {region} (model vs model)")
            # Directly pass both test and reference data to the plotting function
            _save_native_plot(
                parameter,
                var_key,
                region,
                uxds_test=uxds_test,
                uxds_ref=uxds_ref,
            )
        else:
            print(f"Processing {var_key} for region {region} (model-only fallback)")
            # If reference data is missing, fall back to model-only mode
            _save_native_plot(
                parameter,
                var_key,
                region,
                uxds_test=uxds_test,
                uxds_ref=None,  # Force model-only behavior
            )


def _run_diags_3d(
    parameter: LatLonNativeParameter,
    ds_test: xr.Dataset,
    ds_ref: xr.Dataset,
    ds_land_sea_mask: xr.Dataset,
    season: str,
    regions: List[str],
    var_key: str,
    ref_name: str,
    uxds_test: ux.dataset.UxDataset = None,
    uxds_ref: ux.dataset.UxDataset = None,
):
    """Run diagnostics on a 3D variable using native grid.

    This function gets the variable's metrics by region, then saves the
    metrics, metric plots, and data (optional, `CoreParameter.save_netcdf`).

    Parameters
    ----------
    parameter : LatLonNativeParameter
        The parameter object.
    ds_test : xr.Dataset
        The dataset containing the test variable.
    ds_ref : xr.Dataset
        The dataset containing the ref variable.
    ds_land_sea_mask : xr.Dataset
        The land sea mask dataset, which is only used for masking if the region
        is "land" or "ocean".
    season : str
        The season.
    regions : List[str]
        The list of regions.
    var_key : str
        The key of the variable.
    ref_name : str
        The reference name.
    uxds_test : ux.dataset.UxDataset, optional
        The uxarray dataset containing the test native grid information.
    uxds_ref : ux.dataset.UxDataset, optional
        The uxarray dataset containing the reference native grid information.
    """
    plev = parameter.plevs
    logger.info("Selected pressure level(s): {}".format(plev))

    ds_test_rg = regrid_z_axis_to_plevs(ds_test, var_key, parameter.plevs)
    ds_ref_rg = regrid_z_axis_to_plevs(ds_ref, var_key, parameter.plevs)

    for ilev in plev:
        z_axis_key = get_z_axis(ds_test_rg[var_key]).name
        ds_test_ilev = ds_test_rg.sel({z_axis_key: ilev})
        ds_ref_ilev = ds_ref_rg.sel({z_axis_key: ilev})

        for region in regions:
            (
                ds_test_region,
                ds_test_region_regrid,
                ds_ref_region,
                ds_ref_region_regrid,
                ds_diff_region,
            ) = _subset_and_align_native_datasets(
                parameter,
                ds_test_ilev,
                ds_ref_ilev,
                ds_land_sea_mask,
                var_key,
                region,
                uxds_test,
                uxds_ref,
            )

            metrics_dict = _create_metrics_dict(
                var_key,
                ds_test_region,
                ds_test_region_regrid,
                ds_ref_region,
                ds_ref_region_regrid,
                ds_diff_region,
                uxds_test=uxds_test,
                uxds_ref=uxds_ref,
            )

            parameter._set_param_output_attrs(var_key, season, region, ref_name, ilev)
            print(f"Processing {var_key} for region {region} at level {ilev} (model vs obs)")

            # Use the native grid plotting function
            _save_native_plot(
                parameter,
                var_key,
                region,
                uxds_test=uxds_test,
                uxds_ref=uxds_ref,
            )


def _get_ref_climo_dataset(
    dataset: Dataset, var_key: str, season: ClimoFreq
) -> xr.Dataset | None:
    """Get the reference climatology dataset for the variable and season.

    If the reference climatatology does not exist or could not be found, it
    will be considered a model-only run and return `None`.

    This function also stores the reference data file path in the parameter object
    for native grid visualization.

    Parameters
    ----------
    dataset : Dataset
        The dataset object.
    var_key : str
        The key of the variable.
    season : CLIMO_FREQ
        The climatology frequency.

    Returns
    -------
    xr.Dataset | None
        The reference climatology if it exists or None if it does not.
        None indicates a model-only run.

    Raises
    ------
    RuntimeError
        If `self.data_type` is not "ref".
    """
    if dataset.data_type == "ref":
        try:
            # Get the reference climatology dataset
            ds_ref = dataset.get_climo_dataset(var_key, season)

            # The ref_data_file_path should already be set in get_climo_dataset,
            # but we log it here for debugging
            if hasattr(dataset.parameter, 'ref_data_file_path'):
                logger.info(f"Reference data file path: {dataset.parameter.ref_data_file_path}")
            else:
                logger.warning("ref_data_file_path not set in parameter")

            # Additional ways to extract the file path for debugging
            file_path = None
            if hasattr(ds_ref, 'file_path'):
                file_path = ds_ref.file_path
                logger.info(f"Reference data file path from ds_ref.file_path: {file_path}")
            elif hasattr(ds_ref, 'filepath'):
                file_path = ds_ref.filepath
                logger.info(f"Reference data file path from ds_ref.filepath: {file_path}")
            elif hasattr(ds_ref, '_file_obj') and hasattr(ds_ref._file_obj, 'name'):
                file_path = ds_ref._file_obj.name
                logger.info(f"Reference data file path from ds_ref._file_obj.name: {file_path}")

            # Store additional path in parameter if found through other methods
            if file_path and not hasattr(dataset.parameter, 'ref_data_file_path'):
                dataset.parameter.ref_data_file_path = file_path

        except (RuntimeError, IOError) as e:
            ds_ref = None
            logger.info(f"Cannot process reference data due to error: {e}")
            logger.info("Analyzing test data only.")
    else:
        raise RuntimeError(
            "`Dataset._get_ref_dataset` only works with "
            f"`parameter.data_type == 'ref'`, not {dataset.data_type}."
        )

    return ds_ref


def _process_test_dataset(
    parameter: LatLonNativeParameter,
    ds_test: xr.Dataset,
    ds_land_sea_mask: xr.Dataset,
    var_key: str,
    region: str,
) -> xr.Dataset:
    """Process the test dataset for the given region.

    Parameters
    ----------
    parameter : LatLonNativeParameter
        The core parameter object.
    ds_test : xr.Dataset
        The test dataset.
    ds_land_sea_mask : xr.Dataset
        The dataset for land-sea mask.
    var_key : str
        The variable key.
    region : str
        The region to process.

    Returns:
    --------
    xr.Dataset
        The processed test dataset for the region.
    """
    ds_test_region = ds_test.copy()

    if "land" in region or "ocean" in region:
        ds_test_region = _apply_land_sea_mask(
            ds_test,
            ds_land_sea_mask,
            var_key,
            region,  # type: ignore
            parameter.regrid_tool,
            parameter.regrid_method,
        )

    if "global" not in region:
        ds_test_region = _subset_on_region(ds_test_region, var_key, region)

    return ds_test_region


def _subset_and_align_native_datasets(
    parameter: LatLonNativeParameter,
    ds_test: xr.Dataset,
    ds_ref: xr.Dataset,
    ds_land_sea_mask: xr.Dataset,
    var_key: str,
    region: str,
    uxds: ux.dataset.UxDataset = None,
):
    """Subset and align datasets using native grid data.

    For native grid data, this is an adapted version that preserves the native grid
    structure while still computing metrics and difference fields.

    Parameters
    ----------
    parameter : LatLonNativeParameter
        The parameter object.
    ds_test : xr.Dataset
        The dataset containing the test variable.
    ds_ref : xr.Dataset
        The dataset containing the reference variable.
    ds_land_sea_mask : xr.Dataset
        The land sea mask dataset.
    var_key : str
        The variable key.
    region : str
        The region.
    uxds : ux.dataset.UxDataset, optional
        The uxarray dataset containing the native grid information.

    Returns
    -------
    Tuple[xr.Dataset, xr.Dataset, xr.Dataset, xr.Dataset, xr.Dataset]
        A tuple containing the processed datasets:
        - ds_test_region: The test dataset subsetted to the region.
        - ds_test_region_regrid: The test dataset regridded (same as ds_test_region for native).
        - ds_ref_region: The reference dataset subsetted to the region.
        - ds_ref_region_regrid: The reference dataset regridded to match the test data resolution.
        - ds_diff_region: The difference between the test and reference datasets.
    """
    # First, process the test and reference datasets for the given region
    ds_test_region = _process_test_dataset(
        parameter, ds_test, ds_land_sea_mask, var_key, region
    )

    ds_ref_region = _process_test_dataset(
        parameter, ds_ref, ds_land_sea_mask, var_key, region
    )

    # For native grid, we don't regrid the test data, we keep it on its native grid
    ds_test_region_regrid = ds_test_region

    # For model vs obs comparison, we still need to handle the reference data
    # For now, we'll use a simple approach - the reference data stays on its grid
    # but metrics calculations will handle the irregular grid
    ds_ref_region_regrid = ds_ref_region

    # For difference calculation, this is more complex with unstructured grids
    # In a first implementation, we'll compute a simple difference without regridding
    # In future versions, we might want to implement more sophisticated methods
    # Create a copy of the test dataset for the difference
    ds_diff_region = ds_test_region_regrid.copy()

    # For basic implementation, we'll just subtract the means for the difference plot
    test_mean = spatial_avg(ds_test_region_regrid, var_key)
    ref_mean = spatial_avg(ds_ref_region_regrid, var_key)
    diff = test_mean - ref_mean

    # Create a new variable in the difference dataset with a constant value
    # This is a simplification for the first implementation
    # In future versions, we may want to implement proper regridding between native and reference grids
    ds_diff_region[var_key].values = ds_diff_region[var_key].values * 0 + diff

    return (
        ds_test_region,
        ds_test_region_regrid,
        ds_ref_region,
        ds_ref_region_regrid,
        ds_diff_region,
    )


def _create_metrics_dict(
    var_key: str,
    ds_test: xr.Dataset,
    ds_test_regrid: xr.Dataset | None,
    ds_ref: xr.Dataset | None,
    ds_ref_regrid: xr.Dataset | None,
    ds_diff: xr.Dataset | None,
    uxds: ux.dataset.UxDataset = None,
) -> MetricsDict:
    """Calculate metrics using the variable in the datasets for native grid visualization.

    Metrics include min value, max value, spatial average (mean), standard
    deviation, correlation (pearson_r), and RMSE. The default value for
    optional metrics is None.

    Parameters
    ----------
    var_key : str
        The variable key.
    ds_test : xr.Dataset
        The test dataset.
    ds_test_regrid : xr.Dataset | None
        The regridded test Dataset. This arg will be None if a model only
        run is performed.
    ds_ref : xr.Dataset | None
        The optional reference dataset. This arg will be None if a model only
        run is performed.
    ds_ref_regrid : xr.Dataset | None
        The optional regridded reference dataset. This arg will be None if a
        model only run is performed.
    ds_diff : xr.Dataset | None
        The difference between ``ds_test_regrid`` and ``ds_ref_regrid`` if both
        exist. This arg will be None if a model only run is performed.
    uxds : ux.dataset.UxDataset, optional
        The uxarray dataset containing the native grid information.

    Returns
    -------
    Metrics
        A dictionary with the key being a string and the value being either
        a sub-dictionary (key is metric and value is float) or a string
        ("unit").
    """
    # Extract these variables for reuse.
    var_test = ds_test[var_key]

    # xarray.DataArray.min() and max() returns a `np.ndarray` with a single
    # int/float element. Using `.item()` returns that single element.
    metrics_dict: MetricsDict = {
        "test": {
            "min": var_test.min().item(),
            "max": var_test.max().item(),
            "mean": spatial_avg(ds_test, var_key),  # type: ignore
            "std": std(ds_test, var_key),
        },
        "unit": ds_test[var_key].attrs["units"],
    }
    metrics_dict = _set_default_metric_values(metrics_dict)

    # Additional native grid information
    if uxds is not None:
        try:
            grid_info = {
                "grid_type": getattr(uxds, "grid_type", "unknown"),
                "ne": getattr(uxds, "ne", ""),
                "npe": getattr(uxds, "npe", ""),
                "element_count": len(uxds.face) if hasattr(uxds, "face") else 0,
            }
            metrics_dict["grid_info"] = grid_info
        except Exception as e:
            logger.warning(f"Error adding grid info to metrics: {e}")

    if ds_ref is not None:
        var_ref = ds_ref[var_key]

        metrics_dict["ref"] = {
            "min": var_ref.min().item(),
            "max": var_ref.max().item(),
            "mean": spatial_avg(ds_ref, var_key),  # type: ignore
            "std": std(ds_ref, var_key),
        }

    if ds_test_regrid is not None and ds_ref_regrid is not None:
        var_test_regrid = ds_test_regrid[var_key]
        metrics_dict["test_regrid"] = {
            "min": var_test_regrid.min().item(),
            "max": var_test_regrid.max().item(),
            "mean": spatial_avg(ds_test_regrid, var_key),  # type: ignore
            "std": std(ds_test_regrid, var_key),
        }

        var_ref_regrid = ds_ref_regrid[var_key]
        metrics_dict["ref_regrid"] = {
            "min": var_ref_regrid.min().item(),
            "max": var_ref_regrid.max().item(),
            "mean": spatial_avg(ds_ref_regrid, var_key),  # type: ignore
            "std": std(ds_ref_regrid, var_key),
        }

        # For native grid, the correlation and RMSE calculations might need to be adapted
        # In the first implementation, we'll use global means for simplicity
        metrics_dict["misc"] = {
            "rmse": abs(
                metrics_dict["test_regrid"]["mean"] - metrics_dict["ref_regrid"]["mean"]
            ),
            "corr": 0.0,  # Placeholder - proper correlation would require regridding
        }

    # for model-only run
    if ds_test is not None and ds_ref_regrid is None:
        metrics_dict["test_regrid"] = metrics_dict["test"]

    if ds_diff is not None:
        var_diff = ds_diff[var_key]

        metrics_dict["diff"] = {
            "min": var_diff.min().item(),
            "max": var_diff.max().item(),
            "mean": spatial_avg(ds_diff, var_key),  # type: ignore
        }

    return metrics_dict


def _set_default_metric_values(metrics_dict: MetricsDict) -> MetricsDict:
    """Set the default values for the metrics in case they are not calculated.

    Parameters
    ----------
    metrics_dict : MetricsDict
        The dictionary containing the metrics.

    Returns
    -------
    MetricsDict
        The dictionary containing the metrics with default values.
    """
    var_keys = ["test_regrid", "ref", "ref_regrid", "diff"]
    metric_keys = ["min", "max", "mean", "std"]
    for var_key in var_keys:
        metrics_dict[var_key] = {
            metric_key: METRICS_DEFAULT_VALUE for metric_key in metric_keys
        }

    metrics_dict["misc"] = {
        "rmse": METRICS_DEFAULT_VALUE,
        "corr": METRICS_DEFAULT_VALUE,
    }

    return metrics_dict

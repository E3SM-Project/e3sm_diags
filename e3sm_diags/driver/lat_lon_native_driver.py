from __future__ import annotations

from typing import TYPE_CHECKING, List, Tuple

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
from e3sm_diags.plot.lat_lon_native_plot import plot as plot_func

logger = _setup_child_logger(__name__)

if TYPE_CHECKING:
    from e3sm_diags.parameter.lat_lon_native_parameter import LatLonNativeParameter

# The default value for metrics if it is not calculated. This value was
# preserved from the legacy CDAT codebase because the viewer expects this
# value for metrics that aren't calculated.
# TODO: Support for 3d variables is not tested and need furture development
METRICS_DEFAULT_VALUE = 999.999


def run_diag(parameter: LatLonNativeParameter) -> LatLonNativeParameter:  # noqa: C901
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
            if hasattr(ds_test, "file_path"):
                logger.info(f"Dataset file path: {ds_test.file_path}")
            if hasattr(ds_test, "filepath"):
                logger.info(f"Dataset filepath: {ds_test.filepath}")
            if hasattr(ds_test, "_file_obj") and hasattr(ds_test._file_obj, "name"):
                logger.info(f"Dataset file object name: {ds_test._file_obj.name}")
            try:
                for name, var in ds_test.variables.items():
                    if hasattr(var, "file"):
                        logger.info(f"Variable {name} file: {var.file}")
                        break
            except Exception:
                pass

            # Helper function for variable derivation and transformation
            from e3sm_diags.derivations.derivations import (
                DERIVED_VARIABLES,
                FUNC_NEEDS_TARGET_VAR,
            )

            def process_variable_derivations(dataset, variable_key, dataset_name=""):  # noqa: C901
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
                                    original_units = dataset[variable_key].attrs.get(
                                        "units", "unknown"
                                    )
                                    result = func(dataset[variable_key])
                                    if isinstance(result, xr.DataArray):
                                        dataset[variable_key] = result
                                        new_units = dataset[variable_key].attrs.get(
                                            "units", "unknown"
                                        )
                                        logger.info(
                                            f"Applied conversion to {variable_key}{name_suffix}: {original_units} -> {new_units}"
                                        )
                                        break
                                except Exception as e:
                                    logger.warning(
                                        f"Failed to apply conversion to {variable_key}{name_suffix}: {e}"
                                    )
                else:
                    # If not found directly, attempt to derive it
                    logger.info(
                        f"Variable {variable_key} not found directly{name_suffix}, attempting derivation"
                    )
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
                                        result = func(
                                            *[dataset[v] for v in vars_tuple],
                                            variable_key,
                                        )
                                    else:
                                        result = func(*[dataset[v] for v in vars_tuple])

                                    # Add the derived variable to the dataset
                                    if isinstance(result, xr.DataArray):
                                        dataset[variable_key] = result
                                        logger.info(
                                            f"Successfully derived {variable_key} from {vars_tuple}{name_suffix}"
                                        )
                                        derived = True
                                        break
                                except Exception as e:
                                    logger.warning(
                                        f"Failed to derive {variable_key} from {vars_tuple}{name_suffix}: {e}"
                                    )

                        if not derived:
                            logger.warning(
                                f"Could not derive {variable_key}{name_suffix} - required source variables not available"
                            )
                    else:
                        logger.warning(
                            f"{variable_key} is not a recognized derivable variable"
                        )

                # Verify if variable exists and log possible matches if not
                if variable_key not in dataset:
                    logger.warning(f"Variable {variable_key} not found{name_suffix}!")
                    logger.info(
                        f"Available variables{name_suffix}: {list(dataset.data_vars)}"
                    )

                    # Try to find the variable with a different name or naming convention
                    possible_matches = []
                    for data_var in dataset.data_vars:
                        if (
                            variable_key.lower() in data_var.lower()
                            or data_var.lower() in variable_key.lower()
                        ):
                            possible_matches.append(data_var)

                    if possible_matches:
                        logger.info(
                            f"Possible variable matches{name_suffix}: {possible_matches}"
                        )

                    return False
                return True

            # Load the native grid information for test data
            uxds_test = None
            if parameter.test_grid_file:
                try:
                    logger.info(
                        f"Loading test native grid from: {parameter.test_grid_file}"
                    )
                    # When loading the dataset, include the test data to map it onto the grid
                    uxds_test = ux.open_dataset(
                        parameter.test_grid_file, parameter.test_data_file_path
                    )
                    logger.info(
                        "Successfully loaded test native grid data with uxarray"
                    )

                    # Process variable derivations for test dataset
                    process_variable_derivations(uxds_test, var_key, "test")

                except Exception as e:
                    logger.error(f"Failed to load test native grid: {e}")
                    import traceback

                    logger.error(traceback.format_exc())
                    uxds_test = None

            # Load the native grid information for reference data
            uxds_ref = None
            if ds_ref is not None:
                logger.info("Loading reference data")

                if (
                    not hasattr(parameter, "ref_grid_file")
                    or parameter.ref_grid_file is None
                ):
                    logger.warning(
                        "No ref_grid_file specified in parameter. This is required for model_vs_model native grid visualization."
                    )
                    logger.warning(
                        "Make sure your parameter configuration includes 'ref_grid_file'"
                    )
                    # Try to use test_grid_file as a fallback for models with the same grid
                    if (
                        hasattr(parameter, "test_grid_file")
                        and parameter.test_grid_file
                    ):
                        logger.warning(
                            f"Attempting to use test_grid_file as fallback for reference: {parameter.test_grid_file}"
                        )
                        try:
                            # Attempt to use ref_data_file_path if available
                            if (
                                hasattr(parameter, "ref_data_file_path")
                                and parameter.ref_data_file_path
                            ):
                                logger.info(
                                    f"Using ref_data_file_path with test grid: {parameter.ref_data_file_path}"
                                )
                                uxds_ref = ux.open_dataset(
                                    parameter.test_grid_file,
                                    parameter.ref_data_file_path,
                                )
                            else:
                                # Otherwise fall back to ds_ref
                                uxds_ref = ux.open_dataset(
                                    parameter.test_grid_file, ds_ref
                                )

                            logger.info(
                                "Successfully loaded reference data using test grid file as fallback"
                            )
                            logger.info(f"uxds_ref is now {type(uxds_ref)}")
                            logger.info(
                                f"uxds_ref variables: {list(uxds_ref.data_vars)}"
                            )
                            process_variable_derivations(uxds_ref, var_key, "reference")
                        except Exception as e:
                            logger.error(
                                f"Failed to load reference using test grid as fallback: {e}"
                            )
                            import traceback

                            logger.error(traceback.format_exc())
                            uxds_ref = None
                elif parameter.ref_grid_file:
                    try:
                        logger.info(
                            f"Loading reference native grid from: {parameter.ref_grid_file}"
                        )

                        # Try to use ref_data_file_path if available
                        if (
                            hasattr(parameter, "ref_data_file_path")
                            and parameter.ref_data_file_path
                        ):
                            logger.info(
                                f"Using ref_data_file_path: {parameter.ref_data_file_path}"
                            )
                            try:
                                uxds_ref = ux.open_dataset(
                                    parameter.ref_grid_file,
                                    parameter.ref_data_file_path,
                                )
                                logger.info(
                                    "Successfully loaded reference grid with ref_data_file_path"
                                )
                            except Exception as e:
                                logger.warning(
                                    f"Failed to load reference with ref_data_file_path: {e}"
                                )
                                # Fall back to using ds_ref
                                uxds_ref = ux.open_dataset(
                                    parameter.ref_grid_file, ds_ref
                                )
                                logger.info(
                                    "Fallback: Successfully loaded reference using ds_ref"
                                )
                        else:
                            # When loading the dataset, include the reference data to map it onto the grid
                            uxds_ref = ux.open_dataset(parameter.ref_grid_file, ds_ref)
                            logger.info(
                                "Successfully loaded reference native grid data with uxarray"
                            )
                        logger.info(f"uxds_ref is now {type(uxds_ref)}")

                        # Log the variables in the loaded reference dataset
                        logger.info(
                            f"Reference dataset variables: {list(uxds_ref.data_vars)}"
                        )
                        if var_key in uxds_ref:
                            logger.info(f"uxds_ref[{var_key}] already exists")
                        else:
                            logger.info(
                                f"uxds_ref does not contain {var_key} yet, will try derivation"
                            )

                        # Process variable derivations for reference dataset
                        if process_variable_derivations(uxds_ref, var_key, "reference"):
                            logger.info(
                                f"Successfully derived/found {var_key} in reference dataset"
                            )
                            logger.info(
                                f"uxds_ref[{var_key}] shape: {uxds_ref[var_key].shape}"
                            )
                            logger.info(
                                f"uxds_ref[{var_key}] data type: {uxds_ref[var_key].dtype}"
                            )
                            if hasattr(uxds_ref[var_key], "units"):
                                logger.info(
                                    f"uxds_ref[{var_key}] units: {uxds_ref[var_key].units}"
                                )
                        else:
                            logger.warning(
                                f"Unable to find or derive {var_key} in reference dataset"
                            )

                    except Exception as e:
                        logger.error(f"Failed to load reference native grid: {e}")
                        import traceback

                        logger.error(traceback.format_exc())
                        uxds_ref = None

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
                    logger.warning(
                        "Skipping native grid diagnostics: uxds_test is None"
                    )
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
                        season,
                        regions,
                        var_key,
                        ref_name,
                        uxds_test,
                        uxds_ref,
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


def compute_diff_between_grids(
    uxds_test: ux.dataset.UxDataset, uxds_ref: ux.dataset.UxDataset, var_key: str
) -> ux.dataset.UxDataset:
    """Compute the difference between two native grid datasets.

    This function handles the remapping between different grids if needed,
    and computes the difference between test and reference data.

    Parameters
    ----------
    uxds_test : ux.dataset.UxDataset
        The test dataset on native grid
    uxds_ref : ux.dataset.UxDataset
        The reference dataset on native grid
    var_key : str
        The variable key to compute difference for

    Returns
    -------
    ux.dataset.UxDataset or None
        A dataset containing the difference data, or None if computation fails
    """
    try:
        # Check if variables exist in both datasets
        if var_key not in uxds_test or var_key not in uxds_ref:
            if var_key not in uxds_test:
                logger.error(f"Variable {var_key} not found in test dataset")
            if var_key not in uxds_ref:
                logger.error(f"Variable {var_key} not found in reference dataset")
            return None

        logger.info("Computing difference between test and reference data")

        # Check if both grids are identical by comparing grid properties
        same_grid = False
        try:
            # Get grid sizes using uxgrid.sizes
            test_sizes = uxds_test.uxgrid.sizes
            ref_sizes = uxds_ref.uxgrid.sizes

            # Compare face counts for grid similarity
            test_face_count = test_sizes.get("face", 0)
            ref_face_count = ref_sizes.get("face", 0)

            if test_face_count == ref_face_count and test_face_count > 0:
                same_grid = True
                logger.info(f"Same grid detected with {test_face_count} faces")
            else:
                logger.info(
                    f"Different grids detected: test grid ({test_face_count} faces), reference grid ({ref_face_count} faces)"
                )
        except Exception as e:
            logger.warning(f"Error comparing grids: {e}")

        # Create a diff dataset
        uxds_diff = None

        if same_grid:
            # If grids are the same, directly compute difference
            logger.info("Same grid detected - computing difference directly")
            try:
                # Extract the variable data arrays
                test_var = uxds_test[var_key].squeeze()
                ref_var = uxds_ref[var_key].squeeze()

                # Create a copy of the test dataset to store the difference
                uxds_diff = uxds_test.copy()

                # Compute the difference
                uxds_diff[var_key] = test_var - ref_var
                logger.info("Difference computed successfully using direct subtraction")
            except Exception as e:
                logger.error(f"Error computing direct difference: {e}")
                import traceback

                logger.error(traceback.format_exc())
                return None
        else:
            # Determine which grid has lower resolution using uxgrid.sizes
            # Lower resolution grid typically has fewer faces
            test_sizes = uxds_test.uxgrid.sizes
            ref_sizes = uxds_ref.uxgrid.sizes

            test_face_count = test_sizes.get("face", 0)
            ref_face_count = ref_sizes.get("face", 0)

            # By default, use the test grid as target
            target_is_test = True
            if ref_face_count < test_face_count:
                # Reference grid has lower resolution, use it as target
                target_is_test = False
                logger.info(
                    f"Using reference grid as target (lower resolution: {ref_face_count} vs {test_face_count} faces)"
                )
            else:
                logger.info(
                    f"Using test grid as target (lower resolution: {test_face_count} vs {ref_face_count} faces)"
                )

            try:
                if target_is_test:
                    # Regrid reference data to test grid
                    logger.info("Remapping reference data to test grid")
                    # Use nearest neighbor for simplicity and robustness
                    # Extract variable to remap
                    ref_var = uxds_ref[var_key].squeeze()
                    # Perform remapping to test grid
                    uxds_diff = uxds_test.copy()
                    # Get remapped data using reference variable's remap method
                    logger.info(
                        "Remapping reference data to test grid using nearest_neighbor"
                    )
                    ref_remapped = ref_var.remap.nearest_neighbor(
                        uxds_test.uxgrid, remap_to="face centers"
                    )
                    logger.info(
                        f"Successfully remapped reference data with shape: {ref_remapped.shape}"
                    )

                    # Compute difference
                    uxds_diff[var_key] = uxds_test[var_key].squeeze() - ref_remapped
                else:
                    # Regrid test data to reference grid
                    logger.info("Remapping test data to reference grid")
                    # Extract variable to remap
                    test_var = uxds_test[var_key].squeeze()
                    # Perform remapping to reference grid
                    uxds_diff = uxds_ref.copy()
                    # Get remapped data using test variable's remap method
                    logger.info(
                        "Remapping test data to reference grid using nearest_neighbor"
                    )
                    test_remapped = test_var.remap.nearest_neighbor(
                        uxds_ref.uxgrid, remap_to="face centers"
                    )
                    logger.info(
                        f"Successfully remapped test data with shape: {test_remapped.shape}"
                    )

                    # Compute difference
                    uxds_diff[var_key] = test_remapped - uxds_ref[var_key].squeeze()

                logger.info(
                    "Remapping and difference computation completed successfully"
                )
            except Exception as e:
                logger.error(f"Error during remapping and difference computation: {e}")
                import traceback

                logger.error(traceback.format_exc())
                return None

        # Copy attributes from test variable
        if var_key in uxds_diff and var_key in uxds_test:
            for attr, value in uxds_test[var_key].attrs.items():
                uxds_diff[var_key].attrs[attr] = value

            # For proper visualization, add a metadata attribute indicating this is a difference field
            uxds_diff[var_key].attrs["long_name"] = (
                f"Difference in {uxds_diff[var_key].attrs.get('long_name', var_key)}"
            )

        return uxds_diff

    except Exception as e:
        logger.error(f"Error in compute_diff_between_grids: {e}")
        import traceback

        logger.error(traceback.format_exc())
        return None


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
    # Check if the variable exists in the uxarray dataset
    if var_key not in uxds_test:
        logger.warning(f"Variable {var_key} not found in uxarray dataset!")
        logger.info(f"Available variables: {list(uxds_test.data_vars)}")
        return

    # Check if all values are identical - this would result in a single-color plot
    var_data = uxds_test[var_key].values
    n_unique = len(np.unique(var_data))
    if n_unique <= 1:
        logger.warning(
            f"Variable {var_key} contains only {n_unique} unique value! Will appear as solid color."
        )

    # Check for missing units
    if "units" not in uxds_test[var_key].attrs:
        logger.warning(f"Variable {var_key} has no units attribute")

    for region in regions:
        parameter._set_param_output_attrs(var_key, season, region, ref_name, ilev=None)
        print(f"Processing {var_key} for region {region}")

        # Use the native grid plotting function with explicit uxarray dataset
        plot_func(
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
            _process_test_dataset(
                parameter, ds_test_ilev, ds_land_sea_mask, var_key, region
            )
            parameter._set_param_output_attrs(var_key, season, region, ref_name, ilev)
            print(f"Processing {var_key} for region {region} at level {ilev}")

            # Use the native grid plotting function
            plot_func(
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
    season: str,
    regions: List[str],
    var_key: str,
    ref_name: str,
    uxds_test: ux.dataset.UxDataset = None,
    uxds_ref: ux.dataset.UxDataset = None,
):
    """Run diagnostics on a 2D variable using native grid.

    This function creates plots for each region using the native grid datasets.

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
    uxds_test : ux.dataset.UxDataset, optional
        The uxarray dataset containing the test native grid information.
    uxds_ref : ux.dataset.UxDataset, optional
        The uxarray dataset containing the reference native grid information.
    """
    # Check for valid test and reference data
    if uxds_test is None:
        logger.error(
            "No test uxarray dataset - cannot proceed with native grid visualization"
        )
        return
    elif var_key not in uxds_test:
        logger.error(f"Variable {var_key} not found in test uxarray dataset")
        logger.info(f"Available variables: {list(uxds_test.data_vars)}")
        return

    # Check if we have valid reference data
    has_valid_ref = uxds_ref is not None and var_key in uxds_ref
    if not has_valid_ref:
        logger.warning(
            "Reference data not available, will fall back to model-only mode"
        )

    for region in regions:
        parameter._set_param_output_attrs(var_key, season, region, ref_name, ilev=None)

        if has_valid_ref:
            print(f"Processing {var_key} for region {region} (model vs model)")

            # Calculate the difference data before calling the plotting function
            uxds_diff = compute_diff_between_grids(uxds_test, uxds_ref, var_key)

            # Pass test, reference, and the pre-calculated difference data to the plotting function
            plot_func(
                parameter,
                var_key,
                region,
                uxds_test=uxds_test,
                uxds_ref=uxds_ref,
                uxds_diff=uxds_diff,
            )
        else:
            print(f"Processing {var_key} for region {region} (model-only fallback)")
            # If reference data is missing, fall back to model-only mode
            plot_func(
                parameter,
                var_key,
                region,
                uxds_test=uxds_test,
                uxds_ref=None,  # Force model-only behavior
                uxds_diff=None,
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
            # Call the functions without storing their return values
            # since they're not used in this function
            _subset_and_align_native_datasets(
                parameter,
                ds_test_ilev,
                ds_ref_ilev,
                ds_land_sea_mask,
                var_key,
                region,
                uxds_test,
                uxds_ref,
            )

            parameter._set_param_output_attrs(var_key, season, region, ref_name, ilev)
            print(
                f"Processing {var_key} for region {region} at level {ilev} (model vs obs)"
            )

            # Use the native grid plotting function
            plot_func(
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
            if hasattr(dataset.parameter, "ref_data_file_path"):
                logger.info(
                    f"Reference data file path: {dataset.parameter.ref_data_file_path}"
                )
            else:
                logger.warning("ref_data_file_path not set in parameter")

            # Additional ways to extract the file path for debugging
            file_path = None
            if hasattr(ds_ref, "file_path"):
                file_path = ds_ref.file_path
                logger.info(
                    f"Reference data file path from ds_ref.file_path: {file_path}"
                )
            elif hasattr(ds_ref, "filepath"):
                file_path = ds_ref.filepath
                logger.info(
                    f"Reference data file path from ds_ref.filepath: {file_path}"
                )
            elif hasattr(ds_ref, "_file_obj") and hasattr(ds_ref._file_obj, "name"):
                file_path = ds_ref._file_obj.name
                logger.info(
                    f"Reference data file path from ds_ref._file_obj.name: {file_path}"
                )

            # Store additional path in parameter if found through other methods
            if file_path and not hasattr(dataset.parameter, "ref_data_file_path"):
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

from __future__ import annotations

import traceback
from typing import TYPE_CHECKING, Sequence

import uxarray as ux

from e3sm_diags.derivations.default_regions_xr import REGION_SPECS
from e3sm_diags.driver import METRICS_DEFAULT_VALUE
from e3sm_diags.driver.utils.arithmetic import subtract_dataarrays
from e3sm_diags.driver.utils.dataset_native import NativeDataset
from e3sm_diags.driver.utils.type_annotations import MetricsDict
from e3sm_diags.logger import _setup_child_logger
from e3sm_diags.metrics.metrics import native_correlation, native_rmse
from e3sm_diags.plot.lat_lon_native_plot import plot as plot_func

logger = _setup_child_logger(__name__)

if TYPE_CHECKING:
    from e3sm_diags.driver.utils.type_annotations import TimeSelection
    from e3sm_diags.parameter.lat_lon_native_parameter import LatLonNativeParameter


def run_diag(parameter: LatLonNativeParameter) -> LatLonNativeParameter:  # noqa: C901
    """Get metrics for the lat_lon_native diagnostic set.

    This function loops over each variable, season/time_slice, pressure level (if 3-D),
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
    ref_name = getattr(parameter, "ref_name", "")
    regions = parameter.regions

    # Determine whether to use seasons or time_slices
    if len(parameter.time_slices) > 0:
        time_periods: Sequence["TimeSelection"] = parameter.time_slices
        use_time_slices = True
        logger.info(f"Using time_slices: {time_periods}")
    else:
        time_periods = parameter.seasons
        use_time_slices = False
        logger.info(f"Using seasons: {time_periods}")

    test_ds = NativeDataset(parameter, data_type="test")
    ref_ds = NativeDataset(parameter, data_type="ref")

    for var_key in variables:
        logger.info("Variable: {}".format(var_key))
        parameter.var_id = var_key

        for time_period in time_periods:
            if use_time_slices:
                logger.info(f"Processing time slice: {time_period}")
                parameter._set_time_slice_name_yrs_attrs(
                    test_ds.dataset, ref_ds.dataset, time_period
                )
            else:
                logger.info(f"Processing season: {time_period}")
                parameter._set_name_yrs_attrs(
                    test_ds.dataset, ref_ds.dataset, time_period
                )

            ds_xr_test = test_ds.get_native_dataset(
                var_key, time_period, use_time_slices
            )
            ds_xr_ref = ref_ds.get_native_dataset(
                var_key, time_period, use_time_slices, allow_missing=True
            )

            # Log basic dataset info
            if ds_xr_test is not None:
                logger.debug(f"Test dataset variables: {list(ds_xr_test.variables)}")

            uxds_test_grid = None
            if parameter.test_grid_file:
                logger.info(f"Loading test native grid: {parameter.test_grid_file}")

                uxds_test_grid = test_ds.get_grid_dataset()

                # Apply variable derivations if needed
                test_ds._process_variable_derivations(var_key)

            if ds_xr_ref is not None:
                uxds_ref = ref_ds.get_grid_dataset()

                if uxds_ref is not None:
                    # Apply variable derivations if needed
                    ref_ds._process_variable_derivations(var_key)

                    logger.debug(
                        f"Reference dataset variables: {list(uxds_ref.data_vars)}"
                    )

            if ds_xr_ref is None:
                if uxds_test_grid is not None:
                    _run_diags_2d_model_only(
                        parameter,
                        time_period,
                        regions,
                        var_key,
                        ref_name,
                        uxds_test_grid,
                    )
                else:
                    logger.warning(
                        "Skipping native grid diagnostics: uxds_test is None"
                    )
            else:
                _run_diags_2d(
                    parameter,
                    time_period,
                    regions,
                    var_key,
                    ref_name,
                    uxds_test_grid,
                    uxds_ref,
                )

    return parameter


def _run_diags_2d_model_only(
    parameter: LatLonNativeParameter,
    season: str,
    regions: list[str],
    var_key: str,
    ref_name: str,
    uxds_test: ux.UxDataset,
):
    """Run a model-only diagnostics on a 2D variable using native grid.

    This function plots the native grid data directly using uxarray dataset.

    Parameters
    ----------
    parameter : LatLonNativeParameter
        The parameter object.
    season : str
        The season.
    regions : list[str]
        The list of regions.
    var_key : str
        The key of the variable.
    ref_name : str
        The reference name.
    uxds_test : ux.UxDataset
        The uxarray dataset containing the test native grid information.
    """
    # Process each region
    for region in regions:
        parameter._set_param_output_attrs(var_key, season, region, ref_name, ilev=None)
        logger.info(f"Processing {var_key} for region {region}")

        # Apply regional subsetting before metrics calculation
        uxds_test_subset = _apply_regional_subsetting(uxds_test, var_key, region)

        # Create a minimal metrics_dict for model-only mode
        # Only need test metrics, no reference or diff metrics
        parameter.metrics_dict = _create_metrics_dict(
            var_key,
            uxds_test_subset,
            uxds_ref=None,
            uxds_test_remapped=uxds_test_subset,  # For model-only, test_regrid = test
            uxds_ref_remapped=None,
            uxds_diff=None,
        )

        # Create plot with model-only mode
        plot_func(
            parameter,
            var_key,
            region,
            uxds_test=uxds_test,
            uxds_ref=None,
            uxds_diff=None,
        )


def _run_diags_2d(
    parameter: LatLonNativeParameter,
    season: str,
    regions: list[str],
    var_key: str,
    ref_name: str,
    uxds_test: ux.UxDataset = None,
    uxds_ref: ux.UxDataset = None,
):
    """Run diagnostics on a 2D variable using native grid.

    This function creates plots for each region using the native grid datasets.

    Parameters
    ----------
    parameter : LatLonNativeParameter
        The parameter object.
    season : str
        The season.
    regions : list[str]
        The list of regions.
    var_key : str
        The key of the variable.
    ref_name : str
        The reference name.
    uxds_test : ux.UxDataset, optional
        The uxarray dataset containing the test native grid information.
    uxds_ref : ux.UxDataset, optional
        The uxarray dataset containing the reference native grid information.
    """
    # Check if we have valid reference data
    has_valid_ref = uxds_ref is not None and var_key in uxds_ref

    for region in regions:
        parameter._set_param_output_attrs(var_key, season, region, ref_name, ilev=None)

        if has_valid_ref:
            logger.info(f"Processing {var_key} for region {region} (model vs model)")

            uxds_diff, uxds_test_remapped, uxds_ref_remapped = (
                _compute_diff_between_grids(uxds_test, uxds_ref, var_key)
            )

            # Apply regional subsetting to all datasets before metrics calculation
            uxds_test_subset = _apply_regional_subsetting(uxds_test, var_key, region)
            uxds_ref_subset = _apply_regional_subsetting(uxds_ref, var_key, region)
            uxds_test_remapped_subset = _apply_regional_subsetting(
                uxds_test_remapped, var_key, region
            )
            uxds_ref_remapped_subset = _apply_regional_subsetting(
                uxds_ref_remapped, var_key, region
            )
            uxds_diff_subset = _apply_regional_subsetting(uxds_diff, var_key, region)

            # Create metrics dictionary using regionally subsetted datasets
            metrics_dict = _create_metrics_dict(
                var_key,
                uxds_test_subset,
                uxds_ref_subset,
                uxds_test_remapped_subset,
                uxds_ref_remapped_subset,
                uxds_diff_subset,
            )

            # Store metrics in parameter for plot function to access
            parameter.metrics_dict = metrics_dict

            parameter._set_param_output_attrs(
                var_key, season, region, ref_name, ilev=None
            )

            # Call plot function with original datasets for visualization
            plot_func(
                parameter,
                var_key,
                region,
                uxds_test=uxds_test,
                uxds_ref=uxds_ref,
                uxds_diff=uxds_diff,
            )
        else:
            logger.info(f"Processing {var_key} for region {region} (model-only)")

            # Apply regional subsetting to test dataset before metrics calculation
            uxds_test_subset = _apply_regional_subsetting(uxds_test, var_key, region)

            # Create metrics dictionary for model-only run using regionally subsetted dataset
            metrics_dict = _create_metrics_dict(
                var_key,
                uxds_test_subset,
                None,  # No reference dataset
                None,  # No remapped test dataset (not needed for model-only)
                None,  # No remapped reference dataset
                None,  # No difference dataset
            )

            # Store metrics in parameter for plot function to access
            parameter.metrics_dict = metrics_dict

            parameter._set_param_output_attrs(
                var_key, season, region, ref_name, ilev=None
            )

            # Call plot function with original dataset for visualization
            plot_func(
                parameter,
                var_key,
                region,
                uxds_test=uxds_test,
                uxds_ref=None,
                uxds_diff=None,
            )


def _compute_diff_between_grids(
    uxds_test: ux.UxDataset, uxds_ref: ux.UxDataset, var_key: str
) -> tuple[ux.UxDataset | None, ux.UxDataset, ux.UxDataset]:
    """Compute the difference between two native grid datasets.

    This function handles the remapping between different grids if needed,
    and computes the difference between test and reference data.

    FIXME: This function has too many nested blocks and should be refactored.
    The broad exception handling may hide bugs and makes debugging difficult.

    Parameters
    ----------
    uxds_test : ux.UxDataset
        The test dataset on native grid
    uxds_ref : ux.UxDataset
        The reference dataset on native grid
    var_key : str
        The variable key to compute difference for

    Returns
    -------
    tuple[ux.UxDataset | None, ux.UxDataset, ux.UxDataset]
        A tuple containing (difference_dataset, remapped_test, remapped_ref).
        The difference dataset can be None if computation fails.
    """
    try:
        # Check if variables exist in both datasets
        if var_key not in uxds_test or var_key not in uxds_ref:
            if var_key not in uxds_test:
                logger.error(f"Variable {var_key} not found in test dataset")
            if var_key not in uxds_ref:
                logger.error(f"Variable {var_key} not found in reference dataset")

            return None, uxds_test, uxds_ref

        # Determine if both grids are identical by comparing properties and
        # create a difference dataset accordingly. Otherwise return None.
        same_grid, test_face_count, ref_face_count = _compare_grids(uxds_test, uxds_ref)

        if same_grid:
            uxds_diff = _compute_direct_difference(uxds_test, uxds_ref, var_key)
            # For same grid, no remapping needed
            remapped_test = uxds_test
            remapped_ref = uxds_ref
        else:
            # Determine which grid to use as target (prefer lower resolution grid)
            target_is_test = ref_face_count >= test_face_count

            uxds_diff, remapped_test, remapped_ref = _compute_remapped_difference(
                uxds_test, uxds_ref, var_key, target_is_test
            )

        if uxds_diff is None:
            return None, uxds_test, uxds_ref

        # Copy attributes and add diff metadata
        if var_key in uxds_diff and var_key in uxds_test:
            for attr, value in uxds_test[var_key].attrs.items():
                uxds_diff[var_key].attrs[attr] = value

            # Add metadata indicating this is a difference field
            uxds_diff[var_key].attrs["long_name"] = (
                f"Difference in {uxds_diff[var_key].attrs.get('long_name', var_key)}"
            )

        return uxds_diff, remapped_test, remapped_ref

    except Exception as e:
        logger.error(f"Error in compute_diff_between_grids: {e}")

        return None, uxds_test, uxds_ref


def _compare_grids(
    uxds_test: ux.UxDataset, uxds_ref: ux.UxDataset
) -> tuple[bool, int, int]:
    """Compare two grids to determine if they're identical.

    This function compares the grid properties of the test and reference datasets
    to determine if they are on the same grid.

    Parameters
    ----------
    uxds_test : ux.UxDataset
        The test dataset on native grid.
    uxds_ref : ux.UxDataset
        The reference dataset on native grid.
    Returns
    -------
    tuple[bool, int, int]
        A tuple containing (same_grid, test_face_count, ref_face_count).
    """
    test_sizes = uxds_test.uxgrid.sizes
    ref_sizes = uxds_ref.uxgrid.sizes

    test_face_count = test_sizes.get("face", 0)
    ref_face_count = ref_sizes.get("face", 0)

    same_grid = test_face_count == ref_face_count and test_face_count > 0

    if same_grid:
        logger.debug(f"Same grid detected with {test_face_count} faces")
    else:
        logger.debug(
            f"Different grids: test ({test_face_count} faces), ref ({ref_face_count} faces)"
        )

    return same_grid, test_face_count, ref_face_count


def _compute_direct_difference(
    uxds_test: ux.UxDataset, uxds_ref: ux.UxDataset, var_key: str
) -> ux.UxDataset | None:
    """Compute direct difference when grids are identical.

    This function computes the difference directly without remapping.

    FIXME: This function has too many nested blocks and should be refactored.
    The broad exception handling may hide bugs and makes debugging difficult.

    Parameters
    ----------
    uxds_test : ux.UxDataset
        The test dataset on native grid.
    uxds_ref : ux.UxDataset
        The reference dataset on native grid.
    var_key : str
        The variable key to compute difference for.

    Returns
    -------
    ux.UxDataset or None
        A dataset containing the difference data, or None if computation fails.
    """
    try:
        # Extract the variable data arrays and handle time dimension if present
        test_var = uxds_test[var_key]
        ref_var = uxds_ref[var_key]

        # Handle multiple time points in test data
        if "time" in test_var.dims and test_var.sizes["time"] > 1:
            logger.info(
                f"Test variable {var_key} has multiple time points. Using first time point for difference calculation."
            )
            test_var = test_var.isel(time=0)

        # Handle multiple time points in reference data
        if "time" in ref_var.dims and ref_var.sizes["time"] > 1:
            logger.info(
                f"Reference variable {var_key} has multiple time points. Using first time point for difference calculation."
            )
            ref_var = ref_var.isel(time=0)

        # Squeeze any remaining singleton dimensions
        test_var = test_var.squeeze()
        ref_var = ref_var.squeeze()

        # Create a copy of the test dataset to store the difference
        uxds_diff = uxds_test.copy()

        # Compute the difference
        uxds_diff[var_key] = subtract_dataarrays(test_var, ref_var)
        logger.debug("Difference computed using direct subtraction")
        return uxds_diff

    except Exception as e:
        logger.error(f"Error computing direct difference: {e}")
        logger.debug(
            f"Test var shape: {uxds_test[var_key].shape}, dims: {uxds_test[var_key].dims}"
        )
        logger.debug(
            f"Ref var shape: {uxds_ref[var_key].shape}, dims: {uxds_ref[var_key].dims}"
        )

        return None


def _compute_remapped_difference(
    uxds_test: ux.UxDataset, uxds_ref: ux.UxDataset, var_key: str, target_is_test: bool
) -> tuple[ux.UxDataset | None, ux.UxDataset, ux.UxDataset]:
    """Compute difference with remapping for different grids.

    FIXME: This function has too many nested blocks and should be refactored.
    The broad exception handling may hide bugs and makes debugging difficult.

    Parameters
    ----------
    uxds_test : ux.UxDataset
        The test dataset on native grid.
    uxds_ref : ux.UxDataset
        The reference dataset on native grid.
    var_key : str
        The variable key to compute difference for.
    target_is_test : bool
        If True, remap reference to test grid; otherwise remap test to reference
        grid.

    Returns
    -------
    ux.UxDataset or None
        A dataset containing the difference data, or None if computation fails.
    """
    try:
        # Extract variables and handle time dimension
        test_var = uxds_test[var_key]
        ref_var = uxds_ref[var_key]

        # Handle multiple time points in test data
        if "time" in test_var.dims and test_var.sizes["time"] > 1:
            logger.info(
                f"Test variable {var_key} has multiple time points. Using first time point for remapping."
            )
            test_var = test_var.isel(time=0)

        # Handle multiple time points in reference data
        if "time" in ref_var.dims and ref_var.sizes["time"] > 1:
            logger.info(
                f"Reference variable {var_key} has multiple time points. Using first time point for remapping."
            )
            ref_var = ref_var.isel(time=0)

        # Squeeze any remaining singleton dimensions
        test_var = test_var.squeeze()
        ref_var = ref_var.squeeze()

        if target_is_test:
            # Remap reference to test grid
            logger.info("Remapping reference data to test grid")
            uxds_diff = uxds_test.copy()
            remapped_test = uxds_test

            ref_remapped = ref_var.remap.nearest_neighbor(
                uxds_test.uxgrid, remap_to="face centers"
            )
            uxds_diff[var_key] = test_var - ref_remapped

            # Create remapped reference dataset
            remapped_ref = uxds_test.copy()
            remapped_ref[var_key] = ref_remapped

        else:
            # Remap test to reference grid
            logger.info("Remapping test data to reference grid")
            uxds_diff = uxds_ref.copy()
            remapped_ref = uxds_ref

            test_remapped = test_var.remap.nearest_neighbor(
                uxds_ref.uxgrid, remap_to="face centers"
            )
            uxds_diff[var_key] = test_remapped - ref_var

            # Create remapped test dataset
            remapped_test = uxds_ref.copy()
            remapped_test[var_key] = test_remapped

        return uxds_diff, remapped_test, remapped_ref

    except Exception as e:
        logger.error(f"Error during remapping and difference computation: {e}")
        logger.debug(
            f"Test var shape: {uxds_test[var_key].shape}, dims: {uxds_test[var_key].dims}"
        )
        logger.debug(
            f"Ref var shape: {uxds_ref[var_key].shape}, dims: {uxds_ref[var_key].dims}"
        )
        logger.debug(traceback.format_exc())

        return None, uxds_test, uxds_ref


def _create_metrics_dict(
    var_key: str,
    uxds_test: ux.UxDataset,
    uxds_ref: ux.UxDataset | None,
    uxds_test_remapped: ux.UxDataset | None,
    uxds_ref_remapped: ux.UxDataset | None,
    uxds_diff: ux.UxDataset | None,
) -> MetricsDict:
    """Create a metrics dictionary for native grid datasets.

    This function follows the same pattern as lat_lon_driver._create_metrics_dict
    but uses uxarray datasets and native grid operations.

    Parameters
    ----------
    var_key : str
        The variable key.
    uxds_test : ux.UxDataset
        The original test uxarray dataset.
    uxds_ref : ux.UxDataset | None
        The original reference uxarray dataset.
    uxds_test_remapped : ux.UxDataset | None
        The remapped test uxarray dataset.
    uxds_ref_remapped : ux.UxDataset | None
        The remapped reference uxarray dataset.
    uxds_diff : ux.UxDataset | None
        The difference uxarray dataset.

    Returns
    -------
    MetricsDict
        The metrics dictionary.
    """
    # Basic test metrics using original dataset
    var_test = uxds_test[var_key]
    metrics_dict: MetricsDict = {
        "test": {
            "min": [var_test.min().item()],
            "max": [var_test.max().item()],
            "mean": [var_test.weighted_mean().item()],
            "std": METRICS_DEFAULT_VALUE,  # Not implemented yet for native grids
        },
        "unit": uxds_test[var_key].attrs.get("units", ""),
    }

    # Set default values for all optional metrics
    metrics_dict = _set_default_metric_values(metrics_dict)

    # Add reference metrics if available (using original dataset)
    if uxds_ref is not None and var_key in uxds_ref:
        var_ref = uxds_ref[var_key]
        metrics_dict["ref"] = {
            "min": [var_ref.min().item()],
            "max": [var_ref.max().item()],
            "mean": [var_ref.weighted_mean().item()],
            "std": METRICS_DEFAULT_VALUE,  # Not implemented yet for native grids
        }

    # Add remapped test metrics if available
    if uxds_test_remapped is not None and var_key in uxds_test_remapped:
        var_test_remapped = uxds_test_remapped[var_key]
        metrics_dict["test_regrid"] = {
            "min": [var_test_remapped.min().item()],
            "max": [var_test_remapped.max().item()],
            "mean": [var_test_remapped.weighted_mean().item()],
            "std": METRICS_DEFAULT_VALUE,  # Not implemented yet for native grids
        }

    # Add remapped reference metrics if available
    if uxds_ref_remapped is not None and var_key in uxds_ref_remapped:
        var_ref_remapped = uxds_ref_remapped[var_key]
        metrics_dict["ref_regrid"] = {
            "min": [var_ref_remapped.min().item()],
            "max": [var_ref_remapped.max().item()],
            "mean": [var_ref_remapped.weighted_mean().item()],
            "std": METRICS_DEFAULT_VALUE,  # Not implemented yet for native grids
        }

    # Calculate RMSE and correlation on remapped datasets (following lat_lon pattern)
    if uxds_test_remapped is not None and uxds_ref_remapped is not None:
        try:
            rmse_val = native_rmse(uxds_test_remapped, uxds_ref_remapped, var_key)
            corr_val = native_correlation(
                uxds_test_remapped, uxds_ref_remapped, var_key
            )

            metrics_dict["misc"] = {
                "rmse": [rmse_val],
                "corr": [corr_val],
            }
        except Exception as e:
            logger.warning(f"Failed to calculate RMSE/correlation: {e}")
            # Keep default NaN values for misc metrics

    # For model-only run, copy test metrics to test_regrid
    if uxds_test is not None and uxds_ref_remapped is None:
        metrics_dict["test_regrid"] = metrics_dict["test"]

    # Add difference metrics if available
    if uxds_diff is not None and var_key in uxds_diff:
        var_diff = uxds_diff[var_key]
        metrics_dict["diff"] = {
            "min": [var_diff.min().item()],
            "max": [var_diff.max().item()],
            "mean": [var_diff.weighted_mean().item()],
            "std": METRICS_DEFAULT_VALUE,  # Not implemented yet for native grids
        }

    return metrics_dict


def _set_default_metric_values(metrics_dict: MetricsDict) -> MetricsDict:
    """Set default values for optional metrics in the dictionary.

    This function follows the same pattern as lat_lon_driver._set_default_metric_values.
    """
    var_keys = ["test_regrid", "ref", "ref_regrid", "diff"]
    metric_keys = ["min", "max", "mean", "std"]

    for var_key in var_keys:
        if var_key not in metrics_dict:
            metrics_dict[var_key] = {
                metric_key: METRICS_DEFAULT_VALUE for metric_key in metric_keys
            }

    if "misc" not in metrics_dict:
        metrics_dict["misc"] = {
            "rmse": METRICS_DEFAULT_VALUE,
            "corr": METRICS_DEFAULT_VALUE,
        }

    return metrics_dict


def _apply_regional_subsetting(
    uxds: ux.UxDataset | None, var_key: str, region: str
) -> ux.UxDataset | None:
    """Apply regional subsetting to a uxarray dataset based on region specification.

    This function follows the same pattern as the regional subsetting in
    lat_lon_native_plot.py but moves it to the driver for consistency.

    Parameters
    ----------
    uxds : ux.UxDataset or None
        The uxarray dataset to subset.
    var_key : str
        The variable key to subset.
    region : str
        The region specification (e.g., "global", "CONUS", etc.).

    Returns
    -------
    ux.UxDataset or None
        The regionally subsetted dataset, or None if input was None.
    """
    if uxds is None:
        return uxds

    # Get region specs (same logic as in plot function)
    region_specs = REGION_SPECS.get(region, None)

    if region_specs is None:
        # Unknown region, return original dataset
        logger.warning(
            f"Region '{region}' not found in REGION_SPECS. Using global dataset."
        )
        return uxds

    # Get bounds (same logic as in plot function)
    lat_bounds = region_specs.get("lat", (-90, 90))  # type: ignore
    lon_bounds = region_specs.get("lon", (0, 360))  # type: ignore
    is_global_domain = lat_bounds == (-90, 90) and lon_bounds == (0, 360)

    if is_global_domain:
        # Global domain, no subsetting needed
        return uxds

    try:
        # Check if target variable exists
        if var_key not in uxds.data_vars:
            logger.warning(
                f"Variable '{var_key}' not found in dataset. Available vars: {list(uxds.data_vars)}"
            )
            return uxds

        # Apply subsetting to the specific variable
        var_subset = uxds[var_key].subset.bounding_box(lon_bounds, lat_bounds)

        # Create new dataset from subsetted variable
        uxds_subset = var_subset.to_dataset()
        uxds_subset.attrs.update(uxds.attrs)
        uxds_subset[var_key].attrs.update(uxds[var_key].attrs)
        return uxds_subset

    except Exception as e:
        logger.warning(
            f"Failed to apply regional subsetting for region '{region}': {e}"
        )
        logger.warning("Using global dataset instead.")
        return uxds

from __future__ import annotations

import traceback
from typing import TYPE_CHECKING, Sequence

import uxarray as ux

from e3sm_diags.driver.utils.dataset_native import NativeDataset
from e3sm_diags.logger import _setup_child_logger
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
                parameter._set_time_slice_attrs(
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

                # TODO: What is this supposed to do? Just check if the variables
                # can be derived? It returns a bool but the return value is never
                # used downstream.
                test_ds._process_variable_derivations(var_key)

            if ds_xr_ref is not None:
                uxds_ref = ref_ds.get_grid_dataset()

                if uxds_ref is not None:
                    # TODO: What is this supposed to do? Just check if the variables
                    # can be derived? It returns a bool but the return value is never
                    # used downstream.
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

            uxds_diff = _compute_diff_between_grids(uxds_test, uxds_ref, var_key)

            # Plot with comparison mode
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

            # Plot with model-only mode.
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
) -> ux.UxDataset:
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
    ux.UxDataset or None
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

        # Determine if both grids are identical by comparing properties and
        # create a difference dataset accordingly. Otherwise return None.
        same_grid, test_face_count, ref_face_count = _compare_grids(uxds_test, uxds_ref)

        if same_grid:
            uxds_diff = _compute_direct_difference(uxds_test, uxds_ref, var_key)
        else:
            # Determine which grid to use as target (prefer lower resolution grid)
            target_is_test = ref_face_count >= test_face_count

            uxds_diff = _compute_remapped_difference(
                uxds_test, uxds_ref, var_key, target_is_test
            )

        if uxds_diff is None:
            return None

        # Copy attributes and add diff metadata
        if var_key in uxds_diff and var_key in uxds_test:
            for attr, value in uxds_test[var_key].attrs.items():
                uxds_diff[var_key].attrs[attr] = value

            # Add metadata indicating this is a difference field
            uxds_diff[var_key].attrs["long_name"] = (
                f"Difference in {uxds_diff[var_key].attrs.get('long_name', var_key)}"
            )

        return uxds_diff

    except Exception as e:
        logger.error(f"Error in compute_diff_between_grids: {e}")

        return None


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
        uxds_diff[var_key] = test_var - ref_var
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
) -> ux.UxDataset | None:
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

            ref_remapped = ref_var.remap.nearest_neighbor(
                uxds_test.uxgrid, remap_to="face centers"
            )
            uxds_diff[var_key] = test_var - ref_remapped

        else:
            # Remap test to reference grid
            logger.info("Remapping test data to reference grid")
            uxds_diff = uxds_ref.copy()

            test_remapped = test_var.remap.nearest_neighbor(
                uxds_ref.uxgrid, remap_to="face centers"
            )
            uxds_diff[var_key] = test_remapped - ref_var

        return uxds_diff

    except Exception as e:
        logger.error(f"Error during remapping and difference computation: {e}")
        logger.debug(
            f"Test var shape: {uxds_test[var_key].shape}, dims: {uxds_test[var_key].dims}"
        )
        logger.debug(
            f"Ref var shape: {uxds_ref[var_key].shape}, dims: {uxds_ref[var_key].dims}"
        )
        logger.debug(traceback.format_exc())

        return None

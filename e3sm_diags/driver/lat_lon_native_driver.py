from __future__ import annotations

from typing import TYPE_CHECKING, List

import uxarray as ux
import xarray as xr

from e3sm_diags.driver.utils.climo_xr import ClimoFreq
from e3sm_diags.driver.utils.dataset_xr import Dataset
from e3sm_diags.driver.utils.regrid import (
    _apply_land_sea_mask,
    _subset_on_region,
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

            # Log basic dataset info
            logger.debug(f"Test dataset variables: {list(ds_test.variables)}")

            # Helper function for variable derivation and transformation
            from e3sm_diags.derivations.derivations import (
                DERIVED_VARIABLES,
                FUNC_NEEDS_TARGET_VAR,
            )

            def _get_matching_src_vars(dataset, target_var_map):
                """Get matching source variables following dataset_xr pattern."""
                for src_vars, func in target_var_map.items():
                    if all(v in dataset for v in src_vars):
                        return {src_vars: func}
                return None

            def _apply_derivation_func(dataset, func, src_var_keys, target_var_key):
                """Apply derivation function following dataset_xr pattern."""
                func_args = [dataset[var] for var in src_var_keys]

                if func in FUNC_NEEDS_TARGET_VAR:
                    func_args = [target_var_key] + func_args

                derived_var = func(*func_args)
                dataset[target_var_key] = derived_var
                return dataset

            def process_variable_derivations(dataset, var_key, dataset_name=""):
                """Process variable derivations following dataset_xr approach."""
                name_suffix = f" in {dataset_name} dataset" if dataset_name else ""

                # Follow dataset_xr._get_climo_dataset logic:
                # 1. If var is in derived_vars_map, try to derive it
                if var_key in DERIVED_VARIABLES:
                    target_var_map = DERIVED_VARIABLES[var_key]
                    matching_target_var_map = _get_matching_src_vars(
                        dataset, target_var_map
                    )

                    if matching_target_var_map is not None:
                        # Get derivation function and source variables
                        derivation_func = list(matching_target_var_map.values())[0]
                        src_var_keys = list(matching_target_var_map.keys())[0]

                        logger.info(
                            f"Deriving {var_key}{name_suffix} using source variables: {src_var_keys}"
                        )

                        try:
                            _apply_derivation_func(
                                dataset, derivation_func, src_var_keys, var_key
                            )
                            return True
                        except Exception as e:
                            logger.warning(
                                f"Failed to derive {var_key}{name_suffix}: {e}"
                            )
                            # Fall through to check if variable exists directly

                # 2. Check if variable exists directly in dataset
                if var_key in dataset.data_vars:
                    return True

                # 3. Variable not found and couldn't be derived
                logger.warning(
                    f"Variable {var_key} not found{name_suffix} and could not be derived"
                )
                return False

            # Load the native grid information for test data
            uxds_test = None
            if parameter.test_grid_file:
                try:
                    logger.info(f"Loading test native grid: {parameter.test_grid_file}")
                    uxds_test = ux.open_dataset(
                        parameter.test_grid_file, parameter.test_data_file_path
                    )

                    # Process variable derivations for test dataset
                    process_variable_derivations(uxds_test, var_key, "test")

                except Exception as e:
                    logger.error(f"Failed to load test native grid: {e}")
                    import traceback

                    logger.debug(traceback.format_exc())
                    uxds_test = None

            # Load the native grid information for reference data
            uxds_ref = None
            if ds_ref is not None:
                try:
                    has_ref_grid = (
                        hasattr(parameter, "ref_grid_file")
                        and parameter.ref_grid_file is not None
                    )

                    if not has_ref_grid:
                        logger.info(
                            "No ref_grid_file specified. Skipping reference grid loading."
                        )
                        uxds_ref = None
                        continue
                    else:
                        grid_file = parameter.ref_grid_file

                    # Use ref_data_file_path if available, otherwise use ds_ref
                    if (
                        hasattr(parameter, "ref_data_file_path")
                        and parameter.ref_data_file_path
                    ):
                        data_source = parameter.ref_data_file_path
                    else:
                        data_source = ds_ref  # type: ignore

                    # Load the dataset with uxarray
                    uxds_ref = ux.open_dataset(grid_file, data_source)

                    # Process variable derivations
                    process_variable_derivations(uxds_ref, var_key, "reference")

                    # Debug variable info
                    logger.debug(
                        f"Reference dataset variables: {list(uxds_ref.data_vars)}"
                    )

                except Exception as e:
                    logger.error(f"Failed to load reference native grid: {e}")
                    import traceback

                    logger.debug(traceback.format_exc())
                    uxds_ref = None

            if ds_ref is None:
                # Only handle 2D variables for now
                if uxds_test is not None:
                    _run_diags_2d_model_only(
                        parameter,
                        season,
                        regions,
                        var_key,
                        ref_name,
                        uxds_test,
                    )
                else:
                    logger.warning(
                        "Skipping native grid diagnostics: uxds_test is None"
                    )
            else:
                # Only handle 2D variables for now
                _run_diags_2d(
                    parameter,
                    season,
                    regions,
                    var_key,
                    ref_name,
                    uxds_test,
                    uxds_ref,
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

        # Determine if both grids are identical by comparing properties
        same_grid, test_face_count, ref_face_count = _compare_grids(uxds_test, uxds_ref)

        # Create difference dataset
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


def _compare_grids(uxds_test, uxds_ref):
    """Compare two grids to determine if they're identical."""
    try:
        # Get grid sizes using uxgrid.sizes
        test_sizes = uxds_test.uxgrid.sizes
        ref_sizes = uxds_ref.uxgrid.sizes

        # Compare face counts for grid similarity
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

    except Exception as e:
        logger.warning(f"Error comparing grids: {e}")
        return False, 0, 0


def _compute_direct_difference(uxds_test, uxds_ref, var_key):
    """Compute direct difference when grids are identical."""
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


def _compute_remapped_difference(uxds_test, uxds_ref, var_key, target_is_test):
    """Compute difference with remapping for different grids."""
    import traceback

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
    # Check if we have valid reference data
    has_valid_ref = uxds_ref is not None and var_key in uxds_ref

    for region in regions:
        parameter._set_param_output_attrs(var_key, season, region, ref_name, ilev=None)

        if has_valid_ref:
            logger.info(f"Processing {var_key} for region {region} (model vs model)")

            # Calculate the difference data before calling the plotting function
            uxds_diff = compute_diff_between_grids(uxds_test, uxds_ref, var_key)

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

            # Plot with model-only mode
            plot_func(
                parameter,
                var_key,
                region,
                uxds_test=uxds_test,
                uxds_ref=None,
                uxds_diff=None,
            )


def _get_ref_climo_dataset(
    dataset: Dataset, var_key: str, season: ClimoFreq
) -> xr.Dataset | None:
    """Get the reference climatology dataset for the variable and season.

    If the reference climatology does not exist or could not be found, it
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

            # Try to get file_path from different possible sources and store it in parameter
            file_path = None
            if hasattr(ds_ref, "file_path"):
                file_path = ds_ref.file_path
            elif hasattr(ds_ref, "filepath"):
                file_path = ds_ref.filepath
            elif hasattr(ds_ref, "_file_obj") and hasattr(ds_ref._file_obj, "name"):
                file_path = ds_ref._file_obj.name

            # Store path in parameter if found and not already set
            if file_path and not hasattr(dataset.parameter, "ref_data_file_path"):
                dataset.parameter.ref_data_file_path = file_path

            return ds_ref

        except (RuntimeError, IOError) as e:
            logger.info(f"Cannot process reference data: {e}. Using model-only mode.")
            return None
    else:
        raise RuntimeError(
            "`_get_ref_climo_dataset` only works with "
            f"`dataset.data_type == 'ref'`, not {dataset.data_type}."
        )


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
            metrics_dict["grid_info"] = grid_info  # type: ignore
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
                metrics_dict["test_regrid"]["mean"] - metrics_dict["ref_regrid"]["mean"]  # type: ignore
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

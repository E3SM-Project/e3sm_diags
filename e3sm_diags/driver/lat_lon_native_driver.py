from __future__ import annotations

from typing import TYPE_CHECKING, List, Sequence

import uxarray as ux
import xarray as xr

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
    from e3sm_diags.parameter.lat_lon_native_parameter import (
        LatLonNativeParameter,
        TimeSelection,
    )

# The default value for metrics if it is not calculated. This value was
# preserved from the legacy CDAT codebase because the viewer expects this
# value for metrics that aren't calculated.
# TODO: Support for 3d variables is not tested and need furture development
METRICS_DEFAULT_VALUE = 999.999


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
        using_time_slices = True
        logger.info(f"Using time_slices: {time_periods}")
    else:
        time_periods = parameter.seasons
        using_time_slices = False
        logger.info(f"Using seasons: {time_periods}")

    # Variables storing xarray `Dataset` objects start with `ds_` and
    # variables storing e3sm_diags `Dataset` objects end with `_ds`. This
    # is to help distinguish both objects from each other.
    test_ds = Dataset(parameter, data_type="test")
    ref_ds = Dataset(parameter, data_type="ref")

    for var_key in variables:
        logger.info("Variable: {}".format(var_key))
        parameter.var_id = var_key

        for time_period in time_periods:
            if using_time_slices:
                # For time_slices, we need to pass the slice info differently
                logger.info(f"Processing time slice: {time_period}")
                # Set up period-specific attributes for file naming
                parameter._set_time_slice_attrs(test_ds, ref_ds, time_period)
            else:
                # For seasons, use existing logic
                logger.info(f"Processing season: {time_period}")
                parameter._set_name_yrs_attrs(test_ds, ref_ds, time_period)

            # Use the same function for both cases - it will handle the logic internally
            ds_test = _get_native_dataset(
                test_ds, var_key, time_period, using_time_slices
            )
            ds_ref = _get_native_dataset(
                ref_ds, var_key, time_period, using_time_slices, allow_missing=True
            )

            # Log basic dataset info
            if ds_test is not None:
                logger.debug(f"Test dataset variables: {list(ds_test.variables)}")

            # Load the native grid information for test data
            uxds_test = None
            if parameter.test_grid_file:
                try:
                    logger.info(f"Loading test native grid: {parameter.test_grid_file}")
                    uxds_test = ux.open_dataset(
                        parameter.test_grid_file, parameter.test_data_file_path
                    )

                    # Process variable derivations for test dataset
                    _process_variable_derivations(uxds_test, var_key, "test")

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
                    _process_variable_derivations(uxds_ref, var_key, "reference")

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
                        time_period,
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
                    time_period,
                    regions,
                    var_key,
                    ref_name,
                    uxds_test,
                    uxds_ref,
                )

    return parameter


def _get_native_dataset(
    dataset: Dataset,
    var_key: str,
    season: "TimeSelection",
    is_time_slice: bool = False,
    allow_missing: bool = False,
) -> xr.Dataset | None:
    """Get the climatology dataset for the variable and season for native grid processing.

    This function handles both test and reference datasets. For reference datasets,
    if the data cannot be found and allow_missing=True, it will return None to
    enable model-only runs.

    This function also stores the data file path in the parameter object
    for native grid visualization.

    Parameters
    ----------
    dataset : Dataset
        The dataset object (test or reference).
    var_key : str
        The key of the variable.
    season : TimeSelection
        The climatology frequency or time slice string.
    is_time_slice : bool, optional
        If True, treat season as a time slice string rather than climatology frequency.
        Default is False.
    allow_missing : bool, optional
        If True, return None when dataset cannot be loaded instead of raising
        an exception. This enables model-only runs when reference data is missing.
        Default is False.

    Returns
    -------
    xr.Dataset | None
        The climatology dataset if it exists, or None if allow_missing=True
        and the dataset cannot be loaded.

    Raises
    ------
    RuntimeError, IOError
        If the dataset cannot be loaded and allow_missing=False.
    """
    try:
        if is_time_slice:
            # For time slices, get the full dataset without averaging
            ds = _get_full_native_dataset(dataset, var_key)
            # Apply the time slice
            ds = _apply_time_slice(ds, season)
        else:
            # Standard climatology processing
            from e3sm_diags.driver.utils.climo_xr import CLIMO_FREQS

            if season in CLIMO_FREQS:
                ds = dataset.get_climo_dataset(var_key, season)  # type: ignore
            else:
                raise ValueError(f"Invalid season for climatology: {season}")

        # Store file path in parameter for native grid processing
        if is_time_slice:
            # For time slices, we know the exact file path we used
            filepath = dataset._get_climo_filepath_with_params()
            if filepath:
                if dataset.data_type == "test":
                    dataset.parameter.test_data_file_path = filepath
                elif dataset.data_type == "ref":
                    dataset.parameter.ref_data_file_path = filepath
        # Note: For climatology case, get_climo_dataset() already handles file path storage

        return ds

    except (RuntimeError, IOError) as e:
        if allow_missing:
            logger.info(
                f"Cannot process {dataset.data_type} data: {e}. Using model-only mode."
            )
            return None
        else:
            raise


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
        return False


def _get_full_native_dataset(dataset: Dataset, var_key: str) -> xr.Dataset:
    """Get the full native dataset without any time averaging.

    This function uses the dataset's file path parameters to directly open
    the raw data file for time slicing operations.

    Parameters
    ----------
    dataset : Dataset
        The dataset object (test or reference).
    var_key : str
        The key of the variable.

    Returns
    -------
    xr.Dataset
        The full dataset with all time steps.

    Raises
    ------
    RuntimeError
        If unable to get the full dataset.
    """
    import os

    import xarray as xr

    # Get the file path using the parameter-based method
    filepath = dataset._get_climo_filepath_with_params()

    if filepath is None:
        raise RuntimeError(
            f"Unable to get file path for {dataset.data_type} dataset. "
            f"For time slicing, please ensure that "
            f"{'ref_file' if dataset.data_type == 'ref' else 'test_file'} parameter is set."
        )

    if not os.path.exists(filepath):
        raise RuntimeError(f"File not found: {filepath}")

    logger.info(f"Opening full native dataset from: {filepath}")

    try:
        # Open the dataset directly without any averaging
        ds = xr.open_dataset(filepath, decode_times=True)
        logger.info(
            f"Successfully opened dataset with time dimension size: {ds.sizes.get('time', 'N/A')}"
        )
        return ds

    except Exception as e:
        raise RuntimeError(f"Failed to open dataset {filepath}: {e}") from e


def _apply_time_slice(ds: xr.Dataset, time_slice: str) -> xr.Dataset:
    """Apply time slice selection to a dataset.

    Parameters
    ----------
    ds : xr.Dataset
        The input dataset with time dimension.
    time_slice : str
        The time slice specification (e.g., "0:10:2", "5:15", "7").

    Returns
    -------
    xr.Dataset
        The dataset with time slice applied.
    """

    # Parse the time slice string
    time_dim = None
    for dim in ds.dims:
        if str(dim).lower() in ["time", "t"]:
            time_dim = dim
            break

    if time_dim is None:
        logger.warning(
            "No time dimension found in dataset. Returning original dataset."
        )
        return ds

    # Parse slice notation
    if ":" in time_slice:
        # Handle slice notation like "0:10:2" or "5:15" or ":10" or "5:" or "::2"
        parts = time_slice.split(":")

        start = int(parts[0]) if parts[0] else None
        end = int(parts[1]) if len(parts) > 1 and parts[1] else None
        step = int(parts[2]) if len(parts) > 2 and parts[2] else None

        # Apply the slice
        ds_sliced = ds.isel({time_dim: slice(start, end, step)})
    else:
        # Single index
        index = int(time_slice)
        ds_sliced = ds.isel({time_dim: index})

    logger.info(
        f"Applied time slice '{time_slice}' to dataset. "
        f"Original time length: {ds.sizes[time_dim]}, "
        f"Sliced time length: {ds_sliced.sizes.get(time_dim, 1)}"
    )

    return ds_sliced


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


def _get_matching_src_vars(dataset, target_var_map):
    """Get matching source variables following dataset_xr pattern."""
    for src_vars, func in target_var_map.items():
        if all(v in dataset for v in src_vars):
            return {src_vars: func}
    return None


def _apply_derivation_func(dataset, func, src_var_keys, target_var_key):
    """Apply derivation function following dataset_xr pattern."""
    from e3sm_diags.derivations.derivations import FUNC_NEEDS_TARGET_VAR

    func_args = [dataset[var] for var in src_var_keys]

    if func in FUNC_NEEDS_TARGET_VAR:
        func_args = [target_var_key] + func_args

    derived_var = func(*func_args)
    dataset[target_var_key] = derived_var
    return dataset


def _process_variable_derivations(dataset, var_key, dataset_name=""):
    """Process variable derivations following dataset_xr approach."""
    from e3sm_diags.derivations.derivations import DERIVED_VARIABLES

    name_suffix = f" in {dataset_name} dataset" if dataset_name else ""

    # Follow dataset_xr._get_climo_dataset logic:
    # 1. If var is in derived_vars_map, try to derive it
    if var_key in DERIVED_VARIABLES:
        target_var_map = DERIVED_VARIABLES[var_key]
        matching_target_var_map = _get_matching_src_vars(dataset, target_var_map)

        if matching_target_var_map is not None:
            # Get derivation function and source variables
            derivation_func = list(matching_target_var_map.values())[0]
            src_var_keys = list(matching_target_var_map.keys())[0]

            logger.info(
                f"Deriving {var_key}{name_suffix} using source variables: {src_var_keys}"
            )

            try:
                _apply_derivation_func(dataset, derivation_func, src_var_keys, var_key)
                return True
            except Exception as e:
                logger.warning(f"Failed to derive {var_key}{name_suffix}: {e}")
                # Fall through to check if variable exists directly

    # 2. Check if variable exists directly in dataset
    if var_key in dataset.data_vars:
        return True

    # 3. Variable not found and couldn't be derived
    logger.warning(
        f"Variable {var_key} not found{name_suffix} and could not be derived"
    )
    return False

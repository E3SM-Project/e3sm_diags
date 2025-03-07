from __future__ import annotations

from typing import TYPE_CHECKING, List, Tuple

import xarray as xr

from e3sm_diags.driver.utils.climo_xr import ClimoFreq
from e3sm_diags.driver.utils.dataset_xr import Dataset
from e3sm_diags.driver.utils.io import _save_data_metrics_and_plots
from e3sm_diags.driver.utils.regrid import (
    _apply_land_sea_mask,
    _subset_on_region,
    get_z_axis,
    has_z_axis,
    regrid_z_axis_to_plevs,
    subset_and_align_datasets,
)
from e3sm_diags.driver.utils.type_annotations import MetricsDict
from e3sm_diags.logger import custom_logger
from e3sm_diags.metrics.metrics import correlation, rmse, spatial_avg, std
from e3sm_diags.plot.lat_lon_plot import plot as plot_func

logger = custom_logger(__name__)

if TYPE_CHECKING:
    from e3sm_diags.parameter.core_parameter import CoreParameter

# The default value for metrics if it is not calculated. This value was
# preserved from the legacy CDAT codebase because the viewer expects this
# value for metrics that aren't calculated.
# TODO: Update `lat_lon_viewer.py` to handle missing metrics with None value.
METRICS_DEFAULT_VALUE = 999.999


def run_diag(parameter: CoreParameter) -> CoreParameter:
    """Get metrics for the lat_lon diagnostic set.

    This function loops over each variable, season, pressure level (if 3-D),
    and region.

    Parameters
    ----------
    parameter : CoreParameter
        The parameter for the diagnostic.

    Returns
    -------
    CoreParameter
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

            if ds_ref is None:
                is_vars_3d = has_z_axis(ds_test[var_key])

                if not is_vars_3d:
                    _run_diags_2d_model_only(
                        parameter,
                        ds_test,
                        ds_land_sea_mask,
                        season,
                        regions,
                        var_key,
                        ref_name,
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
                    )

    return parameter


def _run_diags_2d_model_only(
    parameter: CoreParameter,
    ds_test: xr.Dataset,
    ds_land_sea_mask: xr.Dataset,
    season: str,
    regions: List[str],
    var_key: str,
    ref_name: str,
):
    """Run a model-only diagnostics on a 2D variable.

    This function gets the variable's metrics by region, then saves the
    metrics, metric plots, and data (optional, `CoreParameter.save_netcdf`).

    Parameters
    ----------
    parameter : CoreParameter
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
    """
    for region in regions:
        ds_test_region = _process_test_dataset(
            parameter, ds_test, ds_land_sea_mask, var_key, region
        )

        metrics_dict = _create_metrics_dict(
            var_key,
            ds_test_region,
            None,
            None,
            None,
            None,
        )

        parameter._set_param_output_attrs(var_key, season, region, ref_name, ilev=None)
        _save_data_metrics_and_plots(
            parameter,
            plot_func,
            var_key,
            ds_test_region,
            None,
            None,
            metrics_dict,
        )


def _run_diags_3d_model_only(
    parameter: CoreParameter,
    ds_test: xr.Dataset,
    ds_land_sea_mask: xr.Dataset,
    season: str,
    regions: List[str],
    var_key: str,
    ref_name: str,
):
    """Run a model-only diagnostics on a 3D variable.

    This function gets the variable's metrics by region, then saves the
    metrics, metric plots, and data (optional, `CoreParameter.save_netcdf`).

    Parameters
    ----------
    parameter : CoreParameter
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
            metrics_dict = _create_metrics_dict(
                var_key,
                ds_test_region,
                None,
                None,
                None,
                None,
            )
            parameter._set_param_output_attrs(var_key, season, region, ref_name, ilev)
            _save_data_metrics_and_plots(
                parameter,
                plot_func,
                var_key,
                ds_test_region,
                None,
                None,
                metrics_dict,
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
    parameter: CoreParameter,
    ds_test: xr.Dataset,
    ds_ref: xr.Dataset,
    ds_land_sea_mask: xr.Dataset,
    season: str,
    regions: List[str],
    var_key: str,
    ref_name: str,
):
    """Run diagnostics on a 2D variable.

    This function gets the variable's metrics by region, then saves the
    metrics, metric plots, and data (optional, `CoreParameter.save_netcdf`).

    Parameters
    ----------
    parameter : CoreParameter
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
    """
    for region in regions:
        (
            ds_test_region,
            ds_test_region_regrid,
            ds_ref_region,
            ds_ref_region_regrid,
            ds_diff_region,
        ) = subset_and_align_datasets(
            parameter,
            ds_test,
            ds_ref,
            ds_land_sea_mask,
            var_key,
            region,
        )

        metrics_dict = _create_metrics_dict(
            var_key,
            ds_test_region,
            ds_test_region_regrid,
            ds_ref_region,
            ds_ref_region_regrid,
            ds_diff_region,
        )

        parameter._set_param_output_attrs(var_key, season, region, ref_name, ilev=None)
        _save_data_metrics_and_plots(
            parameter,
            plot_func,
            var_key,
            ds_test_region,
            ds_ref_region,
            ds_diff_region,
            metrics_dict,
        )


def _run_diags_3d(
    parameter: CoreParameter,
    ds_test: xr.Dataset,
    ds_ref: xr.Dataset,
    ds_land_sea_mask: xr.Dataset,
    season: str,
    regions: List[str],
    var_key: str,
    ref_name: str,
):
    """Run diagnostics on a 3D variable.

    This function gets the variable's metrics by region, then saves the
    metrics, metric plots, and data (optional, `CoreParameter.save_netcdf`).

    Parameters
    ----------
    parameter : CoreParameter
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
            ) = subset_and_align_datasets(
                parameter,
                ds_test_ilev,
                ds_ref_ilev,
                ds_land_sea_mask,
                var_key,
                region,
            )

            metrics_dict = _create_metrics_dict(
                var_key,
                ds_test_region,
                ds_test_region_regrid,
                ds_ref_region,
                ds_ref_region_regrid,
                ds_diff_region,
            )

            parameter._set_param_output_attrs(var_key, season, region, ref_name, ilev)
            _save_data_metrics_and_plots(
                parameter,
                plot_func,
                var_key,
                ds_test_region,
                ds_ref_region,
                ds_diff_region,
                metrics_dict,
            )


def _get_ref_climo_dataset(
    dataset: Dataset, var_key: str, season: ClimoFreq
) -> xr.Dataset | None:
    """Get the reference climatology dataset for the variable and season.

    If the reference climatatology does not exist or could not be found, it
    will be considered a model-only run and return `None`.

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
            ds_ref = dataset.get_climo_dataset(var_key, season)
        except (RuntimeError, IOError):
            ds_ref = None

            logger.info("Cannot process reference data, analyzing test data only.")
    else:
        raise RuntimeError(
            "`Dataset._get_ref_dataset` only works with "
            f"`parameter.data_type == 'ref'`, not {dataset.data_type}."
        )

    return ds_ref


def _process_test_dataset(
    parameter: CoreParameter,
    ds_test: xr.Dataset,
    ds_land_sea_mask: xr.Dataset,
    var_key: str,
    region: str,
) -> xr.Dataset:
    """Process the test dataset for the given region.

    Parameters
    ----------
    parameter : CoreParameter
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
) -> MetricsDict:
    """Calculate metrics using the variable in the datasets.

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

        metrics_dict["misc"] = {
            "rmse": rmse(ds_test_regrid, ds_ref_regrid, var_key),
            "corr": correlation(ds_test_regrid, ds_ref_regrid, var_key),
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

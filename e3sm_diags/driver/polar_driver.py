from __future__ import annotations

from typing import TYPE_CHECKING

import xarray as xr

from e3sm_diags.driver.utils.dataset_xr import Dataset
from e3sm_diags.driver.utils.io import _save_data_metrics_and_plots
from e3sm_diags.driver.utils.regrid import (
    get_z_axis,
    has_z_axis,
    regrid_z_axis_to_plevs,
    subset_and_align_datasets,
)
from e3sm_diags.driver.utils.type_annotations import MetricsDict
from e3sm_diags.logger import _setup_child_logger
from e3sm_diags.metrics.metrics import correlation, rmse, spatial_avg
from e3sm_diags.plot.polar_plot import plot as plot_func

logger = _setup_child_logger(__name__)

if TYPE_CHECKING:
    from e3sm_diags.parameter.core_parameter import CoreParameter


def run_diag(parameter: CoreParameter) -> CoreParameter:
    variables = parameter.variables
    seasons = parameter.seasons
    ref_name = getattr(parameter, "ref_name", "")
    regions = parameter.regions

    test_ds = Dataset(parameter, data_type="test")
    ref_ds = Dataset(parameter, data_type="ref")

    for var_key in variables:
        logger.info("Variable: {}".format(var_key))
        parameter.var_id = var_key

        for season in seasons:
            parameter._set_name_yrs_attrs(test_ds, ref_ds, season)

            # Get land/ocean fraction for masking.
            ds_land_sea_mask: xr.Dataset = test_ds._get_land_sea_mask(season)

            ds_test = test_ds.get_climo_dataset(var_key, season)
            ds_ref = ref_ds.get_climo_dataset(var_key, season)

            # Store the variable's DataArray objects for reuse.
            dv_test = ds_test[var_key]
            dv_ref = ds_ref[var_key]

            is_vars_3d = has_z_axis(dv_test) and has_z_axis(dv_ref)
            is_dims_diff = has_z_axis(dv_test) != has_z_axis(dv_ref)

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


def _run_diags_2d(
    parameter: CoreParameter,
    ds_test: xr.Dataset,
    ds_ref: xr.Dataset,
    ds_land_sea_mask: xr.Dataset,
    season: str,
    regions: list[str],
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
        The dataset containing the ref variable. If this is a model-only run
        then it will be the same dataset as ``ds_test``.
    ds_land_sea_mask : xr.Dataset
        The land sea mask dataset, which is only used for masking if the region
        is "land" or "ocean".
    season : str
        The season.
    regions : list[str]
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
            ds_test_regridded=ds_test_region_regrid,
            ds_ref_regridded=ds_ref_region_regrid,
        )


def _run_diags_3d(
    parameter: CoreParameter,
    ds_test: xr.Dataset,
    ds_ref: xr.Dataset,
    ds_land_sea_mask: xr.Dataset,
    season: str,
    regions: list[str],
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
        The dataset containing the ref variable. If this is a model-only run
        then it will be the same dataset as ``ds_test``.
    ds_land_sea_mask : xr.Dataset
        The land sea mask dataset, which is only used for masking if the region
        is "land" or "ocean".
    season : str
        The season.
    regions : list[str]
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
                ds_test_regridded=ds_test_region_regrid,
                ds_ref_regridded=ds_ref_region_regrid,
            )


def _create_metrics_dict(
    var_key: str,
    ds_test: xr.Dataset,
    ds_test_regrid: xr.Dataset,
    ds_ref: xr.Dataset,
    ds_ref_regrid: xr.Dataset,
    ds_diff: xr.Dataset,
) -> MetricsDict:
    """Calculate metrics using the variable in the datasets.

    Metrics include min value, max value, spatial average (mean), standard
    deviation, correlation (pearson_r), and RMSE.

    Parameters
    ----------
    var_key : str
        The variable key.
    ds_test : xr.Dataset
        The test dataset.
    ds_test_regrid : xr.Dataset
        The regridded test Dataset.
    ds_ref : xr.Dataset
        The reference dataset.
    ds_ref_regrid : xr.Dataset
        The regridded reference dataset.
    ds_diff : xr. Dataset
        The difference between ``ds_test_regrid`` and ``ds_ref_regrid``.

    Returns
    -------
    Metrics
        A dictionary with the key being a string and the value being either
        a sub-dictionary (key is metric and value is float) or a string
        ("unit").
    """
    metrics_dict = {}

    metrics_dict["units"] = ds_test[var_key].attrs["units"]
    metrics_dict["ref"] = {
        "min": ds_ref[var_key].min().item(),
        "max": ds_ref[var_key].max().item(),
        "mean": spatial_avg(ds_ref, var_key, axis=["X", "Y"]),
    }
    metrics_dict["test"] = {
        "min": ds_test[var_key].min().item(),
        "max": ds_test[var_key].max().item(),
        "mean": spatial_avg(ds_test, var_key, axis=["X", "Y"]),
    }

    metrics_dict["diff"] = {
        "min": ds_diff[var_key].min().item(),
        "max": ds_diff[var_key].max().item(),
        "mean": spatial_avg(ds_diff, var_key, axis=["X", "Y"]),
    }

    metrics_dict["misc"] = {
        "rmse": rmse(ds_test_regrid, ds_ref_regrid, var_key, axis=["X", "Y"]),
        "corr": correlation(ds_test_regrid, ds_ref_regrid, var_key, axis=["X", "Y"]),
    }

    return metrics_dict

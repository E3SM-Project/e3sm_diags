import copy
from typing import List, Tuple

import xarray as xr

from e3sm_diags.driver.utils.dataset_xr import Dataset
from e3sm_diags.driver.utils.io import _save_data_metrics_and_plots
from e3sm_diags.driver.utils.regrid import (
    align_grids_to_lower_res,
    has_z_axis,
    regrid_z_axis_to_plevs,
    subset_and_align_datasets,
)
from e3sm_diags.driver.utils.type_annotations import MetricsDict
from e3sm_diags.logger import custom_logger
from e3sm_diags.metrics.metrics import correlation, rmse
from e3sm_diags.parameter.zonal_mean_2d_parameter import (
    DEFAULT_PLEVS,
    ZonalMean2dParameter,
)
from e3sm_diags.plot.cartopy.zonal_mean_2d_plot import plot as plot_func

logger = custom_logger(__name__)

DEFAULT_PLEVS = copy.deepcopy(DEFAULT_PLEVS)


def run_diag(
    parameter: ZonalMean2dParameter, default_plevs=DEFAULT_PLEVS
) -> ZonalMean2dParameter:
    variables = parameter.variables
    seasons = parameter.seasons
    ref_name = getattr(parameter, "ref_name", "")
    regions = parameter.regions

    if not parameter._is_plevs_set():
        parameter.plevs = default_plevs

    test_ds = Dataset(parameter, data_type="test")
    ref_ds = Dataset(parameter, data_type="ref")

    for var_key in variables:
        logger.info("Variable: {}".format(var_key))
        parameter.var_id = var_key

        for season in seasons:
            parameter._set_name_yrs_attrs(test_ds, ref_ds, season)

            ds_test = test_ds.get_climo_dataset(var_key, season)
            ds_ref = ref_ds.get_ref_climo_dataset(var_key, season, ds_test)

            # Store the variable's DataArray objects for reuse.
            dv_test = ds_test[var_key]
            dv_ref = ds_ref[var_key]

            is_vars_3d = has_z_axis(dv_test) and has_z_axis(dv_ref)
            is_dims_diff = has_z_axis(dv_test) != has_z_axis(dv_ref)

            if not is_vars_3d:
                ds_land_sea_mask: xr.Dataset = test_ds._get_land_sea_mask(season)

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
                    season,
                    var_key,
                    ref_name,
                )
            elif is_dims_diff:
                raise RuntimeError(
                    "Dimensions of the two variables are different. Aborting."
                )

    return parameter


def _run_diags_2d(
    parameter: ZonalMean2dParameter,
    ds_test: xr.Dataset,
    ds_ref: xr.Dataset,
    ds_land_sea_mask: xr.Dataset,
    season: str,
    regions: List[str],
    var_key: str,
    ref_name: str,
):
    for region in regions:
        parameter._set_param_output_attrs(var_key, season, region, ref_name, ilev=None)

        (
            ds_test_region,
            ds_ref_region,
            ds_test_region_regrid,
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
    parameter: ZonalMean2dParameter,
    ds_test: xr.Dataset,
    ds_ref: xr.Dataset,
    season: str,
    var_key: str,
    ref_name: str,
):
    plevs = parameter.plevs
    logger.info(f"Selected pressure level: {plevs}")

    ds_t_plevs = regrid_z_axis_to_plevs(ds_test, var_key, plevs)
    ds_r_plevs = regrid_z_axis_to_plevs(ds_ref, var_key, plevs)

    ds_t_plevs = _convert_g_kg_to_ppm_units(parameter, ds_t_plevs, var_key)
    ds_r_plevs = _convert_g_kg_to_ppm_units(parameter, ds_r_plevs, var_key)

    # Calculate the spatial average of the variables on their original pressure
    # level grids. These variables are used for the visualization on original
    # horizonal grids for and for metrics output.
    ds_t_plevs_avg = ds_t_plevs.spatial.average(var_key, axis="X")
    ds_r_plevs_avg = ds_r_plevs.spatial.average(var_key, axis="X")

    # Align the grids for the variables and calculate their spatial averages
    # and the difference between their spatial averages. These variables are
    # used for calculating rmse and correlation.
    (
        ds_t_plevs_rg_avg,
        ds_r_plevs_rg_avg,
        ds_diff_avg,
    ) = _get_avg_for_regridded_datasets(parameter, ds_t_plevs, ds_r_plevs, var_key)

    metrics_dict = _create_metrics_dict(
        var_key,
        ds_t_plevs_avg,
        ds_t_plevs_rg_avg,
        ds_r_plevs_avg,
        ds_r_plevs_rg_avg,
        ds_diff_avg,
    )

    # Set parameter attributes for output files.
    parameter.var_region = "global"
    parameter.output_file = "-".join([ref_name, var_key, season, parameter.regions[0]])
    parameter.main_title = str(" ".join([var_key, season]))

    # NOTE: Taken from all plot function.
    if parameter.diff_type == "relative":
        ds_t_plevs[var_key].units = "%"
        ds_r_plevs[var_key].units = "%"

    _save_data_metrics_and_plots(
        parameter,
        plot_func,
        var_key,
        ds_t_plevs_avg,
        ds_r_plevs_avg,
        ds_diff_avg,
        metrics_dict,
    )


def _convert_g_kg_to_ppm_units(
    parameter: ZonalMean2dParameter, ds: xr.Dataset, var_key: str
) -> xr.Dataset:
    """Adjust the units for a variable to handle special cases.

    This is a special case to handle small values of stratosphere specific
    humidity. The general derived variable process converts specific humidity to
    units [g/kg]. This function converts from "g/kg" to "ppm".

    Parameters
    ----------
    parameter : ZonalMean2dParameter
        The parameter object.
    ds : xr.Dataset
        The dataset.
    var_key : str
        They key of the variable.

    Returns
    -------
    xr.Dataset
        The dataset with units converted.
    """
    ds_new = ds.copy()

    with xr.set_options(keep_attrs=True):
        if (
            parameter.current_set == "zonal_mean_2d_stratosphere"
            and parameter.var_id == "Q"
        ):
            ds_new[var_key] = ds_new[var_key] * 1000.0
            ds_new[var_key].attrs["units"] = "ppm"

    return ds_new


def _get_avg_for_regridded_datasets(
    parameter: ZonalMean2dParameter,
    ds_test_plevs: xr.Dataset,
    ds_ref_plevs: xr.Dataset,
    var_key: str,
) -> Tuple[xr.Dataset, xr.Dataset, xr.Dataset]:
    """Get the average and difference between averages for the plevs datasets.

    Parameters
    ----------
    parameter : ZonalMean2dParameter
        The parameter object.
    ds_test_plevs : xr.Dataset
        The test dataset on pressure level coordinates.
    ds_ref_plevs : xr.Dataset
        The reference dataset on pressure level coordinates.
    var_key : str
        The key of the variable.

    Returns
    -------
    Tuple[xr.Dataset, xr.Dataset, xr.Dataset]
        A tuple consisting of the average of the test dataset, the
        average of the ref dataset, and the difference between averages.
    """
    ds_test_rg, ds_ref_rg = align_grids_to_lower_res(
        ds_test_plevs,
        ds_ref_plevs,
        var_key,
        parameter.regrid_tool,
        parameter.regrid_method,
    )

    # Get the difference between the regridded variables and use it to
    # make sure the regridded variables have the same mask.
    ds_diff_rg = ds_test_rg - ds_ref_rg

    ds_test_rg = ds_ref_rg + ds_diff_rg
    ds_ref_rg = ds_test_rg - ds_diff_rg

    # Calculate the spatial averages for the masked variables.
    ds_test_rg_avg = ds_test_rg.spatial.average(var_key, axis=["X"])
    ds_ref_rg_avg = ds_ref_rg.spatial.average(var_key, axis=["X"])

    # Calculate the spatial average for the differences
    ds_diff_avg = ds_diff_rg.spatial.average(var_key, axis="X")
    if parameter.diff_type == "relative":
        ds_diff_avg = ds_diff_avg / ds_ref_rg_avg * 100.0

    return ds_test_rg_avg, ds_test_rg_avg, ds_diff_avg


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
    deviation, correlation (pearson_r), and RMSE. The default value for
    optional metrics is None.

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
        The difference between ``ds_test_regrid`` and ``ds_ref_regrid`` if both
        exist. This arg will be None if a model only run is performed.

    Returns
    -------
    Metrics
        A dictionary with the key being a string and the value being either
        a sub-dictionary (key is metric and value is float) or a string
        ("unit").
    """
    # TODO: Make sure bounds are set for lev on all datasets.
    metrics_dict: MetricsDict = {}

    metrics_dict["ref"] = {
        "min": ds_ref[var_key].min().item(),
        "max": ds_test[var_key].max().item(),
        # TODO: Axes is "yz", xCDAT spatial average does not support "Z".
        # "mean": mean(ds_ref, axis="yz"),
    }
    metrics_dict["test"] = {
        "min": ds_test[var_key].min().item(),
        "max": ds_test[var_key].max().item(),
        # TODO: Axes is "yz", xCDAT spatial average does not support "Z".
        # "mean": mean(ds_test, axis="yz"),
    }

    metrics_dict["diff"] = {
        "min": ds_diff[var_key].min().item(),
        "max": ds_diff[var_key].max().item(),
        # TODO: Axes is "yz", xCDAT spatial average does not support "Z".
        # "mean": mean(ds_diff, axis="yz"),
    }

    # TODO: Test if these work on "yz" axes.
    metrics_dict["misc"] = {
        "rmse": rmse(ds_test_regrid, ds_ref_regrid, var_key, axis=["Y", "Z"]),
        "corr": correlation(ds_test_regrid, ds_ref_regrid, var_key, axis=["Y", "Z"]),
    }

    return metrics_dict

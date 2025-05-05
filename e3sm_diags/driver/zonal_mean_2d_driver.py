import copy
from typing import Tuple

import xarray as xr
import xcdat as xc  # noqa: F401

from e3sm_diags.driver.utils.dataset_xr import Dataset
from e3sm_diags.driver.utils.io import _save_data_metrics_and_plots
from e3sm_diags.driver.utils.regrid import (
    align_grids_to_lower_res,
    has_z_axis,
    regrid_z_axis_to_plevs,
)
from e3sm_diags.driver.utils.type_annotations import MetricsDict
from e3sm_diags.logger import _setup_child_logger
from e3sm_diags.metrics.metrics import correlation, rmse, spatial_avg
from e3sm_diags.parameter.zonal_mean_2d_parameter import (
    DEFAULT_PLEVS,
    ZonalMean2dParameter,
)
from e3sm_diags.plot.zonal_mean_2d_plot import plot as plot_func

logger = _setup_child_logger(__name__)

DEFAULT_PLEVS = copy.deepcopy(DEFAULT_PLEVS)


def run_diag(
    parameter: ZonalMean2dParameter, default_plevs=DEFAULT_PLEVS
) -> ZonalMean2dParameter:
    variables = parameter.variables
    seasons = parameter.seasons
    ref_name = getattr(parameter, "ref_name", "")

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
            ds_ref = ref_ds.get_climo_dataset(var_key, season)

            # Store the variable's DataArray objects for reuse.
            dv_test = ds_test[var_key]
            dv_ref = ds_ref[var_key]

            is_vars_3d = has_z_axis(dv_test) and has_z_axis(dv_ref)
            is_dims_diff = has_z_axis(dv_test) != has_z_axis(dv_ref)

            if is_dims_diff:
                raise RuntimeError(
                    "The dimensions of the test and reference variables are different, "
                    f"({dv_test.dims} vs. {dv_ref.dims})."
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
            else:
                raise RuntimeError("Dimensions of the two variables are different.")

    return parameter


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
    # horizonal grids for metrics output.
    ds_t_plevs_avg = ds_t_plevs.spatial.average(var_key, axis="X")
    ds_r_plevs_avg = ds_r_plevs.spatial.average(var_key, axis="X")

    # Align the grids for the variables and calculate their spatial averages
    # and the difference between their spatial averages. These variables are
    # used for calculating rmse and correlation.
    (
        ds_t_plevs_rg_avg,
        ds_r_plevs_rg_avg,
        ds_diff_rg_avg,
    ) = _get_avg_for_regridded_datasets(parameter, ds_t_plevs, ds_r_plevs, var_key)

    metrics_dict = _create_metrics_dict(
        var_key,
        ds_t_plevs_avg,
        ds_t_plevs_rg_avg,
        ds_r_plevs_avg,
        ds_r_plevs_rg_avg,
        ds_diff_rg_avg,
    )

    # Set parameter attributes for output files.
    parameter._set_param_output_attrs(
        var_key, season, parameter.regions[0], ref_name, ilev=None
    )
    _save_data_metrics_and_plots(
        parameter,
        plot_func,
        var_key,
        ds_t_plevs_avg,
        ds_r_plevs_avg,
        ds_diff_rg_avg,
        metrics_dict,
        ds_test_regridded=ds_t_plevs_rg_avg,
        ds_ref_regridded=ds_r_plevs_rg_avg,
    )


def _convert_g_kg_to_ppm_units(
    parameter: ZonalMean2dParameter, ds: xr.Dataset, var_key: str
) -> xr.Dataset:
    """Adjust the units for a variable to handle special cases.

    This is a special case to handle small values of stratosphere specific
    humidity. The general derived variable process converts specific humidity to
    units [g/kg]. This function converts from "g/kg" to "ppm" by volume.

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
        if parameter.current_set == "zonal_mean_2d_stratosphere" and (
            parameter.var_id == "Q" or parameter.var_id == "H2OLNZ"
        ):
            ds_new[var_key] = ds_new[var_key] * 28.97 / 18.0 * 1000.0
            ds_new[var_key].attrs["units"] = "ppmv"

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

    with xr.set_options(keep_attrs=True):
        ds_diff_rg = ds_test_rg.copy()
        ds_diff_rg[var_key] = ds_test_rg[var_key] - ds_ref_rg[var_key]

        ds_test_rg[var_key] = ds_ref_rg[var_key] + ds_diff_rg[var_key]
        ds_ref_rg[var_key] = ds_test_rg[var_key] - ds_diff_rg[var_key]

    # Calculate the spatial averages for the masked variables.
    ds_test_rg_avg = ds_test_rg.spatial.average(var_key, axis=["X"])
    ds_ref_rg_avg = ds_ref_rg.spatial.average(var_key, axis=["X"])

    # Calculate the spatial average for the differences
    ds_diff_rg_avg = ds_diff_rg.spatial.average(var_key, axis=["X"])

    if parameter.diff_type == "relative":
        ds_diff_rg_avg[var_key] = (
            ds_diff_rg_avg[var_key] / ds_ref_rg_avg[var_key] * 100.0
        )
        ds_diff_rg_avg[var_key].attrs["units"] = "%"

    return ds_test_rg_avg, ds_ref_rg_avg, ds_diff_rg_avg


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
        "mean": spatial_avg(ds_ref, var_key, axis=["Y", "Z"]),
    }
    metrics_dict["test"] = {
        "min": ds_test[var_key].min().item(),
        "max": ds_test[var_key].max().item(),
        "mean": spatial_avg(ds_test, var_key, axis=["Y", "Z"]),
    }

    metrics_dict["diff"] = {
        "min": ds_diff[var_key].min().item(),
        "max": ds_diff[var_key].max().item(),
        "mean": spatial_avg(ds_diff, var_key, axis=["Y", "Z"]),
    }

    metrics_dict["misc"] = {
        "rmse": rmse(ds_test_regrid, ds_ref_regrid, var_key, axis=["Y", "Z"]),
        "corr": correlation(ds_test_regrid, ds_ref_regrid, var_key, axis=["Y", "Z"]),
    }

    return metrics_dict

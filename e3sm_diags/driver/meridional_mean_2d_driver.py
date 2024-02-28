from __future__ import annotations

from copy import deepcopy
from typing import TYPE_CHECKING

import xarray as xr

from e3sm_diags.driver.utils.dataset_xr import Dataset
from e3sm_diags.driver.utils.io import _save_data_metrics_and_plots
from e3sm_diags.driver.utils.regrid import (
    align_grids_to_lower_res,
    has_z_axis,
    regrid_z_axis_to_plevs,
)
from e3sm_diags.driver.utils.type_annotations import MetricsDict
from e3sm_diags.logger import custom_logger
from e3sm_diags.metrics.metrics import correlation, rmse, spatial_avg
from e3sm_diags.parameter.zonal_mean_2d_parameter import DEFAULT_PLEVS

# TODO: Update this ref after moving the plotter down a dir level.
from e3sm_diags.plot.cartopy.meridional_mean_2d_plot import plot as plot_func

logger = custom_logger(__name__)

if TYPE_CHECKING:
    from e3sm_diags.parameter.meridional_mean_2d_parameter import (
        MeridionalMean2dParameter,
    )

DEFAULT_PLEVS = deepcopy(DEFAULT_PLEVS)


def run_diag(parameter: MeridionalMean2dParameter) -> MeridionalMean2dParameter:
    """Run the meridional_mean_2d diagnostics.

    Parameters
    ----------
    parameter : MeridionalMean2dParameter
        The parameter for the diagnostic.

    Returns
    -------
    MeridionalMean2dParameter
        The parameter for the diagnostic with the result (completed or failed)

    Raises
    ------
    RuntimeError
        If the dimensions of the test and ref variables differ.
    RuntimeError
        If the test or ref variables do are not 3-D (no Z-axis).
    """
    variables = parameter.variables
    seasons = parameter.seasons
    ref_name = getattr(parameter, "ref_name", "")

    test_ds = Dataset(parameter, data_type="test")
    ref_ds = Dataset(parameter, data_type="ref")

    for var_key in variables:
        logger.info("Variable: {}".format(var_key))
        parameter.var_id = var_key

        for season in seasons:
            parameter._set_name_yrs_attrs(test_ds, ref_ds, season)

            ds_test = test_ds.get_climo_dataset(var_key, season)
            ds_ref = ref_ds.get_ref_climo_dataset(var_key, season, ds_test)

            dv_test = ds_test[var_key]
            dv_ref = ds_ref[var_key]

            is_vars_3d = has_z_axis(dv_test) and has_z_axis(dv_ref)
            is_dims_diff = has_z_axis(dv_test) != has_z_axis(dv_ref)

            if is_dims_diff:
                raise RuntimeError(
                    "Dimensions of the test and ref variables are different."
                )
            elif not is_vars_3d:
                raise RuntimeError(
                    "The test and/or ref variables are not 3-D (no Z axis)."
                )
            elif is_vars_3d:
                # Since the default is now stored in MeridionalMean2dParameter,
                # we must get it from there if the plevs param is blank.
                if not parameter._is_plevs_set():
                    parameter.plevs = DEFAULT_PLEVS

                _run_diags_3d(parameter, ds_test, ds_ref, season, var_key, ref_name)

    return parameter


def _run_diags_3d(
    parameter: MeridionalMean2dParameter,
    ds_test: xr.Dataset,
    ds_ref: xr.Dataset,
    season: str,
    var_key: str,
    ref_name: str,
):
    plevs = parameter.plevs

    ds_t_plevs = regrid_z_axis_to_plevs(ds_test, var_key, plevs)
    ds_r_plevs = regrid_z_axis_to_plevs(ds_ref, var_key, plevs)

    ds_t_plevs_avg = ds_t_plevs.spatial.average(var_key, axis=["Y"])
    ds_r_plevs_avg = ds_r_plevs.spatial.average(var_key, axis=["Y"])

    # TODO: Make sure this the logic and output matches the old CDAT version
    # of this code. Can we use xESMF instead here?
    ds_t_plevs_rg_avg, ds_r_plevs_rg_avg = align_grids_to_lower_res(
        ds_t_plevs_avg, ds_r_plevs_avg, var_key, tool="regrid2", method="conservative"
    )

    # Get the difference between final regridded variables.
    with xr.set_options(keep_attrs=True):
        ds_diff_plevs_rg_avg = ds_t_plevs_rg_avg.copy()
        ds_diff_plevs_rg_avg[var_key] = (
            ds_t_plevs_rg_avg[var_key] - ds_r_plevs_rg_avg[var_key]
        )

    metrics_dict = _create_metrics_dict(
        var_key,
        ds_t_plevs_avg,
        ds_t_plevs_rg_avg,
        ds_r_plevs_avg,
        ds_r_plevs_rg_avg,
        ds_diff_plevs_rg_avg,
    )

    # TODO: This section can be turned into a CoreParameter method since it is
    # repeated in various drivers.
    # Set parameter attributes for output files.
    parameter.var_region = "global"
    parameter.output_file = "-".join([ref_name, var_key, season, parameter.regions[0]])
    parameter.main_title = str(" ".join([var_key, season]))

    _save_data_metrics_and_plots(
        parameter,
        plot_func,
        var_key,
        ds_t_plevs_avg,
        ds_r_plevs_avg,
        ds_diff_plevs_rg_avg,
        metrics_dict,
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
        "max": ds_test[var_key].max().item(),
        "mean": spatial_avg(ds_ref, var_key, axis=["X", "Z"]),
    }
    metrics_dict["test"] = {
        "min": ds_test[var_key].min().item(),
        "max": ds_test[var_key].max().item(),
        "mean": spatial_avg(ds_test, var_key, axis=["X", "Z"]),
    }

    metrics_dict["diff"] = {
        "min": ds_diff[var_key].min().item(),
        "max": ds_diff[var_key].max().item(),
        "mean": spatial_avg(ds_diff, var_key, axis=["X", "Z"]),
    }

    metrics_dict["misc"] = {
        "rmse": rmse(ds_test_regrid, ds_ref_regrid, var_key, axis=["X", "Z"]),
        "corr": correlation(ds_test_regrid, ds_ref_regrid, var_key, axis=["X", "Z"]),
    }
    return metrics_dict

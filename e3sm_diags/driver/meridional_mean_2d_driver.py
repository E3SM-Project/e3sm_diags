from __future__ import annotations

from copy import deepcopy
from typing import TYPE_CHECKING

import xarray as xr
import xcdat as xc

from e3sm_diags.driver.utils.dataset_xr import Dataset
from e3sm_diags.driver.utils.io import (
    _get_xarray_datasets,
    _save_data_metrics_and_plots,
)
from e3sm_diags.driver.utils.regrid import (
    align_grids_to_lower_res,
    has_z_axis,
    regrid_z_axis_to_plevs,
)
from e3sm_diags.logger import _setup_child_logger
from e3sm_diags.metrics.metrics import correlation, rmse, spatial_avg
from e3sm_diags.parameter.zonal_mean_2d_parameter import DEFAULT_PLEVS
from e3sm_diags.plot.meridional_mean_2d_plot import plot as plot_func

logger = _setup_child_logger(__name__)

if TYPE_CHECKING:
    from e3sm_diags.driver.utils.type_annotations import MetricsDict, TimeSelection
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
    ref_name = getattr(parameter, "ref_name", "")

    test_ds = Dataset(parameter, data_type="test")
    ref_ds = Dataset(parameter, data_type="ref")

    time_selection_type, time_selections = parameter._get_time_selection_to_use()

    for var_key in variables:
        logger.info("Variable: {}".format(var_key))
        parameter.var_id = var_key

        for time_selection in time_selections:
            ds_test, ds_ref, _ = _get_xarray_datasets(
                test_ds, ref_ds, var_key, time_selection_type, time_selection
            )

            # Set name_yrs after loading data because time sliced datasets
            # have the required attributes only after loading the data.
            parameter._set_name_yrs_attrs(test_ds, ref_ds, time_selection)

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
                if not parameter._is_plevs_set():
                    parameter.plevs = DEFAULT_PLEVS

                _run_diags_3d(
                    parameter, ds_test, ds_ref, time_selection, var_key, ref_name
                )

    return parameter


def _run_diags_3d(
    parameter: MeridionalMean2dParameter,
    ds_test: xr.Dataset,
    ds_ref: xr.Dataset,
    time_selection: TimeSelection,
    var_key: str,
    ref_name: str,
):
    plevs = parameter.plevs

    ds_t_plevs = regrid_z_axis_to_plevs(ds_test, var_key, plevs)
    ds_r_plevs = regrid_z_axis_to_plevs(ds_ref, var_key, plevs)

    ds_t_plevs_avg = ds_t_plevs.spatial.average(var_key, axis=["Y"])
    ds_r_plevs_avg = ds_r_plevs.spatial.average(var_key, axis=["Y"])

    # A placeholder Y axis must be added back to the averaged variables to
    # align grids via horizontal regridding. Afterwards, the Y axis is dropped.
    ds_t_plevs_avg = _add_placeholder_y_axis(ds_test, ds_t_plevs_avg, var_key)
    ds_r_plevs_avg = _add_placeholder_y_axis(ds_ref, ds_r_plevs_avg, var_key)

    ds_t_plevs_rg_avg, ds_r_plevs_rg_avg = align_grids_to_lower_res(
        ds_t_plevs_avg,
        ds_r_plevs_avg,
        var_key,
        tool="xesmf",
        method="conservative_normed",
        axis_to_compare="X",
    )

    # After regridding, squeeze the placeholder Y axis.
    ds_t_plevs_avg = ds_t_plevs_avg.squeeze()
    ds_r_plevs_avg = ds_r_plevs_avg.squeeze()
    ds_t_plevs_rg_avg = ds_t_plevs_rg_avg.squeeze()
    ds_r_plevs_rg_avg = ds_r_plevs_rg_avg.squeeze()

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

    parameter._set_param_output_attrs(
        var_key, time_selection, parameter.regions[0], ref_name, ilev=None
    )
    _save_data_metrics_and_plots(
        parameter,
        plot_func,
        var_key,
        ds_t_plevs_avg,
        ds_r_plevs_avg,
        ds_diff_plevs_rg_avg,
        metrics_dict,
        ds_test_regridded=ds_t_plevs_rg_avg,
        ds_ref_regridded=ds_r_plevs_rg_avg,
    )


def _add_placeholder_y_axis(ds_original: xr.Dataset, ds_avg: xr.Dataset, var_key: str):
    lat_original = xc.get_dim_coords(ds_original, axis="Y")

    # Make sure to drop the old Y axis before adding the new Y axis. Otherwise,
    # the new Y axis can't be added.
    lat_key = lat_original.name
    ds_avg_new = ds_avg.drop_dims(lat_key)

    ds_avg_new[var_key] = ds_avg_new[var_key].expand_dims({lat_key: [0]})
    ds_avg_new[lat_key].attrs = lat_original.attrs.copy()

    lat_bnds_key = lat_original.attrs["bounds"]
    lat_bnds_original = ds_original[lat_bnds_key].copy()
    ds_avg_new[lat_bnds_key] = xr.DataArray(
        name=lat_bnds_key,
        dims=lat_bnds_original.dims,
        data=[[-1, 1]],
        attrs=lat_bnds_original.attrs.copy(),
    )

    return ds_avg_new


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

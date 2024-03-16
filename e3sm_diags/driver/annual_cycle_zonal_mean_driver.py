from __future__ import annotations

from typing import TYPE_CHECKING, List, Tuple

import xarray as xr
import xcdat as xc
from scipy import interpolate

from e3sm_diags.driver.utils.dataset_xr import Dataset
from e3sm_diags.driver.utils.io import _save_data_metrics_and_plots
from e3sm_diags.driver.utils.regrid import (
    get_z_axis,
    has_z_axis,
    regrid_z_axis_to_plevs,
)
from e3sm_diags.logger import custom_logger
from e3sm_diags.metrics.metrics import spatial_avg
from e3sm_diags.plot.zonal_mean_xy_plot import plot as plot_func

logger = custom_logger(__name__)

if TYPE_CHECKING:
    from e3sm_diags.parameter.core_parameter import CoreParameter


def run_diag(parameter: CoreParameter) -> CoreParameter:
    """Get annual cycle zonal mean results for the annual_cycle_zonal_mean diagnostic set.

    This function loops over each 2D variable

    NOTE: zonal_mean_xy set only supports "global" region.

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
    ref_name = getattr(parameter, "ref_name", "")
    regions = parameter.regions

    for region in regions:
        if region != "global":
            raise RuntimeError(
                f"Region ({region}) is not supported. Only global region is currently "
                "supported for the annual_cycle_zonal set."
            )

    # Variables storing xarray `Dataset` objects start with `ds_` and
    # variables storing e3sm_diags `Dataset` objects end with `_ds`. This
    # is to help distinguish both objects from each other.
    test_ds = Dataset(parameter, data_type="test")
    ref_ds = Dataset(parameter, data_type="ref")

    for var_key in variables:
        logger.info("Variable: {}".format(var_key))
        parameter.var_id = var_key

        parameter._set_name_yrs_attrs(test_ds, ref_ds, "01")

        ds_test = test_ds.get_climo_dataset(var_key, "01")
        # TODO consider to refactor the behavior of get_ref_climo_dataset
        ds_ref = ref_ds.get_ref_climo_dataset(var_key, "01", ds_test)

        # Store the variable's DataArray objects for reuse.
        dv_test = ds_test[var_key]
        dv_ref = ds_ref[var_key]

        is_vars_3d = has_z_axis(dv_test) and has_z_axis(dv_ref)
        is_dims_diff = has_z_axis(dv_test) != has_z_axis(dv_ref)

        if not is_vars_3d:
            _run_diags_2d(
                parameter,
                ds_test,
                ds_ref,
                "AC",
                regions,
                var_key,
                ref_name,
            )
        elif is_vars_3d:
            raise RuntimeError(
                "3 D variables are not supported in annual cycle zonal mean set. Aborting."
            )
        elif is_dims_diff:
            raise RuntimeError(
                "Dimensions of the two variables are different. Aborting."
            )

    return parameter


def _run_diags_2d(
    parameter: CoreParameter,
    ds_test: xr.Dataset,
    ds_ref: xr.Dataset,
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
        The dataset containing the ref variable. If this is a model-only run
        then it will be the same dataset as ``ds_test``.
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
        logger.info(f"Selected region: {region}")

        parameter._set_param_output_attrs(var_key, season, region, ref_name, ilev=None)
        # Calculate annual cycle
        # Regridding
        # Calculate zonal mean
        da_test_1d, da_ref_1d = _calc_zonal_mean(ds_test, ds_ref, var_key)
        da_diff_1d = _get_diff_of_zonal_means(da_test_1d, da_ref_1d)

        _save_data_metrics_and_plots(
            parameter,
            plot_func,
            var_key,
            da_test_1d.to_dataset(),
            da_ref_1d.to_dataset(),
            da_diff_1d.to_dataset(),
            metrics_dict=None,
        )


def _get_annual_cycle(
    ds_test: xr.Dataset,
    var_key: str,
) -> xr.DataArray:
    """Get annual cycle.

    # TODO: Write unit tests for this function.

    Parameters
    ----------
    ds_test : xr.Dataset
        The dataset containing the test variable.
    var_key : str
        The key of the variable.

    Returns
    -------
    xr.DataArray
        xarray DatAarray
    """
    months = range(1, 13)
    month_list = [f"{x:02}" for x in list(months)]

    for index, month in enumerate(month_list):
        ds_test_mon = ds_test.get_climo_dataset(var_key, month)
        print(ds_test_mon)


def _calc_zonal_mean(
    ds_test: xr.Dataset,
    ds_ref: xr.Dataset,
    var_key: str,
) -> Tuple[xr.DataArray, xr.DataArray]:
    """Calculate zonal mean metrics.

    # TODO: Write unit tests for this function.

    Parameters
    ----------
    ds_test : xr.Dataset
        The dataset containing the test variable.
    ds_ref : xr.Dataset
        The dataset containing the ref variable. If this is a model-only run
        then it will be the same dataset as ``ds_test``.
    var_key : str
        The key of the variable.

    Returns
    -------
    Tuple[xr.DataArray, xr.DataArray]
        A Tuple containing the zonal mean for the test variable and the ref
        variable.
    """
    da_test = spatial_avg(ds_test, var_key, axis=["X"], as_list=False)
    da_ref = spatial_avg(ds_ref, var_key, axis=["X"], as_list=False)

    return da_test, da_ref  # type: ignore

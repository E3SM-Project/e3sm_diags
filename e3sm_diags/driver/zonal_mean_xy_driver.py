from typing import TYPE_CHECKING

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
from e3sm_diags.logger import _setup_child_logger
from e3sm_diags.metrics.metrics import spatial_avg
from e3sm_diags.plot.zonal_mean_xy_plot import plot as plot_func

logger = _setup_child_logger(__name__)

if TYPE_CHECKING:
    from e3sm_diags.parameter.core_parameter import CoreParameter


def run_diag(parameter: CoreParameter) -> CoreParameter:
    """Get zonal mean results for the zonal_mean_xy diagnostic set.

    This function loops over each variable, season, pressure level (if 3-D).

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
    seasons = parameter.seasons
    ref_name = getattr(parameter, "ref_name", "")
    regions = parameter.regions

    for region in regions:
        if region != "global":
            raise RuntimeError(
                f"Region ({region}) is not supported. Only global region is currently "
                "supported for the zonal_mean_xy set."
            )

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
            ds_ref = ref_ds.get_climo_dataset(var_key, season)

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
                    regions,
                    var_key,
                    ref_name,
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
        logger.info(f"Selected region: {region}")

        da_test_1d, da_ref_1d = _calc_zonal_mean(ds_test, ds_ref, var_key)
        da_diff_1d = _get_diff_of_zonal_means(da_test_1d, da_ref_1d)

        parameter._set_param_output_attrs(var_key, season, region, ref_name, ilev=None)
        _save_data_metrics_and_plots(
            parameter,
            plot_func,
            var_key,
            da_test_1d.to_dataset(),
            da_ref_1d.to_dataset(),
            da_diff_1d.to_dataset(),
            metrics_dict=None,
        )


def _run_diags_3d(
    parameter: CoreParameter,
    ds_test: xr.Dataset,
    ds_ref: xr.Dataset,
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
            logger.info(f"Selected region: {region}")

            da_test_1d, da_ref_1d = _calc_zonal_mean(ds_test_ilev, ds_ref_ilev, var_key)
            da_diff_1d = _get_diff_of_zonal_means(da_test_1d, da_ref_1d)

            parameter._set_param_output_attrs(var_key, season, region, ref_name, ilev)
            _save_data_metrics_and_plots(
                parameter,
                plot_func,
                var_key,
                da_test_1d.to_dataset(),
                da_ref_1d.to_dataset(),
                da_diff_1d.to_dataset(),
                metrics_dict=None,
            )


def _calc_zonal_mean(
    ds_test: xr.Dataset,
    ds_ref: xr.Dataset,
    var_key: str,
) -> tuple[xr.DataArray, xr.DataArray]:
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
    tuple[xr.DataArray, xr.DataArray]
        A Tuple containing the zonal mean for the test variable and the ref
        variable.
    """
    da_test_1d = spatial_avg(ds_test, var_key, axis=["X"], as_list=False)
    da_ref_1d = spatial_avg(ds_ref, var_key, axis=["X"], as_list=False)

    return da_test_1d, da_ref_1d  # type: ignore


def _get_diff_of_zonal_means(da_a: xr.DataArray, da_b: xr.DataArray) -> xr.DataArray:
    """Get the difference between the zonal means of two variables.

    Both variables are aligned to the same grid (lower resolution of the two)
    and the difference is calculated.

    # TODO: Write unit tests for this function

    Parameters
    ----------
    da_a : xr.DataArray
        The first variable.
    da_b : xr.DataArray
        The second variable.

    Returns
    -------
    xr.DataArray
        The difference between the zonal means of two variables.
    """

    with xr.set_options(keep_attrs=True):
        lat_a = xc.get_dim_coords(da_a, axis="Y")
        lat_b = xc.get_dim_coords(da_b, axis="Y")
        if len(lat_a) > len(lat_b):
            interpf = interpolate.interp1d(lat_a, da_a.values, bounds_error=False)
            da_a_new = da_b.copy(data=interpf(lat_b))

            return da_a_new - da_b.copy()
        else:
            interpf = interpolate.interp1d(lat_b, da_b.values, bounds_error=False)
            da_b_new = da_a.copy(data=interpf(lat_a))

            return da_a.copy() - da_b_new

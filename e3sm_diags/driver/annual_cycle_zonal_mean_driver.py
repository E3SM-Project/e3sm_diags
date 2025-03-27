from __future__ import annotations

from typing import TYPE_CHECKING, List

import xarray as xr
import xcdat as xc

from e3sm_diags.driver.utils.dataset_xr import Dataset
from e3sm_diags.driver.utils.io import _save_data_metrics_and_plots
from e3sm_diags.driver.utils.regrid import align_grids_to_lower_res, has_z_axis
from e3sm_diags.logger import _setup_child_logger
from e3sm_diags.metrics.metrics import spatial_avg
from e3sm_diags.plot.annual_cycle_zonal_mean_plot import plot as plot_func

logger = _setup_child_logger(__name__)

if TYPE_CHECKING:
    from e3sm_diags.parameter.core_parameter import CoreParameter


def run_diag(parameter: CoreParameter) -> CoreParameter:
    """Get annual cycle zonal mean results for the annual_cycle_zonal_mean diagnostic set.

    This function loops over each 2D variable

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

    test_ds = Dataset(parameter, data_type="test")
    ref_ds = Dataset(parameter, data_type="ref")

    for var_key in variables:
        logger.info("Variable: {}".format(var_key))
        parameter.var_id = var_key

        parameter._set_name_yrs_attrs(test_ds, ref_ds, "01")

        ds_test = test_ds.get_climo_dataset(var_key, "ANNUALCYCLE")
        ds_ref = ref_ds.get_climo_dataset(var_key, "ANNUALCYCLE")

        # Encode the time coordinate to month integer (1 to 12). This is
        # necessary to avoid cases where "months since ..." units are used
        # and the calendar attribute is not set to "360_day". It also
        # replicates the CDAT code behavior.
        ds_test = _encode_time_coords(ds_test)
        ds_ref = _encode_time_coords(ds_ref)

        dv_test = ds_test[var_key]
        dv_ref = ds_ref[var_key]

        is_vars_3d = has_z_axis(dv_test) and has_z_axis(dv_ref)
        is_dims_diff = has_z_axis(dv_test) != has_z_axis(dv_ref)

        if is_dims_diff:
            raise RuntimeError(
                "The dimensions of the test and reference variables are different, "
                f"({dv_test.dims} vs. {dv_ref.dims})."
            )
        elif not is_vars_3d:
            _run_diags_annual_cycle(
                parameter,
                ds_test,
                ds_ref,
                regions,
                var_key,
                ref_name,
            )
        elif is_vars_3d:
            raise RuntimeError(
                "3-D variables are not supported in annual cycle zonal mean set."
            )
        else:
            raise RuntimeError("Dimensions of the two variables are different.")

    return parameter


def _run_diags_annual_cycle(
    parameter: CoreParameter,
    ds_test: xr.Dataset,
    ds_ref: xr.Dataset,
    regions: List[str],
    var_key: str,
    ref_name: str,
):
    """Run annual cycle zonal run diagnostics.

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
    regions : List[str]
        The list of regions.
    var_key : str
        The key of the variable.
    ref_name : str
        The reference name.
    """
    for region in regions:
        logger.info(f"Selected region: {region}")

        ds_test_reg, ds_ref_reg = align_grids_to_lower_res(
            ds_test,
            ds_ref,
            var_key,
            parameter.regrid_tool,
            parameter.regrid_method,
        )

        test_zonal_mean: xr.DataArray = spatial_avg(
            ds_test, var_key, axis=["X"], as_list=False
        )  # type: ignore
        test_reg_zonal_mean: xr.DataArray = spatial_avg(
            ds_test_reg,
            var_key,
            axis=["X"],
            as_list=False,  # type: ignore
        )

        if (
            parameter.ref_name == "OMI-MLS"
        ):  # SCO from OMI-MLS only available as (time, lat)
            test_zonal_mean = test_zonal_mean.sel(lat=slice(-60, 60))
            test_reg_zonal_mean = test_reg_zonal_mean.sel(lat=slice(-60, 60))

            if var_key == "SCO":
                ref_zonal_mean = ds_ref[var_key]
                ref_reg_zonal_mean = ds_ref[var_key]
            else:
                ref_zonal_mean = spatial_avg(ds_ref, var_key, axis=["X"], as_list=False)  # type: ignore
                ref_reg_zonal_mean = spatial_avg(
                    ds_ref_reg,
                    var_key,
                    axis=["X"],
                    as_list=False,  # type: ignore
                )

        else:
            ref_zonal_mean = spatial_avg(ds_ref, var_key, axis=["X"], as_list=False)  # type: ignore
            ref_reg_zonal_mean = spatial_avg(
                ds_ref_reg,
                var_key,
                axis=["X"],
                as_list=False,  # type: ignore
            )

        # Make a copy of dataset to preserve time dimension
        with xr.set_options(keep_attrs=True):
            diff = test_reg_zonal_mean - ref_reg_zonal_mean

        parameter._set_param_output_attrs(
            var_key, "ANNUALCYCLE", region, ref_name, ilev=None
        )
        _save_data_metrics_and_plots(
            parameter,
            plot_func,
            var_key,
            test_zonal_mean.to_dataset(),
            ref_zonal_mean.to_dataset(),
            diff.to_dataset(),
            metrics_dict=None,
            ds_test_regridded=test_reg_zonal_mean.to_dataset(),
            ds_ref_regridded=ref_reg_zonal_mean.to_dataset(),
        )


def _encode_time_coords(ds: xr.Dataset) -> xr.Dataset:
    """Encode the time coordinates to month integers (1 to 12).

    Parameters
    ----------
    ds : xr.Dataset
        The dataset with decoded time coordinates in `cftime`.

    Returns
    -------
    xr.Dataset
        The dataset with decoded time coordinates as month integers.
    """
    ds_new = ds.copy()

    time_coords = xc.get_dim_coords(ds_new, axis="T")
    dim_key = time_coords.name

    encoded_time, _ = xc.create_axis(dim_key, list(range(1, 13)))
    encoded_time.attrs = time_coords.attrs

    ds_new = ds_new.assign_coords({dim_key: encoded_time})

    return ds_new

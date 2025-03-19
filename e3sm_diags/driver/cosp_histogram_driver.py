from __future__ import annotations

from typing import TYPE_CHECKING

import xarray as xr

from e3sm_diags.driver.utils.dataset_xr import Dataset
from e3sm_diags.driver.utils.io import _save_data_metrics_and_plots
from e3sm_diags.driver.utils.regrid import _subset_on_region
from e3sm_diags.logger import _setup_child_logger
from e3sm_diags.metrics.metrics import spatial_avg
from e3sm_diags.plot.cosp_histogram_plot import plot as plot_func

if TYPE_CHECKING:
    from e3sm_diags.parameter.core_parameter import CoreParameter

logger = _setup_child_logger(__name__)


def run_diag(parameter: CoreParameter) -> CoreParameter:
    """Get metrics for the cosp_histogram diagnostic set.

    This function loops over each variable, season, pressure level, and region.

    It subsets the test and reference variables on the selected region, then
    calculates the spatial average for both variables. The difference between
    the test and reference spatial averages is calculated. Afterwards, the
    spatial averages for the test, ref, and differences are plotted.

    Parameters
    ----------
    parameter : CoreParameter
        The parameter for the diagnostic.

    Returns
    -------
    CoreParameter
        The parameter for the diagnostic with the result (completed or failed).
    """
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

            ds_test = test_ds.get_climo_dataset(var_key, season)
            ds_ref = ref_ds.get_climo_dataset(var_key, season)

            for region in regions:
                logger.info("Selected region: {}".format(region))

                ds_test_region = _subset_on_region(ds_test, var_key, region)
                ds_ref_region = _subset_on_region(ds_ref, var_key, region)

                # Make a copy of the regional datasets to overwrite the existing
                # variable with its spatial average.
                ds_test_avg = ds_test.copy()
                ds_ref_avg = ds_test.copy()
                ds_test_avg[var_key] = spatial_avg(
                    ds_test_region, var_key, as_list=False
                )
                ds_ref_avg[var_key] = spatial_avg(ds_ref_region, var_key, as_list=False)

                # The dimension names of both variables must be aligned to
                # perform arithmetic with Xarray. Sometimes the dimension names
                # might differ based on the derived variable (e.g.,
                # "cosp_htmisr" vs. "misr_cth").
                ds_test_avg[var_key] = _align_test_to_ref_dims(
                    ds_test_avg[var_key], ds_ref_avg[var_key]
                )
                ds_diff_avg = _get_diff_of_avg(var_key, ds_test_avg, ds_ref_avg)

                parameter._set_param_output_attrs(
                    var_key, season, region, ref_name, ilev=None
                )
                _save_data_metrics_and_plots(
                    parameter,
                    plot_func,
                    var_key,
                    ds_test_avg,
                    ds_ref_avg,
                    ds_diff_avg,
                    metrics_dict=None,
                )

    return parameter


def _get_diff_of_avg(
    var_key: str, ds_test_avg: xr.Dataset, ds_ref_avg: xr.Dataset
) -> xr.Dataset:
    # Use the test dataset as the base dataset to subtract the reference dataset
    # from.
    ds_diff_avg = ds_test_avg.copy()

    # There are case where the axes of the test and ref datasets aren't in the
    # same units. We avoid label-based Xarray arithmetic which expect coordinates
    # with the same units and will produce np.nan results by subtracting.
    # Instead, we subtract using the reference xr.DataArray's `np.array`
    # (`.values`).
    ds_diff_avg[var_key] = ds_diff_avg[var_key] - ds_ref_avg[var_key].values

    return ds_diff_avg


def _align_test_to_ref_dims(
    da_test: xr.DataArray, da_ref: xr.DataArray
) -> xr.DataArray:
    """Align the dimensions of the test data to the ref data.

    This is useful for situations where label-based arithmetic needs to be
    performed using Xarray.

    Parameters
    ----------
    da_test : xr.DataArray
        The test dataarray.
    da_ref : xr.DataArray
        The ref dataarray.

    Returns
    -------
    xr.DataArray
        The test dataarray with dimensions aligned to the ref dataarray.
    """
    da_test_new = da_test.copy()

    # NOTE: This logic assumes that prs and tau are in the same order for
    # the test and ref variables. If they are not, then this will break or
    # perform incorrect arithmetic.
    # FIXME: B905: zip() without an explicit strict= parameter
    da_test_new = da_test_new.rename(
        {dim1: dim2 for dim1, dim2 in zip(da_test.dims, da_ref.dims)}
    )

    return da_test_new

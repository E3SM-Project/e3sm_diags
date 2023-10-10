from __future__ import annotations

from typing import TYPE_CHECKING

from e3sm_diags.driver.utils.dataset_xr import Dataset
from e3sm_diags.driver.utils.io import _save_data_metrics_and_plots
from e3sm_diags.driver.utils.regrid import _subset_on_region
from e3sm_diags.logger import custom_logger
from e3sm_diags.metrics.metrics import spatial_avg
from e3sm_diags.plot.cosp_histogram_plot import plot as plot_func

if TYPE_CHECKING:
    from e3sm_diags.parameter.core_parameter import CoreParameter

logger = custom_logger(__name__)


def run_diag(parameter: CoreParameter) -> CoreParameter:
    """Get metrics for the cosp_histogram diagnostic set.

    This funciton loops over each variable, season, pressure level, and region.

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
            ds_ref = ref_ds.get_ref_climo_dataset(var_key, season, ds_test)

            for region in regions:
                logger.info("Selected region: {}".format(region))

                ds_test_region = _subset_on_region(ds_test, var_key, region)
                ds_ref_region = _subset_on_region(ds_ref, var_key, region)

                parameter._set_param_output_attrs(
                    var_key, season, region, ref_name, ilev=None
                )

                # Make a copy of the regional datasets to overwrite the existing
                # variable with its spatial average.
                ds_test_region_avg = ds_test_region.copy()
                ds_ref_region_avg = ds_ref_region.copy()
                ds_test_region_avg[var_key] = spatial_avg(
                    ds_test_region, var_key, as_list=False
                )
                ds_ref_region_avg[var_key] = spatial_avg(
                    ds_ref_region, var_key, as_list=False
                )
                ds_diff_region_avg = ds_test_region_avg - ds_ref_region_avg

                # TODO: Need to update this function to use cosp_histogram_plot.py
                _save_data_metrics_and_plots(
                    parameter,
                    plot_func,
                    var_key,
                    ds_test_region_avg,
                    ds_ref_region_avg,
                    ds_diff_region_avg,
                    metrics_dict=None,
                )

    return parameter

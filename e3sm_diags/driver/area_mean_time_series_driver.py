from __future__ import annotations

import collections
import json
import os
from typing import TYPE_CHECKING

import cdutil
import xarray as xr
import xcdat as xc

from e3sm_diags.driver import LAND_OCEAN_MASK_PATH, utils
from e3sm_diags.driver.utils.dataset_xr import Dataset, squeeze_time_dim
from e3sm_diags.driver.utils.regrid import _apply_land_sea_mask
from e3sm_diags.logger import custom_logger
from e3sm_diags.metrics.metrics import spatial_avg
from e3sm_diags.plot.cartopy import area_mean_time_series_plot

if TYPE_CHECKING:
    from e3sm_diags.parameter.area_mean_time_series_parameter import (
        AreaMeanTimeSeriesParameter,
    )


logger = custom_logger(__name__)

RefsTestMetrics = collections.namedtuple("RefsTestMetrics", ["refs", "test", "metrics"])


def run_diag(parameter: AreaMeanTimeSeriesParameter) -> AreaMeanTimeSeriesParameter:
    """Run the diagnostics for area_mean_time_series.

    Parameters
    ----------
    parameter : AreaMeanTimeSeriesParameter
        The parameter for area_mean_time_series.

    Returns
    -------
    AreaMeanTimeSeriesParameter
        The parameter for area_mean_time_series with the results of the
        diagnostic run.
    """
    variables = parameter.variables
    regions = parameter.regions
    ref_names = parameter.ref_names

    for var in variables:
        logger.info("Variable: {}".format(var))

        metrics_dict = {}
        save_data = {}

        test_ds = Dataset(parameter, data_type="test")
        ds_mask = _get_default_land_sea_mask()

        for region in regions:
            logger.info("Selected region: {}".format(region))

            # NOTE:  Test dataset portion
            # ------------------------------------------------------------------
            ds_test = test_ds.get_time_series_dataset(var)
            time_coords = xc.get_dim_coords(ds_test, axis="T").values

            logger.info(
                "Start and end time for selected time slices for test data: "
                f"{time_coords[0]} "
                f"{time_coords[1]}",
            )

            parameter.viewer_descr[var] = getattr(ds_test, "long_name", var)
            parameter.test_name_yrs = test_ds.get_name_yrs_attr()
            ds_test_region = _apply_land_sea_mask(
                ds_test,
                ds_mask,
                var,
                region,  # type: ignore
                parameter.regrid_tool,
                parameter.regrid_method,
            )

            # Average over selected region, and average over months to get the
            # yearly mean.
            ds_test_region_avg: xr.DataArray = spatial_avg(  # type: ignore
                ds_test_region, var, axis=["X", "Y"]
            )
            ds_test_region_avg = ds_test_region_avg.bounds.add_time_bounds(
                "freq", freq="month"
            )

            ds_test_region_avg_yr = cdutil.YEAR(ds_test_region)

            # TODO: I don't think setting attributes again is needed
            ds_test_region_avg_yr.long_name = ds_test.long_name
            ds_test_region_avg_yr.units = ds_test.units
            save_data[parameter.test_name_yrs] = ds_test_region_avg_yr.asma().tolist()

            # NOTE:  Reference dataset portion
            # ------------------------------------------------------------------
            refs = []
            for ref_name in ref_names:
                parameter.ref_name = ref_name
                ref_ds = Dataset(parameter, data_type="ref")

                parameter.ref_name_yrs = ref_ds.get_name_yrs_attr()
                try:
                    ds_ref = ref_ds.get_time_series_dataset(var)

                    ref_domain = _apply_land_sea_mask(
                        ds_ref,
                        ds_mask,
                        var,
                        region,  # type:ignore
                        parameter.regrid_tool,
                        parameter.regrid_method,
                    )

                    ref_domain = cdutil.averager(ref_domain, axis="xy")
                    cdutil.setTimeBoundsMonthly(ref_domain)
                    logger.info(
                        (
                            "Start and end time for selected time slices for ref data: "
                            f"{ref_domain.getTime().asComponentTime()[0]} "
                            f"{ref_domain.getTime().asComponentTime()[-1]}"
                        )
                    )

                    ref_domain_year = cdutil.YEAR(ref_domain)
                    ref_domain_year.ref_name = ref_name
                    save_data[ref_name] = ref_domain_year.asma().tolist()

                    refs.append(ref_domain_year)
                except Exception:
                    logger.exception(
                        "No valid value for reference datasets available for the specified time range"
                    )

            # NOTE: I/O and plotting portion
            # ------------------------------------------------------------------
            parameter.output_file = "-".join([var, region])
            fnm = os.path.join(
                utils.general.get_output_dir(parameter.current_set, parameter),
                parameter.output_file + ".json",
            )

            with open(fnm, "w") as outfile:
                json.dump(save_data, outfile)

            metrics_dict[region] = RefsTestMetrics(
                test=ds_test_region_avg_yr, refs=refs, metrics=[]
            )

        area_mean_time_series_plot.plot(var, metrics_dict, parameter)

    return parameter


def _get_default_land_sea_mask() -> xr.Dataset:
    """Get the e3sm_diags default land sea mask.

    Returns
    -------
    xr.Dataset
        The land sea mask dataset object.
    """
    ds_mask = xr.open_dataset(LAND_OCEAN_MASK_PATH)
    ds_mask = squeeze_time_dim(ds_mask)

    return ds_mask

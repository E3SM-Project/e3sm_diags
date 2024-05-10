from __future__ import annotations

import collections
import json
import os
from typing import TYPE_CHECKING

import xarray as xr
import xcdat as xc

from e3sm_diags.driver import LAND_OCEAN_MASK_PATH, utils
from e3sm_diags.driver.utils.dataset_xr import Dataset, squeeze_time_dim
from e3sm_diags.driver.utils.regrid import _apply_land_sea_mask, _subset_on_region
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

    test_ds = Dataset(parameter, data_type="test")
    parameter.test_name_yrs = test_ds.get_name_yrs_attr()

    ds_mask = _get_default_land_sea_mask()

    for var in variables:
        logger.info("Variable: {}".format(var))

        metrics_dict = {}
        save_data = {}

        for region in regions:
            logger.info("Selected region: {}".format(region))

            ds_test = test_ds.get_time_series_dataset(var)
            _log_time(ds_test, "test")

            ds_test_domain_avg = _get_test_region_yearly_avg(
                parameter, ds_test, ds_mask, var, region
            )
            save_data[parameter.test_name_yrs] = ds_test_domain_avg[var].values.tolist()

            # NOTE:  Reference dataset portion
            # ------------------------------------------------------------------
            refs = []

            for ref_name in ref_names:
                parameter.ref_name = ref_name

                ref_ds = Dataset(parameter, data_type="ref")
                parameter.ref_name_yrs = ref_ds.get_name_yrs_attr()

                ds_ref_domain_avg_yr = _get_ref_region_yearly_avg(
                    parameter, ref_ds, ds_mask, var, region
                )

                if ds_ref_domain_avg_yr is not None:
                    save_data[ref_name] = ds_ref_domain_avg_yr.asma().tolist()
                    refs.append(ds_ref_domain_avg_yr)

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
                test=ds_test_domain_avg, refs=refs, metrics=[]
            )

        parameter.viewer_descr[var] = getattr(ds_test, "long_name", var)
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


def _get_test_region_yearly_avg(
    parameter: AreaMeanTimeSeriesParameter,
    ds_test: xr.Dataset,
    ds_mask: xr.Dataset,
    var: str,
    region: str,
) -> xr.Dataset:
    ds_test_new = ds_test.copy()

    if "land" in region or "ocean" in region:
        ds_test_new = _apply_land_sea_mask(
            ds_test_new,
            ds_mask,
            var,
            region,  # type: ignore
            parameter.regrid_tool,
            parameter.regrid_method,
        )

    if "global" not in region:
        ds_test_new = _subset_on_region(ds_test_new, var, region)

        # Average over selected region, and average over months to get the
        # yearly mean.
    da_test_region_avg: xr.DataArray = spatial_avg(  # type: ignore
        ds_test_new, var, axis=["X", "Y"], as_list=False
    )
    ds_test_region_avg = da_test_region_avg.to_dataset()
    ds_test_region_avg = ds_test_region_avg.bounds.add_time_bounds("freq", freq="month")
    ds_test_region_avg_yr = ds_test_region_avg.temporal.average(var)

    return ds_test_region_avg_yr


def _get_ref_region_yearly_avg(
    parameter: AreaMeanTimeSeriesParameter,
    ref_ds: Dataset,
    ds_mask: xr.Dataset,
    var: str,
    region: str,
) -> xr.Dataset | None:
    ds_ref_domain_avg = None

    try:
        ds_ref = ref_ds.get_time_series_dataset(var)

        if "land" in region or "ocean" in region:
            ds_ref_domain = _apply_land_sea_mask(
                ds_ref,
                ds_mask,
                var,
                region,  # type:ignore
                parameter.regrid_tool,
                parameter.regrid_method,
            )

            # Get the spatial average of the X and Y axes.
        ds_ref_domain_avg = spatial_avg(  # type: ignore
            ds_ref_domain, var, axis=["X", "Y"], as_list=False
        ).to_dataset()
        _log_time(ds_ref_domain_avg, "ref")

        # Add time bounds and calculate the annual average.
        ds_ref_domain_avg = ds_ref_domain_avg.bounds.add_time_bounds(
            "freq", freq="month"
        )
        ds_ref_domain_avg_yr = ds_ref_domain_avg.temporal.average(var)
        # FIXME: Should this be ValueError>?
    except Exception:
        logger.exception(
            "No valid value for reference datasets available for the specified time range"
        )

    return ds_ref_domain_avg_yr


def _log_time(ds: xr.Dataset, data_type: str):
    time_coords = xc.get_dim_coords(ds, axis="T").values

    logger.info(
        f"Start and end time for selected time slices for {data_type} data: "
        f"{time_coords[0]} "
        f"{time_coords[1]}",
    )

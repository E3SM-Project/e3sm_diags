from __future__ import annotations

import json
import os
from typing import TYPE_CHECKING, Literal, NamedTuple

import xarray as xr
import xcdat as xc

from e3sm_diags.driver import LAND_OCEAN_MASK_PATH
from e3sm_diags.driver.utils.dataset_xr import Dataset, squeeze_time_dim
from e3sm_diags.driver.utils.io import _get_output_dir, _write_to_netcdf
from e3sm_diags.driver.utils.regrid import _apply_land_sea_mask, _subset_on_region
from e3sm_diags.logger import _setup_child_logger
from e3sm_diags.metrics.metrics import spatial_avg
from e3sm_diags.plot import area_mean_time_series_plot

if TYPE_CHECKING:
    from e3sm_diags.parameter.area_mean_time_series_parameter import (
        AreaMeanTimeSeriesParameter,
    )


logger = _setup_child_logger(__name__)


class RefsTestMetrics(NamedTuple):
    test: xr.Dataset
    refs: list
    metrics: list


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

            # Test variable metrics.
            # ------------------------------------------------------------------
            ds_test = test_ds.get_time_series_dataset(var)
            _log_time(ds_test, "test")

            ds_test_domain_avg = _get_region_yearly_avg(
                parameter, ds_test, ds_mask, var, region, data_type="test"
            )
            save_data[parameter.test_name_yrs] = ds_test_domain_avg[var].values.tolist()

            # Calculate reference metrics (optional, if valid time range).
            # ------------------------------------------------------------------
            refs = []

            for ref_name in ref_names:
                parameter.ref_name = ref_name

                ref_ds = Dataset(parameter, data_type="ref")
                parameter.ref_name_yrs = ref_ds.get_name_yrs_attr()

                ds_ref_domain_avg_yr = None
                try:
                    ds_ref = ref_ds.get_time_series_dataset(var)
                except (IOError, ValueError):
                    logger.exception(
                        "No valid value for reference datasets available for the specified time range"
                    )
                else:
                    ds_ref_domain_avg_yr = _get_region_yearly_avg(
                        parameter, ds_ref, ds_mask, var, region, data_type="ref"
                    )

                if ds_ref_domain_avg_yr is not None:
                    save_data[ref_name] = ds_ref_domain_avg_yr.asma().tolist()

                    # Set the ref name attribute on the ref variable for the
                    # plot label.
                    ds_ref_domain_avg_yr[var].attrs["ref_name"] = ref_name
                    refs.append(ds_ref_domain_avg_yr)

            # I/O and plotting.
            # ------------------------------------------------------------------
            parameter.output_file = "-".join([var, region])
            fnm = os.path.join(
                _get_output_dir(parameter),
                parameter.output_file + ".json",
            )

            with open(fnm, "w") as outfile:
                json.dump(save_data, outfile)

            metrics_dict[region] = RefsTestMetrics(
                test=ds_test_domain_avg, refs=refs, metrics=[]
            )

            if parameter.save_netcdf:
                _write_to_netcdf(parameter, ds_test_domain_avg[var], var, "test")

        parameter.viewer_descr[var] = ds_test[var].attrs.get("long_name", var)
        area_mean_time_series_plot.plot(var, parameter, metrics_dict)

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


def _get_region_yearly_avg(
    parameter: AreaMeanTimeSeriesParameter,
    ds: xr.Dataset,
    ds_mask: xr.Dataset,
    var: str,
    region: str,
    data_type: Literal["test", "ref"],
) -> xr.Dataset:
    ds_new = ds.copy()

    if "land" in region or "ocean" in region:
        ds_new = _apply_land_sea_mask(
            ds_new,
            ds_mask,
            var,
            region,  # type: ignore
            parameter.regrid_tool,
            parameter.regrid_method,
        )

    if "global" not in region and data_type == "test":
        ds_new = _subset_on_region(ds_new, var, region)

    da_region_avg: xr.DataArray = spatial_avg(  # type: ignore
        ds_new, var, axis=["X", "Y"], as_list=False
    )

    # Convert the DataArray back to a Dataset to add time bounds and
    # to calculate the yearly means.
    ds_region_avg = da_region_avg.to_dataset()
    ds_region_avg = ds_region_avg.bounds.add_time_bounds("freq", freq="month")
    if data_type == "ref":
        _log_time(ds_region_avg, "ref")

    ds_region_avg_yr = ds_region_avg.temporal.group_average(var, freq="year")

    return ds_region_avg_yr


def _log_time(ds: xr.Dataset, data_type: str):
    time_coords = xc.get_dim_coords(ds, axis="T").values

    logger.info(
        f"Start and end time for selected time slices for {data_type} data: "
        f"{time_coords[0]} "
        f"{time_coords[1]}",
    )

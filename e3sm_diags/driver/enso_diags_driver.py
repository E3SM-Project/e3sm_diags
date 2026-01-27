from __future__ import annotations

import math
import os
from typing import TYPE_CHECKING, Literal, TypedDict

import numpy as np
import xarray as xr
import xcdat as xc
import xskillscore as xs

from e3sm_diags import INSTALL_PATH
from e3sm_diags.driver.utils.dataset_xr import Dataset
from e3sm_diags.driver.utils.io import (
    _save_data_metrics_and_plots,
    _write_vars_to_netcdf,
)
from e3sm_diags.driver.utils.regrid import _subset_on_region, align_grids_to_lower_res
from e3sm_diags.logger import _setup_child_logger
from e3sm_diags.metrics.metrics import correlation, rmse, spatial_avg, std
from e3sm_diags.plot.enso_diags_plot import plot_map, plot_scatter

if TYPE_CHECKING:
    from e3sm_diags.parameter.enso_diags_parameter import EnsoDiagsParameter


logger = _setup_child_logger(__name__)


class MetricsDictScatter(TypedDict):
    """A TypedDict class representing metrics for the scatter driver"""

    var: str
    units: str
    region: str
    test: xr.DataArray
    ref: xr.DataArray


class MetricsSubDict(TypedDict):
    """A TypedDict class representing the metrics sub-dictionary."""

    min: float
    max: float
    mean: list[float]
    std: list[float]
    rmse: float | None
    corr: float | None


# A type annotation representing the metrics dictionary.
MetricsDictMap = dict[str, MetricsSubDict | str]


def run_diag(parameter: EnsoDiagsParameter) -> EnsoDiagsParameter:
    if parameter.plot_type == "map":
        return run_diag_map(parameter)
    elif parameter.plot_type == "scatter":
        return run_diag_scatter(parameter)
    else:
        raise Exception("Invalid plot_type={}".format(parameter.plot_type))


def run_diag_map(parameter: EnsoDiagsParameter) -> EnsoDiagsParameter:
    variables = parameter.variables
    seasons = parameter.seasons
    regions = parameter.regions
    nino_region_str = parameter.nino_region
    run_type = parameter.run_type

    logger.info("run_type: {}".format(run_type))

    test_ds = Dataset(parameter, data_type="test")
    ref_ds = Dataset(parameter, data_type="ref")

    if run_type == "model_vs_model":
        da_test_nino = calculate_nino_index_model(test_ds, parameter, nino_region_str)
        da_ref_nino = calculate_nino_index_model(ref_ds, parameter, nino_region_str)
    elif run_type == "model_vs_obs":
        da_test_nino = calculate_nino_index_model(test_ds, parameter, nino_region_str)
        da_ref_nino = calculate_nino_index_obs(
            parameter, nino_region_str, data_type="ref"
        )
    else:
        raise Exception("Invalid run_type={}".format(run_type))

    for season in seasons:
        logger.info("Season: {}".format(season))
        parameter._set_name_yrs_attrs(test_ds, ref_ds, season)

        for var_key in variables:
            logger.info("Variable: {}".format(var_key))
            parameter.var_id = "{}-regression-over-nino".format(var_key)

            ds_test = test_ds.get_time_series_dataset(var_key)
            ds_ref = ref_ds.get_time_series_dataset(var_key)

            for region in regions:
                logger.info("Selected region: {}".format(region))

                ds_test_reg_coe, da_test_conf_lvls = calc_linear_regression(
                    ds_test, da_test_nino, var_key, region
                )
                ds_ref_reg_coe, da_ref_conf_lvls = calc_linear_regression(
                    ds_ref, da_ref_nino, var_key, region
                )

                (
                    ds_test_reg_coe_regrid,
                    ds_ref_reg_coe_regrid,
                ) = align_grids_to_lower_res(
                    ds_test_reg_coe,
                    ds_ref_reg_coe,
                    var_key,
                    parameter.regrid_tool,
                    parameter.regrid_method,
                )

                ds_diff_reg_coe = ds_test_reg_coe_regrid.copy()
                ds_diff_reg_coe[var_key] = (
                    ds_diff_reg_coe[var_key] - ds_ref_reg_coe_regrid[var_key]
                )

                metrics_dict = _create_metrics_dict(
                    ds_test_reg_coe,
                    ds_test_reg_coe_regrid,
                    ds_ref_reg_coe,
                    ds_ref_reg_coe_regrid,
                    ds_diff_reg_coe,
                    var_key,
                )

                if not parameter.contour_levels or not parameter.diff_levels:
                    contour_levels = _get_contour_levels(metrics_dict)
                    parameter.contour_levels = parameter.diff_levels = contour_levels

                parameter.main_title = (
                    f"Regression coefficient, {var_key} anomaly to {nino_region_str}"
                )
                parameter.output_file = (
                    f"regression-coefficient-{var_key.lower()}-over-"
                    f"{nino_region_str.lower()}"
                )
                parameter.var_region = region

                _save_data_metrics_and_plots(
                    parameter,
                    plot_map,
                    var_key,
                    ds_test_reg_coe,
                    ds_ref_reg_coe,
                    ds_diff_reg_coe,
                    metrics_dict,  # type: ignore
                    plot_kwargs={
                        "da_test_conf_lvls": da_test_conf_lvls,
                        "da_ref_conf_lvls": da_ref_conf_lvls,
                    },
                    viewer_descr=", ".join([parameter.main_title, region]),
                )

    return parameter


def run_diag_scatter(parameter: EnsoDiagsParameter) -> EnsoDiagsParameter:
    variables = parameter.variables
    run_type = parameter.run_type

    logger.info("run_type: {}".format(run_type))

    metrics_dict: MetricsDictScatter = {
        "var": "TS",
        "units": "degC",
        "region": "NINO3",
        "test": xr.DataArray(),
        "ref": xr.DataArray(),
    }

    test_ds = Dataset(parameter, data_type="test")
    ref_ds = Dataset(parameter, data_type="ref")

    parameter._set_name_yrs_attrs(test_ds, ref_ds, None)

    if run_type == "model_vs_model":
        metrics_dict["test"] = calculate_nino_index_model(
            test_ds, parameter, metrics_dict["region"]
        )
        metrics_dict["ref"] = calculate_nino_index_model(
            ref_ds, parameter, metrics_dict["region"]
        )
    elif run_type == "model_vs_obs":
        metrics_dict["test"] = calculate_nino_index_model(
            test_ds, parameter, metrics_dict["region"]
        )
        metrics_dict["ref"] = calculate_nino_index_obs(
            parameter, metrics_dict["region"], "ref"
        )
    else:
        raise Exception("Invalid run_type={}".format(run_type))

    for var_key in variables:
        logger.info("Variable: {}".format(var_key))
        if var_key == "TAUX":
            regions = ["NINO4"]
        else:
            regions = ["NINO3"]

        for region in regions:
            y: MetricsDictScatter = {
                "var": var_key,
                "region": region,
                "units": "",
                "test": xr.DataArray(),
                "ref": xr.DataArray(),
            }

            ds_test = test_ds.get_time_series_dataset(var_key)
            ds_ref = ref_ds.get_time_series_dataset(var_key)

            ds_test_region = _subset_on_region(ds_test, var_key, region)
            ds_ref_region = _subset_on_region(ds_ref, var_key, region)

            # Domain average
            ds_test_avg = ds_test_region.spatial.average(var_key)
            ds_ref_avg = ds_ref_region.spatial.average(var_key)

            # Get anomaly from annual cycle climatology
            y["test"] = ds_test_avg.temporal.departures(var_key, freq="month")[var_key]
            y["ref"] = ds_ref_avg.temporal.departures(var_key, freq="month")[var_key]
            y["units"] = y["test"].attrs["units"]

            if var_key == "TAUX":
                y["test"] *= 1000
                y["ref"] *= 1000
                y["units"] = "10^3 {}".format(y["units"])

            parameter.var_id = f"{y['var']}-feedback"

            title_tuple = (
                y["var"],
                y["region"],
                metrics_dict["var"],
                metrics_dict["region"],
            )
            parameter.main_title = "{} anomaly ({}) vs. {} anomaly ({})".format(
                *title_tuple
            )
            parameter.viewer_descr[y["var"]] = parameter.main_title
            parameter.output_file = "feedback-{}-{}-{}-{}".format(*title_tuple)

            plot_scatter(parameter, metrics_dict, y)

            if parameter.save_netcdf:
                _write_vars_to_netcdf(
                    parameter,
                    var_key,
                    y["test"].to_dataset(),
                    y["ref"].to_dataset(),
                    None,
                )

    return parameter


def calculate_nino_index_model(
    ds_obj: Dataset,
    parameter: EnsoDiagsParameter,
    nino_region_str: str,
) -> xr.DataArray:
    """Calculate nino index based on model output SST or TS.

    If neither of these is available, then use the built-in HadISST nino index.

    Parameters
    ----------
    ds_obj : Dataset
        The dataset object.
    parameter : EnsoDiagsParameter
        The parameter object.
    nino_region_str : str
        The nino region.

    Returns
    -------
    xr.DataArray
        The nino index based on the model output.

    Raises
    ------
    RunTimeError
        If the "SST" or "TS" variables are not found.
    """
    try:
        try:
            # Try sea surface temperature first.
            sst = ds_obj.get_time_series_dataset("SST")
            nino_var_key = "SST"
        except IOError as e1:
            if str(e1).startswith("No files found for target variable SST"):
                logger.info(
                    "Handling the following exception by looking for surface "
                    f"temperature: {e1}",
                )
                # Try surface temperature.
                sst = ds_obj.get_time_series_dataset("TS")
                nino_var_key = "TS"

                logger.info(
                    "Simulated sea surface temperature not found, using surface "
                    "temperature instead."
                )
            else:
                raise e1

        sst_nino = _subset_on_region(sst, nino_var_key, nino_region_str)

        # Domain average
        sst_avg = sst_nino.spatial.average(nino_var_key, axis=["X", "Y"])

        # Get anomaly from annual cycle climatology
        sst_avg_anomaly = sst_avg.temporal.departures(nino_var_key, freq="month")
        da_nino = sst_avg_anomaly[nino_var_key]

    except IOError as e2:
        logger.info(
            "Handling the following exception by trying built-in HadISST nino index "
            f"time series: {e2}"
        )
        da_nino = calculate_nino_index_obs(parameter, nino_region_str, ds_obj.data_type)
        logger.info(
            "Simulated surface temperature not found, using built-in HadISST nino "
            "index time series instead."
        )

    return da_nino


def calculate_nino_index_obs(
    parameter: EnsoDiagsParameter,
    nino_region_str: str,
    data_type: Literal["test", "ref"],
) -> xr.DataArray:
    """Calculate the nino index using default observational datasets.

    This function uses the default HadISST nino index time series from
    http://www.esrl.noaa.gov/psd/gcos_wgsp/Timeseries/ for observation datasets
    and when model output is not available.

    Relevant files are e3sm_diags/driver/default_diags/enso_NINO{3, 34, 4}.long.data

    Parameters
    ----------
    parameter : EnsoDiagsParameter
        The parameter object.
    nino_region_str : str
        The nino region.
    data_type : {"test", "ref"}
        The data type, either "test" or "ref".

    Returns
    -------
    xr.DataArray
        The nino index.

    Raises
    ------
    RuntimeError
        If the requested years are outside the SST observational records.
    """
    data_file = "".join(["enso_", nino_region_str, ".long.data"])
    nino_index_path = os.path.join(INSTALL_PATH, "enso_diags", data_file)

    # Load data up to year 2018 from 1870
    sst_orig = np.loadtxt(nino_index_path, skiprows=1, max_rows=149)

    if data_type == "test":
        start = int(parameter.test_start_yr)
        end = int(parameter.test_end_yr)
    elif data_type == "ref":
        start = int(parameter.ref_start_yr)
        end = int(parameter.ref_end_yr)

    sst_years = sst_orig[:, 0].astype(int)

    try:
        start_ind = np.where(sst_years == start)[0][0]
        end_ind = np.where(sst_years == end)[0][0]
    except Exception as err:
        msg = "Requested years are outside of available sst obs records."
        raise RuntimeError(msg) from err

    sst = sst_orig[start_ind : end_ind + 1, 1:]

    # Get anomaly from annual cycle climatology
    annual_cycle = np.mean(sst, axis=0)
    sst_anomaly = np.empty_like(sst)
    num_years = end - start + 1

    for iyr in range(num_years):
        sst_anomaly[iyr, :] = sst[iyr, :] - annual_cycle

    len_months = num_years * 12
    nino_index = np.reshape(sst_anomaly, len_months)

    da_nino_index = xr.DataArray(
        name="SST", data=nino_index, dims="time", coords={"time": np.arange(len_months)}
    )

    return da_nino_index


def calc_linear_regression(
    ds: xr.Dataset,
    da_nino: xr.DataArray,
    var_key: str,
    region: str,
) -> tuple[xr.Dataset, xr.DataArray]:
    """Calculate the linear regression between the variable and the nino index.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset containing the variable.
    da_nino : xr.DataArray
        The nino index.
    var_key : str
        The key of the variable.
    region : str
        The nino region.

    Returns
    -------
    tuple[xr.Dataset, xr.DataArray]
        A tuple containing the regression coefficient dataset and the
        confidence levels dataarray.
    """
    # Average over selected region, and average over months to get the yearly mean.
    domain = _subset_on_region(ds, var_key, region)

    # Get anomaly from annual cycle climatology
    anomaly = domain.temporal.departures(var_key, freq="month")
    anomaly_var = anomaly[var_key].copy()

    # Align the time coordinates to enable Xarray alignment before calculating
    # linear slope. This is necessary when the reference nino index is
    # from observation data.
    independent_var = da_nino.copy()
    independent_var = _align_time_coords(anomaly_var, independent_var)

    reg_coe = xs.linslope(independent_var, anomaly_var, keep_attrs=True)

    # Set confidence level to 1 if significant and 0 if not.
    # p-value < 5%, implies significance at 95% confidence level.
    conf_lvls = xs.pearson_r_p_value(independent_var, anomaly_var, keep_attrs=True)
    conf_lvls_masked = xr.where(conf_lvls < 0.05, 1, 0, keep_attrs=True)

    sst_units = "degC"
    reg_coe.attrs["units"] = "{}/{}".format(ds[var_key].units, sst_units)

    reg_coe.name = var_key
    conf_lvls_masked.name = var_key

    ds_reg_coe = reg_coe.to_dataset()
    ds_reg_coe = ds_reg_coe.bounds.add_missing_bounds(axes=["X", "Y", "T"])

    return ds_reg_coe, conf_lvls_masked


def _align_time_coords(da: xr.DataArray, da_nino_index: xr.DataArray) -> xr.DataArray:
    """Align the nino index time coordinates to the variable.

    The reference data is calculated using a nino index from observation data.
    The time coordinates of the nino index is updated to align with the
    reference data for downstream operations using Xarray broadcasting
    and alignment.

    Parameters
    ----------
    da_test : xr.DataArray
        The variable.
    da_nino_index : xr.DataArray
        The nino index variable.

    Returns
    -------
    xr.DataArray
        The nino index variable with time coordinates aligned to the variable.
    """
    time_coords = xc.get_dim_coords(da, axis="T")

    da_ref_new = da_nino_index.assign_coords({"time": time_coords})

    return da_ref_new


def _create_metrics_dict(
    ds_test: xr.Dataset,
    ds_test_regrid: xr.Dataset,
    ds_ref: xr.Dataset,
    ds_ref_regrid: xr.Dataset,
    ds_diff: xr.Dataset,
    var_key: str,
) -> MetricsDictMap:
    """Calculate metrics using the variable in the datasets.

    Metrics include min value, max value, spatial average (mean), and standard
    deviation.

    Parameters
    ----------
    var_key : str
        The variable key.
    ds_test : xr.Dataset
        The test dataset.
    ds_test_regrid : xr.Dataset
        The regridded test dataset.
    ds_ref : xr.Dataset
        The reference dataset.
    ds_ref_regrid : xr.Dataset
        The regridded reference dataset.
    ds_diff : xr.Dataset | None
        The difference between ``ds_test_regrid`` and ``ds_ref_regrid``.

    Returns
    -------
    MetricsDict
        A dictionary with the key being a string and the value being a
        sub-dictionary (key is metric and value is float).
    """
    metrics_dict: MetricsDictMap = {}

    metrics_dict["ref"] = get_metrics_subdict(ds_ref, var_key)
    metrics_dict["ref_regrid"] = get_metrics_subdict(ds_ref_regrid, var_key)
    metrics_dict["test"] = get_metrics_subdict(ds_test, var_key)
    metrics_dict["test_regrid"] = get_metrics_subdict(ds_test_regrid, var_key)
    metrics_dict["diff"] = get_metrics_subdict(ds_diff, var_key)
    metrics_dict["diff"]["rmse"] = rmse(ds_test_regrid, ds_ref_regrid, var_key)  # type: ignore
    metrics_dict["diff"]["corr"] = correlation(ds_test_regrid, ds_ref_regrid, var_key)  # type: ignore

    metrics_dict["unit"] = ds_test[var_key].units

    return metrics_dict


def get_metrics_subdict(ds: xr.Dataset, var_key: str) -> MetricsSubDict:
    """Get the metrics sub-dictionary.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset object.
    var_key : str
        The key of the variable to calculate metrics for.

    Returns
    -------
    MetricsSubDict
        A dictionary of metrics.
    """
    metrics_subdict: MetricsSubDict = {
        "min": ds[var_key].min().item(),
        "max": ds[var_key].max().item(),
        "mean": spatial_avg(ds, var_key, axis=["X", "Y"], as_list=True),  # type: ignore
        "std": std(ds, var_key),
    }
    return metrics_subdict


def _get_contour_levels(metrics_dict: MetricsDictMap) -> list[float]:
    """Get the contour levels for the map plot.

    The non-regridded data ("test" and "ref") are used for their min and max
    values.

    Parameters
    ----------
    metrics_dict : MetricsDictMap
        The metrics dictionary.

    Returns
    -------
    list[float]
        A list of floats representing contour levels.
    """
    min_contour_level = math.floor(
        min(
            metrics_dict["ref"]["min"],  # type: ignore
            metrics_dict["test"]["min"],  # type: ignore
        )
    )
    max_contour_level = math.ceil(
        max(
            metrics_dict["ref"]["max"],  # type: ignore
            metrics_dict["test"]["max"],  # type: ignore
        )
    )

    # Center on zero.
    bound = max(abs(max_contour_level), abs(min_contour_level))
    lower_bound = -bound
    upper_bound = bound
    contour_level_range = 2 * bound

    step_size = contour_level_range / 10
    if step_size > 1:
        step_size = int(step_size)
    elif step_size == 0:
        step_size = 1 / 10

    contour_levels = [
        float(x) for x in np.arange(lower_bound, upper_bound + 1, step_size)
    ]

    return contour_levels

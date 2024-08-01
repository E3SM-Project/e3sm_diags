from __future__ import annotations

import json
import math
import os
from typing import TYPE_CHECKING

import cdutil
import numpy as np
import scipy.stats
import xarray as xr
import xcdat as xc
import xskillscore as xs

import e3sm_diags
from e3sm_diags.derivations import default_regions
from e3sm_diags.driver import utils
from e3sm_diags.driver.utils.dataset_xr import Dataset
from e3sm_diags.driver.utils.regrid import _subset_on_region
from e3sm_diags.logger import custom_logger
from e3sm_diags.metrics import corr, max_cdms, mean, min_cdms, rmse, std
from e3sm_diags.plot.cartopy.enso_diags_plot import plot_map, plot_scatter

if TYPE_CHECKING:
    from e3sm_diags.parameter.enso_diags_parameter import EnsoDiagsParameter


logger = custom_logger(__name__)


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

        # FIXME: Time coordinates should align with da_test_nino.
        da_ref_nino = calculate_nino_index(nino_region_str, parameter, ref=True)
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

                test_reg_coe, test_conf_lvls = perform_regression(
                    ds_test, da_test_nino, var_key, region
                )

                # FIXME: These do not align with main.
                ref_reg_coe, ref_conf_lvls = perform_regression(
                    ds_ref, da_ref_nino, var_key, region
                )

                (
                    test_reg_coe_regrid,
                    ref_reg_coe_regrid,
                ) = utils.general.regrid_to_lower_res(
                    test_reg_coe,
                    ref_reg_coe,
                    parameter.regrid_tool,
                    parameter.regrid_method,
                )
                diff = test_reg_coe_regrid - ref_reg_coe_regrid

                # Metrics
                metrics_dict = create_metrics(
                    ref_reg_coe,
                    test_reg_coe,
                    ref_reg_coe_regrid,
                    test_reg_coe_regrid,
                    diff,
                )
                # If not defined, determine contour_levels
                if not parameter.contour_levels or not parameter.diff_levels:
                    # We want contour levels for the plot,
                    # which uses original (non-regridded) test and ref,
                    # so we use those min and max values.
                    min_contour_level = math.floor(
                        min(
                            metrics_dict["ref"]["min"],
                            metrics_dict["test"]["min"],
                        )
                    )
                    max_contour_level = math.ceil(
                        max(
                            metrics_dict["ref"]["max"],
                            metrics_dict["test"]["max"],
                        )
                    )
                    CENTER_ON_ZERO = True

                    if CENTER_ON_ZERO:
                        bound = max(abs(max_contour_level), abs(min_contour_level))
                        lower_bound = -bound
                        upper_bound = bound
                        contour_level_range = 2 * bound
                    else:
                        lower_bound = min_contour_level
                        upper_bound = max_contour_level
                        contour_level_range = max_contour_level - min_contour_level
                    step_size = contour_level_range / 10
                    if step_size > 1:
                        step_size = int(step_size)
                    elif step_size == 0:
                        step_size = 1 / 10
                    contour_levels = list(
                        np.arange(lower_bound, upper_bound + 1, step_size)
                    )
                    if not parameter.contour_levels:
                        parameter.contour_levels = contour_levels
                    if not parameter.diff_levels:
                        parameter.diff_levels = contour_levels

                parameter.main_title = (
                    "Regression coefficient, {} anomaly to {}".format(
                        var_key, nino_region_str
                    )
                )
                parameter.viewer_descr[var_key] = ", ".join(
                    [parameter.main_title, region]
                )
                parameter.output_file = "regression-coefficient-{}-over-{}".format(
                    var_key.lower(), nino_region_str.lower()
                )

                # Saving the metrics as a json.
                metrics_dict["unit"] = test_reg_coe_regrid.units
                metrics_output_file_name = os.path.join(
                    utils.general.get_output_dir(parameter.current_set, parameter),
                    parameter.output_file + ".json",
                )
                with open(metrics_output_file_name, "w") as outfile:
                    json.dump(metrics_dict, outfile)
                # Get the file name that the user has passed in and display that.
                metrics_output_file_name = os.path.join(
                    utils.general.get_output_dir(parameter.current_set, parameter),
                    parameter.output_file + ".json",
                )
                logger.info("Metrics saved in: {}".format(metrics_output_file_name))

                # Plot
                parameter.var_region = region
                # Plot original ref and test, not regridded versions.
                plot_map(
                    ref_reg_coe,
                    test_reg_coe,
                    diff,
                    metrics_dict,
                    ref_conf_lvls,
                    test_conf_lvls,
                    parameter,
                )
                utils.general.save_ncfiles(
                    parameter.current_set,
                    test_reg_coe,
                    ref_reg_coe,
                    diff,
                    parameter,
                )

    return parameter


def run_diag_scatter(parameter: EnsoDiagsParameter) -> EnsoDiagsParameter:
    variables = parameter.variables
    run_type = parameter.run_type

    # We will always use the same regions, so we don't do the following:
    # x['region'] = parameter.nino_region
    # y['region'] = parameter.regions[0]
    x = {"var": "TS", "units": "degC", "region": "NINO3"}
    test_data = utils.dataset.Dataset(parameter, test=True)
    ref_data = utils.dataset.Dataset(parameter, ref=True)

    if parameter.print_statements:
        logger.info("run_type: {}".format(run_type))
    if run_type == "model_vs_model":
        x["test"] = calculate_nino_index_model(test_data, parameter, x["region"])
        x["ref"] = calculate_nino_index_model(ref_data, parameter, x["region"])
    elif run_type == "model_vs_obs":
        x["test"] = calculate_nino_index_model(test_data, parameter, x["region"])
        x["ref"] = calculate_nino_index(x["region"], parameter, ref=True)
    else:
        raise Exception("Invalid run_type={}".format(run_type))

    parameter.test_name_yrs = utils.general.get_name_and_yrs(parameter, test_data)
    parameter.ref_name_yrs = utils.general.get_name_and_yrs(parameter, ref_data)

    for y_var in variables:
        if y_var == "TAUX":
            regions = ["NINO4"]
        else:
            regions = ["NINO3"]
        for region in regions:
            y = {"var": y_var, "region": region}
            test_data_ts = test_data.get_timeseries_variable(y_var)
            ref_data_ts = ref_data.get_timeseries_variable(y_var)
            y_region = default_regions.regions_specs[region]["domain"]  # type: ignore
            test_data_ts_regional = test_data_ts(y_region)
            ref_data_ts_regional = ref_data_ts(y_region)
            # Domain average
            test_avg = cdutil.averager(test_data_ts_regional, axis="xy")
            ref_avg = cdutil.averager(ref_data_ts_regional, axis="xy")
            # Get anomaly from annual cycle climatology
            y["test"] = cdutil.ANNUALCYCLE.departures(test_avg)
            y["ref"] = cdutil.ANNUALCYCLE.departures(ref_avg)
            y["units"] = test_avg.units
            if y_var == "TAUX":
                y["test"] *= 1000
                y["ref"] *= 1000
                y["units"] = "10^3 {}".format(y["units"])

            parameter.var_id = "{}-feedback".format(y["var"])
            title_tuple = (y["var"], y["region"], x["var"], x["region"])
            parameter.main_title = "{} anomaly ({}) vs. {} anomaly ({})".format(
                *title_tuple
            )
            parameter.viewer_descr[y["var"]] = parameter.main_title
            parameter.output_file = "feedback-{}-{}-{}-{}".format(*title_tuple)

            plot_scatter(x, y, parameter)

    return parameter


def calculate_nino_index(
    nino_region_str: str, parameter: EnsoDiagsParameter, test=False, ref=False
) -> xr.DataArray:
    """
    Use the built-in HadISST nino index time series from http://www.esrl.noaa.gov/psd/gcos_wgsp/Timeseries/
    for observation datasets and when model output is not available.

    Relevant files are e3sm_diags/driver/default_diags/enso_NINO{3, 34, 4}.long.data
    """
    data_file = "".join(["enso_", nino_region_str, ".long.data"])
    nino_index_path = os.path.join(e3sm_diags.INSTALL_PATH, "enso_diags", data_file)

    # load data up to year 2018 from 1870
    sst_orig = np.loadtxt(nino_index_path, skiprows=1, max_rows=149)

    if test:
        start = int(parameter.test_start_yr)
        end = int(parameter.test_end_yr)
    elif ref:
        start = int(parameter.ref_start_yr)
        end = int(parameter.ref_end_yr)
    else:
        start = int(parameter.start_yr)
        end = int(parameter.end_yr)
    sst_years = sst_orig[:, 0].astype(int)

    try:
        start_ind = np.where(sst_years == start)[0][0]
        end_ind = np.where(sst_years == end)[0][0]
    except Exception:
        msg = "Requested years are outside of available sst obs records."
        raise RuntimeError(msg)

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


def calculate_nino_index_model(
    ds_obj: Dataset,
    parameter: EnsoDiagsParameter,
    nino_region_str: str,
) -> xr.DataArray:
    """
    Calculate nino index based on model output SST or TS. If neither of these is available,
    then use the built-in HadISST nino index.
    """

    try:
        try:
            # Try sea surface temperature first.
            sst = ds_obj.get_time_series_dataset("SST")
            nino_var_key = "SST"
        except RuntimeError as e1:
            if str(e1).startswith("Neither does SST nor the variables in"):
                logger.info(
                    "Handling the following exception by looking for surface "
                    f"temperature: {e1}",
                )
                # Try surface temperature.
                sst = ds_obj.get_time_series_dataset("TS")
                nino_var_key = "TS"

                logger.info(
                    "Simulated sea surface temperature not found, using surface temperature instead."
                )
            else:
                raise e1

        sst_nino = _subset_on_region(sst, nino_var_key, nino_region_str)

        # Domain average
        sst_avg = sst_nino.spatial.average(nino_var_key, axis=["X", "Y"])

        # Get anomaly from annual cycle climatology
        sst_avg_anomaly = sst_avg.temporal.departures(nino_var_key, freq="month")
        da_nino = sst_avg_anomaly[nino_var_key]

    except RuntimeError as e2:
        logger.info(
            "Handling the following exception by trying built-in HadISST nino index "
            f"time series: {e2}"
        )
        test = ds_obj.test
        ref = ds_obj.ref
        da_nino = calculate_nino_index(nino_region_str, parameter, test=test, ref=ref)
        logger.info(
            "Simulated surface temperature not found, using built-in HadISST nino index time series instead."
        )

    return da_nino


def perform_regression(
    ds: xr.Dataset,
    da_nino_index: xr.DataArray,
    var_key: str,
    region: str,
):
    # Average over selected region, and average over months to get the yearly mean.
    domain = _subset_on_region(ds, var_key, region)

    # Get anomaly from annual cycle climatology
    anomaly = domain.temporal.departures(var_key, freq="month")
    anomaly_var = anomaly[var_key].copy()

    time_dim = xc.get_dim_keys(anomaly_var, axis="T")
    reg_coe = anomaly_var.isel({time_dim: 0}).drop(time_dim).copy()
    confidence_levels = anomaly_var.isel({time_dim: 0}).drop(time_dim).copy()

    lat = xc.get_dim_coords(anomaly_var, axis="Y")
    lon = xc.get_dim_coords(anomaly_var, axis="X")
    nlat = len(lat)
    nlon = len(lon)

    independent_var = da_nino_index.copy()

    for ilat in range(nlat):
        for ilon in range(nlon):
            dependent_var = anomaly_var[:, ilat, ilon]

            slope, _, _, pvalue, _ = scipy.stats.linregress(
                independent_var, dependent_var
            )
            reg_coe[ilat, ilon] = slope

            # Set confidence level to 1 if significant and 0 if not.
            if pvalue < 0.05:
                # p-value < 5%, implies significance at 95% confidence level.
                confidence_levels[ilat, ilon] = 1
            else:
                confidence_levels[ilat, ilon] = 0

    sst_units = "degC"
    reg_coe.attrs["units"] = "{}/{}".format(ds[var_key].units, sst_units)

    return reg_coe, confidence_levels


def perform_regression_new(
    ds: xr.Dataset,
    ds_nino: xr.DataArray,
    var_key: str,
    region: str,
):
    # Average over selected region, and average over months to get the yearly mean.
    domain = _subset_on_region(ds, var_key, region)

    # Get anomaly from annual cycle climatology
    anomaly = domain.temporal.departures(var_key, freq="month")
    anomaly_var = anomaly[var_key].copy()

    independent_var = ds_nino.copy()
    reg_coe = xs.linslope(independent_var, anomaly_var, keep_attrs=True)

    # Set confidence level to 1 if significant and 0 if not.
    # p-value < 5%, implies significance at 95% confidence level.
    conf_lvls = xs.pearson_r_p_value(independent_var, anomaly_var, keep_attrs=True)
    conf_lvls_masked = xr.where(conf_lvls < 0.05, 1, 0, keep_attrs=True)

    return reg_coe, conf_lvls_masked


def create_single_metrics_dict(values):
    d = {
        "min": float(min_cdms(values)),
        "max": float(max_cdms(values)),
        "mean": float(mean(values)),
        "std": float(std(values)),
    }
    return d


def create_metrics(ref, test, ref_regrid, test_regrid, diff):
    """Creates the relevant metrics in a dictionary"""
    metrics_dict = {}

    metrics_dict["ref"] = create_single_metrics_dict(ref)
    metrics_dict["ref_regrid"] = create_single_metrics_dict(ref_regrid)
    metrics_dict["test"] = create_single_metrics_dict(test)
    metrics_dict["test_regrid"] = create_single_metrics_dict(test_regrid)
    metrics_dict["diff"] = create_single_metrics_dict(diff)
    d = metrics_dict["diff"]
    d["rmse"] = float(rmse(test_regrid, ref_regrid))
    d["corr"] = float(corr(test_regrid, ref_regrid))

    return metrics_dict

from __future__ import annotations

import json
import math
import os
from typing import TYPE_CHECKING, Literal, TypedDict

import numpy as np
import xarray as xr
import xcdat as xc
import xskillscore as xs

from e3sm_diags import INSTALL_PATH
from e3sm_diags.driver.utils.arithmetic import subtract_dataarrays
from e3sm_diags.driver.utils.climo_xr import climo
from e3sm_diags.driver.utils.dataset_xr import Dataset
from e3sm_diags.driver.utils.io import (
    _get_output_dir,
    _get_output_filename_filepath,
    _save_data_metrics_and_plots,
    _write_vars_to_netcdf,
)
from e3sm_diags.driver.utils.regrid import _subset_on_region, align_grids_to_lower_res
from e3sm_diags.logger import _setup_child_logger
from e3sm_diags.metrics.metrics import correlation, rmse, spatial_avg, std
from e3sm_diags.plot.enso_diags_plot import (
    plot_equatorial_soi,
    plot_interannual_variability,
    plot_lead_lag,
    plot_map,
    plot_nino_index_timeseries,
    plot_scatter,
    plot_seasonality,
)

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
    elif parameter.plot_type == "nino_index_timeseries":
        return run_diag_nino_index_timeseries(parameter)
    elif parameter.plot_type == "seasonality":
        return run_diag_seasonality(parameter)
    elif parameter.plot_type == "interannual_variability":
        return run_diag_interannual_variability(parameter)
    elif parameter.plot_type == "equatorial_soi":
        return run_diag_equatorial_soi(parameter)
    elif parameter.plot_type == "lead_lag":
        return run_diag_lead_lag(parameter)
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
                ds_diff_reg_coe[var_key] = subtract_dataarrays(
                    ds_diff_reg_coe[var_key], ds_ref_reg_coe_regrid[var_key]
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


def run_diag_nino_index_timeseries(parameter: EnsoDiagsParameter) -> EnsoDiagsParameter:
    """Plot the monthly nino index anomaly time series for each nino region.

    The nino regions (e.g. NINO3, NINO34, NINO4) are stacked as subplot rows,
    with the test case in the left column and the reference in the right column.
    This is a port of a-prime's ``plot_multiple_index`` diagnostic.
    """
    regions = parameter.nino_regions
    run_type = parameter.run_type

    logger.info("run_type: {}".format(run_type))

    test_ds = Dataset(parameter, data_type="test")
    ref_ds = Dataset(parameter, data_type="ref")

    parameter._set_name_yrs_attrs(test_ds, ref_ds, None)

    test_indices: dict[str, xr.DataArray] = {}
    ref_indices: dict[str, xr.DataArray] = {}

    for region in regions:
        logger.info("Nino region: {}".format(region))

        test_indices[region] = calculate_nino_index_model(test_ds, parameter, region)

        if run_type == "model_vs_model":
            ref_indices[region] = calculate_nino_index_model(ref_ds, parameter, region)
        elif run_type == "model_vs_obs":
            ref_indices[region] = calculate_nino_index_obs(parameter, region, "ref")
        else:
            raise Exception("Invalid run_type={}".format(run_type))

    parameter.var_id = "NINO-index"
    parameter.main_title = "Nino index time series"
    parameter.output_file = "nino-index-timeseries"

    # The viewer adds one row per variable; the index is SST/TS based, so use a
    # single representative key for the description.
    var_key = parameter.variables[0] if parameter.variables else "TS"
    parameter.viewer_descr[var_key] = parameter.main_title

    plot_nino_index_timeseries(parameter, test_indices, ref_indices)

    if parameter.save_netcdf:
        ds_test_out = xr.Dataset({region: test_indices[region] for region in regions})
        ds_ref_out = xr.Dataset({region: ref_indices[region] for region in regions})

        _, test_filepath = _get_output_filename_filepath(parameter, "test")
        _, ref_filepath = _get_output_filename_filepath(parameter, "ref")

        ds_test_out.to_netcdf(test_filepath)
        ds_ref_out.to_netcdf(ref_filepath)

        logger.info(f"Nino index test output saved in: {test_filepath}")
        logger.info(f"Nino index ref output saved in: {ref_filepath}")

    return parameter


def run_diag_seasonality(parameter: EnsoDiagsParameter) -> EnsoDiagsParameter:
    """Plot the seasonality (per-calendar-month std dev) of the nino index.

    The nino regions (e.g. NINO3, NINO34, NINO4) are stacked as subplot rows,
    with the test and reference standard deviations overlaid in each panel. This
    is a port of a-prime's ``plot_multiple_index_seasonality`` diagnostic.
    """
    regions = parameter.nino_regions
    run_type = parameter.run_type

    logger.info("run_type: {}".format(run_type))

    test_ds = Dataset(parameter, data_type="test")
    ref_ds = Dataset(parameter, data_type="ref")

    parameter._set_name_yrs_attrs(test_ds, ref_ds, None)

    test_seasonality: dict[str, np.ndarray] = {}
    ref_seasonality: dict[str, np.ndarray] = {}

    for region in regions:
        logger.info("Nino region: {}".format(region))

        da_test = calculate_nino_index_model(test_ds, parameter, region)

        if run_type == "model_vs_model":
            da_ref = calculate_nino_index_model(ref_ds, parameter, region)
        elif run_type == "model_vs_obs":
            da_ref = calculate_nino_index_obs(parameter, region, "ref")
        else:
            raise Exception("Invalid run_type={}".format(run_type))

        test_seasonality[region] = _calculate_seasonality(da_test)
        ref_seasonality[region] = _calculate_seasonality(da_ref)

    parameter.var_id = "NINO-index"
    parameter.main_title = "ENSO Seasonality: NINO index"
    parameter.output_file = "nino-index-seasonality"

    var_key = parameter.variables[0] if parameter.variables else "TS"
    parameter.viewer_descr[var_key] = parameter.main_title

    plot_seasonality(parameter, test_seasonality, ref_seasonality)

    if parameter.save_netcdf:
        months = np.arange(1, 13)
        ds_test_out = xr.Dataset(
            {region: ("month", test_seasonality[region]) for region in regions},
            coords={"month": months},
        )
        ds_ref_out = xr.Dataset(
            {region: ("month", ref_seasonality[region]) for region in regions},
            coords={"month": months},
        )

        _, test_filepath = _get_output_filename_filepath(parameter, "test")
        _, ref_filepath = _get_output_filename_filepath(parameter, "ref")

        ds_test_out.to_netcdf(test_filepath)
        ds_ref_out.to_netcdf(ref_filepath)

        logger.info(f"Nino index seasonality test output saved in: {test_filepath}")
        logger.info(f"Nino index seasonality ref output saved in: {ref_filepath}")

    return parameter


def _calculate_seasonality(da_index: xr.DataArray) -> np.ndarray:
    """Compute the per-calendar-month standard deviation of a nino index.

    Returns a 12-element array (January through December) of the standard
    deviation across years for each calendar month, matching a-prime's
    seasonality computation.

    Parameters
    ----------
    da_index : xr.DataArray
        The monthly nino index anomaly time series. The model index carries a
        datetime time coordinate; the built-in observational index uses an
        integer month counter that starts in January.

    Returns
    -------
    np.ndarray
        The 12 monthly standard deviations (Jan..Dec).
    """
    time_coords = np.asarray(xc.get_dim_coords(da_index, axis="T"))

    if np.issubdtype(time_coords.dtype, np.number):
        # Observational index: a plain month counter starting in January.
        values = da_index.values
        n_years = values.shape[0] // 12
        return np.std(values[: n_years * 12].reshape(n_years, 12), axis=0)

    # Model index: group by calendar month using the datetime coordinate.
    std_by_month = da_index.groupby("time.month").std("time")
    return std_by_month.sortby("month").values


# The two sea level pressure regions whose standardized anomalies are
# differenced to form the Equatorial SOI (eastern equatorial Pacific minus the
# equatorial Indian Ocean), matching a-prime's eqsoi diagnostic.
EQSOI_REGIONS = ("EPAC", "INDO")

# The nino index region overlaid with the Equatorial SOI, and the label used
# for it in the figure.
EQSOI_NINO_REGION = "NINO34"
EQSOI_NINO_LABEL = "NINO3.4"


def run_diag_equatorial_soi(parameter: EnsoDiagsParameter) -> EnsoDiagsParameter:
    """Plot the Equatorial SOI and Nino3.4 monthly index time series.

    Two indices are stacked as subplot rows -- the Equatorial SOI (EQSOI),
    computed from sea level pressure, and the Nino3.4 SST anomaly index -- with
    the test case in the left column and the reference in the right column. The
    correlation between the two indices is annotated in each panel. This is a
    port of a-prime's ``plot_multiple_index_same_plot`` (the ``SOI_Nino`` set).
    """
    run_type = parameter.run_type
    logger.info("run_type: {}".format(run_type))

    test_ds = Dataset(parameter, data_type="test")
    ref_ds = Dataset(parameter, data_type="ref")

    parameter._set_name_yrs_attrs(test_ds, ref_ds, None)

    # EQSOI from PSL; Nino3.4 as the raw SST anomaly (not standardized), matching
    # a-prime's plot-time ``--stdize 0 0``.
    test_eqsoi = calculate_eqsoi_index(test_ds, parameter)
    test_nino = calculate_nino_index_model(test_ds, parameter, EQSOI_NINO_REGION)

    if run_type == "model_vs_model":
        ref_eqsoi = calculate_eqsoi_index(ref_ds, parameter)
        ref_nino = calculate_nino_index_model(ref_ds, parameter, EQSOI_NINO_REGION)
    elif run_type == "model_vs_obs":
        ref_eqsoi = calculate_eqsoi_index(ref_ds, parameter)
        ref_nino = calculate_nino_index_obs(parameter, EQSOI_NINO_REGION, "ref")
    else:
        raise Exception("Invalid run_type={}".format(run_type))

    # Ensure the nino index carries a unit for the y-axis label.
    test_nino.attrs.setdefault("units", "degC")
    ref_nino.attrs.setdefault("units", "degC")

    test_indices = {"EQSOI": test_eqsoi, EQSOI_NINO_LABEL: test_nino}
    ref_indices = {"EQSOI": ref_eqsoi, EQSOI_NINO_LABEL: ref_nino}

    correlations = {
        "test": _index_correlation(test_eqsoi, test_nino),
        "ref": _index_correlation(ref_eqsoi, ref_nino),
    }

    parameter.var_id = "EQSOI"
    parameter.main_title = "Equatorial SOI and Nino3.4 monthly indices"
    parameter.output_file = "equatorial-soi"

    var_key = parameter.variables[0] if parameter.variables else "PSL"
    parameter.viewer_descr[var_key] = parameter.main_title

    plot_equatorial_soi(parameter, test_indices, ref_indices, correlations)

    if parameter.save_netcdf:
        ds_test_out = xr.Dataset(
            {
                label: ("time", np.asarray(test_indices[label].values))
                for label in test_indices
            }
        )
        ds_ref_out = xr.Dataset(
            {
                label: ("time", np.asarray(ref_indices[label].values))
                for label in ref_indices
            }
        )

        _, test_filepath = _get_output_filename_filepath(parameter, "test")
        _, ref_filepath = _get_output_filename_filepath(parameter, "ref")

        ds_test_out.to_netcdf(test_filepath)
        ds_ref_out.to_netcdf(ref_filepath)

        logger.info(f"Equatorial SOI test output saved in: {test_filepath}")
        logger.info(f"Equatorial SOI ref output saved in: {ref_filepath}")

    return parameter


def calculate_eqsoi_index(
    ds_obj: Dataset, parameter: EnsoDiagsParameter
) -> xr.DataArray:
    """Calculate the Equatorial Southern Oscillation Index (EQSOI) from PSL.

    The EQSOI is the difference between the standardized monthly sea level
    pressure anomalies of the eastern equatorial Pacific (EPAC) and the
    equatorial Indian Ocean / Indonesia (INDO). For each region the area-average
    monthly anomaly is computed (annual cycle removed), then standardized
    (subtract the mean and divide by the standard deviation). This is a port of
    a-prime's ``compute_diff_index`` applied to PSL over the EPAC and INDO
    regions.

    Parameters
    ----------
    ds_obj : Dataset
        The dataset object.
    parameter : EnsoDiagsParameter
        The parameter object.

    Returns
    -------
    xr.DataArray
        The (unitless) monthly EQSOI time series.
    """
    ds_psl = ds_obj.get_time_series_dataset("PSL")

    standardized: list[xr.DataArray] = []
    for region in EQSOI_REGIONS:
        psl_region = _subset_on_region(ds_psl, "PSL", region)

        # Area average over the region.
        psl_avg = psl_region.spatial.average("PSL", axis=["X", "Y"])

        # Anomaly from the annual cycle climatology.
        psl_anomaly = psl_avg.temporal.departures("PSL", freq="month")
        da = psl_anomaly["PSL"]

        # Standardize the anomaly time series (population std, matching a-prime).
        standardized.append((da - da.mean()) / da.std())

    index = standardized[0] - standardized[1]
    index.name = "EQSOI"
    index.attrs["units"] = "unitless"

    return index


def _index_correlation(da_a: xr.DataArray, da_b: xr.DataArray) -> float:
    """Pearson correlation between two index time series.

    The shorter length is used defensively in case the two indices differ in
    length (e.g. mismatched observational records).
    """
    a = np.asarray(da_a.values)
    b = np.asarray(da_b.values)
    n = min(a.shape[0], b.shape[0])

    return float(np.corrcoef(a[:n], b[:n])[0, 1])


def run_diag_lead_lag(parameter: EnsoDiagsParameter) -> EnsoDiagsParameter:
    """Plot the lead-lag regression and correlation of a field on the nino index.

    For each lag (in months) the field anomaly is regressed on (and correlated
    with) the nino index, producing tiled maps with the lags as rows and the
    test case, reference, and their difference as columns. Two figures are made
    per field/region: the regression coefficients (with significance hatching on
    the test and reference columns) and the correlations. Positive lags indicate
    the nino index leading the field. This is a port of a-prime's
    ``plot_regress_lead_lag_index_field`` ("ENSO Evolution") diagnostic.
    """
    variables = parameter.variables
    regions = parameter.regions
    run_type = parameter.run_type
    nino_region = parameter.nino_region
    lags = parameter.lead_lag_months

    logger.info("run_type: {}".format(run_type))

    test_ds = Dataset(parameter, data_type="test")
    ref_ds = Dataset(parameter, data_type="ref")

    parameter._set_name_yrs_attrs(test_ds, ref_ds, None)

    da_test_nino = calculate_nino_index_model(test_ds, parameter, nino_region)
    if run_type == "model_vs_model":
        da_ref_nino = calculate_nino_index_model(ref_ds, parameter, nino_region)
    elif run_type == "model_vs_obs":
        da_ref_nino = calculate_nino_index_obs(parameter, nino_region, "ref")
    else:
        raise Exception("Invalid run_type={}".format(run_type))

    for var_key in variables:
        logger.info("Variable: {}".format(var_key))

        ds_test = test_ds.get_time_series_dataset(var_key)
        ds_ref = ref_ds.get_time_series_dataset(var_key)

        for region in regions:
            logger.info("Selected region: {}".format(region))
            parameter.var_region = region

            regr_panels: list[dict] = []
            corr_panels: list[dict] = []
            for lag in lags:
                logger.info("Lag (months): {}".format(lag))

                test_regr, test_corr, test_conf = _calc_lead_lag_regression(
                    ds_test, da_test_nino, var_key, region, lag
                )
                ref_regr, ref_corr, ref_conf = _calc_lead_lag_regression(
                    ds_ref, da_ref_nino, var_key, region, lag
                )

                regr_panels.append(
                    {
                        "lag": lag,
                        "test": test_regr[var_key],
                        "ref": ref_regr[var_key],
                        "diff": _regrid_and_diff(
                            test_regr, ref_regr, var_key, parameter
                        ),
                        "test_conf": test_conf,
                        "ref_conf": ref_conf,
                    }
                )
                corr_panels.append(
                    {
                        "lag": lag,
                        "test": test_corr[var_key],
                        "ref": ref_corr[var_key],
                        "diff": _regrid_and_diff(
                            test_corr, ref_corr, var_key, parameter
                        ),
                    }
                )

            parameter.var_id = var_key

            for kind, kind_panels in [
                ("regression", regr_panels),
                ("correlation", corr_panels),
            ]:
                parameter.main_title = (
                    f"ENSO Evolution: lead-lag {kind}, {var_key} on {nino_region} "
                    f"({region})"
                )
                parameter.output_file = "lead-lag-{}-{}-{}".format(
                    kind, var_key.lower(), region.lower()
                )
                plot_lead_lag(parameter, var_key, kind_panels, kind)
                _save_lead_lag_outputs(parameter, var_key, lags, kind_panels)

                # Record the figure so the viewer can link each one (a single
                # run produces both the regression and correlation figures). The
                # row label must be unique per figure, otherwise the viewer
                # collapses them into one row and links only the last one.
                parameter.lead_lag_entries.append(
                    {
                        "row": "{} ({})".format(parameter.case_id, kind),
                        "output_file": parameter.output_file,
                        "descr": parameter.main_title,
                    }
                )

    return parameter


def _regrid_and_diff(
    ds_test: xr.Dataset,
    ds_ref: xr.Dataset,
    var_key: str,
    parameter: EnsoDiagsParameter,
) -> xr.DataArray:
    """Regrid the test and reference fields to a common grid and difference them.

    a-prime differences the regression/correlation matrices directly because
    both cases were interpolated to a common grid; e3sm_diags regrids the test
    and reference to the lower-resolution grid first.
    """
    ds_test_regrid, ds_ref_regrid = align_grids_to_lower_res(
        ds_test,
        ds_ref,
        var_key,
        parameter.regrid_tool,
        parameter.regrid_method,
    )

    return subtract_dataarrays(ds_test_regrid[var_key], ds_ref_regrid[var_key])


def _calc_lead_lag_regression(
    ds: xr.Dataset,
    da_nino: xr.DataArray,
    var_key: str,
    region: str,
    lag: int,
) -> tuple[xr.Dataset, xr.Dataset, xr.DataArray]:
    """Regress and correlate a field anomaly on the nino index at a given lag.

    The field anomaly (annual cycle removed) and the nino index are sliced
    relative to each other by ``lag`` months: for ``lag >= 0`` the field at
    ``t + lag`` is regressed on the index at ``t`` (the index leads); for
    ``lag < 0`` the field leads. This mirrors a-prime's ``regress_index_field``.

    Returns the regression-coefficient dataset, the correlation dataset, and the
    significance mask (1 where p-value < 0.05).
    """
    domain = _subset_on_region(ds, var_key, region)

    # Anomaly from the annual cycle climatology.
    anomaly = domain.temporal.departures(var_key, freq="month")
    anomaly_var = anomaly[var_key].copy()

    time_dim = xc.get_dim_keys(anomaly_var, axis="T")

    # Align the nino index time coordinates to the field for Xarray alignment.
    index = _align_time_coords(anomaly_var, da_nino.copy())

    nt = anomaly_var.sizes[time_dim]
    if lag >= 0:
        x = index.isel({time_dim: slice(0, nt - lag)})
        y = anomaly_var.isel({time_dim: slice(lag, nt)})
    else:
        x = index.isel({time_dim: slice(-lag, nt)})
        y = anomaly_var.isel({time_dim: slice(0, nt + lag)})

    # The two windows cover different absolute times; reset to a common counter
    # so xskillscore aligns them positionally.
    n = x.sizes[time_dim]
    x = x.assign_coords({time_dim: np.arange(n)})
    y = y.assign_coords({time_dim: np.arange(n)})

    reg_coe = xs.linslope(x, y, dim=time_dim, keep_attrs=True)
    corr = xs.pearson_r(x, y, dim=time_dim, keep_attrs=True)
    p_val = xs.pearson_r_p_value(x, y, dim=time_dim)
    conf = xr.where(p_val < 0.05, 1, 0, keep_attrs=False)

    reg_coe.attrs["units"] = "{}/{}".format(ds[var_key].units, "degC")
    corr.attrs["units"] = "unitless"

    reg_coe.name = var_key
    corr.name = var_key
    conf.name = var_key

    ds_reg_coe = reg_coe.to_dataset().bounds.add_missing_bounds(axes=["X", "Y"])
    ds_corr = corr.to_dataset().bounds.add_missing_bounds(axes=["X", "Y"])

    return ds_reg_coe, ds_corr, conf


def _save_lead_lag_outputs(
    parameter: EnsoDiagsParameter,
    var_key: str,
    lags: list[int],
    panels: list[dict],
) -> None:
    """Save the test / ref / diff lead-lag fields (stacked over lag) to netCDF."""
    if not parameter.save_netcdf:
        return

    for data_type in ("test", "ref", "diff"):
        da = xr.concat([panel[data_type] for panel in panels], dim="lag")
        da = da.assign_coords(lag=lags)

        _, filepath = _get_output_filename_filepath(parameter, data_type)
        da.to_dataset(name=var_key).to_netcdf(filepath)
        logger.info(f"Lead-lag {data_type} output saved in: {filepath}")


# Number of standard deviations and contour levels used to derive the
# climatology / interannual-std map contour ranges, matching the defaults of
# a-prime's plot_stddev ``compute_contour_levels``.
INTERANNUAL_N_STDDEV = 5
INTERANNUAL_NUM_LEVELS = 11


def run_diag_interannual_variability(
    parameter: EnsoDiagsParameter,
) -> EnsoDiagsParameter:
    """Plot the climatology and interannual standard deviation maps.

    For each variable and region the six panels of a-prime's ``plot_stddev``
    are produced as two three-panel map figures: the climatological mean
    (test/ref/diff) and the interannual standard deviation (test/ref/diff). This
    is a port of a-prime's ``plot_stddev`` /
    ``compute_reg_seasonal_climo_and_stddev`` diagnostic.
    """
    variables = parameter.variables
    regions = parameter.regions
    run_type = parameter.run_type

    logger.info("run_type: {}".format(run_type))

    test_ds = Dataset(parameter, data_type="test")
    ref_ds = Dataset(parameter, data_type="ref")

    parameter._set_name_yrs_attrs(test_ds, ref_ds, None)

    for var_key in variables:
        logger.info("Variable: {}".format(var_key))

        ds_test = test_ds.get_time_series_dataset(var_key)
        ds_ref = ref_ds.get_time_series_dataset(var_key)

        for region in regions:
            logger.info("Selected region: {}".format(region))
            parameter.var_region = region

            test_stats = _compute_climo_and_interannual_std(ds_test, var_key, region)
            ref_stats = _compute_climo_and_interannual_std(ds_ref, var_key, region)

            # Build the test/ref/diff fields and contour levels for each of the
            # two columns (climatological mean and interannual std dev).
            panels: dict[str, dict] = {}
            metrics: dict[str, MetricsDictMap] = {}
            for stat in ("mean", "std"):
                ds_test_stat = test_stats[stat]
                ds_ref_stat = ref_stats[stat]

                ds_test_regrid, ds_ref_regrid = align_grids_to_lower_res(
                    ds_test_stat,
                    ds_ref_stat,
                    var_key,
                    parameter.regrid_tool,
                    parameter.regrid_method,
                )

                ds_diff = ds_test_regrid.copy()
                ds_diff[var_key] = subtract_dataarrays(
                    ds_test_regrid[var_key], ds_ref_regrid[var_key]
                )

                metrics[stat] = _create_metrics_dict(
                    ds_test_stat,
                    ds_test_regrid,
                    ds_ref_stat,
                    ds_ref_regrid,
                    ds_diff,
                    var_key,
                )

                panels[stat] = {
                    "test": ds_test_stat[var_key],
                    "ref": ds_ref_stat[var_key],
                    "diff": ds_diff[var_key],
                    "metrics": metrics[stat],
                    # Derive the contour levels from the reference field so the
                    # test and reference panels share the same scale, matching
                    # a-prime.
                    "levels": _compute_variability_contour_levels(
                        ds_ref_stat[var_key].values
                    ),
                    "diff_levels": _compute_diff_contour_levels(
                        ds_ref_stat[var_key].values
                    ),
                }

            parameter.var_id = var_key
            parameter.main_title = (
                f"Climatology and Interannual Std. Dev., {var_key} ({region})"
            )
            parameter.output_file = "interannual-variability-{}-{}".format(
                var_key.lower(), region.lower()
            )
            parameter.viewer_descr[var_key] = parameter.main_title

            plot_interannual_variability(parameter, var_key, panels)

            _save_interannual_variability_outputs(parameter, var_key, panels, metrics)

    return parameter


def _save_interannual_variability_outputs(
    parameter: EnsoDiagsParameter,
    var_key: str,
    panels: dict[str, dict],
    metrics: dict[str, MetricsDictMap],
) -> None:
    """Save the metrics JSON and (optionally) netCDF for the variability figure.

    The mean and std fields are written as separate variables (``<var>_mean``
    and ``<var>_std``) in the standard ``_test`` / ``_ref`` / ``_diff`` netCDF
    files so the viewer's data links resolve.
    """
    output_dir = _get_output_dir(parameter)
    json_path = os.path.join(output_dir, f"{parameter.output_file}.json")
    with open(json_path, "w") as outfile:
        json.dump({stat: metrics[stat] for stat in panels}, outfile)
    logger.info(f"Metrics saved in {json_path}")

    if not parameter.save_netcdf:
        return

    for data_type in ("test", "ref", "diff"):
        ds_out = xr.Dataset(
            {f"{var_key}_{stat}": panels[stat][data_type] for stat in panels}
        )
        _, filepath = _get_output_filename_filepath(parameter, data_type)
        ds_out.to_netcdf(filepath)
        logger.info(f"Interannual variability {data_type} output saved in: {filepath}")


def _compute_climo_and_interannual_std(
    ds: xr.Dataset, var_key: str, region: str
) -> dict[str, xr.Dataset]:
    """Compute the climatological mean and interannual std dev maps.

    The monthly time series is subset to ``region``. The climatological mean
    reuses the shared e3sm_diags :func:`climo` routine (day-weighted annual
    climatology). The interannual standard deviation aggregates the series into
    day-weighted annual means, then takes the standard deviation across years.
    This mirrors a-prime's ``compute_reg_seasonal_climo_and_stddev``, which
    ``climo`` alone cannot produce because it collapses the time axis.

    Parameters
    ----------
    ds : xr.Dataset
        The monthly time series dataset.
    var_key : str
        The variable key.
    region : str
        The region to subset on.

    Returns
    -------
    dict[str, xr.Dataset]
        A dictionary with keys ``"mean"`` and ``"std"`` mapping to the
        climatological mean and interannual standard deviation datasets.
    """
    ds_region = _subset_on_region(ds, var_key, region)

    # Climatological mean via the shared e3sm_diags climo routine (day-weighted).
    da_mean = climo(ds_region, var_key, "ANN")

    # Interannual std dev: day-weighted annual means, then std across years.
    ds_annual = ds_region.temporal.group_average(var_key, freq="year", weighted=True)
    da_annual = ds_annual[var_key]
    time_dim = xc.get_dim_keys(da_annual, axis="T")
    da_std = da_annual.std(dim=time_dim, keep_attrs=True)

    units = ds_region[var_key].attrs.get("units", "")

    stats: dict[str, xr.Dataset] = {}
    for stat, da_stat in [("mean", da_mean), ("std", da_std)]:
        da_stat.name = var_key
        da_stat.attrs["units"] = units

        ds_stat = da_stat.to_dataset()
        ds_stat = ds_stat.bounds.add_missing_bounds(axes=["X", "Y"])
        stats[stat] = ds_stat

    return stats


def _round_to_first(x: float) -> float:
    """Round ``x`` to its first significant digit (a-prime ``round_to_first``)."""
    if x == 0:
        return 0.0
    return round(x, -int(math.floor(math.log10(abs(x)))))


def _round_to_first_given_range(x: float, range_x: float) -> float:
    """Round ``x`` to the first significant digit of ``range_x``.

    Port of a-prime's ``round_to_first_given_range``.
    """
    if x == 0 or range_x == 0:
        return 0.0
    return round(x, -int(math.floor(math.log10(abs(range_x)))))


def _compute_variability_contour_levels(
    field: np.ndarray,
    n_stddev: int = INTERANNUAL_N_STDDEV,
    num_levels: int = INTERANNUAL_NUM_LEVELS,
) -> list[float]:
    """Compute contour levels spanning ``mean`` +/- ``n_stddev`` std of ``field``.

    Port of a-prime's ``compute_contour_levels``. Fields that straddle zero are
    given a symmetric range; otherwise the range is clipped to the data extent.
    """
    field = np.asarray(field)
    fmin = float(np.nanmin(field))
    fmax = float(np.nanmax(field))
    fmean = float(np.nanmean(field))
    fstd = float(np.nanstd(field))

    if fmin < 0 and fmax > 0:
        max_plot_temp = fmean + n_stddev * fstd
        range_plot = 2 * max_plot_temp
        max_plot = _round_to_first_given_range(max_plot_temp, range_plot)
        levels = np.linspace(-max_plot, max_plot, num=num_levels)
    else:
        max_plot_temp = min(fmean + n_stddev * fstd, fmax)
        min_plot_temp = max(fmean - n_stddev * fstd, fmin)
        range_plot = max_plot_temp - min_plot_temp
        max_plot = _round_to_first_given_range(max_plot_temp, range_plot)
        min_plot = _round_to_first_given_range(min_plot_temp, range_plot)
        levels = np.linspace(min_plot, max_plot, num=num_levels)

    return [float(x) for x in levels]


def _compute_diff_contour_levels(
    ref_field: np.ndarray, num_levels: int = INTERANNUAL_NUM_LEVELS
) -> list[float]:
    """Compute symmetric difference contour levels from ``ref_field``.

    Matches a-prime's plot_stddev difference panel, which uses a symmetric range
    of +/- 2 standard deviations of the reference field.
    """
    max_plot = _round_to_first(2.0 * float(np.nanstd(np.asarray(ref_field))))
    if max_plot == 0:
        max_plot = 1.0
    return [float(x) for x in np.linspace(-max_plot, max_plot, num=num_levels)]


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
    # NOTE: Ensure attributes from the variable are preserved rather than
    # overriden with the mask value of 1 using `keep_attrs=False`.
    conf_lvls_masked = xr.where(conf_lvls < 0.05, 1, 0, keep_attrs=False)

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

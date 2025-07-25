from __future__ import annotations

import collections
import os
from datetime import datetime, timedelta
from typing import TYPE_CHECKING, Any

import numpy as np
import xarray as xr
from netCDF4 import Dataset as netcdffile

import e3sm_diags
from e3sm_diags.plot import tc_analysis_plot

if TYPE_CHECKING:
    from numpy.ma.core import MaskedArray

    from e3sm_diags.parameter.tc_analysis_parameter import TCAnalysisParameter

from e3sm_diags.logger import _setup_child_logger

logger = _setup_child_logger(__name__)

# Years include 1979–2018 according to Balaguru et al. 2020
OBS_START_YR = 1979
OBS_END_YR = 2018
OBS_YEARS = np.arange(OBS_START_YR, OBS_END_YR + 1)

# (basin name, E bound, W bound, S bound, N bound, observed hurricane number per year)
BasinInfo = tuple[str, float, float, float, float, float]
BASIN_DICT: dict[str, BasinInfo] = {
    "NA": ("North Atlantic", 270, 360, 0, 45, 8.6),
    "WP": ("Northwest Pacific", 100, 180, 0, 45, 26.7),
    "EP": ("Eastern Pacific", 180, 270, 0, 45, 18.1),
    "NI": ("North Indian", 30, 100, 0, 45, 4.6),
    "SI": ("South Indian", 20, 135, -45, 0, 15.4),
    "SP": ("South Pacific", 135, 270, -45, 0, 10.4),
}


def run_diag(parameter: TCAnalysisParameter) -> TCAnalysisParameter:
    """Runs the tropical cyclone analysis diagnostic.

    :param parameter: Parameters for the run
    :type parameter: CoreParameter
    :raises ValueError: Invalid run type
    :return: Parameters for the run
    :rtype: CoreParameter
    """
    test_data_path = parameter.test_data_path
    reference_data_path = parameter.reference_data_path
    test_name = parameter.test_name
    run_type = parameter.run_type
    test_start_yr = parameter.test_start_yr
    test_end_yr = parameter.test_end_yr

    # Use the basin location info from [Balaguru et al. 2020]
    # obs for TC frequency of each basins

    # https://www.jstage.jst.go.jp/article/jmsj/84/2/84_2_259/_pdf/-char/ja

    # ns_obs = 26.7 + 10.4 + 18.1 + 8.6 + 4.6 + 15.4
    # pdf_wp_obs = 26.7 / ns_obs
    # pdf_sp_obs = 10.4 / ns_obs
    # pdf_ep_obs = 18.1 / ns_obs
    # pdf_at_obs = 8.6 / ns_obs
    # pdf_ni_obs = 4.6 / ns_obs
    # pdf_si_obs = 15.4 / ns_obs

    test_te_file = os.path.join(
        test_data_path,
        "cyclones_stitch_{}_{}_{}.dat".format(test_name, test_start_yr, test_end_yr),
    )
    test_cyclones_file = os.path.join(
        test_data_path,
        "cyclones_hist_{}_{}_{}.nc".format(test_name, test_start_yr, test_end_yr),
    )
    test_cyclones_hist = xr.open_dataset(test_cyclones_file).sel(lat=slice(-60, 60))[
        "density"
    ]
    test_aew_file = os.path.join(
        test_data_path,
        "aew_hist_{}_{}_{}.nc".format(test_name, test_start_yr, test_end_yr),
    )
    test_aew_hist = xr.open_dataset(test_aew_file).sel()["density"]

    test_data = collections.OrderedDict()
    ref_data = collections.OrderedDict()

    test_data["metrics"] = generate_tc_metrics_from_te_stitch_file(test_te_file)
    test_data["cyclone_density"] = test_cyclones_hist
    test_data["aew_density"] = test_aew_hist
    test_num_years = int(test_end_yr) - int(test_start_yr) + 1
    test_data["aew_num_years"] = test_num_years  # type: ignore
    test_data["cyclone_num_years"] = test_num_years  # type: ignore
    parameter.test_title = "{} ({}-{})".format(
        parameter.test_name, test_start_yr, test_end_yr
    )

    if run_type == "model_vs_model":
        ref_name = parameter.ref_name
        ref_start_yr = parameter.ref_start_yr
        ref_end_yr = parameter.ref_end_yr
        ref_te_file = os.path.join(
            reference_data_path,
            "cyclones_stitch_{}_{}_{}.dat".format(ref_name, ref_start_yr, ref_end_yr),
        )
        ref_cyclones_file = os.path.join(
            reference_data_path,
            "cyclones_hist_{}_{}_{}.nc".format(ref_name, ref_start_yr, ref_end_yr),
        )

        ref_cyclones_hist = xr.open_dataset(ref_cyclones_file).sel(lat=slice(-60, 60))[
            "density"
        ]

        ref_aew_file = os.path.join(
            reference_data_path,
            "aew_hist_{}_{}_{}.nc".format(ref_name, ref_start_yr, ref_end_yr),
        )
        # Note the refactor included subset that was missed in original implementation
        ref_aew_hist = xr.open_dataset(ref_aew_file).sel()["density"]
        ref_data["metrics"] = generate_tc_metrics_from_te_stitch_file(ref_te_file)
        ref_data["cyclone_density"] = ref_cyclones_hist
        ref_data["aew_density"] = ref_aew_hist
        ref_num_years = int(ref_end_yr) - int(ref_start_yr) + 1
        ref_data["aew_num_years"] = ref_num_years  # type: ignore
        ref_data["cyclone_num_years"] = ref_num_years  # type: ignore
        parameter.ref_title = "{} ({}-{})".format(
            parameter.ref_name, ref_start_yr, ref_end_yr
        )
    elif run_type == "model_vs_obs":
        ref_data["metrics"] = generate_tc_metrics_from_obs_files(reference_data_path)
        ref_cyclones_file = os.path.join(
            reference_data_path, "cyclones_hist_IBTrACS_1979_2018.nc"
        )
        ref_cyclones_hist = xr.open_dataset(ref_cyclones_file).sel(lat=slice(-60, 60))[
            "density"
        ]

        ref_aew_file = os.path.join(reference_data_path, "aew_hist_ERA5_2010_2014.nc")
        ref_aew_hist = xr.open_dataset(ref_aew_file).sel(
            lat=slice(0, 35), lon=slice(180, 360)
        )["density"]
        ref_data["cyclone_density"] = ref_cyclones_hist
        ref_data["cyclone_num_years"] = 40  # type: ignore
        ref_data["aew_density"] = ref_aew_hist
        # Question, should the num_years = 5?
        ref_data["aew_num_years"] = 1  # type: ignore
        parameter.ref_name = "Observation"
        parameter.ref_title = "Observation"
    else:
        raise ValueError("Invalid run_type={}".format(run_type))

    tc_analysis_plot.plot(test_data, ref_data, parameter, BASIN_DICT)

    return parameter


def generate_tc_metrics_from_te_stitch_file(te_stitch_file: str) -> dict[str, Any]:
    """Generates tropical cyclone metrics from TE stitch file.

    :param te_stitch_file: TE stitch file path
    :type te_stitch_file: str
    :return: Tropical cyclone metrics
    :rtype: dict[str, Any]

    # TODO: Add tests to cover this function
    """
    logger.info("\nGenerating TC Metrics from TE Stitch Files")
    logger.info("============================================")
    if not os.path.exists(te_stitch_file):
        raise FileNotFoundError(f"The file {te_stitch_file} does not exist.")

    with open(te_stitch_file) as f:
        lines_orig = f.readlines()

    if not lines_orig:
        raise ValueError(f"The file {te_stitch_file} is empty.")

    data_start_year = int(te_stitch_file.split(".")[-2].split("_")[-2])
    data_end_year = int(te_stitch_file.split(".")[-2].split("_")[-1])

    lines = _filter_lines_within_year_bounds(lines_orig, data_end_year)

    if not lines:
        raise ValueError(f"The file {te_stitch_file} is empty.")

    # Calculate number of storms and max length
    num_storms, max_len = _calc_num_storms_and_max_len(lines)
    # Parse variables from TE stitch file
    te_stitch_vars = _get_vars_from_te_stitch(lines, max_len, num_storms)
    # Add year info
    te_stitch_vars["year_start"] = data_start_year
    te_stitch_vars["year_end"] = data_end_year
    te_stitch_vars["num_years"] = data_end_year - data_start_year + 1
    logger.info(
        f"TE Start Year: {te_stitch_vars['year_start']}, TE End Year: {te_stitch_vars['year_end']}, Total Years: {te_stitch_vars['num_years']}"
    )

    # Use E3SM land-sea mask
    mask_path = os.path.join(e3sm_diags.INSTALL_PATH, "acme_ne30_ocean_land_mask.nc")
    ocnfrac = xr.open_dataset(mask_path)["OCNFRAC"].squeeze(dim="time", drop=True)

    # From model data, this dict stores a tuple for each basin.
    # (mean ace, tc_intensity_dist, seasonal_cycle, # storms, # of storms over the ocean)
    result_mod: dict[str, Any] = {}
    result_mod["num_years"] = te_stitch_vars["num_years"]

    for basin, basin_info in BASIN_DICT.items():
        mod_vars = _derive_metrics_per_basin(
            num_storms,
            te_stitch_vars,
            ocnfrac,
            basin_info,
        )

        pdf_mod_intensity = _calc_ts_intensity_dist(mod_vars["mod_wnd"])
        pdf_mod_seasonal_cycle = _calc_seasonal_cycle(mod_vars["mod_mon"])
        storms_overall_per_yr = mod_vars["mod_num"] / te_stitch_vars["num_years"]
        storms_over_ocean_per_yr = mod_vars["mod_num_ocn"] / te_stitch_vars["num_years"]

        result_mod[basin] = [
            mod_vars["mod_ace_mean"],
            pdf_mod_intensity,
            pdf_mod_seasonal_cycle,
            storms_overall_per_yr,
            storms_over_ocean_per_yr,
        ]

    return result_mod


def _filter_lines_within_year_bounds(
    lines_orig: list[str], data_end_year: int
) -> list[str]:
    """Filters lines within the specified year bounds.

    This function processes a list of strings, each representing a line of data.
    It filters out lines based on a year extracted from each line, ensuring that
    only lines with years less than or equal to `data_end_year` are retained.
    Additionally, it removes excessive time points crossing year bounds from
    6-hourly data.

    Parameters
    ----------
    lines_orig : list[str]
        A list of strings where each string represents a line of data.
    data_end_year : int
        The end year for filtering lines. Only lines with years less than or
        equal to this value will be retained.
    Returns
    -------
    list[str]
        A list of strings filtered based on the specified year bounds.
    """
    line_ind = []
    for i in range(0, np.size(lines_orig)):
        if lines_orig[i][0] == "s":
            year = int(lines_orig[i].split("\t")[2])

            if year <= data_end_year:
                line_ind.append(i)

    end_ind = line_ind[-1]

    new_lines = lines_orig[0:end_ind]
    return new_lines


def _calc_num_storms_and_max_len(lines: list[str]) -> tuple[int, int]:
    """Calculate number of storms and max length using lines from a TE stitch file.

    :param lines: Lines from TE stitch file
    :type lines: list[str]
    :return: Number of storms and max storm length
    :rtype: tuple[int, int]
    """
    num_storms = 0
    max_len = 0

    for line in lines:
        if line[0] == "s":
            num_storms = num_storms + 1
            num_points = 0
        else:
            num_points = num_points + 1
            max_len = max(max_len, num_points)

    logger.info(f"Number of storms: {num_storms}")
    logger.info(f"Max length of storms: {max_len}")
    return num_storms, max_len


def _get_vars_from_te_stitch(
    lines: list[str], max_len: int, num_storms: int
) -> dict[str, Any]:
    """Extracts variables from lines of a TE stitch file.

    :param lines: Lines from a TE stitch file
    :type lines: list[str]
    :param max_len: Max length of storms
    :type max_len: int
    :param num_storms: Number of storms
    :type num_storms: int
    :return: Dictionary of variables from TE stitch file
    :rtype: dict[str, Any]
    """
    keys = ("longmc", "latmc", "vsmc", "yearmc", "monthmc")
    vars_dict = {k: np.empty((max_len, num_storms)) * np.nan for k in keys}

    index = 0

    for line in lines:
        line_split = line.split("\t")
        if line[0] == "s":
            index = index + 1
            k = 0
        else:
            k = k + 1
            vars_dict["longmc"][k - 1, index - 1] = float(line_split[2])
            vars_dict["latmc"][k - 1, index - 1] = float(line_split[3])
            # Convert wind speed from units m/s to knot by multiplying 1.94
            vars_dict["vsmc"][k - 1, index - 1] = float(line_split[5]) * 1.94
            vars_dict["yearmc"][k - 1, index - 1] = float(line_split[6])
            vars_dict["monthmc"][k - 1, index - 1] = float(line_split[7])

    return vars_dict


def _derive_metrics_per_basin(
    num_storms: int,
    vars: dict[str, Any],
    ocnfrac: xr.DataArray,
    basin_info: BasinInfo,
) -> dict[str, Any]:
    """Derives metrics for each basin using TE stitch variables and other information.

    :param num_storms: Number of storms
    :type num_storms: int
    :param vars: TE stitch variables
    :type vars: dict[str, Any]
    :param ocnfrac: Ocnfrac xarray dataarray variable
    :type ocnfrac: xarray.DataArray
    :param basin_info: Basin information
    :type basin_info: BasinInfo
    :return: A dictionary containing mod variables
    :rtype: Dist[str, Any]

    # TODO: Add tests to cover this function
    # TODO: Refactor this function to avoid using mod vars and dict separately
    """
    mod_mon = []
    mod_wnd = []
    mod_num = 0
    mod_num_ocn = 0

    years = np.asarray(np.arange(vars["year_start"], vars["year_end"] + 1))
    mod_ace = np.zeros((np.size(years)))

    latmc = vars["latmc"]
    longmc = vars["longmc"]
    vsmc = vars["vsmc"]
    monthmc = vars["monthmc"]
    yearmc = vars["yearmc"]

    for k in range(0, num_storms):
        if (
            latmc[0, k] > basin_info[3]
            and latmc[0, k] < basin_info[4]
            and longmc[0, k] > basin_info[1]
            and longmc[0, k] < basin_info[2]
        ):
            mod_num = mod_num + 1

        lat = latmc[:, k][~np.isnan(latmc[:, k])]
        lon = longmc[:, k][~np.isnan(latmc[:, k])]
        wind = vsmc[:, k][~np.isnan(latmc[:, k])]
        mon = monthmc[:, k][~np.isnan(latmc[:, k])]
        yer = yearmc[:, k][~np.isnan(latmc[:, k])]

        # Get the nearest location on land-sea mask to the first point of a TC Track
        p = np.abs(ocnfrac.lat.values - lat[0])
        loc_y = int(np.argmin(p))
        p = np.abs(ocnfrac.lon.values - lon[0])
        loc_x = int(np.argmin(p))
        ocn_frac_0 = ocnfrac[loc_y, loc_x]

        if (
            ocn_frac_0 > 0
            and lon[0] > basin_info[1]
            and lon[0] < basin_info[2]
            and lat[0] > basin_info[3]
            and lat[0] < basin_info[4]
        ):
            mod_num_ocn = mod_num_ocn + 1
            mod_mon.append(mon[0])
            mod_wnd.append(np.max(wind))

            wind_new = wind[wind > 35]
            mod_ace[np.where(years - yer[0] == 0)] = (
                mod_ace[np.where(years - yer[0] == 0)] + np.sum(wind_new**2) / 1e4
            )

    mod_vars = {
        "mod_mon": mod_mon,
        "mod_wnd": mod_wnd,
        "mod_num": mod_num,
        "mod_num_ocn": mod_num_ocn,
        "mod_ace_mean": np.mean(mod_ace),
    }
    return mod_vars


def generate_tc_metrics_from_obs_files(reference_data_path: str) -> dict[str, Any]:
    """Generates tropical cyclone metrics from observation files.

    :param reference_data_path: Reference data path
    :type reference_data_path: str
    :return: TC Metrics
    :rtype: dict[str, Any]

    # TODO: Add tests to cover this function
    """
    logger.info("\nGenerating TC Metrics from Obs Files")
    logger.info("======================================")
    # Using IBTrACS data, store a tuple for each basin
    # (mean ace, tc_intensity_dist, seasonal_cycle, # observed hurricane per year,
    # # of storms over the ocean)
    result_obs: dict[str, Any] = {}
    result_obs["num_years"] = OBS_YEARS.size

    for basin, basin_info in BASIN_DICT.items():
        nc = netcdffile(
            os.path.join(reference_data_path, "IBTrACS.{}.v04r00.nc".format(basin))
        )

        # Extract and parse variables from .nc file
        vsmc = np.squeeze(nc["wmo_wind"][:, :])
        time = np.squeeze(nc["time"][:, :])
        monthmc, yearic = _get_monthmc_yearic(time)
        num_rows = time.shape[0]

        # Calculate mean ace
        mean_ace = _calc_mean_ace(vsmc, yearic, num_rows)

        # Calculate TS intensity distribution and seasonal cycle
        mon, wnd = _get_mon_wind(vsmc, monthmc, yearic, num_rows)
        tc_intensity_dist = _calc_ts_intensity_dist(wnd)
        pdf_obs_seasonal_cycle = _calc_seasonal_cycle(mon)

        # Number of hurricanes and storms over the ocean (placeholder of 1)
        num_hurricanes_per_yr = basin_info[5]
        num_storms_over_ocean_per_yr = 1

        # Final result
        result_obs[basin] = (
            mean_ace,
            tc_intensity_dist,
            pdf_obs_seasonal_cycle,
            num_hurricanes_per_yr,
            num_storms_over_ocean_per_yr,
        )

    return result_obs


def _get_monthmc_yearic(time: "MaskedArray") -> tuple[np.ndarray, np.ndarray]:
    """Extracts the monthmc and yearic by parsing the time variable.

    All masked (missing) data points are ignored.

    :param time: Time variable, units are days since 1858-11-17 00:00:00
    :type time: MaskedArray
    :return: Arrays for the months and years based on the day of a hurricane
    :rtype: tuple[np.ndarray, np.ndarray]
    """
    monthmc = np.zeros(time.shape)
    yearic = np.zeros(time.shape[0])

    # The units attribute of the time variable
    days_since = datetime(1858, 11, 17, 0, 0, 0)

    for i in range(0, time.shape[0]):
        for j in range(0, time.shape[1]):
            if time[i, j] is not np.ma.masked:
                day_hurr = days_since + timedelta(time[i, j])
                monthmc[i, j] = day_hurr.month
                yearic[i] = day_hurr.year

    return monthmc, yearic


def _get_mon_wind(
    vsmc: "MaskedArray",
    monthmc: np.ndarray,
    yearic: np.ndarray,
    num_rows: int,
) -> tuple[list[int], list[int]]:
    """Extracts the months and max wind speeds.

    :param vsmc: Maximum sustained wind speed from official WMO agency.
    :type vsmc: MaskedArray
    :param monthmc: Months
    :type monthmc: np.ndarray
    :param yearic: Years
    :type yearic: np.ndarray
    :param num_rows: Number of rows in the data matrix
    :type num_rows: int
    :return: Array of months and max wind speeds
    :rtype: tuple[list[int], list[int]]
    """
    mon = []
    wnd = []

    for i in range(0, num_rows):
        if OBS_START_YR <= yearic[i] <= OBS_END_YR:
            mon.append(monthmc[i, 0])
            wnd.append(np.max(vsmc[i, :]))

    return mon, wnd


def _calc_mean_ace(vsmc: "MaskedArray", yearic: np.ndarray, num_rows: int) -> float:
    """Calculates the mean ace using wind speeds.

    :param vsmc: Maximum sustained wind speed from official WMO agency.
    :type vsmc: MaskedArray
    :param yearic: Years
    :type yearic: np.ndarray
    :param num_rows: Number of rows in the data matrix
    :type num_rows: int
    :return: Mean ace
    :rtype: float
    """
    num_years = OBS_YEARS.size
    ace = np.zeros(num_years)

    for i in range(0, num_years):
        for j in range(0, num_rows):
            if yearic[j] == OBS_YEARS[i]:
                wind = vsmc[j, :]
                wind_ts = wind[wind >= 35]
                ace[i] = ace[i] + np.sum(wind_ts**2) / 1e4

    return np.mean(ace)  # type: ignore


def _calc_ts_intensity_dist(wind_speeds: list[int]) -> np.ndarray:
    """Calculate TC intensity distribution based on wind speed.

    :param wind_speeds: Wind speeds
    :type wind_speeds: list[int]
    :return: Number of storms in each hurricane category (tc intensity distribution)
    :rtype: np.ndarray
    """
    num_categories = 6
    dist = np.zeros(num_categories)

    for speed in wind_speeds:
        if speed <= 34:
            continue
        elif speed > 34 and speed <= 63:
            dist[0] = dist[0] + 1
        elif speed > 63 and speed <= 82:
            dist[1] = dist[1] + 1
        elif speed > 82 and speed <= 95:
            dist[2] = dist[2] + 1
        elif speed > 95 and speed <= 112:
            dist[3] = dist[3] + 1
        elif speed > 113 and speed <= 136:
            dist[4] = dist[4] + 1
        elif speed > 136:
            dist[5] = dist[5] + 1

    return dist


def _calc_seasonal_cycle(mon: list[int]) -> np.ndarray:
    """Calculates the seasonal cycle using a list of months.

    :param mon: List of months
    :type mon: list[int]
    :return: Seasonal cycle
    :rtype: np.ndarray
    """
    seasonal_cycle = np.zeros(12)
    for i in range(0, len(mon)):
        k = int(mon[i])
        seasonal_cycle[k - 1] = seasonal_cycle[k - 1] + 1

    return seasonal_cycle / np.sum(seasonal_cycle)

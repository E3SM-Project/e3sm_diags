from __future__ import print_function

import collections
import os
from datetime import datetime, timedelta

import cdms2
import numpy as np
from netCDF4 import Dataset as netcdffile

import acme_diags
from acme_diags.plot.cartopy import tc_analysis_plot


def tc_intensity(wind):
    # Calculate TC intensity distribution based on wind speed
    tc_intensity_distr = np.zeros((6))

    for i in range(0, len(wind)):

        if wind[i] > 34 and wind[i] <= 63:

            tc_intensity_distr[0] = tc_intensity_distr[0] + 1

        elif wind[i] > 63 and wind[i] <= 82:

            tc_intensity_distr[1] = tc_intensity_distr[1] + 1

        elif wind[i] > 82 and wind[i] <= 95:

            tc_intensity_distr[2] = tc_intensity_distr[2] + 1

        elif wind[i] > 95 and wind[i] <= 112:

            tc_intensity_distr[3] = tc_intensity_distr[3] + 1

        elif wind[i] > 113 and wind[i] <= 136:

            tc_intensity_distr[4] = tc_intensity_distr[4] + 1

        elif wind[i] > 136:

            tc_intensity_distr[5] = tc_intensity_distr[5] + 1

    # tc_intensity_distr = tc_intensity_distr / np.sum(tc_intensity_distr)
    return tc_intensity_distr


def generate_tc_metrics_from_obs_files(reference_data_path, basin_dict):
    # result_obs saves [mean ace, tc_intensity, seasonal_cycle] from IBTrACS data for each basins
    result_obs = {}
    for basin, basin_info in basin_dict.items():
        nc = netcdffile(
            os.path.join(reference_data_path, "IBTrACS.{}.v04r00.nc".format(basin))
        )
        latmc = np.squeeze(nc["lat"][:, :])
        longmc = np.squeeze(nc["lon"][:, :])
        vsmc = np.squeeze(nc["wmo_wind"][:, :])
        time = np.squeeze(nc["time"][:, :])
        monthmc = np.zeros((np.shape(latmc)[0], np.shape(latmc)[1]))
        daymc = np.zeros((np.shape(latmc)[0], np.shape(latmc)[1]))
        hourmc = np.zeros((np.shape(latmc)[0], np.shape(latmc)[1]))
        yearic = np.zeros((np.shape(latmc)[0]))

        for i in range(0, np.shape(time)[0]):
            for j in range(0, np.shape(time)[1]):
                if time[i, j] is not np.ma.masked:
                    day_hurr = datetime(1858, 11, 17, 0, 0, 0) + timedelta(time[i, j])
                    yearic[i] = day_hurr.year
                    monthmc[i, j] = day_hurr.month
                    daymc[i, j] = day_hurr.day
                    hourmc[i, j] = day_hurr.hour

        nummc = np.shape(latmc)[0]
        lat = []
        lon = []
        mon = []
        wnd = []
        num = 0

        # Years include 1979–2018 according to Balaguru et al. 2020
        obs_start = 1979
        obs_end = 2018
        years = np.arange(obs_start, obs_end + 1)
        for i in range(0, nummc):
            # if yearic[i] > 1970 and yearic[i] < 2016:
            if yearic[i] >= obs_start and yearic[i] <= obs_end:

                lat.append(latmc[i, 0])
                lon.append(longmc[i, 0])
                mon.append(monthmc[i, 0])
                wnd.append(np.max(vsmc[i, :]))
                num = num + 1

        # years = np.arange(1971, 2016)
        ace = np.zeros((np.size(years)))

        for i in range(0, np.size(years)):

            for j in range(0, nummc):

                if yearic[j] == years[i]:

                    wind = vsmc[j, :]
                    wind_ts = wind[wind >= 35]
                    ace[i] = ace[i] + np.sum(wind_ts ** 2) / 1e4

        print(basin_info[0])

        ##############################################
        # Calculate number of storms in each hurricane category
        tc_intensity_distr = tc_intensity(wnd)
        # Calculate seasonal cycle
        pdf_obs_seasonal_cycle = np.zeros((12))

        for i in range(0, len(mon)):
            k = int(mon[i])
            pdf_obs_seasonal_cycle[k - 1] = pdf_obs_seasonal_cycle[k - 1] + 1

        pdf_obs_seasonal_cycle = pdf_obs_seasonal_cycle / np.sum(pdf_obs_seasonal_cycle)
        result_obs[basin] = [
            np.mean(ace),
            tc_intensity_distr,
            pdf_obs_seasonal_cycle,
            basin_info[5],
            1,
        ]

    result_obs["num_years"] = 40

    return result_obs


def generate_tc_metrics_from_te_stitch_files(te_stitch_file, basin_dict):

    # with open(os.path.join(test_data_path, "cyclones_all_stitch.txt")) as f:
    with open(te_stitch_file) as f:
        lines = f.readlines()

    num_storms = 0
    max_len = 0
    for i in range(0, np.size(lines)):
        if lines[i][0] == "s":
            num_storms = num_storms + 1
            num_points = 0
        else:
            num_points = num_points + 1
            max_len = max(max_len, num_points)
    print("number of storms", num_storms)
    print("max length of storms", max_len)

    latmc = np.empty((max_len, num_storms)) * np.nan
    longmc = np.empty((max_len, num_storms)) * np.nan
    vsmc = np.empty((max_len, num_storms)) * np.nan
    pcmc = np.empty((max_len, num_storms)) * np.nan
    monthmc = np.empty((max_len, num_storms)) * np.nan
    yearmc = np.empty((max_len, num_storms)) * np.nan

    index = 0
    year_start = int(lines[0].split("\t")[2])
    year_end = year_start

    for i in range(0, np.size(lines)):

        if lines[i][0] == "s":
            index = index + 1
            year = int(lines[i].split("\t")[2])
            year_end = max(year, year_start)
            k = 0
        else:
            k = k + 1
            longmc[k - 1, index - 1] = float(lines[i].split("\t")[2])
            latmc[k - 1, index - 1] = float(lines[i].split("\t")[3])
            pcmc[k - 1, index - 1] = float(lines[i].split("\t")[4])
            # Convert wind speed from units m/s to knot by multiplying 1.94
            vsmc[k - 1, index - 1] = float(lines[i].split("\t")[5]) * 1.94
            # yearmc[k - 1, index - 1] = float(lines[i].split("\t")[7])
            # monthmc[k - 1, index - 1] = float(lines[i].split("\t")[8])
            yearmc[k - 1, index - 1] = float(lines[i].split("\t")[6])
            monthmc[k - 1, index - 1] = float(lines[i].split("\t")[7])

    print(year_start, year_end)
    num_year = year_end - year_start + 1
    print(num_year)

    # Use Navy bathymetry as land sea mask to takes in accounts only ocean grids fall in the basin lat-lon range. esp. for intensity analysis
    # nc = netcdffile("./bath_data.nc")
    # lat_bath = nc["Y"][:]
    # lon_bath = nc["X"][:]
    # bath = nc["bath"][:]
    # Use E3SM land-sea mask instead

    mask_path = os.path.join(acme_diags.INSTALL_PATH, "acme_ne30_ocean_land_mask.nc")
    ocnfrac = cdms2.open(mask_path)("OCNFRAC", squeeze=1)

    # result_mod dictionary saves [mean ace, tc_intensity, seasonal_cycle, number of storms, number of storms over ocean] from model data
    result_mod = {}
    for basin, basin_info in basin_dict.items():

        mod_mon = []
        mod_wnd = []
        mod_num = 0
        mod_num_ocn = 0
        years = np.asarray(np.arange(year_start, year_end + 1))
        mod_ace = np.zeros((np.size(years)))

        # for k in range(0, np.shape(longmc)[1]):
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
            p = np.abs(ocnfrac.getLatitude()[:] - lat[0])
            loc_y = int(np.argmin(p))
            p = np.abs(ocnfrac.getLongitude()[:] - lon[0])
            loc_x = int(np.argmin(p))
            ocn_frac_0 = ocnfrac[loc_y, loc_x]

            # Use Navy bathymetry alternatively
            # p = np.abs(lat_bath - lat[0])
            # loc_y = np.argmin(p)
            # p = np.abs(lon_bath - lon[0])
            # loc_x = np.argmin(p)
            #
            # depth = bath[loc_y,loc_x]

            # A more time consuming method to get the nearest point
            # ocn_frac_0 = ocnfrac(latitude=(lat[0], lat[0], 'cob'),longitude=(lon[0], lon[0], 'cob'),squeeze=1)

            if (
                # depth < 0
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
                    mod_ace[np.where(years - yer[0] == 0)] + np.sum(wind_new ** 2) / 1e4
                )

        # Calculate storm intensity
        pdf_mod_intensity = tc_intensity(mod_wnd)

        # Calculate seasonal cycle
        pdf_mod_seasonal_cycle = np.zeros((12))

        for i in range(0, len(mod_mon)):
            k = int(mod_mon[i])
            pdf_mod_seasonal_cycle[k - 1] = pdf_mod_seasonal_cycle[k - 1] + 1

        pdf_mod_seasonal_cycle = pdf_mod_seasonal_cycle / np.sum(pdf_mod_seasonal_cycle)

        result_mod[basin] = [
            np.mean(mod_ace),
            pdf_mod_intensity,
            pdf_mod_seasonal_cycle,
            mod_num / num_year,
            mod_num_ocn,
        ]
    result_mod["num_years"] = num_year

    return result_mod


def run_diag(parameter):
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
    # Each key in basin_dict associated a list including [basin name, E bound, W bound, S bound, N bound, observed hurricane number per year]
    basin_dict = {
        "NA": ["North Atlantic", 270, 360, 0, 45, 8.6],
        "WP": ["Northwest Pacific", 100, 180, 0, 45, 26.7],
        "EP": ["Eastern Pacific", 180, 270, 0, 45, 18.1],
        "NI": ["North Indian", 30, 100, 0, 45, 4.6],
        "SI": ["South Indian", 20, 135, -45, 0, 15.4],
        "SP": ["South Pacific", 135, 270, -45, 0, 10.4],
    }

    test_te_file = os.path.join(
        test_data_path,
        "cyclones_stitch_{}_{}_{}.dat".format(test_name, test_start_yr, test_end_yr),
    )
    test_cyclones_file = os.path.join(
        test_data_path,
        "cyclones_hist_{}_{}_{}.nc".format(test_name, test_start_yr, test_end_yr),
    )
    test_cyclones_hist = cdms2.open(test_cyclones_file)(
        "density", lat=(-60, 60, "ccb"), squeeze=1
    )
    test_aew_file = os.path.join(
        test_data_path,
        "aew_hist_{}_{}_{}.nc".format(test_name, test_start_yr, test_end_yr),
    )
    test_aew_hist = cdms2.open(test_aew_file)(
        "density", lat=(0, 35, "ccb"), lon=(-180, 0, "ccb"), squeeze=1
    )

    test_data = collections.OrderedDict()
    ref_data = collections.OrderedDict()

    test_data["metrics"] = generate_tc_metrics_from_te_stitch_files(
        test_te_file, basin_dict
    )
    test_data["cyclone_density"] = test_cyclones_hist
    test_data["aew_density"] = test_aew_hist
    test_num_years = int(test_end_yr) - int(test_start_yr) + 1
    test_data["aew_num_years"] = test_num_years
    test_data["cyclone_num_years"] = test_num_years
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
            "cyclones_hist_{}_{}_{}.nc".format(test_name, test_start_yr, test_end_yr),
        )
        ref_cyclones_hist = cdms2.open(ref_cyclones_file)("density", squeeze=1)
        ref_aew_file = os.path.join(
            reference_data_path,
            "aew_hist_{}_{}_{}.nc".format(test_name, test_start_yr, test_end_yr),
        )
        ref_aew_hist = cdms2.open(ref_aew_file)("density", squeeze=1)
        ref_data["metrics"] = generate_tc_metrics_from_te_stitch_files(
            ref_te_file, basin_dict
        )
        ref_data["cyclone_density"] = ref_cyclones_hist
        ref_data["aew_density"] = ref_aew_hist
        ref_num_years = int(ref_end_yr) - int(ref_start_yr) + 1
        ref_data["aew_num_years"] = ref_num_years
        ref_data["cyclone_num_years"] = ref_num_years
        parameter.ref_title = "{} ({}-{})".format(
            parameter.ref_name, ref_start_yr, ref_end_yr
        )
    elif run_type == "model_vs_obs":
        ref_data["metrics"] = generate_tc_metrics_from_obs_files(
            reference_data_path, basin_dict
        )
        ref_cyclones_file = os.path.join(
            reference_data_path, "cyclones_hist_IBTrACS_1979_2018.nc"
        )
        ref_cyclones_hist = cdms2.open(ref_cyclones_file)(
            "density", lat=(-60, 60, "ccb"), squeeze=1
        )
        ref_aew_file = os.path.join(reference_data_path, "aew_hist_ERA5_2010_2014.nc")
        ref_aew_hist = cdms2.open(ref_aew_file)(
            "density", lat=(0, 35, "ccb"), lon=(180, 360, "ccb"), squeeze=1
        )
        ref_data["cyclone_density"] = ref_cyclones_hist
        ref_data["cyclone_num_years"] = 40
        ref_data["aew_density"] = ref_aew_hist
        ref_data["aew_num_years"] = 1
        parameter.ref_name = "Observation"
        parameter.ref_title = "Observation"
    else:
        raise Exception("Invalid run_type={}".format(run_type))

    tc_analysis_plot.plot(test_data, ref_data, parameter, basin_dict)

    return parameter

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mar 17, 2021

@author: Karthik Balaguru, refactored by Jill Zhang
"""

import os
from datetime import datetime, timedelta

import cartopy.crs as ccrs
import cdms2
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset as netcdffile


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

    tc_intensity_distr = tc_intensity_distr / np.sum(tc_intensity_distr)
    return tc_intensity_distr


# Start driver
test_data_path = "/Users/zhang40/Documents/ACME/e3sm_tc_diags/original"
reference_data_path = "/Users/zhang40/Documents/ACME/e3sm_tc_diags/original"

# Use the basin location info from [Balaguru et al. 2020]
basin_dict = {
    "NA": ["North Atlantic", 270, 360, 0, 45],  # , 270, 360, 0, 90],
    "WP": ["Northwest Pacific", 100, 180, 0, 45],  # , 100, 180, 0, 90],
    "EP": ["Eastern Pacific", 180, 270, 0, 45],  # , 180, 270, 0, 90],
    "NI": ["North Indian", 30, 100, 0, 45],  # , 40, 100, 0, 90],
    "SI": ["South Indian", 20, 135, -45, 0],  # , 20, 140, -90, 0],
    "SP": ["South Pacific", 135, 270, -45, 0],  # , 140, 260, -90, 0],
}

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
    years = np.arange(1971, obs_end + 1)
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
    result_obs[basin] = [np.mean(ace), tc_intensity_distr, pdf_obs_seasonal_cycle]
###################################################
# Start to process .txt files generated from Tempest-extremes

with open(os.path.join(test_data_path, "cyclones_all_stitch.txt")) as f:
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
        line_split = lines[i].split("\t")
        longmc[k - 1, index - 1] = float(lines[i].split("\t")[2])
        latmc[k - 1, index - 1] = float(lines[i].split("\t")[3])
        pcmc[k - 1, index - 1] = float(lines[i].split("\t")[4])
        # Convert wind speed from units m/s to knot by multiplying 1.94
        vsmc[k - 1, index - 1] = float(lines[i].split("\t")[5]) * 1.94
        yearmc[k - 1, index - 1] = float(lines[i].split("\t")[7])
        monthmc[k - 1, index - 1] = float(lines[i].split("\t")[8])

print(year_start, year_end)
num_year = year_end - year_start + 1
print(num_year)


# Use Navy bathymetry as land sea mask to takes in accounts only ocean grids fall in the basin lat-lon range. esp. for intensity analysis
# nc = netcdffile("./bath_data.nc")
# lat_bath = nc["Y"][:]
# lon_bath = nc["X"][:]
# bath = nc["bath"][:]
# Use E3SM land-sea mask instead
ocnfrac = cdms2.open(os.path.join(test_data_path, "acme_ne30_ocean_land_mask.nc"))[
    "OCNFRAC"
](squeeze=1)

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
        mod_num,
        mod_num_ocn,
    ]


fig, axes = plt.subplots(2, 3, figsize=(12, 6), sharex=True, sharey=True)
fig.subplots_adjust(hspace=0.4, wspace=0.15)
axes = axes.ravel()

ace_obs = []
ace_mod = []
num_mod = []
for ifig, (basin, basin_info) in enumerate(basin_dict.items()):
    ace_obs.append(result_obs[basin][0])
    ace_mod.append(result_mod[basin][0])
    num_mod.append(result_mod[basin][3])
    ax = axes[ifig]
    ax.plot(np.arange(1, 7), result_obs[basin][1], "k", linewidth=2, label="Obs")
    ax.plot(np.arange(1, 7), result_mod[basin][1], "--k", linewidth=2, label="E3SM-HR")
    ax.legend()
    ax.set_title(basin_info[0])
    ax.set_ylabel("Fraction")
    plt.xticks([1, 2, 3, 4, 5, 6], ["Cat0", "Cat1", "Cat2", "Cat3", "Cat4", "Cat5"])
    ax.xaxis.set_tick_params(labelbottom=True)

plt.ylim(-0.05, 0.9)


# obs for TC frequency of each basins

# https://www.jstage.jst.go.jp/article/jmsj/84/2/84_2_259/_pdf/-char/ja

ns_obs = 26.7 + 10.4 + 18.1 + 8.6 + 4.6 + 15.4
pdf_wp_obs = 26.7 / ns_obs
pdf_sp_obs = 10.4 / ns_obs
pdf_ep_obs = 18.1 / ns_obs
pdf_at_obs = 8.6 / ns_obs
pdf_ni_obs = 4.6 / ns_obs
pdf_si_obs = 15.4 / ns_obs

fig0 = plt.figure(figsize=(12, 6))
ax = fig0.add_subplot(111)

N = 6
ind = np.arange(N)  # the x locations for the groups
width = 0.27

obs_vals = [pdf_at_obs, pdf_wp_obs, pdf_ep_obs, pdf_ni_obs, pdf_si_obs, pdf_sp_obs]
rects1 = ax.bar(ind + width / 2, obs_vals, width, color="darkgrey")
mod_vals = num_mod / np.sum(num_mod)
rects2 = ax.bar(ind - width / 2, mod_vals, width, color="black")
print(
    "total number of storms",
    num_storms,
    "total number based on 6 basins",
    np.sum(num_mod),
    num_mod,
)

ax.set_xticks(ind)
ax.set_xticklabels(
    (
        "North Atlantic",
        "Northwest Pacific",
        "Eastern Pacific",
        "North Indian",
        "South Indian",
        "South Pacific",
    )
)
ax.legend(
    (rects2[0], rects1[0]),
    (
        "Observations (~{})".format(int(ns_obs)),
        "E3SM-HR (~{})".format(int(np.sum(num_mod) / num_year)),
    ),
)
ax.set_ylabel("Fraction")
ax.set_title("Relative frequency of TCs for each ocean basins")


# plt.ylim(-.05,.9)
fig1 = plt.figure(figsize=(12, 6))
ax = fig1.add_subplot(111)

N = 6
ind = np.arange(N)  # the x locations for the groups
width = 0.27

obs_vals = ace_obs / np.sum(ace_obs)
print(obs_vals, ace_mod, "ace")
rects2 = ax.bar(ind - width / 2, obs_vals, width, color="black")
mod_vals = ace_mod / np.sum(ace_mod)
rects1 = ax.bar(ind + width / 2, mod_vals, width, color="darkgrey")

ax.set_xticks(ind)
ax.set_xticklabels(
    (
        "North Atlantic",
        "Northwest Pacific",
        "Eastern Pacific",
        "North Indian",
        "South Indian",
        "South Pacific",
    )
)
ax.legend((rects2[0], rects1[0]), ("Obs", "E3SM-HR"))
ax.set_ylabel("Fraction")
ax.set_title(
    "Distribution of accumulated cyclone energy (ACE) among various ocean basins"
)


fig2, axes = plt.subplots(2, 3, figsize=(12, 6), sharex=True, sharey=True)
fig2.subplots_adjust(hspace=0.4, wspace=0.15)
axes = axes.ravel()


for ifig, (basin, basin_info) in enumerate(basin_dict.items()):
    ax = axes[ifig]
    ax.plot(np.arange(1, 13), result_obs[basin][2], "k", linewidth=2, label="Obs")
    ax.plot(np.arange(1, 13), result_mod[basin][2], "--k", linewidth=2, label="E3SM-HR")
    ax.legend()
    ax.set_title(basin_info[0])
    ax.set_ylabel("Fraction")
    plt.xticks(
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
        ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"],
    )
    ax.xaxis.set_tick_params(labelbottom=True)


##########################################################
# Plotting tracks density
file_name = "cyclones_all_obs.nc"
fin = cdms2.open(os.path.join(test_data_path, file_name))
var = fin("density", squeeze=1) / 38

clevs = np.arange(0, 1.1, 0.1)

fig, ax = plt.subplots(
    1,
    1,
    figsize=(10, 5),
    subplot_kw={"projection": ccrs.PlateCarree(central_longitude=180)},
)
ax.coastlines()
p1 = ax.contourf(
    var.getLongitude(),
    var.getLatitude(),
    var,
    transform=ccrs.PlateCarree(),
    extend="both",
    levels=clevs,
    cmap="jet",
)

cbax = fig.add_axes([0.92, 0.1, 0.01, 0.7])
cbar = fig.colorbar(p1, cax=cbax)


file_name = "cyclones_all_hires.nc"
fin = cdms2.open(os.path.join(test_data_path, file_name))
var = fin("density", squeeze=1) / 50

fig, ax = plt.subplots(
    1,
    1,
    figsize=(10, 5),
    subplot_kw={"projection": ccrs.PlateCarree(central_longitude=180)},
)
ax.coastlines()
p1 = ax.contourf(
    var.getLongitude(),
    var.getLatitude(),
    var,
    transform=ccrs.PlateCarree(),
    levels=clevs,
    extend="both",
    cmap="jet",
)


cbax = fig.add_axes([0.92, 0.1, 0.01, 0.7])
cbar = fig.colorbar(p1, cax=cbax)

##########################################################
# African Easterly Waves (AEWs) counts
file_name = "aew_era5_stitch_hist_5e-6.nc"
fin = cdms2.open(os.path.join(reference_data_path, file_name))
var = fin("density", squeeze=1)

clevs = np.arange(0, 15.1, 1)
fig, ax = plt.subplots(
    1,
    1,
    figsize=(10, 5),
    subplot_kw={"projection": ccrs.PlateCarree(central_longitude=180)},
)
ax.coastlines()
p1 = ax.contourf(
    var.getLongitude(),
    var.getLatitude(),
    var,
    transform=ccrs.PlateCarree(),
    extend="both",
    levels=clevs,
    cmap="jet",
)

cbax = fig.add_axes([0.92, 0.1, 0.01, 0.7])
cbar = fig.colorbar(p1, cax=cbax)


file_name = "aew_e3sm_lowres_stitch_hist_5e-6.nc"
fin = cdms2.open(os.path.join(test_data_path, file_name))
var = fin("density", squeeze=1)

fig, ax = plt.subplots(
    1,
    1,
    figsize=(10, 5),
    subplot_kw={"projection": ccrs.PlateCarree(central_longitude=180)},
)
ax.coastlines()
p1 = ax.contourf(
    var.getLongitude(),
    var.getLatitude(),
    var,
    transform=ccrs.PlateCarree(),
    levels=clevs,
    extend="both",
    cmap="jet",
)


cbax = fig.add_axes([0.92, 0.1, 0.01, 0.7])
cbar = fig.colorbar(p1, cax=cbax)
plt.show()

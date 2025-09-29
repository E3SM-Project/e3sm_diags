import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib
import numpy as np
import xcdat as xc
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter

from e3sm_diags.logger import _setup_child_logger
from e3sm_diags.plot.utils import MAIN_TITLE_FONTSIZE, _save_main_plot

matplotlib.use("agg")
import matplotlib.pyplot as plt  # isort:skip  # noqa: E402

logger = _setup_child_logger(__name__)


# Position and sizes of subplot axes in page coordinates (0 to 1)
PANEL_CFG = [
    (0.1691, 0.55, 0.6465, 0.2758),
    (0.1691, 0.27, 0.6465, 0.2758),
]


# Each key gives a list with ax extent, x ticks , y ticks, title, clevs, reference and time resolution ratio (convert 3hrly to 6hrly data, density needs to be devided by 2)
# TODO flexible to apply to 3hrly model output when compare track density.
PLOT_INFO = {
    "aew": {
        "ax_extent": [182, 359, 0, 35],
        "x_ticks": [240, 300],
        "y_ticks": [0, 15, 30],
        "title": "African Easterly Wave Density",
        "clevs": np.arange(0, 16.1, 2),
        "reference": "EAR5 (2000-2014)",
        "time_resolution_ratio": 1,
    },
    "cyclone": {
        "ax_extent": [0, 359, -60, 60],
        "x_ticks": [0, 60, 120, 180, 240, 300, 359.99],
        "y_ticks": [-60, -30, 0, 30, 60],
        "title": "TC Tracks Density",
        "clevs": [0, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25],
        "reference": "IBTrACS (1979-2018)",
        "time_resolution_ratio": 2,
    },
}


def plot(test, ref, parameter, basin_dict):
    test_metrics = test["metrics"]
    ref_metrics = ref["metrics"]

    test_num_year = test_metrics["num_years"]
    ref_num_year = ref_metrics["num_years"]

    if parameter.short_test_name:
        test_name = parameter.short_test_name
    else:
        test_name = parameter.test_name
    ref_name = parameter.ref_name

    # TC intensity of each basins
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), sharex=True, sharey=True)
    fig.subplots_adjust(hspace=0.4, wspace=0.15)
    axes = axes.ravel()

    ace_ref = []
    ace_test = []
    num_ref = []
    num_test = []
    for ifig, (basin, basin_info) in enumerate(basin_dict.items()):
        ace_ref.append(ref_metrics[basin][0])
        ace_test.append(test_metrics[basin][0])
        num_ref.append(ref_metrics[basin][3])
        num_test.append(test_metrics[basin][3])
        ax = axes[ifig]
        ax.plot(
            np.arange(1, 7),
            test_metrics[basin][1] / test_num_year,
            "--k",
            linewidth=2,
            label="Test",
        )
        ax.plot(
            np.arange(1, 7),
            ref_metrics[basin][1] / ref_num_year,
            "k",
            linewidth=2,
            label="Ref",
        )
        ax.legend()
        ax.set_title(basin_info[0])
        ax.set_ylabel("TCs per year")
        plt.xticks([1, 2, 3, 4, 5, 6], ["Cat0", "Cat1", "Cat2", "Cat3", "Cat4", "Cat5"])
        ax.xaxis.set_tick_params(labelbottom=True)

    plt.suptitle(
        "Test: {}".format(test_name) + "\n" + "Ref: {}".format(ref_name),
        ha="left",
        x=0.1,
        y=0.99,
    )

    parameter.output_file = "tc-intensity"
    _save_main_plot(parameter)

    plt.close(fig)

    # TC frequency of each basins
    fig = plt.figure(figsize=(12, 7))
    ax = fig.add_subplot(111)

    N = 6
    ind = np.arange(N)  # the x locations for the groups
    width = 0.27

    ref_vals = num_ref / np.sum(num_ref)
    rects2 = ax.bar(ind - width / 2, ref_vals, width, color="black")
    test_vals = num_test / np.sum(num_test)
    rects1 = ax.bar(ind + width / 2, test_vals, width, color="darkgrey")
    logger.info("total number based on 6 basins")

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
            "{}: (~{})storms per year".format(ref_name, int(np.sum(num_ref))),
            "{}: (~{}) storms per year".format(test_name, int(np.sum(num_test))),
        ),
    )
    ax.set_ylabel("Fraction")
    ax.set_title("Relative frequency of TCs for each ocean basins")

    parameter.output_file = "tc-frequency"
    _save_main_plot(parameter)

    plt.close(fig)

    fig1 = plt.figure(figsize=(12, 6))
    ax = fig1.add_subplot(111)

    N = 6
    ind = np.arange(N)  # the x locations for the groups
    width = 0.27

    ref_vals = ace_ref / np.sum(ace_ref)
    rects2 = ax.bar(ind - width / 2, ref_vals, width, color="black")
    test_vals = ace_test / np.sum(ace_test)
    rects1 = ax.bar(ind + width / 2, test_vals, width, color="darkgrey")

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

    ax.legend((rects2[0], rects1[0]), (ref_name, test_name))
    ax.set_ylabel("Fraction")
    ax.set_title(
        "Distribution of accumulated cyclone energy (ACE) among various ocean basins"
    )

    parameter.output_file = "ace-distribution"
    _save_main_plot(parameter)

    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(12, 6), sharex=True, sharey=True)
    fig.subplots_adjust(hspace=0.4, wspace=0.15)
    axes = axes.ravel()

    for ifig, (basin, basin_info) in enumerate(basin_dict.items()):
        ax = axes[ifig]
        ax.plot(np.arange(1, 13), ref_metrics[basin][2], "k", linewidth=2, label="Test")
        ax.plot(
            np.arange(1, 13),
            test_metrics[basin][2],
            "--k",
            linewidth=2,
            label="Ref",
        )
        ax.legend()
        ax.set_title(basin_info[0])
        ax.set_ylabel("Fraction")
        plt.xticks(
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
            ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"],
        )
        ax.xaxis.set_tick_params(labelbottom=True)

    plt.suptitle(
        "Test: {}".format(test_name) + "\n" + "Ref: {}".format(ref_name),
        ha="left",
        x=0.1,
        y=0.99,
    )

    parameter.output_file = "tc-frequency-annual-cycle"
    _save_main_plot(parameter)

    plt.close(fig)

    ##########################################################
    # Plot TC tracks density
    plot_map(test, ref, "aew", parameter)

    # Plot AEW density
    plot_map(test, ref, "cyclone", parameter)


def plot_map(test_data, ref_data, region, parameter):
    """Create figure, projection for maps"""

    test = test_data["{}_density".format(region)]
    test_num_years = test_data["{}_num_years".format(region)]

    ref = ref_data["{}_density".format(region)]
    ref_num_years = ref_data["{}_num_years".format(region)]

    fig = plt.figure(figsize=(8.5, 8.5), dpi=parameter.dpi)
    proj = ccrs.PlateCarree(central_longitude=180)

    # First panel
    plot_panel(
        0,
        fig,
        proj,
        test,
        test_num_years,
        region,
        parameter.test_title,
    )

    # Second panel
    plot_panel(
        1,
        fig,
        proj,
        ref,
        ref_num_years,
        region,
        parameter.ref_title,
    )

    # Figure title
    fig.suptitle(PLOT_INFO[region]["title"], x=0.5, y=0.9, fontsize=14)
    parameter.output_file = "{}-density-map".format(region)

    _save_main_plot(parameter)
    plt.close(fig)


def plot_panel(n, fig, proj, var, var_num_years, region, title):
    ax = fig.add_axes(PANEL_CFG[n], projection=proj)
    ax.set_extent(PLOT_INFO[region]["ax_extent"], ccrs.PlateCarree())

    clevs = PLOT_INFO[region]["clevs"]

    lat = xc.get_dim_coords(var, axis="Y")
    lon = xc.get_dim_coords(var, axis="X")

    var = var.squeeze()

    p1 = ax.contourf(
        lon,
        lat,
        var / var_num_years / PLOT_INFO[region]["time_resolution_ratio"],
        transform=ccrs.PlateCarree(),
        levels=clevs,
        extend="both",
        cmap="jet",
    )
    ax.coastlines(lw=0.3)
    ax.add_feature(cfeature.LAND, zorder=100, edgecolor="k")

    if title != "Observation":
        ax.set_title("{}".format(title), fontdict={"fontsize": MAIN_TITLE_FONTSIZE})
    else:
        ax.set_title(
            "{}".format(PLOT_INFO[region]["reference"]),
            fontdict={"fontsize": MAIN_TITLE_FONTSIZE},
        )
    ax.set_xticks(PLOT_INFO[region]["x_ticks"], crs=ccrs.PlateCarree())
    ax.set_yticks(PLOT_INFO[region]["y_ticks"], crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True, number_format=".0f")
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.tick_params(labelsize=8.0, direction="out", width=1)
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")

    cbax = fig.add_axes(
        (PANEL_CFG[n][0] + 0.6635, PANEL_CFG[n][1] + 0.0415, 0.0326, 0.1792)
    )
    cbar = fig.colorbar(p1, cax=cbax)
    cbar.set_ticks(clevs)
    cbar.ax.tick_params(labelsize=9.0, length=0)
    return

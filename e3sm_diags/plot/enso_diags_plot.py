from __future__ import annotations

import os
from typing import TYPE_CHECKING, List, Tuple

import cartopy.crs as ccrs
import matplotlib
import numpy as np
import xarray as xr
import xcdat as xc
from numpy.polynomial.polynomial import polyfit

from e3sm_diags.derivations.default_regions_xr import REGION_SPECS
from e3sm_diags.driver.utils.general import get_output_dir
from e3sm_diags.logger import custom_logger
from e3sm_diags.parameter.enso_diags_parameter import EnsoDiagsParameter
from e3sm_diags.plot import get_colormap
from e3sm_diags.plot.utils import (
    _configure_titles,
    _configure_x_and_y_axes,
    _get_c_levels_and_norm,
    _get_x_ticks,
    _get_y_ticks,
    _make_lon_cyclic,
    _save_plot,
)

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # isort:skip  # noqa: E402

if TYPE_CHECKING:
    from e3sm_diags.driver.enso_diags_driver import MetricsDict, MetricsSubDict

logger = custom_logger(__name__)

plotTitle = {"fontsize": 11.5}
plotSideTitle = {"fontsize": 9.5}

# Position and sizes of subplot axes in page coordinates (0 to 1)
PANEL_CFG = [
    (0.1691, 0.6810, 0.6465, 0.2258),
    (0.1691, 0.3961, 0.6465, 0.2258),
    (0.1691, 0.1112, 0.6465, 0.2258),
]

# Border padding relative to subplot axes for saving individual panels
# (left, bottom, right, top) in page coordinates
BORDER_PADDING = (-0.06, -0.03, 0.13, 0.03)

# Use 179.99 as central longitude due to https://github.com/SciTools/cartopy/issues/946
PROJECTION = ccrs.PlateCarree(central_longitude=179.99)


def plot_map(
    parameter: EnsoDiagsParameter,
    da_test: xr.DataArray,
    da_ref: xr.DataArray,
    da_diff: xr.DataArray,
    metrics_dict: MetricsDict,
    da_test_conf_lvls: xr.DataArray,
    da_ref_conf_lvls: xr.DataArray,
):
    # Create figure
    fig = plt.figure(figsize=parameter.figsize, dpi=parameter.dpi)
    fig.suptitle(parameter.main_title, x=0.5, y=0.97, fontsize=15)

    # Use non-regridded test and ref for stats, so we have the original stats
    # displayed
    _add_colormap(
        0,
        da_test,
        fig,
        parameter,
        parameter.test_colormap,
        parameter.contour_levels,
        (parameter.test_name_yrs, parameter.test_title, da_test.units),
        metrics_dict["test"],  # type: ignore
        conf=da_test_conf_lvls,
    )
    _add_colormap(
        1,
        da_ref,
        fig,
        parameter,
        parameter.reference_colormap,
        parameter.contour_levels,
        (parameter.ref_name_yrs, parameter.reference_title, da_ref.units),
        metrics_dict["ref"],  # type: ignore
        conf=da_ref_conf_lvls,
    )
    _add_colormap(
        2,
        da_diff,
        fig,
        parameter,
        parameter.diff_colormap,
        parameter.diff_levels,
        (None, parameter.diff_title, da_test.units),
        metrics_dict["diff"],  # type: ignore
    )

    _save_plot(fig, parameter)

    plt.close()


def plot_scatter(x, y, parameter):
    fig = plt.figure(figsize=parameter.figsize, dpi=parameter.dpi)
    test_color = "black"
    ref_color = "red"
    test_title = "Test" if parameter.test_title == "" else parameter.test_title
    if parameter.test_name_yrs:
        test_title += " : {}".format(parameter.test_name_yrs)
    ref_title = (
        "Reference" if parameter.reference_title == "" else parameter.reference_title
    )
    if parameter.ref_name_yrs:
        ref_title += " : {}".format(parameter.ref_name_yrs)
    # https://stackoverflow.com/questions/14827650/pyplot-scatter-plot-marker-size
    plt.scatter(
        x["test"],
        y["test"],
        label=test_title,
        color=test_color,
        marker="s",
        s=8,
    )
    plt.scatter(x["ref"], y["ref"], label=ref_title, color=ref_color, marker="o", s=8)
    for value_type in ["test", "ref"]:
        if value_type == "test":
            type_str = "Test"
            type_color = test_color
            x_range = (min(x["test"]), max(x["test"]))
        else:
            type_str = "Reference"
            type_color = ref_color
            x_range = (min(x["ref"]), max(x["ref"]))
        # https://stackoverflow.com/questions/35091879/merge-2-arrays-vertical-to-tuple-numpy
        # Two parallel lists of values
        values = np.array((x[value_type], y[value_type]))
        # Zip the values together (i.e., list of (x,y) instead of (list of x, list of y)
        values = values.T
        if y["var"] == "TAUX":
            value_strs = [""]
        else:
            value_strs = ["positive ", "negative "]
        for value_str in value_strs:
            # https://stackoverflow.com/questions/24219841/numpy-where-on-a-2d-matrix
            if value_str == "positive ":
                # For all rows (x,y), choose the rows where x > 0
                rows = np.where(values[:, 0] > 0)
                smaller_x_range = (0, x_range[1])
                linestyle = "-"
            elif value_str == "negative ":
                # For all rows (x,y), choose the rows where x < 0
                rows = np.where(values[:, 0] < 0)
                smaller_x_range = (x_range[0], 0)
                linestyle = "--"
            elif value_str == "":
                rows = None
                smaller_x_range = x_range
                linestyle = "-"
            if rows:
                # Get the filtered zip (list of (x,y) where x > 0 or x < 0)
                filtered_values = values[rows]
            else:
                filtered_values = values
            # Get the filtered parallel lists (i.e., (list of x, list of y))
            filtered_values = filtered_values.T
            # https://stackoverflow.com/questions/19068862/how-to-overplot-a-line-on-a-scatter-plot-in-python
            b, m = polyfit(filtered_values[0], filtered_values[1], 1)
            label = "Linear fit for %sTS anomalies: %s (slope = %.2f)" % (
                value_str,
                type_str,
                m,
            )
            ys = [b + m * x for x in smaller_x_range]
            plt.plot(
                smaller_x_range,
                ys,
                label=label,
                color=type_color,
                linestyle=linestyle,
            )
    max_test = max(abs(min(y["test"])), abs(max(y["test"])))
    max_ref = max(abs(min(y["ref"])), abs(max(y["ref"])))
    max_value = max(max_test, max_ref) + 1
    plt.ylim(-max_value, max_value)
    plt.xlabel("{} anomaly ({})".format(x["var"], x["units"]))
    plt.ylabel("{} anomaly ({})".format(y["var"], y["units"]))
    plt.legend()

    # Prepare to save figure
    # get_output_dir => {parameter.results_dir}/{set_name}/{parameter.case_id}
    # => {parameter.results_dir}/enso_diags/{parameter.case_id}
    output_dir = get_output_dir(parameter.current_set, parameter)
    if parameter.print_statements:
        logger.info("Output dir: {}".format(output_dir))
    # get_output_dir => {parameter.orig_results_dir}/{set_name}/{parameter.case_id}
    # => {parameter.orig_results_dir}/enso_diags/{parameter.case_id}
    original_output_dir = get_output_dir(parameter.current_set, parameter)
    if parameter.print_statements:
        logger.info("Original output dir: {}".format(original_output_dir))
    # parameter.output_file is defined in e3sm_diags/driver/enso_diags_driver.py
    # {parameter.results_dir}/enso_diags/{parameter.case_id}/{parameter.output_file}
    file_path = os.path.join(output_dir, parameter.output_file)
    # {parameter.orig_results_dir}/enso_diags/{parameter.case_id}/{parameter.output_file}
    original_file_path = os.path.join(original_output_dir, parameter.output_file)

    # Figure title
    fig.suptitle(parameter.main_title, x=0.5, y=0.93, fontsize=15)

    # Save figure
    for f in parameter.output_format:
        f = f.lower().split(".")[-1]
        plot_suffix = "." + f
        plot_file_path = file_path + plot_suffix
        plt.savefig(plot_file_path)
        # Get the filename that the user has passed in and display that.

        original_plot_file_path = original_file_path + plot_suffix
        logger.info(f"Plot saved in: {original_plot_file_path}")

    plt.close()


def _add_colormap(
    subplot_num: int,
    var: xr.DataArray,
    fig: plt.Figure,
    parameter: EnsoDiagsParameter,
    color_map: str,
    contour_levels: List[float],
    title: Tuple[str | None, str, str],
    stats: MetricsSubDict,
    conf: xr.DataArray | None = None,
):
    var = _make_lon_cyclic(var)
    lat = xc.get_dim_coords(var, axis="Y")
    lon = xc.get_dim_coords(var, axis="X")

    var = var.squeeze()

    # Configure contour levels and boundary norm.
    # --------------------------------------------------------------------------
    c_levels, norm = _get_c_levels_and_norm(contour_levels)

    ax = fig.add_axes(PANEL_CFG[subplot_num], projection=PROJECTION)

    # Get region info and X and Y plot ticks.
    # --------------------------------------------------------------------------
    region_key = parameter.regions[0]
    region_specs = REGION_SPECS[region_key]

    # Get the region's domain slices for latitude and longitude if set, or
    # use the default value. If both are not set, then the region type is
    # considered "global".
    lat_slice = region_specs.get("lat", (-90, 90))  # type: ignore
    lon_slice = region_specs.get("lon", (0, 360))  # type: ignore

    # Boolean flags for configuring plots.
    is_global_domain = lat_slice == (-90, 90) and lon_slice == (0, 360)
    is_lon_full = lon_slice == (0, 360)

    # Determine X and Y ticks using longitude and latitude domains respectively.
    lon_west, lon_east = lon_slice
    x_ticks = _get_x_ticks(lon_west, lon_east, is_global_domain, is_lon_full)

    lat_south, lat_north = lat_slice
    y_ticks = _get_y_ticks(lat_south, lat_north)

    # Get the figure Axes object using the projection above.
    # --------------------------------------------------------------------------
    ax.set_extent([lon_west, lon_east, lat_south, lat_north], crs=PROJECTION)
    color_map = get_colormap(color_map, parameter)
    contours = ax.contourf(
        lon,
        lat,
        var,
        transform=ccrs.PlateCarree(),
        norm=norm,
        levels=c_levels,
        cmap=color_map,
        extend="both",
    )

    if conf is not None:
        conf = _make_lon_cyclic(conf)
        conf = conf.squeeze()  # type: ignore

        # Values in conf will be either 0 or 1. Thus, there are only two levels -
        # represented by the no-hatching and hatching levels.
        ax.contourf(
            lon,
            lat,
            conf,
            2,
            transform=ccrs.PlateCarree(),
            norm=norm,
            colors="none",
            extend="both",
            hatches=[None, "//"],
        )

    # Full world would be aspect 360/(2*180) = 1
    ax.set_aspect((lon_east - lon_west) / (2 * (lat_north - lat_south)))
    ax.coastlines(lw=0.3)

    # Configure the titles, x and y axes
    # --------------------------------------------------------------------------
    _configure_titles(ax, title)
    _configure_x_and_y_axes(
        ax, x_ticks, y_ticks, ccrs.PlateCarree(), parameter.current_set
    )

    # Place a vertical line in the middle of the plot - i.e. 180 degrees
    ax.axvline(x=0.5, color="k", linewidth=0.5)

    # Add the colorbar.
    # --------------------------------------------------------------------------
    cbax = fig.add_axes(
        (
            PANEL_CFG[subplot_num][0] + 0.6635,
            PANEL_CFG[subplot_num][1] + 0.0115,
            0.0326,
            0.1792,
        )
    )
    cbar = fig.colorbar(contours, cax=cbax)

    if c_levels is None:
        cbar.ax.tick_params(labelsize=9.0, length=0)
    else:
        maxval = np.amax(np.absolute(c_levels[1:-1]))
        if maxval < 1.0:
            fmt = "%5.3f"
            pad = 30
        elif maxval < 10.0:
            fmt = "%5.2f"
            pad = 25
        elif maxval < 100.0:
            fmt = "%5.1f"
            pad = 25
        else:
            fmt = "%6.1f"
            pad = 30
        cbar.set_ticks(c_levels[1:-1])
        labels = [fmt % level for level in c_levels[1:-1]]
        cbar.ax.set_yticklabels(labels, ha="right")
        cbar.ax.tick_params(labelsize=9.0, pad=pad, length=0)

    # Add metrics text to the figure.
    # --------------------------------------------------------------------------
    top_stats = (stats["max"], stats["min"], stats["mean"], stats["std"])
    top_text = "Max\nMin\nMean\nSTD"
    fig.text(
        PANEL_CFG[subplot_num][0] + 0.6635,
        PANEL_CFG[subplot_num][1] + 0.2107,
        top_text,
        ha="left",
        fontdict=plotSideTitle,
    )
    fig.text(
        PANEL_CFG[subplot_num][0] + 0.7635,
        PANEL_CFG[subplot_num][1] + 0.2107,
        "%.2f\n%.2f\n%.2f\n%.2f" % top_stats,  # type: ignore
        ha="right",
        fontdict=plotSideTitle,
    )

    # Hatch text
    if conf is not None:
        hatch_text = "Hatched when pvalue < 0.05"
        fig.text(
            PANEL_CFG[subplot_num][0] + 0.25,
            PANEL_CFG[subplot_num][1] - 0.0355,
            hatch_text,
            ha="right",
            fontdict=plotSideTitle,
        )


def add_cyclic(var):
    lon = var.getLongitude()
    return var(longitude=(lon[0], lon[0] + 360.0, "coe"))


def get_ax_size(fig, ax):
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height
    width *= fig.dpi
    height *= fig.dpi
    return width, height


def determine_tick_step(degrees_covered):
    if degrees_covered > 180:
        return 60
    if degrees_covered > 60:
        return 30
    elif degrees_covered > 20:
        return 10
    else:
        return 1

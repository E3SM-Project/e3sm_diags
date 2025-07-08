from typing import TYPE_CHECKING

import cartopy.crs as ccrs
import matplotlib
import numpy as np
import xarray as xr
import xcdat as xc
from numpy.polynomial.polynomial import polyfit

from e3sm_diags.derivations.default_regions_xr import REGION_SPECS
from e3sm_diags.logger import _setup_child_logger
from e3sm_diags.parameter.enso_diags_parameter import EnsoDiagsParameter
from e3sm_diags.plot.utils import (
    DEFAULT_PANEL_CFG,
    SECONDARY_TITLE_FONTSIZE,
    _add_colorbar,
    _add_contour_plot,
    _configure_titles,
    _configure_x_and_y_axes,
    _get_c_levels_and_norm,
    _get_x_ticks,
    _get_y_ticks,
    _make_lon_cyclic,
    _save_main_plot,
    _save_plot,
    _save_single_subplot,
)

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # isort:skip  # noqa: E402

if TYPE_CHECKING:
    from e3sm_diags.driver.enso_diags_driver import (
        MetricsDictMap,
        MetricsDictScatter,
        MetricsSubDict,
    )

logger = _setup_child_logger(__name__)

# Use 179.99 as central longitude due to https://github.com/SciTools/cartopy/issues/946
PROJECTION = ccrs.PlateCarree(central_longitude=179.99)

# Border padding relative to subplot axes for saving individual panels
# (left, bottom, right, top) in page coordinates
ENSO_BORDER_PADDING_MAP = (-0.07, -0.025, 0.17, 0.022)


def plot_scatter(
    parameter: EnsoDiagsParameter, x: MetricsDictScatter, y: MetricsDictScatter
):
    fig = plt.figure(figsize=parameter.figsize, dpi=parameter.dpi)
    fig.suptitle(parameter.main_title, x=0.5, y=0.93, fontsize=15)

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

        values = np.array((x[value_type], y[value_type]))  # type: ignore
        values = values.T

        if y["var"] == "TAUX":
            value_strs = [""]
        else:
            value_strs = ["positive ", "negative "]

        for value_str in value_strs:
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
                filtered_values = values[rows]
            else:
                filtered_values = values

            filtered_values = filtered_values.T

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

    _save_plot_scatter(fig, parameter)
    plt.close(fig)


def plot_map(
    parameter: EnsoDiagsParameter,
    da_test: xr.DataArray,
    da_ref: xr.DataArray,
    da_diff: xr.DataArray,
    metrics_dict: MetricsDictMap,
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
    _plot_diff_rmse_and_corr(fig, metrics_dict["diff"])  # type: ignore

    _save_plot(fig, parameter, DEFAULT_PANEL_CFG, ENSO_BORDER_PADDING_MAP)

    plt.close(fig)


def _add_colormap(
    subplot_num: int,
    var: xr.DataArray,
    fig: plt.Figure,
    parameter: EnsoDiagsParameter,
    color_map: str,
    contour_levels: list[float],
    title: tuple[str | None, str, str],
    metrics: MetricsSubDict,
    conf: xr.DataArray | None = None,
):
    var = _make_lon_cyclic(var)
    lat = xc.get_dim_coords(var, axis="Y")
    lon = xc.get_dim_coords(var, axis="X")

    var = var.squeeze()

    # Configure contour levels and boundary norm.
    # --------------------------------------------------------------------------
    c_levels, norm = _get_c_levels_and_norm(contour_levels)

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
    x_ticks = _get_x_ticks(
        lon_west,
        lon_east,
        is_global_domain,
        is_lon_full,
        tick_step_func=_determine_tick_step,
    )

    lat_south, lat_north = lat_slice
    y_ticks = _get_y_ticks(lat_south, lat_north, tick_step_func=_determine_tick_step)

    # Get the figure Axes object using the projection above.
    # --------------------------------------------------------------------------
    ax = fig.add_axes(DEFAULT_PANEL_CFG[subplot_num], projection=PROJECTION)
    ax.set_extent([lon_west, lon_east, lat_south, lat_north], crs=PROJECTION)
    contour_plot = _add_contour_plot(
        ax, var, lon, lat, color_map, ccrs.PlateCarree(), norm, c_levels
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

    # Place a vertical line in the middle of the plot - i.e. 180 degrees
    ax.axvline(x=0.5, color="k", linewidth=0.5)

    # Configure the titles, x and y axes, and colorbar.
    # --------------------------------------------------------------------------
    _configure_titles(ax, title)
    _configure_x_and_y_axes(
        ax, x_ticks, y_ticks, ccrs.PlateCarree(), parameter.current_set
    )
    _add_colorbar(
        fig,
        subplot_num,
        DEFAULT_PANEL_CFG,
        contour_plot,
        c_levels=c_levels,
        rect=None,
        c_label_fmt_and_pad_func=_get_contour_label_format_and_pad,
    )

    # Add metrics text to the figure.
    # --------------------------------------------------------------------------
    metrics_values = (metrics["max"], metrics["min"], metrics["mean"], metrics["std"])
    top_text = "Max\nMin\nMean\nSTD"
    fig.text(
        DEFAULT_PANEL_CFG[subplot_num][0] + 0.6635,
        DEFAULT_PANEL_CFG[subplot_num][1] + 0.2,
        top_text,
        ha="left",
        fontdict={"fontsize": 9},
    )
    fig.text(
        DEFAULT_PANEL_CFG[subplot_num][0] + 0.7635,
        DEFAULT_PANEL_CFG[subplot_num][1] + 0.2,
        "%.2f\n%.2f\n%.2f\n%.2f" % metrics_values,  # type: ignore
        ha="right",
        fontdict={"fontsize": 9},
    )

    # Hatch text
    if conf is not None:
        hatch_text = "Hatched when pvalue < 0.05"
        fig.text(
            DEFAULT_PANEL_CFG[subplot_num][0] + 0.25,
            DEFAULT_PANEL_CFG[subplot_num][1] - 0.0355,
            hatch_text,
            ha="right",
            fontdict={"fontsize": SECONDARY_TITLE_FONTSIZE},
        )


def _plot_diff_rmse_and_corr(fig: plt.Figure, metrics_dict: MetricsSubDict):
    bottom_stats = (metrics_dict["rmse"], metrics_dict["corr"])
    bottom_text = "RMSE\nCORR"

    fig.text(
        DEFAULT_PANEL_CFG[2][0] + 0.6635,
        DEFAULT_PANEL_CFG[2][1] - 0.0205,
        bottom_text,
        ha="left",
        fontdict={"fontsize": SECONDARY_TITLE_FONTSIZE},
    )
    fig.text(
        DEFAULT_PANEL_CFG[2][0] + 0.7635,
        DEFAULT_PANEL_CFG[2][1] - 0.0205,
        "%.2f\n%.2f" % bottom_stats,  # type: ignore
        ha="right",
        fontdict={"fontsize": SECONDARY_TITLE_FONTSIZE},
    )


def _determine_tick_step(degrees_covered: float) -> int:
    """Determine the number of tick steps based on the degrees covered by the axis.

    Parameters
    ----------
    degrees_covered : float
        The degrees covered by the axis.

    Returns
    -------
    int
        The number of tick steps.
    """
    if degrees_covered > 180:
        return 60
    if degrees_covered > 60:
        return 30
    elif degrees_covered > 20:
        return 10
    else:
        return 1


def _get_contour_label_format_and_pad(c_levels: list[float]) -> tuple[str, int]:
    """Get the label format and padding for each contour level.

    Parameters
    ----------
    c_levels : list[float]
        The contour levels.

    Returns
    -------
    tuple[str, int]
        A tuple for the label format and padding.
    """
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

    return fmt, pad


def _save_plot_scatter(fig: plt.Figure, parameter: EnsoDiagsParameter):
    """Save the scatter plot using the shared _save_single_subplot function."""
    _save_main_plot(parameter)

    # Save the single subplot using shared helper (panel_config=None for full figure)
    if parameter.output_format_subplot:
        _save_single_subplot(fig, parameter, 0, None, None)

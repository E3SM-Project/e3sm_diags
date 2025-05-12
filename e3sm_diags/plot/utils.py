from __future__ import annotations

import os
from typing import Callable, List, Literal, Tuple

import cartopy.crs as ccrs
import matplotlib
import matplotlib.contour as mcontour
import numpy as np
import xarray as xr
import xcdat as xc
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
from matplotlib.transforms import Bbox

from e3sm_diags import INSTALL_PATH
from e3sm_diags.driver.utils.io import _get_output_dir
from e3sm_diags.logger import _setup_child_logger
from e3sm_diags.parameter.core_parameter import CoreParameter

matplotlib.use("Agg")
from matplotlib import colors  # isort:skip  # noqa: E402
import matplotlib.pyplot as plt  # isort:skip  # noqa: E402

logger = _setup_child_logger(__name__)

# Plot title and side title configurations.
MAIN_TITLE_FONTSIZE = 11.5
SECONDARY_TITLE_FONTSIZE = 9.5

# Position and sizes of subplot axes in page coordinates (0 to 1)
PanelConfig = List[Tuple[float, float, float, float]]
DEFAULT_PANEL_CFG: PanelConfig = [
    (0.1691, 0.6810, 0.6465, 0.2258),
    (0.1691, 0.3961, 0.6465, 0.2258),
    (0.1691, 0.1112, 0.6465, 0.2258),
]

# Border padding relative to subplot axes for saving individual panels
# (left, bottom, right, top) in page coordinates
BorderPadding = Tuple[float, float, float, float]
DEFAULT_BORDER_PADDING: BorderPadding = (-0.06, -0.03, 0.13, 0.03)

# The type annotation for the rect arg used for creating the color bar axis.
Rect = Tuple[float, float, float, float]

# Sets that use the lat_lon formatter to configure the X and Y axes of the plot.
SETS_USING_LAT_LON_FORMATTER = [
    "lat_lon",
    "lat_lon_land",
    "lat_lon_river",
    "diurnal_cycle",
    "enso_diags",
    "meridional_mean_2d",
    "streamflow",
    "tc_analysis",
    "aerosol_aeronet",
]


def _save_plot(
    fig: plt.Figure,
    parameter: CoreParameter,
    panel_configs: PanelConfig = DEFAULT_PANEL_CFG,
    border_padding: BorderPadding = DEFAULT_BORDER_PADDING,
):
    """Save the plot using the figure object and parameter configs.

    This function creates the output filename to save the plot. It also
    saves each individual subplot if the reference name is an empty string ("").

    Parameters
    ----------
    fig : plt.Figure
        The plot figure.
    parameter : CoreParameter
        The CoreParameter with file configurations.
    panel_configs : PanelConfig
        A list of panel configs consisting of positions and sizes, with each
        element representing a panel. By default, set to ``DEFAULT_PANEL_CFG``.
    border_padding : Tuple[float, float, float, float]
        A tuple of border padding configs (left, bottom, right, top) for each
        panel relative to the subplot axes. By default, set to
        ``DEFAULT_BORDER_PADDING``.
    """
    for f in parameter.output_format:
        f = f.lower().split(".")[-1]
        fnm = os.path.join(
            _get_output_dir(parameter),
            parameter.output_file + "." + f,
        )
        plt.savefig(fnm)
        logger.info(f"Plot saved in: {fnm}")

    # Save individual subplots
    if parameter.ref_name == "":
        panels = [panel_configs[0]]
    else:
        panels = panel_configs

    for f in parameter.output_format_subplot:
        fnm = os.path.join(
            _get_output_dir(parameter),
            parameter.output_file,
        )
        page = fig.get_size_inches()

        for idx, panel in enumerate(panels):
            # Extent of subplot
            subpage = np.array(panel).reshape(2, 2)
            subpage[1, :] = subpage[0, :] + subpage[1, :]
            subpage = subpage + np.array(border_padding).reshape(2, 2)
            subpage = list(((subpage) * page).flatten())  # type: ignore
            extent = Bbox.from_extents(*subpage)

            # Save subplot
            fname = fnm + ".%i." % idx + f
            plt.savefig(fname, bbox_inches=extent)

            orig_fnm = os.path.join(
                _get_output_dir(parameter),
                parameter.output_file,
            )
            fname = orig_fnm + ".%i." % idx + f
            logger.info(f"Sub-plot saved in: {fname}")


def _add_grid_res_info(fig, subplot_num, region_key, lat, lon, panel_configs):
    if subplot_num == 2 and "RRM" in region_key:
        dlat = lat[2] - lat[1]
        dlon = lon[2] - lon[1]
        fig.text(
            panel_configs[subplot_num][0] + 0.4635,
            panel_configs[subplot_num][1] - 0.04,
            "Resolution: {:.2f}x{:.2f}".format(dlat, dlon),
            ha="left",
            fontdict={"fontsize": SECONDARY_TITLE_FONTSIZE},
        )


def _make_lon_cyclic(var: xr.DataArray):
    """Make the longitude axis cyclic by adding a new coordinate point with 360.

    This function appends a new longitude coordinate point by taking the last
    coordinate point and adding 360 to it. It is used for sets such as lat_lon.

    Parameters
    ----------
    var : xr.DataArray
        The variable.

    Returns
    -------
    xr.DataArray
        The variable with a 360 coordinate point.
    """
    coords = xc.get_dim_coords(var, axis="X")
    dim = coords.name

    new_pt = var.isel({f"{dim}": 0})
    new_pt = new_pt.assign_coords({f"{dim}": (new_pt[dim] + 360)})

    new_var = xr.concat([var, new_pt], dim=dim)

    return new_var


def _get_c_levels_and_norm(
    contour_levels: List[float],
) -> Tuple[List[float] | None, colors.BoundaryNorm | None]:
    """Get the contour levels and boundary norm.

    If custom contour_levels are used (> 0 elements), then adjust contour_levels with
    endpoints and add a boundary norm.

    Parameters
    ----------
    contour_levels : List[float]
        The contour levels.

    Returns
    -------
    Tuple[List[float] | None, colors.BoundaryNorm | None]
        A tuple of optional contour levels and boundary norm.
    """
    c_levels = None
    norm = None

    if len(contour_levels) > 0:
        c_levels = [-1.0e8] + contour_levels + [1.0e8]
        norm = colors.BoundaryNorm(boundaries=c_levels, ncolors=256)

    return c_levels, norm


def _add_contour_plot(
    ax: matplotlib.axes.Axes,
    var: xr.DataArray,
    x: xr.DataArray,
    y: xr.DataArray,
    color_map: str,
    projection: ccrs.PlateCarree | None,
    norm: colors.BoundaryNorm | None,
    c_levels: List[float] | None,
) -> mcontour.QuadContourSet:
    """Add the contour plot to the figure axes object.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The figure axes object.
    var : xr.DataArray
        The variable to plot.
    x : xr.DataArray
        The coordinates of the X axis for the plot.
    y : xr.DataArray
        The coordinates of the Y axis for the plot.
    color_map : str
        The color map file path.
    projection : ccrs.PlateCarree | None
        The optional cartopy projection.
    norm : colors.BoundaryNorm | None
        The optional norm boundaries.
    c_levels : List[float] | None
        The optional contour levels.

    Returns
    -------
    mcontour.QuadContourSet
        The contour plot object.
    """
    cmap = _get_colormap(color_map)

    c_plot = ax.contourf(
        x,
        y,
        var,
        cmap=cmap,
        transform=projection,
        norm=norm,
        levels=c_levels,
        extend="both",
    )

    return c_plot


def _get_colormap(color_map: str):
    """Get the colormap (string or mpl colormap object).

    This function retrieves a colormap which can be a predefined matplotlib
    colormap, a colormap defined in a local .rgb file, or a colormap installed
    in a predefined path.

    Parameters
    ----------
    color_map : str
        The name of the colormap or the path to a .rgb file.

    Returns
    -------
    str or matplotlib.colors.LinearSegmentedColormap
        The colormap as a string if it's a predefined colormap, or a
        LinearSegmentedColormap object if it's loaded from a .rgb file.

    Raises
    ------
    IOError
        If the .rgb file is not found in the current working directory or the
        installed colormaps directory.
    """
    color_map = str(color_map)  # unicode don't seem to work well with string.endswith()
    if not color_map.endswith(".rgb"):  # predefined vcs/mpl colormap
        return color_map

    installed_colormap = os.path.join(INSTALL_PATH, "colormaps", color_map)

    if os.path.exists(color_map):
        # colormap is an .rgb in the current directory
        pass
    elif not os.path.exists(color_map) and os.path.exists(installed_colormap):
        # use the colormap from /plot/colormaps
        color_map = installed_colormap
    elif not os.path.exists(color_map) and not os.path.exists(installed_colormap):
        pth = os.path.join(INSTALL_PATH, "colormaps")
        msg = "File {} isn't in the current working directory or installed in {}"
        raise IOError(msg.format(color_map, pth))

    rgb_arr = np.loadtxt(color_map)
    rgb_arr = rgb_arr / 255.0

    cmap = colors.LinearSegmentedColormap.from_list(name=color_map, colors=rgb_arr)
    return cmap


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
    elif degrees_covered > 60:
        return 30
    elif degrees_covered > 30:
        return 10
    elif degrees_covered > 20:
        return 5
    else:
        return 1


def _get_x_ticks(
    lon_west: float,
    lon_east: float,
    is_global_domain: bool,
    is_lon_full: bool,
    axis_orientation: Literal[180, 360] = 180,
    tick_step_func: Callable = _determine_tick_step,
) -> np.ndarray:
    """Get the X axis ticks based on the longitude domain slice.

    Parameters
    ----------
    lon_west : float
        The west point (e.g., 0).
    lon_east : float
        The east point (e.g., 360).
    is_global_domain : bool
        If the domain type is "global".
    is_lon_full : bool
        True if the longitude domain is (0, 360).
    axis_orientation : Literal[180, 360]
        The longitude axis orientation, by default 180.
    tick_step_func : Callable
        A function to determine the tick step, which might vary between sets,
        by default `_determine_tick_step`.

    Returns
    -------
    np.array
        An array of floats representing X axis ticks.
    """
    # NOTE: cartopy does not support region cross dateline yet so longitude
    # needs to be adjusted if > 180.
    # https://github.com/SciTools/cartopy/issues/821.
    # https://github.com/SciTools/cartopy/issues/276
    if lon_west > 180 and lon_east > 180:
        lon_west = lon_west - 360
        lon_east = lon_east - 360

    lon_covered = lon_east - lon_west
    lon_step = tick_step_func(lon_covered)

    x_ticks = np.arange(lon_west, lon_east, lon_step)

    if is_global_domain or is_lon_full:
        # Subtract 0.50 to get 0 W to show up on the right side of the plot.
        # If less than 0.50 is subtracted, then 0 W will overlap 0 E on the
        # left side of the plot.  If a number is added, then the value won't
        # show up at all.
        if axis_orientation == 360:
            x_ticks = np.array([0, 60, 120, 180, 240, 300, 359.99], dtype=float)
        else:
            x_ticks = np.append(x_ticks, lon_east - 0.50)
    else:
        x_ticks = np.append(x_ticks, lon_east)

    return x_ticks


def _get_y_ticks(
    lat_south: float,
    lat_north: float,
    tick_step_func: Callable = _determine_tick_step,
) -> np.ndarray:
    """Get Y axis ticks.

    Parameters
    ----------
    lat_south : float
        The south point (e.g., -180).
    lat_north : float
        The north point (e.g., 180).
    tick_step_func : Callable
        A function to determine the tick step, which might vary between sets,
        by default `_determine_tick_step`.

    Returns
    -------
    np.array
        An array of floats representing Y axis ticks
    """
    lat_covered = lat_north - lat_south

    lat_step = tick_step_func(lat_covered)
    y_ticks = np.arange(lat_south, lat_north, lat_step)
    y_ticks = np.append(y_ticks, lat_north)

    return y_ticks


def _configure_titles(
    ax: matplotlib.axes.Axes,
    title: Tuple[str | None, str | None, str | None],
    main_fontsize: float = MAIN_TITLE_FONTSIZE,
    secondary_fontsize: float = SECONDARY_TITLE_FONTSIZE,
):
    """Configure the axes titles.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The figure axes object.
    title : Tuple[str | None, str | None, str | None]
        A tuple of strings to form the title of the colormap, in the format
        (<optional> years, title, units).
    main_fontsize : float
        The main title font size, by default 11.5.
    secondary_fontsize : float
        The secondary title font sizes, by default 9.5.

    Returns
    -------
    matplotlib.axes.Axes
        The axes objects.
    """
    if title[0] is not None:
        ax.set_title(title[0], loc="left", fontdict={"fontsize": secondary_fontsize})
    if title[1] is not None:
        ax.set_title(title[1], fontdict={"fontsize": main_fontsize})
    if title[2] is not None:
        # NOTE: loc="right"  doesn't work for polar projection
        ax.set_title(title[2], loc="right", fontdict={"fontsize": secondary_fontsize})


def _configure_x_and_y_axes(
    ax: matplotlib.axes.Axes,
    x_ticks: np.ndarray,
    y_ticks: np.ndarray | None,
    projection: ccrs.PlateCarree | None,
    set_name: str,
):
    """Configure the X and Y axes.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The figure axes object.
    x_ticks : np.ndarray
        The array of X ticks.
    y_ticks : np.ndarray | None
        The optional array of Y ticks. Some set plotters pass None to configure
        the Y axis ticks using other ticks such as Z axis plevs instead.
    projection : ccrs.PlateCarree | None
        The optional cartopy projection to use for X and Y ticks.
    set_name : set_name
        The name of the current set which determines whether the latitude and
        longitude major formatters are used.
    """
    # For `ax.set_xticks` and `ax.set_yticks`, `crs` cannot be `None` and we
    # must split up arguments passed to these methods using a conditional
    # statement. Otherwise, this error is raised: `ValueError: Incorrect use of
    # keyword argument 'crs'. Keyword arguments other than 'minor' modify the
    # text labels and can only be used if 'labels' are passed as well.`
    if projection is not None:
        ax.set_xticks(x_ticks, crs=projection)

        if y_ticks is not None:
            ax.set_yticks(y_ticks, crs=projection)
    else:
        ax.set_xticks(x_ticks)

        if y_ticks is not None:
            ax.set_yticks(y_ticks)

    ax.tick_params(labelsize=8.0, direction="out", width=1)

    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")

    if set_name in SETS_USING_LAT_LON_FORMATTER:
        lon_formatter = LongitudeFormatter(
            zero_direction_label=True, number_format=".0f"
        )
        lat_formatter = LatitudeFormatter()
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.yaxis.set_major_formatter(lat_formatter)


def _get_contour_label_format_and_pad(c_levels: List[float]) -> Tuple[str, int]:
    """Get the label format and padding for each contour level.

    Parameters
    ----------
    c_levels : List[float]
        The contour levels.

    Returns
    -------
    Tuple[str, int]
        A tuple for the label format and padding.
    """
    maxval = np.amax(np.absolute(c_levels[1:-1]))

    if maxval < 0.01:
        fmt = "%.1e"
        pad = 35
    elif maxval < 0.2:
        fmt = "%5.3f"
        pad = 28
    elif maxval < 10.0:
        fmt = "%5.2f"
        pad = 25
    elif maxval < 100.0:
        fmt = "%5.1f"
        pad = 25
    elif maxval > 9999.0:
        fmt = "%.0f"
        pad = 40
    else:
        fmt = "%6.1f"
        pad = 30

    return fmt, pad


def _add_colorbar(
    fig: plt.Figure,
    subplot_num: int,
    panel_configs: PanelConfig,
    contour_plot: mcontour.QuadContourSet,
    c_levels: List[float] | None,
    rect: Rect | None = None,
    c_label_fmt_and_pad_func: Callable = _get_contour_label_format_and_pad,
):
    """Configure the colorbar on a colormap.

    Parameters
    ----------
    fig : plt.Figure
        The figure object.
    subplot_num : int
        The subplot number.
    panel_configs : PanelConfig
        A list of panel configs consisting of positions and sizes, with each
        element representing a panel.
    contour_plot : mcontour.QuadContourSet
        The contour plot object.
    c_levels : List[float] | None
        The optional contour levels used to configure the colorbar.
    rect : Rect
        An optional adjustment to the dimensions (left, bottom, width, height)
        of the new `~.axes.Axes`. All quantities are in fractions of figure
        width and height.
    c_label_fmt_and_pad_func : Callable
        An optional function for configuring the contour level label format
        and padding.
    """
    cbax_rect = _get_rect(subplot_num, panel_configs, rect)
    cbax = fig.add_axes(cbax_rect)
    cbar = fig.colorbar(contour_plot, cax=cbax)

    if c_levels is None:
        cbar.ax.tick_params(labelsize=9.0, length=0)
    else:
        cbar.set_ticks(c_levels[1:-1])

        label_format, pad = c_label_fmt_and_pad_func(c_levels)
        labels = [label_format % level for level in c_levels[1:-1]]

        cbar.ax.set_yticklabels(labels, ha="right")
        cbar.ax.tick_params(labelsize=9.0, pad=pad, length=0)


def _get_rect(
    subplot_num: int,
    panel_configs: PanelConfig,
    rect: Rect | None,
) -> Rect:
    """Get the rect arg for the color bar axis.

    Parameters
    ----------
    subplot_num : int
        The subplot number.
    panel_configs : PanelConfig
        A list of panel configs consisting of positions and sizes, with each
        element representing a panel.
    rect : Rect
        An optional adjustment to the dimensions (left, bottom, width, height)
        of the new `~.axes.Axes`. All quantities are in fractions of figure
        width and height.

    Returns
    -------
    Rect
        The rect arg for the color bar axis.
    """
    if rect is None:
        rect = (0.6635, 0.0215, 0.0326, 0.1792)

    return (
        panel_configs[subplot_num][0] + rect[0],
        panel_configs[subplot_num][1] + rect[1],
        rect[2],
        rect[3],
    )


def _add_min_mean_max_text(
    fig: plt.Figure,
    subplot_num: int,
    panel_configs: PanelConfig,
    metrics: Tuple[float, ...],
    set_name: str | None = None,
    fontsize: float = SECONDARY_TITLE_FONTSIZE,
    left_text_pos: Tuple[float, float] | None = None,
    right_text_pos: Tuple[float, float] | None = None,
):
    """Add min, mean, and max text to the figure.

    Parameters
    ----------
    fig : plt.Figure
        The figure object.
    subplot_num : int
        The subplot number.
    panel_configs : PanelConfig
        A list of panel configs consisting of positions and sizes, with each
        element representing a panel.
    metrics : Tuple[float, ...]
        The tuple of metrics, with the first three elements being max, mean,
        and min.
    set_name : str | None
        The optional set name used to determine float format, by default None.
    fontsize : float
        The text font size, by default 9.5.
    left_text_pos: Tuple[float, float] | None
        An optional adjustment to the x, y position of the left text.
    right_text_post: Tuple[float, float] | None
        An optional adjustment to the x, y position of the right text.
    """
    fontdict = {"fontsize": fontsize}

    if left_text_pos is None:
        left_text_pos = (0.6635, 0.2107)

    if right_text_pos is None:
        right_text_pos = (0.7635, 0.2107)

    fig.text(
        panel_configs[subplot_num][0] + left_text_pos[0],
        panel_configs[subplot_num][1] + left_text_pos[1],
        "Max\nMean\nMin",
        ha="left",
        fontdict=fontdict,
    )

    fmt_metrics = _get_float_format(metrics, set_name)

    fig.text(
        panel_configs[subplot_num][0] + right_text_pos[0],
        panel_configs[subplot_num][1] + right_text_pos[1],
        fmt_metrics % metrics[0:3],
        ha="right",
        fontdict=fontdict,
    )


def _get_float_format(metrics: Tuple[float, ...], set_name: str | None) -> str:
    """Get the float format for string text based on decimal places of metrics.

    Parameters
    ----------
    metrics : Tuple[float, ...]
        The tuple of metrics, with the first three elements being max, mean, and
        min.
    set_name : str | None
        The optional name of the set.

    Returns
    -------
    str
        The float format.
    """
    # FIXME: This conditional code was ported over from two plot functions and
    # can be implemented better.
    if set_name in ["zonal_mean_2d", "zonal_mean_2d_stratosphere"]:
        # if positive Max is smaller than 0.01, use scientific notation
        if metrics[0] < 0.01 and metrics[0] > 0:
            float_format = "%.e\n%.e\n%.e"
        else:
            float_format = "%.2f\n%.2f\n%.2f"
    else:
        fmt_m = []

        # Print in scientific notation if value is greater than 10^5
        for i in range(len(metrics[0:3])):
            fs = "1e" if metrics[i] > 100000.0 else "2f"
            fmt_m.append(fs)

        float_format = f"%.{fmt_m[0]}\n%.{fmt_m[1]}\n%.{fmt_m[2]}"

    return float_format


def _add_rmse_corr_text(
    fig: plt.Figure,
    subplot_num: int,
    panel_configs: PanelConfig,
    metrics: Tuple[float, ...],
    fontsize: float = SECONDARY_TITLE_FONTSIZE,
    left_text_pos: Tuple[float, float] | None = None,
    right_text_pos: Tuple[float, float] | None = None,
):
    """Add RMSE and CORR metrics text to the figure.

    Parameters
    ----------
    fig : plt.Figure
        The figure object.
    subplot_num : int
        The subplot number.
    panel_configs : PanelConfig
        A list of panel configs consisting of positions and sizes, with each
        element representing a panel.
    metrics : Tuple[float, ...]
        The tuple of metrics, with the last two elements being RMSE and CORR.
    fontsize : float
        The text font size, by default 9.5.
    left_text_pos: Tuple[float, float] | None
        An optional adjustment to the x, y position of the left text.
    right_text_pos: Tuple[float, float] | None
        An optional adjustment to the x, y position of the right text.
    """
    fontdict = {"fontsize": fontsize}

    if left_text_pos is None:
        left_text_pos = (0.6635, -0.0105)

    if right_text_pos is None:
        right_text_pos = (0.7635, -0.0105)

    fig.text(
        panel_configs[subplot_num][0] + left_text_pos[0],
        panel_configs[subplot_num][1] + left_text_pos[1],
        "RMSE\nCORR",
        ha="left",
        fontdict=fontdict,
    )
    fig.text(
        panel_configs[subplot_num][0] + right_text_pos[0],
        panel_configs[subplot_num][1] + right_text_pos[1],
        "%.2f\n%.2f" % metrics[3:5],
        ha="right",
        fontdict=fontdict,
    )

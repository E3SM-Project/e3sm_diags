from __future__ import annotations

import os
from typing import List, Tuple

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib
import numpy as np
import xarray as xr
import xcdat as xc
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
from matplotlib.transforms import Bbox

from e3sm_diags.derivations.default_regions_xr import REGION_SPECS
from e3sm_diags.driver.utils.general import get_output_dir
from e3sm_diags.logger import custom_logger
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.plot import get_colormap

matplotlib.use("Agg")
from matplotlib import colors  # isort:skip  # noqa: E402
import matplotlib.pyplot as plt  # isort:skip  # noqa: E402

logger = custom_logger(__name__)

# Plot title and side title configurations.
PLOT_TITLE = {"fontsize": 11.5}
PLOT_SIDE_TITLE = {"fontsize": 9.5}

# Position and sizes of subplot axes in page coordinates (0 to 1)
PANEL = [
    (0.1691, 0.6810, 0.6465, 0.2258),
    (0.1691, 0.3961, 0.6465, 0.2258),
    (0.1691, 0.1112, 0.6465, 0.2258),
]

# Border padding relative to subplot axes for saving individual panels
# (left, bottom, right, top) in page coordinates
BORDER_PADDING = (-0.06, -0.03, 0.13, 0.03)


def _save_plot(fig: plt.figure, parameter: CoreParameter):
    """Save the plot using the figure object and parameter configs.

    This function creates the output filename to save the plot. It also
    saves each individual subplot if the reference name is an empty string ("").

    Parameters
    ----------
    fig : plt.figure
        The plot figure.
    parameter : CoreParameter
        The CoreParameter with file configurations.
    """
    for f in parameter.output_format:
        f = f.lower().split(".")[-1]
        fnm = os.path.join(
            get_output_dir(parameter.current_set, parameter),
            parameter.output_file + "." + f,
        )
        plt.savefig(fnm)
        logger.info(f"Plot saved in: {fnm}")

    # Save individual subplots
    if parameter.ref_name == "":
        panels = [PANEL[0]]
    else:
        panels = PANEL

    for f in parameter.output_format_subplot:
        fnm = os.path.join(
            get_output_dir(parameter.current_set, parameter),
            parameter.output_file,
        )
        page = fig.get_size_inches()

        for idx, panel in enumerate(panels):
            # Extent of subplot
            subpage = np.array(panel).reshape(2, 2)
            subpage[1, :] = subpage[0, :] + subpage[1, :]
            subpage = subpage + np.array(BORDER_PADDING).reshape(2, 2)
            subpage = list(((subpage) * page).flatten())  # type: ignore
            extent = Bbox.from_extents(*subpage)

            # Save subplot
            fname = fnm + ".%i." % idx + f
            plt.savefig(fname, bbox_inches=extent)

            orig_fnm = os.path.join(
                get_output_dir(parameter.current_set, parameter),
                parameter.output_file,
            )
            fname = orig_fnm + ".%i." % idx + f
            logger.info(f"Sub-plot saved in: {fname}")


def _add_colormap(
    subplot_num: int,
    var: xr.DataArray,
    fig: plt.figure,
    parameter: CoreParameter,
    color_map: str,
    contour_levels: List[float],
    title: Tuple[str | None, str, str],
    metrics: Tuple[float, ...],
):
    """Adds a colormap containing the variable data and metrics to the figure.

    This function is used by:
      - `lat_lon_plot.py`
      - `aerosol_aeronet_plot.py` (TODO)

    Parameters
    ----------
    subplot_num : int
        The subplot number.
    var : xr.DataArray
        The variable to plot.
    fig : plt.figure
        The figure object to add the subplot to.
    parameter : CoreParameter
        The CoreParameter object containing plot configurations.
    color_map : str
        The colormap styling to use (e.g., "cet_rainbow.rgb").
    contour_levels : List[float]
        The map contour levels.
    title : Tuple[str | None, str, str]
        A tuple of strings to form the title of the colormap, in the format
        (<optional> years, title, units).
    metrics : Tuple[float, ...]
        A tuple of metrics for this subplot.
    """
    var = _make_lon_cyclic(var)
    lat = xc.get_dim_coords(var, axis="Y")
    lon = xc.get_dim_coords(var, axis="X")

    var = var.squeeze()

    # Configure contour levels
    # --------------------------------------------------------------------------
    c_levels = None
    norm = None

    if len(contour_levels) > 0:
        c_levels = [-1.0e8] + contour_levels + [1.0e8]
        norm = colors.BoundaryNorm(boundaries=c_levels, ncolors=256)

    # Configure plot tickets based on longitude and latitude.
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

    # Add the contour plot.
    # --------------------------------------------------------------------------
    projection = ccrs.PlateCarree()
    if is_global_domain or is_lon_full:
        projection = ccrs.PlateCarree(central_longitude=180)

    ax = fig.add_axes(PANEL[subplot_num], projection=projection)
    ax.set_extent([lon_west, lon_east, lat_south, lat_north], crs=projection)
    color_map = get_colormap(color_map, parameter)
    p1 = ax.contourf(
        lon,
        lat,
        var,
        transform=ccrs.PlateCarree(),
        norm=norm,
        levels=c_levels,
        cmap=color_map,
        extend="both",
    )

    # Configure the aspect ratio and coast lines.
    # --------------------------------------------------------------------------
    # Full world would be aspect 360/(2*180) = 1
    ax.set_aspect((lon_east - lon_west) / (2 * (lat_north - lat_south)))
    ax.coastlines(lw=0.3)

    if not is_global_domain and "RRM" in region_key:
        ax.coastlines(resolution="50m", color="black", linewidth=1)
        state_borders = cfeature.NaturalEarthFeature(
            category="cultural",
            name="admin_1_states_provinces_lakes",
            scale="50m",
            facecolor="none",
        )
        ax.add_feature(state_borders, edgecolor="black")

    # Configure the titles.
    # --------------------------------------------------------------------------
    if title[0] is not None:
        ax.set_title(title[0], loc="left", fontdict=PLOT_SIDE_TITLE)
    if title[1] is not None:
        ax.set_title(title[1], fontdict=PLOT_TITLE)
    if title[2] is not None:
        ax.set_title(title[2], loc="right", fontdict=PLOT_SIDE_TITLE)

    # Configure x and y axis.
    # --------------------------------------------------------------------------
    ax.set_xticks(x_ticks, crs=ccrs.PlateCarree())
    ax.set_yticks(y_ticks, crs=ccrs.PlateCarree())

    lon_formatter = LongitudeFormatter(zero_direction_label=True, number_format=".0f")
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)

    ax.tick_params(labelsize=8.0, direction="out", width=1)

    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")

    # Add and configure the color bar.
    # --------------------------------------------------------------------------
    cbax = fig.add_axes(
        (PANEL[subplot_num][0] + 0.6635, PANEL[subplot_num][1] + 0.0215, 0.0326, 0.1792)
    )
    cbar = fig.colorbar(p1, cax=cbax)

    if c_levels is None:
        cbar.ax.tick_params(labelsize=9.0, length=0)
    else:
        cbar.set_ticks(c_levels[1:-1])

        label_format, pad = _get_contour_label_format_and_pad(c_levels)
        labels = [label_format % level for level in c_levels[1:-1]]
        cbar.ax.set_yticklabels(labels, ha="right")
        cbar.ax.tick_params(labelsize=9.0, pad=pad, length=0)

    # Add metrics text.
    # --------------------------------------------------------------------------
    # Min, Mean, Max
    fig.text(
        PANEL[subplot_num][0] + 0.6635,
        PANEL[subplot_num][1] + 0.2107,
        "Max\nMean\nMin",
        ha="left",
        fontdict=PLOT_SIDE_TITLE,
    )

    fmt_m = []

    # Print in scientific notation if value is greater than 10^5
    for i in range(len(metrics[0:3])):
        fs = "1e" if metrics[i] > 100000.0 else "2f"
        fmt_m.append(fs)

    fmt_metrics = f"%.{fmt_m[0]}\n%.{fmt_m[1]}\n%.{fmt_m[2]}"

    fig.text(
        PANEL[subplot_num][0] + 0.7635,
        PANEL[subplot_num][1] + 0.2107,
        # "%.2f\n%.2f\n%.2f" % stats[0:3],
        fmt_metrics % metrics[0:3],
        ha="right",
        fontdict=PLOT_SIDE_TITLE,
    )

    # RMSE, CORR
    if len(metrics) == 5:
        fig.text(
            PANEL[subplot_num][0] + 0.6635,
            PANEL[subplot_num][1] - 0.0105,
            "RMSE\nCORR",
            ha="left",
            fontdict=PLOT_SIDE_TITLE,
        )
        fig.text(
            PANEL[subplot_num][0] + 0.7635,
            PANEL[subplot_num][1] - 0.0105,
            "%.2f\n%.2f" % metrics[3:5],
            ha="right",
            fontdict=PLOT_SIDE_TITLE,
        )

    # Add grid resolution info.
    # --------------------------------------------------------------------------
    if subplot_num == 2 and "RRM" in region_key:
        dlat = lat[2] - lat[1]
        dlon = lon[2] - lon[1]
        fig.text(
            PANEL[subplot_num][0] + 0.4635,
            PANEL[subplot_num][1] - 0.04,
            "Resolution: {:.2f}x{:.2f}".format(dlat, dlon),
            ha="left",
            fontdict=PLOT_SIDE_TITLE,
        )


def _make_lon_cyclic(var: xr.DataArray):
    """Make the longitude axis cyclic by adding a new coordinate point with 360.

    This function appends a new longitude coordinate point by taking the last
    coordinate point and adding 360 to it.

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


def _get_x_ticks(
    lon_west: float, lon_east: float, is_global_domain: bool, is_lon_full: bool
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
    lon_step = _determine_tick_step(lon_covered)

    x_ticks = np.arange(lon_west, lon_east, lon_step)

    if is_global_domain or is_lon_full:
        # Subtract 0.50 to get 0 W to show up on the right side of the plot.
        # If less than 0.50 is subtracted, then 0 W will overlap 0 E on the
        # left side of the plot.  If a number is added, then the value won't
        # show up at all.
        x_ticks = np.append(x_ticks, lon_east - 0.50)
    else:
        x_ticks = np.append(x_ticks, lon_east)

    return x_ticks


def _get_y_ticks(lat_south: float, lat_north: float) -> np.ndarray:
    """Get Y axis ticks.

    Parameters
    ----------
    lat_south : float
        The south point (e.g., -180).
    lat_north : float
        The north point (e.g., 180).

    Returns
    -------
    np.array
        An array of floats representing Y axis ticks
    """
    lat_covered = lat_north - lat_south

    lat_step = _determine_tick_step(lat_covered)
    y_ticks = np.arange(lat_south, lat_north, lat_step)
    y_ticks = np.append(y_ticks, lat_north)

    return y_ticks


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
    elif degrees_covered > 30:
        return 10
    elif degrees_covered > 20:
        return 5
    else:
        return 1


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

    if maxval < 0.2:
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

from __future__ import annotations

import os
from typing import TYPE_CHECKING, List, Optional, Tuple

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib
import numpy as np
import xarray as xr
import xcdat as xc
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter

from e3sm_diags.derivations.default_regions_xr import REGION_SPECS
from e3sm_diags.driver.utils.general import get_output_dir
from e3sm_diags.logger import custom_logger
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.plot import get_colormap

if TYPE_CHECKING:
    from e3sm_diags.driver.lat_lon_driver import Metrics


matplotlib.use("Agg")
import matplotlib.colors as colors  # isort:skip  # noqa: E402
import matplotlib.pyplot as plt  # isort:skip  # noqa: E402

logger = custom_logger(__name__)

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


def plot(
    da_ref: xr.DataArray,
    da_test: xr.DataArray,
    da_diff: xr.DataArray,
    metrics_dict: Metrics,
    parameter: CoreParameter,
):
    # Create figure, projection
    fig = plt.figure(figsize=parameter.figsize, dpi=parameter.dpi)
    proj = ccrs.PlateCarree()

    # Figure title
    fig.suptitle(parameter.main_title, x=0.5, y=0.96, fontsize=18)

    # The variable units.
    units = metrics_dict["unit"]

    # First two panels
    min1 = metrics_dict["test"]["min"]  # type: ignore
    mean1 = metrics_dict["test"]["mean"]  # type: ignore
    max1 = metrics_dict["test"]["max"]  # type: ignore

    plot_panel(
        0,
        fig,
        proj,
        da_test,
        parameter.contour_levels,
        parameter.test_colormap,
        (parameter.test_name_yrs, parameter.test_title, units),  # type: ignore
        parameter,
        stats=(max1, mean1, min1),  # type: ignore
    )

    if not parameter.model_only:
        min2 = metrics_dict["ref"]["min"]  # type: ignore
        mean2 = metrics_dict["ref"]["mean"]  # type: ignore
        max2 = metrics_dict["ref"]["max"]  # type: ignore

        plot_panel(
            1,
            fig,
            proj,
            da_ref,
            parameter.contour_levels,
            parameter.reference_colormap,
            (parameter.ref_name_yrs, parameter.reference_title, units),  # type: ignore
            parameter,
            stats=(max2, mean2, min2),  # type: ignore
        )

        # Third panel
        min3 = metrics_dict["diff"]["min"]  # type: ignore
        mean3 = metrics_dict["diff"]["mean"]  # type: ignore
        max3 = metrics_dict["diff"]["max"]  # type: ignore
        r = metrics_dict["misc"]["rmse"]  # type: ignore
        c = metrics_dict["misc"]["corr"]  # type: ignore

        plot_panel(
            2,
            fig,
            proj,
            da_diff,
            parameter.diff_levels,
            parameter.diff_colormap,
            (None, parameter.diff_title, units),  # type: ignore
            parameter,
            stats=(max3, mean3, min3, r, c),  # type: ignore
        )

    # Save figure
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
        i = 0
        for p in panels:
            # Extent of subplot
            subpage = np.array(p).reshape(2, 2)
            subpage[1, :] = subpage[0, :] + subpage[1, :]
            subpage = subpage + np.array(BORDER_PADDING).reshape(2, 2)
            subpage = list(((subpage) * page).flatten())
            extent = matplotlib.transforms.Bbox.from_extents(*subpage)
            # Save subplot
            fname = fnm + ".%i." % (i) + f
            plt.savefig(fname, bbox_inches=extent)

            orig_fnm = os.path.join(
                get_output_dir(parameter.current_set, parameter),
                parameter.output_file,
            )
            fname = orig_fnm + ".%i." % (i) + f
            logger.info(f"Sub-plot saved in: {fname}")

            i += 1

    plt.close()


def plot_panel(
    n: int,
    fig: plt.figure,
    proj: ccrs.PlateCarree,
    var: xr.DataArray,
    clevels: List[str],
    cmap: str,
    title: Tuple[Optional[str], str, str],
    parameters: CoreParameter,
    stats: Tuple[float, ...],
):
    # TODO: Refactor the function parameters since clevels, cmap, and title all
    # come from the attributes of the CoreParameter object, `parameters`.
    var = _make_lon_cyclic(var)
    lat = xc.get_dim_coords(var, axis="Y")
    lon = xc.get_dim_coords(var, axis="X")

    var = var.squeeze()

    # Configure contour levels
    # --------------------------------------------------------------------------
    levels = None
    norm = None

    if len(clevels) > 0:
        levels = [-1.0e8] + clevels + [1.0e8]
        norm = colors.BoundaryNorm(boundaries=levels, ncolors=256)

    # Configure plot tickets based on longitude and latitude.
    # --------------------------------------------------------------------------
    region_key = parameters.regions[0]
    region_specs = REGION_SPECS[region_key]

    # Get the region's domain slices for latitude and longitude if set, or
    # use the default value. If both are not set, then the region type is
    # considered "global".
    lat_slice = region_specs.get("lat", (-90, 90))  # type: ignore
    lon_slice = region_specs.get("lon", (0, 360))  # type: ignore
    is_lon_full = lon_slice == (0, 360)

    if lat_slice == (-90, 90) and lon_slice == (0, 360):
        domain_type = "global"
    else:
        domain_type = "region"

    lat_south, lat_north = lat_slice
    lon_west, lon_east = lon_slice

    # Determin X ticks using longitude domain
    x_ticks = _get_x_ticks(lon_west, lon_east, domain_type, is_lon_full)

    # Determine Y axis ticks using latitude domain:
    y_ticks = _get_y_ticks(lat_south, lat_north)

    # Add the contour plot.
    # --------------------------------------------------------------------------
    if domain_type == "global" or is_lon_full:
        proj = ccrs.PlateCarree(central_longitude=180)

    ax = fig.add_axes(PANEL[n], projection=proj)
    ax.set_extent([lon_west, lon_east, lat_south, lat_north], crs=proj)
    cmap = get_colormap(cmap, parameters)
    p1 = ax.contourf(
        lon,
        lat,
        var,
        transform=ccrs.PlateCarree(),
        norm=norm,
        levels=levels,
        cmap=cmap,
        extend="both",
    )

    # Configure the aspect ratio and coast lines.
    # --------------------------------------------------------------------------
    # Full world would be aspect 360/(2*180) = 1
    ax.set_aspect((lon_east - lon_west) / (2 * (lat_north - lat_south)))
    ax.coastlines(lw=0.3)

    if domain_type != "global" and "RRM" in region_key:
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
    ax.set_x_ticks(x_ticks, crs=ccrs.PlateCarree())
    ax.set_y_ticks(y_ticks, crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True, number_format=".0f")
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.tick_params(labelsize=8.0, direction="out", width=1)
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")

    # Add and configure the color bar.
    # --------------------------------------------------------------------------
    cbax = fig.add_axes((PANEL[n][0] + 0.6635, PANEL[n][1] + 0.0215, 0.0326, 0.1792))
    cbar = fig.colorbar(p1, cax=cbax)

    if levels is None:
        cbar.ax.tick_params(labelsize=9.0, length=0)

    else:
        maxval = np.amax(np.absolute(levels[1:-1]))
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

        cbar.set_ticks(levels[1:-1])
        labels = [fmt % level for level in levels[1:-1]]
        cbar.ax.set_yticklabels(labels, ha="right")
        cbar.ax.tick_params(labelsize=9.0, pad=pad, length=0)

    # Add metrics text.
    # --------------------------------------------------------------------------
    # Min, Mean, Max
    fig.text(
        PANEL[n][0] + 0.6635,
        PANEL[n][1] + 0.2107,
        "Max\nMean\nMin",
        ha="left",
        fontdict=PLOT_SIDE_TITLE,
    )

    fmt_m = []

    # Print in scientific notation if value is greater than 10^5
    for i in range(len(stats[0:3])):
        fs = "1e" if stats[i] > 100000.0 else "2f"
        fmt_m.append(fs)

    fmt_metrics = f"%.{fmt_m[0]}\n%.{fmt_m[1]}\n%.{fmt_m[2]}"

    fig.text(
        PANEL[n][0] + 0.7635,
        PANEL[n][1] + 0.2107,
        # "%.2f\n%.2f\n%.2f" % stats[0:3],
        fmt_metrics % stats[0:3],
        ha="right",
        fontdict=PLOT_SIDE_TITLE,
    )

    # RMSE, CORR
    if len(stats) == 5:
        fig.text(
            PANEL[n][0] + 0.6635,
            PANEL[n][1] - 0.0105,
            "RMSE\nCORR",
            ha="left",
            fontdict=PLOT_SIDE_TITLE,
        )
        fig.text(
            PANEL[n][0] + 0.7635,
            PANEL[n][1] - 0.0105,
            "%.2f\n%.2f" % stats[3:5],
            ha="right",
            fontdict=PLOT_SIDE_TITLE,
        )

    # Add grid resolution info.
    # --------------------------------------------------------------------------
    if n == 2 and "RRM" in region_key:
        dlat = lat[2] - lat[1]
        dlon = lon[2] - lon[1]
        fig.text(
            PANEL[n][0] + 0.4635,
            PANEL[n][1] - 0.04,
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
    lon_west: float, lon_east: float, domain_type: str, full_lon: bool
) -> np.array:
    """Get the X axis ticks based on the longitude domain slice.

    Parameters
    ----------
    lon_west : float
        The west point (e.g., 0).
    lon_east : float
        The east point (e.g., 360).
    domain_type : str
        The domain type, either "domain" or "region".
    full_lon : bool
        True if the longitude domain is (0, 360).

    Returns
    -------
    np.array
        An array of floats representing X axis ticks.
    """
    lon_covered = lon_east - lon_west
    lon_step = _determine_tick_step(lon_covered)

    x_ticks = np.arange(lon_west, lon_east, lon_step)

    if domain_type == "global" or full_lon:
        # Subtract 0.50 to get 0 W to show up on the right side of the plot.
        # If less than 0.50 is subtracted, then 0 W will overlap 0 E on the
        # left side of the plot.  If a number is added, then the value won't
        # show up at all.
        x_ticks = np.append(x_ticks, lon_east - 0.50)
    else:
        x_ticks = np.append(x_ticks, lon_east)

    return x_ticks


def _get_y_ticks(lat_south: float, lat_north: float) -> np.array:
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


def _determine_tick_step(degrees_covered: float):
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

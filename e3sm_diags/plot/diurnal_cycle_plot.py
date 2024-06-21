from typing import Tuple

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cdutil
import matplotlib
import numpy as np
import xarray as xr
import xcdat as xc
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
from matplotlib.colors import hsv_to_rgb

from e3sm_diags.derivations.default_regions import regions_specs
from e3sm_diags.logger import custom_logger
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.plot.utils import _save_plot

logger = custom_logger(__name__)

matplotlib.use("Agg")  # noqa: E402
import matplotlib.pyplot as plt  # isort:skip  # noqa: E402


MAIN_TITLE_FONTSIZE = {"fontsize": 11.5}
SECONDARY_TITLE_FONTSIZE = {"fontsize": 9.5}

# Position and sizes of subplot axes in page coordinates (0 to 1)
PANEL_CFG = [
    (0.1691, 0.55, 0.6465, 0.2758),
    (0.1691, 0.15, 0.6465, 0.2758),
]
# Border padding relative to subplot axes for saving individual panels
# (left, bottom, right, top) in page coordinates
BORDER_PADDING = (-0.06, -0.03, 0.13, 0.03)


def plot(
    test_tmax: xr.DataArray,
    test_amp: xr.DataArray,
    ref_tmax: xr.DataArray,
    ref_amp: xr.DataArray,
    parameter: CoreParameter,
):
    # Create figure, projection
    fig = plt.figure(figsize=(8.5, 8.5), dpi=parameter.dpi)

    # First panel
    plot_panel(
        0,
        test_tmax,
        test_amp,
        ref_amp,
        fig,
        parameter,
        (parameter.test_name_yrs, parameter.test_title),
    )

    # Second panel
    plot_panel(
        1,
        ref_tmax,
        ref_amp,
        ref_amp,
        fig,
        parameter,
        (parameter.ref_name_yrs, parameter.reference_title),
    )

    # Figure title
    fig.suptitle(parameter.main_title, x=0.5, y=0.9, fontsize=14)
    _save_plot(fig, parameter, PANEL_CFG, BORDER_PADDING)

    plt.close()


def plot_panel(
    n: int,
    var: xr.DataArray,
    amp: xr.DataArray,
    amp_ref: xr.DataArray,
    fig: plt.Figure,
    parameter,
    title: Tuple[str, str],
):
    normalize_test_amp = parameter.normalize_test_amp
    specified_max_amp = parameter.normalize_amp_int

    lat = xc.get_dim_coords(var, axis="Y")
    var = var.squeeze()
    max_amp = round(amp.max())
    max_amp_ref = round(amp_ref.max())
    amp = amp.squeeze()
    amp_ref = amp_ref.squeeze()

    if normalize_test_amp:
        img = np.dstack(
            ((var / 24 - 0.5) % 1, (amp / max_amp_ref) ** 0.5, np.ones_like(amp))
        )
        max_amp = max_amp_ref
        logger.info(
            f"Scale test diurnal cycle amplitude to max of reference ({max_amp_ref}) mm/day"
        )
    else:
        if specified_max_amp != 0:
            max_amp = specified_max_amp

        img = np.dstack(
            ((var / 24 - 0.5) % 1, (amp / max_amp) ** 0.5, np.ones_like(amp))
        )
        logger.info(f"Scale test diurnal cycle amplitude to specified {max_amp} mm/day")
    # Note: hsv_to_rgb would clipping input data to the valid range for imshow with RGB data ([0..1]
    img = hsv_to_rgb(img)

    # Get region info and X and Y plot ticks.
    # --------------------------------------------------------------------------
    # TODO: This section can be refactored using code from lat_lon_plot.py
    # --------------------------------------------------------------------------
    region_str = parameter.regions[0]
    region = regions_specs[region_str]
    global_domain = True
    full_lon = True

    if "domain" in region.keys():  # type: ignore
        # Get domain to plot
        domain = region["domain"]  # type: ignore
        global_domain = False
    else:
        # Assume global domain
        domain = cdutil.region.domain(latitude=(-90.0, 90, "ccb"))
    kargs = domain.components()[0].kargs
    # lon_west, lon_east, lat_south, lat_north = (0, 360, -90, 90)
    lon_west, lon_east, lat_south, lat_north = (-180, 180, -90, 90)
    if "longitude" in kargs:
        full_lon = False
        lon_west, lon_east, _ = kargs["longitude"]
        # Note cartopy Problem with gridlines across the dateline:https://github.com/SciTools/cartopy/issues/821. Region cross dateline is not supported yet.
        if lon_west > 180 and lon_east > 180:
            lon_west = lon_west - 360
            lon_east = lon_east - 360

    if "latitude" in kargs:
        lat_south, lat_north, _ = kargs["latitude"]

    lon_covered = lon_east - lon_west
    lon_step = determine_tick_step(lon_covered)
    xticks = np.arange(lon_west, lon_east, lon_step)

    lat_covered = lat_north - lat_south
    lat_step = determine_tick_step(lat_covered)
    yticks = np.arange(lat_south, lat_north, lat_step)
    yticks = np.append(yticks, lat_north)

    # Get the cartopy projection based on region info.
    # --------------------------------------------------------------------------
    proj = ccrs.PlateCarree()
    if global_domain or full_lon:
        xticks = [0, 60, 120, 180, 240, 300, 359.99]  # type: ignore
    else:
        xticks = np.append(xticks, lon_east)

    # imshow plot
    projection = ccrs.PlateCarree(central_longitude=180)
    ax = fig.add_axes(PANEL_CFG[n], projection)

    # Configure the aspect ratio and coast lines.
    # --------------------------------------------------------------------------
    # Full world would be aspect 360/(2*180) = 1
    ax.set_extent([lon_west, lon_east, lat_south, lat_north])
    ax.coastlines(lw=0.3)

    # Configure the titles, x and y axes, and colorbar.
    # --------------------------------------------------------------------------
    if title[0] is not None:
        ax.set_title(title[0], loc="left", fontdict=SECONDARY_TITLE_FONTSIZE)
    if title[1] is not None:
        ax.set_title(title[1], fontdict=MAIN_TITLE_FONTSIZE)
    ax.set_xticks(xticks, crs=ccrs.PlateCarree())
    # ax.set_xticks([0, 60, 120, 180, 240, 300, 359.99], crs=ccrs.PlateCarree())
    ax.set_yticks(yticks, crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True, number_format=".0f")
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.tick_params(labelsize=8.0, direction="out", width=1)
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")
    # set a margin around the data
    ax.set_xmargin(0.05)
    ax.set_ymargin(0.10)

    ax.set_ylim([yticks[0], yticks[-1]])

    # add the image. Because this image was a tif, the "origin" of the image is in the
    # upper left corner
    img_extent = [lon_west, lon_east, lat_south, lat_north]

    # When the requested region is inconsistent with what the data covered (`global` is secified but TRMM only has 50S-5N), set an arbitrary threhold.
    if global_domain and lat_covered - abs(lat[0] - lat[-1]) > 10:
        img_extent = [lon_west, lon_east, lat[0], lat[-1]]
    # ax.imshow(img, origin='lower', extent=img_extent, transform=ccrs.PlateCarree())
    ax.imshow(img, origin="lower", extent=img_extent, transform=proj)
    if region_str == "CONUS":
        ax.coastlines(resolution="50m", color="black", linewidth=0.3)
        state_borders = cfeature.NaturalEarthFeature(
            category="cultural",
            name="admin_1_states_provinces_lakes",
            scale="50m",
            facecolor="none",
        )
        ax.add_feature(state_borders, edgecolor="black")

    # Color bar
    bar_ax = fig.add_axes(
        (PANEL_CFG[n][0] + 0.67, PANEL_CFG[n][1] + 0.2, 0.07, 0.07), polar=True
    )
    theta, R = np.meshgrid(np.linspace(0, 2 * np.pi, 24), np.linspace(0, 1, 8))
    H, S = np.meshgrid(np.linspace(0, 1, 24), np.linspace(0, 1, 8))
    image = np.dstack(((H - 0.5) % 1, S**0.5, np.ones_like(S)))
    image = hsv_to_rgb(image)
    # bar_ax.set_theta_zero_location('N')
    bar_ax.set_theta_direction(-1)
    bar_ax.set_theta_offset(np.pi / 2)
    bar_ax.set_xticklabels(["0h", "3h", "6h", "9h", "12h", "15h", "18h", "21h"])
    bar_ax.set_yticklabels(["", "", f"{int(max_amp)}"])
    bar_ax.set_rlabel_position(350)
    bar_ax.get_yticklabels()[-2].set_weight("bold")
    # We change the fontsize of minor ticks label
    bar_ax.tick_params(axis="both", labelsize=7, pad=0, length=0)
    bar_ax.text(
        0.2,
        -0.3,
        "Local Time",
        transform=bar_ax.transAxes,
        fontsize=7,
        verticalalignment="center",
    )
    bar_ax.text(
        -0.1,
        -0.9,
        "Max DC amp {:.2f}{}".format(amp.max(), "mm/d"),
        transform=bar_ax.transAxes,
        fontsize=7,
        fontweight="bold",
        verticalalignment="center",
    )
    bar_ax.text(
        -0.1,
        -0.5,
        "DC phase (Hue)",
        transform=bar_ax.transAxes,
        fontsize=7,
        verticalalignment="center",
    )
    bar_ax.text(
        -0.1,
        -0.7,
        "DC amplitude (Saturation)",
        transform=bar_ax.transAxes,
        fontsize=7,
        verticalalignment="center",
    )
    color = image.reshape((image.shape[0] * image.shape[1], image.shape[2]))
    pc = bar_ax.pcolormesh(theta, R, np.zeros_like(R), color=color, shading="auto")
    pc.set_array(None)


def determine_tick_step(degrees_covered):
    if degrees_covered >= 270:
        return 60
    if degrees_covered >= 180:
        return 30
    if degrees_covered >= 90:
        return 25
    if degrees_covered >= 60:
        return 20
    elif degrees_covered >= 30:
        return 10
    elif degrees_covered >= 20:
        return 5
    else:
        return 1

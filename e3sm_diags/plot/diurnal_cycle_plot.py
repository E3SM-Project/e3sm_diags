import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib
import numpy as np
import xarray as xr
import xcdat as xc
from matplotlib.colors import hsv_to_rgb

from e3sm_diags.derivations.default_regions_xr import REGION_SPECS
from e3sm_diags.logger import _setup_child_logger
from e3sm_diags.parameter.diurnal_cycle_parameter import DiurnalCycleParameter
from e3sm_diags.plot.utils import (
    _configure_titles,
    _configure_x_and_y_axes,
    _get_x_ticks,
    _get_y_ticks,
    _save_plot,
)

logger = _setup_child_logger(__name__)

matplotlib.use("Agg")  # noqa: E402
import matplotlib.pyplot as plt  # isort:skip  # noqa: E402


# Position and sizes of subplot axes in page coordinates (0 to 1)
PANEL_CFG = [
    (0.1691, 0.55, 0.6465, 0.2758),
    (0.1691, 0.15, 0.6465, 0.2758),
]


def plot(
    test_maxtime: xr.DataArray,
    test_amplitude: xr.DataArray,
    ref_maxtime: xr.DataArray,
    ref_amplitude: xr.DataArray,
    parameter: DiurnalCycleParameter,
):
    fig = plt.figure(figsize=(8.5, 8.5), dpi=parameter.dpi)
    fig.suptitle(parameter.main_title, x=0.5, y=0.9, fontsize=14)

    _add_colormap(
        0,
        test_maxtime,
        test_amplitude,
        ref_amplitude,
        fig,
        parameter,
        (parameter.test_name_yrs, parameter.test_title, None),
    )

    _add_colormap(
        1,
        ref_maxtime,
        ref_amplitude,
        ref_amplitude,
        fig,
        parameter,
        (parameter.ref_name_yrs, parameter.reference_title, None),
    )

    _save_plot(fig, parameter, PANEL_CFG)

    plt.close(fig)


def _add_colormap(
    subplot_num: int,
    test_maxtime: xr.DataArray,
    test_amplitude: xr.DataArray,
    ref_maxtime: xr.DataArray,
    fig: plt.Figure,
    parameter: DiurnalCycleParameter,
    title: tuple[str, str, None],
):
    """Adds a colormap containing the test max time, test amplitude, and ref amplitude.

    Parameters
    ----------
    subplot_num : int
        The subplot number
    test_maxtime : xr.DataArray
        The test max time.
    test_amplitude : xr.DataArray
        The test amplitude.
    ref_amplitude : xr.DataArray
        The ref amplitude.
    fig : plt.Figure
        The figure object to add the subplot to.
    parameter : DiurnalCycleParameter
        The parameter object containing plot configurations.
    title : tuple[str, str, None]
        A tuple of strings to fomr the title of the color in the format
        (test_name_yrs, reference_title). None is the third value because
        the child function expects a tuple of three optional strings.
    """
    normalize_test_amp = parameter.normalize_test_amp
    specified_max_amp = parameter.normalize_amp_int

    lat = xc.get_dim_coords(test_maxtime, axis="Y")
    test_maxtime = test_maxtime.squeeze()
    max_amp = round(test_amplitude.max().item())
    max_amp_ref = round(ref_maxtime.max().item())
    test_amplitude = test_amplitude.squeeze()
    ref_maxtime = ref_maxtime.squeeze()

    if normalize_test_amp:
        img = np.dstack(
            (
                (test_maxtime / 24 - 0.5) % 1,
                (test_amplitude / max_amp_ref) ** 0.5,
                np.ones_like(test_amplitude),
            )
        )
        max_amp = max_amp_ref
        logger.info(
            f"Scale test diurnal cycle amplitude to max of reference ({max_amp_ref}) mm/day"
        )
    else:
        if specified_max_amp != 0:
            max_amp = specified_max_amp

        img = np.dstack(
            (
                (test_maxtime / 24 - 0.5) % 1,
                (test_amplitude / max_amp) ** 0.5,
                np.ones_like(test_amplitude),
            )
        )
        logger.info(f"Scale test diurnal cycle amplitude to specified {max_amp} mm/day")

    # NOTE: hsv_to_rgb would clipping input data to the valid range for imshow
    # with RGB data ([0..1]
    img = hsv_to_rgb(img)

    # Get region info and X and Y plot ticks.
    # --------------------------------------------------------------------------
    region_key = parameter.regions[0]
    region_specs = REGION_SPECS[region_key]

    # Get the region's domain slices for latitude and longitude if set, or
    # use the default value. If both are not set, then the region type is
    # considered "global".
    lat_slice = region_specs.get("lat", (-90, 90))  # type: ignore
    lon_slice = region_specs.get("lon", (-180, 180))  # type: ignore

    # Boolean flags for configuring plots.
    is_global_domain = lat_slice == (-90, 90) and lon_slice == (-180, 180)
    is_lon_full = lon_slice == (-180, 180)

    # Determine X and Y ticks using longitude and latitude domains respectively.
    lon_west, lon_east = lon_slice
    x_ticks = _get_x_ticks(
        lon_west,
        lon_east,
        is_global_domain,
        is_lon_full,
        axis_orientation=360,
        tick_step_func=_determine_tick_step,
    )

    lat_south, lat_north = lat_slice
    y_ticks = _get_y_ticks(lat_south, lat_north, tick_step_func=_determine_tick_step)

    # Get the figure Axes object using the projection above.
    # --------------------------------------------------------------------------
    ax = fig.add_axes(
        PANEL_CFG[subplot_num], projection=ccrs.PlateCarree(central_longitude=180)
    )

    # Configure the aspect ratio and coast lines.
    # --------------------------------------------------------------------------
    # Full world would be aspect 360/(2*180) = 1
    ax.set_extent([lon_west, lon_east, lat_south, lat_north])
    ax.coastlines(lw=0.3)

    # Configure the titles and x and y axes.
    # --------------------------------------------------------------------------
    _configure_titles(ax, title)

    _configure_x_and_y_axes(
        ax, x_ticks, y_ticks, ccrs.PlateCarree(), parameter.current_set
    )

    ax.set_xmargin(0.05)
    ax.set_ymargin(0.10)
    ax.set_ylim([y_ticks[0], y_ticks[-1]])

    # # Add the image. Because this image was a tif, the "origin" of the image
    # is in the upper left corner
    # --------------------------------------------------------------------------
    img_extent = [lon_west, lon_east, lat_south, lat_north]

    # Get the cartopy projection based on region info.
    # --------------------------------------------------------------------------
    # When the requested region is inconsistent with what the data covered
    # (`global`) is secified but TRMM only has 50S-5N), set an arbitrary
    # threhold.
    lat_covered = lat_north - lat_south
    if is_global_domain and lat_covered - abs(lat[0] - lat[-1]) > 10:
        img_extent = [lon_west, lon_east, lat[0], lat[-1]]

    imshow_projection = ccrs.PlateCarree()
    if is_global_domain or is_lon_full:
        imshow_projection = ccrs.PlateCarree(central_longitude=180)

    ax.imshow(img, origin="lower", extent=img_extent, transform=imshow_projection)

    if region_key == "CONUS":
        ax.coastlines(resolution="50m", color="black", linewidth=0.3)
        state_borders = cfeature.NaturalEarthFeature(
            category="cultural",
            name="admin_1_states_provinces_lakes",
            scale="50m",
            facecolor="none",
        )
        ax.add_feature(state_borders, edgecolor="black")

    # Configure the colorbar
    # --------------------------------------------------------------------------
    bar_ax = fig.add_axes(
        (PANEL_CFG[subplot_num][0] + 0.67, PANEL_CFG[subplot_num][1] + 0.2, 0.07, 0.07),
        polar=True,
    )

    H, S = np.meshgrid(np.linspace(0, 1, 24), np.linspace(0, 1, 8))
    image = np.dstack(((H - 0.5) % 1, S**0.5, np.ones_like(S)))
    image = hsv_to_rgb(image)

    bar_ax.set_theta_direction(-1)
    bar_ax.set_theta_offset(np.pi / 2)
    bar_ax.set_rlabel_position(350)

    bar_ax.set_xticklabels(["0h", "3h", "6h", "9h", "12h", "15h", "18h", "21h"])
    bar_ax.set_yticklabels(["", "", f"{int(max_amp)}"])
    bar_ax.get_yticklabels()[-2].set_weight("bold")

    # Update the fontsize of minor ticks label.
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
        "Max DC amp {:.2f}{}".format(test_amplitude.max(), "mm/d"),
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

    theta, R = np.meshgrid(np.linspace(0, 2 * np.pi, 24), np.linspace(0, 1, 8))
    color = image.reshape((image.shape[0] * image.shape[1], image.shape[2]))
    pc = bar_ax.pcolormesh(theta, R, np.zeros_like(R), color=color, shading="auto")
    pc.set_array(None)


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

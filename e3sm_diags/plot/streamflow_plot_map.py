from typing import List, Tuple, Union

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib
import numpy as np

from e3sm_diags.derivations.default_regions_xr import REGION_SPECS
from e3sm_diags.logger import custom_logger
from e3sm_diags.parameter.streamflow_parameter import StreamflowParameter
from e3sm_diags.plot.utils import (
    _configure_titles,
    _configure_x_and_y_axes,
    _get_x_ticks,
    _get_y_ticks,
    _save_plot,
)

matplotlib.use("Agg")
import matplotlib.colors as colors  # isort:skip  # noqa: E402
import matplotlib.pyplot as plt  # isort:skip  # noqa: E402

logger = custom_logger(__name__)

# Border padding relative to subplot axes for saving individual panels
# (left, bottom, width, height) in page coordinates.
BORDER_PADDING = (-0.14, -0.06, 0.04, 0.08)

PROJECTION = ccrs.PlateCarree(central_longitude=0)
PROJECTION_FUNC = ccrs.PlateCarree

# Position and sizes of subplot axes in page coordinates (0 to 1)
# (left, bottom, width, height) in page coordinates.
ANNUAL_MAP_PANEL_CFG = [
    (0.1691, 0.6810, 0.6465, 0.2258),
    (0.1691, 0.3961, 0.6465, 0.2258),
    (0.1691, 0.1112, 0.6465, 0.2258),
]


def plot_annual_map(parameter: StreamflowParameter, export_data: np.ndarray):
    """Plot the streamflow annual map.

    Parameters
    ----------
    parameter : StreamflowParameter
        The streamflow parameter.
    export_data : np.ndarray
        The export data.
    """
    fig = plt.figure(figsize=parameter.figsize, dpi=parameter.dpi)
    fig.suptitle(parameter.main_title_annual_map, x=0.5, y=0.97, fontsize=15)

    # Bias between test and ref as a percentage.
    # Relative error as a percentage.
    # 100*((annual_mean_test - annual_mean_ref) / annual_mean_ref)
    bias = 100 * ((export_data[:, 1] - export_data[:, 0]) / export_data[:, 0])

    _plot_panel_annual_map(0, fig, parameter, export_data, bias)
    _plot_panel_annual_map(1, fig, parameter, export_data, bias)
    _plot_panel_annual_map(2, fig, parameter, export_data, bias)

    # NOTE: Need to set the output filename to the name of the specific
    # streamflow plot before saving the plot, otherwise the filename will
    # be blank.
    parameter.output_file = parameter.output_file_annual_map
    _save_plot(fig, parameter, border_padding=BORDER_PADDING)

    plt.close()


def _plot_panel_annual_map(
    panel_index: int,
    fig: plt.Figure,
    parameter: StreamflowParameter,
    export_data: np.ndarray,
    bias_array: np.ndarray,
):
    """Plot the panel for each annual map based on the data type.

    Parameters
    ----------
    panel_index : int
        The panel index.
    fig : plt.Figure
        The figure object.
    parameter : StreamflowParameter
        The streamflow parameter.
    export_data : np.ndarray
        The export data.
    bias_array : np.ndarray
        The bias array.
    """
    if panel_index == 0:
        panel_type = "test"
    elif panel_index == 1:
        panel_type = "ref"
    elif panel_index == 2:
        panel_type = "bias"

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
        tick_step_func=_determine_tick_step,
    )

    lat_south, lat_north = lat_slice
    y_ticks = _get_y_ticks(lat_south, lat_north, tick_step_func=_determine_tick_step)

    # Get the figure Axes object using the projection above and configure the
    # aspect ratio, coastlines, and add RIVERS.
    # --------------------------------------------------------------------------
    ax = fig.add_axes(ANNUAL_MAP_PANEL_CFG[panel_index], projection=PROJECTION)
    ax.set_extent([lon_west, lon_east, lat_south, lat_north], crs=PROJECTION)
    ax.set_aspect((lon_east - lon_west) / (2 * (lat_north - lat_south)))
    ax.coastlines(lw=0.3)
    ax.add_feature(cfeature.RIVERS)

    # Plot of streamflow gauges.
    # --------------------------------------------------------------------------
    color_list, value_min, value_max, norm = setup_annual_map(panel_type)
    _plot_gauges(
        ax,
        panel_type,
        export_data,
        bias_array,
        value_min,
        value_max,
        color_list,
    )

    # Configure the titles, x and y axes.
    # --------------------------------------------------------------------------
    if panel_type == "test":
        title = parameter.test_title
    elif panel_type == "ref":
        title = parameter.reference_title
    elif panel_type == "bias":
        title = "Relative Bias"

    _configure_titles(ax, (None, title, None))
    _configure_x_and_y_axes(
        ax, x_ticks, y_ticks, ccrs.PlateCarree(), parameter.current_set
    )

    # Configure the colorbar.
    # --------------------------------------------------------------------------
    cbax = fig.add_axes(
        (
            ANNUAL_MAP_PANEL_CFG[panel_index][0] + 0.6635,
            ANNUAL_MAP_PANEL_CFG[panel_index][1] + 0.0115,
            0.0326,
            0.1792,
        )
    )
    cmap = colors.ListedColormap(color_list)

    if panel_type in ["test", "ref"]:
        cbar_label = "Mean annual discharge ($m^3$/$s$)"
    elif panel_type == "bias":
        cbar_label = "Bias of mean annual discharge (%)\n(test-ref)/ref"

    cbar = fig.colorbar(
        matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm),
        cax=cbax,
        label=cbar_label,
        extend="both",
    )

    if panel_type == "bias":
        step_size = (value_max - value_min) // 5
        ticks = np.arange(int(value_min), int(value_max) + step_size, step_size)
        cbar.ax.tick_params(labelsize=9.0, length=0)
        cbar.ax.set_yticklabels(ticks)


def setup_annual_map(
    panel_type: str,
) -> Tuple[
    List[str],
    float,
    float,
    Union[matplotlib.colors.LogNorm, matplotlib.colors.Normalize],
]:
    """Set up the annual map based on the panel type.

    Parameters
    ----------
    panel_type : str
        The panel type.

    Returns
    -------
    Tuple[ List[str], float, float, matplotlib.colors.LogNorm | matplotlib.colors.Normalize ]
        A tuple for the color list, the minimum value, the maximum value, and
        the color normalization.
    """
    colormap = plt.get_cmap("jet_r")
    color_list = list(map(lambda index: colormap(index)[:3], range(colormap.N)))

    if panel_type in ["test", "ref"]:
        value_min, value_max = 1, 1e4
        norm = matplotlib.colors.LogNorm(vmin=value_min, vmax=value_max)
    elif panel_type == "bias":
        value_min = -100
        value_max = 100
        norm = matplotlib.colors.Normalize()

    return color_list, value_min, value_max, norm


def _plot_gauges(
    ax: plt.Axes,
    panel_type: str,
    export_data: np.ndarray,
    bias_array: np.ndarray,
    value_min: float,
    value_max: float,
    color_list: List[str],
):
    """Plot the streamflow gauges.

    This function plots each each gauge as a single marker point.

    Parameters
    ----------
    ax : plt.Axes
        The matplotlib axes object.
    panel_type : str
        The panel type.
    export_data : np.ndarray
        The export data.
    bias_array : np.ndarray
        The bias array.
    value_min : float
        The minimum value of the map.
    value_max : float
        The maximum value of the map.
    color_list : List[str]
        The list of colors to use for markers.
    """
    for gauge, i in zip(export_data, range(len(export_data)), strict=False):
        if panel_type == "test":
            value = gauge[1]
        elif panel_type == "ref":
            value = gauge[0]
        elif panel_type == "bias":
            value = bias_array[i]

        if np.isnan(value):
            continue

        if value < value_min:
            value = value_min
        elif value > value_max:
            value = value_max

        if panel_type in ["test", "ref"]:
            # Logarithmic Rescale (min-max normalization) to [-1,1] range
            normalized_value = (np.log10(value) - np.log10(value_min)) / (
                np.log10(value_max) - np.log10(value_min)
            )
        elif panel_type == "bias":
            # Rescale (min-max normalization) to [-1,1] range
            normalized_value = (value - value_min) / (value_max - value_min)

        lat = gauge[7]
        lon = gauge[8]

        color = color_list[int(normalized_value * (len(color_list) - 1))]

        plt.plot(
            lon,
            lat,
            marker="o",
            markersize=2,
            color=color,
            transform=PROJECTION_FUNC(),
        )

        # NOTE: the "plt.annotate call" does not have a "transform=" keyword,
        # so for this one we transform the coordinates with a Cartopy call.
        ax.projection.transform_point(lon, lat, src_crs=PROJECTION_FUNC())


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

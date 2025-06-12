import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib
import numpy as np

from e3sm_diags.derivations.default_regions_xr import REGION_SPECS
from e3sm_diags.logger import _setup_child_logger
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
import matplotlib.lines as lines  # isort:skip  # noqa: E402
import matplotlib.pyplot as plt  # isort:skip  # noqa: E402

logger = _setup_child_logger(__name__)

# Border padding relative to subplot axes for saving individual panels
# (left, bottom, width, height) in page coordinates
BORDER_PADDING = (-0.14, -0.06, 0.04, 0.08)

# Position and sizes of subplot axes in page coordinates (0 to 1)
# (left, bottom, width, height) in page coordinates.
PANEL_CFG = [
    (0.0900, 0.5500, 0.7200, 0.3000),
    (0.0900, 0.1300, 0.7200, 0.3000),
]

# Test and ref color lists, slected from 'hsv' colormap.
COLOR_LIST = [
    (0.05, 0.00, 0.99),
    (0.03, 0.30, 0.98),
    (0.12, 0.92, 0.99),
    (0.13, 1.00, 0.65),
    (0.14, 1.00, 0.05),
    (0.98, 0.99, 0.04),
    (0.99, 0.67, 0.04),
    (0.99, 0.34, 0.03),
    (0.99, 0.07, 0.03),
    (0.99, 0.00, 0.53),
    (0.68, 0.00, 1.00),
    (0.29, 0.00, 1.00),
]

# Dictionary mapping seasonality index to marker size.
SEASONALITY_INDEX = {"si_2": 2, "si_4": 3, "si_6": 4, "si_large": 5}
# Legend elements based on the marker size using the seasonality index dict.
LEGEND_ELEMENTS = [
    lines.Line2D(
        [0],
        [0],
        marker="o",
        color="w",
        label="1 <= SI < 2",
        markerfacecolor="black",
        markersize=SEASONALITY_INDEX["si_2"],
    ),
    lines.Line2D(
        [0],
        [0],
        marker="o",
        color="w",
        label="2 <= SI < 4",
        markerfacecolor="black",
        markersize=SEASONALITY_INDEX["si_4"],
    ),
    lines.Line2D(
        [0],
        [0],
        marker="o",
        color="w",
        label="4 <= SI < 6",
        markerfacecolor="black",
        markersize=SEASONALITY_INDEX["si_6"],
    ),
    lines.Line2D(
        [0],
        [0],
        marker="o",
        color="w",
        label="6 <= SI <= 12",
        markerfacecolor="black",
        markersize=SEASONALITY_INDEX["si_large"],
    ),
]

# Projections to use for the seasonality map.
PROJECTION = ccrs.PlateCarree(central_longitude=0)
PROJECTION_FUNC = ccrs.PlateCarree

# Month labels for the Y Axis.
MONTHS_Y_AXIS_LABEL = [
    "Jan",
    "Feb",
    "Mar",
    "Apr",
    "May",
    "Jun",
    "Jul",
    "Aug",
    "Sep",
    "Oct",
    "Nov",
    "Dec",
]


def plot_seasonality_map(parameter: StreamflowParameter, export_data: np.ndarray):
    """Plot the streamflow seasonality map.

    Parameters
    ----------
    parameter : StreamflowParameter
        The streamflow parameter.
    export_data : np.ndarray
        The export data.
    """
    fig = plt.figure(figsize=parameter.figsize, dpi=parameter.dpi)
    fig.suptitle(parameter.main_title_seasonality_map, x=0.5, y=0.97, fontsize=15)

    _plot_panel_seasonality_map(fig, parameter, "test", export_data)
    _plot_panel_seasonality_map(fig, parameter, "ref", export_data)

    # NOTE: Need to set the output filename to the name of the specific
    # streamflow plot before saving the plot, otherwise the filename will
    # be blank.
    parameter.output_file = parameter.output_file_seasonality_map
    _save_plot(fig, parameter, PANEL_CFG, BORDER_PADDING)
    plt.close(fig)


def _plot_panel_seasonality_map(
    fig: plt.Figure,
    parameter: StreamflowParameter,
    plot_type: str,
    export_data: np.ndarray,
):
    """Plot the panel each seasonality map.

    Parameters
    ----------
    fig : plt.Figure
        The figure object.
    parameter : StreamflowParameter
        The parameter.
    plot_type : str
        The plot type.
    export_data : np.ndarray
        The export data.
    """
    if plot_type == "test":
        panel_idx = 0
        seasonality_idx = 5
        peak_month_idx = 6
        title = (None, parameter.test_title, None)
    elif plot_type == "ref":
        panel_idx = 1
        seasonality_idx = 3
        peak_month_idx = 4
        title = (None, parameter.reference_title, None)

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
    ax = fig.add_axes(PANEL_CFG[panel_idx], projection=PROJECTION)
    ax.set_extent([lon_west, lon_east, lat_south, lat_north], crs=PROJECTION)
    ax.set_aspect((lon_east - lon_west) / (2 * (lat_north - lat_south)))
    ax.coastlines(lw=0.3)
    ax.add_feature(cfeature.RIVERS)

    # Plot the streamflow gauges.
    # --------------------------------------------------------------------------
    _plot_gauges(ax, export_data, seasonality_idx, peak_month_idx)

    # Configure legend.
    # --------------------------------------------------------------------------
    plt.legend(handles=LEGEND_ELEMENTS, title="Seasonality (SI)", prop={"size": 8})

    # Configure the titles, x and y axes.
    # --------------------------------------------------------------------------
    _configure_titles(ax, title)
    _configure_x_and_y_axes(
        ax, x_ticks, y_ticks, ccrs.PlateCarree(), parameter.current_set
    )
    # Configure the colorbar.
    # --------------------------------------------------------------------------
    cbax = fig.add_axes(
        (
            PANEL_CFG[panel_idx][0] + 0.7535,
            PANEL_CFG[panel_idx][1] + 0.0515,
            0.0326,
            0.1792,
        )
    )

    cmap = colors.ListedColormap(COLOR_LIST)

    # Set ticks to be in between the bounds
    num_colors = len(COLOR_LIST)
    bounds = list(range(num_colors))
    ticks = list(map(lambda bound: bound + 0.5, bounds))

    # Add one more bound at the bottom of the colorbar.
    # `bounds` should be one longer than `ticks`.
    bounds += [bounds[-1] + 1]
    norm = colors.BoundaryNorm(bounds, cmap.N)
    cbar = fig.colorbar(
        matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm),
        cax=cbax,
        boundaries=bounds,
        ticks=ticks,
        spacing="uniform",
        orientation="vertical",
        label="Peak month",
    )

    cbar.ax.set_yticklabels(MONTHS_Y_AXIS_LABEL)
    cbar.ax.invert_yaxis()

    cbar.ax.tick_params(labelsize=9.0, length=0)


def _plot_gauges(
    ax: plt.axes,
    export_data: np.ndarray,
    seasonality_idx: int,
    peak_month_idx: int,
):
    """Plot the streamflow gauges.

    This function plots each each gauge as a single marker point.

    Parameters
    ----------
    ax : plt.axes
        The matplotlib axes object.
    export : np.ndarray
        An array of gauges, with each gauge having multiple fields (e.g., lat is
        index 7).
    seasonality_idx_export_idx : int
        The index of the seasonality based on the export index, which determines
        the size of the plot marker for each gauge.
    peak_month_export_idx : int
        The index of the peak month export that determines the color of the
        plot marker for each gauge.

    Raises
    ------
    RuntimeError
        Invalid seasonality index found.
    """
    for gauge in export_data:
        lat = gauge[7]
        lon = gauge[8]
        seasonality_index = gauge[seasonality_idx]

        if seasonality_index < 2:
            markersize = SEASONALITY_INDEX["si_2"]
        elif seasonality_index < 4:
            markersize = SEASONALITY_INDEX["si_4"]
        elif seasonality_index < 6:
            markersize = SEASONALITY_INDEX["si_6"]
        elif seasonality_index <= 12:
            markersize = SEASONALITY_INDEX["si_large"]
        else:
            raise RuntimeError(f"Invalid seasonality index={seasonality_index}")

        if seasonality_index == 1:
            color = "black"
        else:
            peak_month = int(gauge[peak_month_idx])
            color = COLOR_LIST[peak_month]  # type: ignore

        plt.plot(
            lon,
            lat,
            marker="o",
            color=color,
            markersize=markersize,
            transform=PROJECTION_FUNC(),
        )

        # NOTE: The "plt.annotate call" does not have a "transform=" keyword,
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

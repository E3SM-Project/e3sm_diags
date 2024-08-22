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
import matplotlib.lines as lines  # isort:skip  # noqa: E402
import matplotlib.pyplot as plt  # isort:skip  # noqa: E402

logger = custom_logger(__name__)

# Border padding relative to subplot axes for saving individual panels
# (left, bottom, width, height) in page coordinates
BORDER_PADDING = (-0.14, -0.06, 0.04, 0.08)

# Position and sizes of subplot axes in page coordinates (0 to 1)
# (left, bottom, width, height) in page coordinates
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

SEASONALITY_INDEX = {"si_2": 2, "si_4": 4, "si_6": 6, "si_large": 5}
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


PROJECTION = ccrs.PlateCarree(central_longitude=0)
PROJECTION_FUNC = ccrs.PlateCarree

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


def plot_seasonality_map(export: np.ndarray, parameter: StreamflowParameter):
    fig = plt.figure(figsize=parameter.figsize, dpi=parameter.dpi)
    fig.suptitle(parameter.main_title_seasonality_map, x=0.5, y=0.97, fontsize=15)

    _plot_panel_seasonality_map(fig, parameter, "test", export)
    _plot_panel_seasonality_map(fig, parameter, "ref", export)

    # Set the output file name before saving the plot.
    parameter.output_file = parameter.output_file_seasonality_map
    _save_plot(fig, parameter, PANEL_CFG, BORDER_PADDING)

    plt.close()


def _plot_panel_seasonality_map(
    fig: plt.Figure, parameter: StreamflowParameter, plot_type: str, export: np.ndarray
):
    if plot_type == "test":
        panel_index = 0
        seasonality_index_export_index = 5
        peak_month_export_index = 6
        title = (None, parameter.test_title, None)
    elif plot_type == "ref":
        panel_index = 1
        seasonality_index_export_index = 3
        peak_month_export_index = 4
        title = (None, parameter.reference_title, None)
    else:
        raise RuntimeError("Invalid plot_type={}".format(plot_type))

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

    # Get the figure Axes object using the projection above.
    # --------------------------------------------------------------------------
    ax = fig.add_axes(PANEL_CFG[panel_index], projection=PROJECTION)

    # Configure the aspect ratio and coast lines.
    # --------------------------------------------------------------------------
    ax.set_extent([lon_west, lon_east, lat_south, lat_north], crs=PROJECTION)
    # Full world would be aspect 360/(2*180) = 1
    ax.set_aspect((lon_east - lon_west) / (2 * (lat_north - lat_south)))
    ax.coastlines(lw=0.3)

    # Plot of streamflow gauges. Color -> peak month, marker size -> seasonality index.
    # --------------------------------------------------------------------------
    # `export` is the array of gauges. Each gauge has multiple fields -- e.g., lat is index 7
    for gauge in export:
        lat = gauge[7]
        lon = gauge[8]
        seasonality_index = gauge[seasonality_index_export_index]
        if seasonality_index < 2:
            markersize = SEASONALITY_INDEX["si_2"]
        elif seasonality_index < 4:
            markersize = SEASONALITY_INDEX["si_4"]
        elif seasonality_index < 6:
            markersize = SEASONALITY_INDEX["si_6"]
        elif seasonality_index <= 12:
            markersize = SEASONALITY_INDEX["si_large"]
        else:
            raise Exception("Invalid seasonality index={}".format(seasonality_index))
        if seasonality_index == 1:
            color = "black"
        else:
            peak_month = int(gauge[peak_month_export_index])
            color = COLOR_LIST[peak_month]  # type: ignore
        # https://scitools.org.uk/iris/docs/v1.9.2/examples/General/projections_and_annotations.html
        # Place a single marker point for each gauge.
        plt.plot(
            lon,
            lat,
            marker="o",
            color=color,
            markersize=markersize,
            transform=PROJECTION_FUNC(),
        )
        # NOTE: the "plt.annotate call" does not have a "transform=" keyword,
        # so for this one we transform the coordinates with a Cartopy call.
        ax.projection.transform_point(lon, lat, src_crs=PROJECTION_FUNC())

    # Configure legend and add RIVERS feature.
    # --------------------------------------------------------------------------
    seasonality_legend_title = "Seasonality (SI)"
    plt.legend(
        handles=LEGEND_ELEMENTS,
        title=seasonality_legend_title,
        prop={"size": 8},
    )

    ax.add_feature(cfeature.RIVERS)

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
            PANEL_CFG[panel_index][0] + 0.7535,
            PANEL_CFG[panel_index][1] + 0.0515,
            0.0326,
            0.1792,
        )
    )
    # https://matplotlib.org/tutorials/colors/colorbar_only.html
    num_colors = len(COLOR_LIST)
    cmap = colors.ListedColormap(COLOR_LIST)
    cbar_label = "Peak month"

    # Set ticks to be in between the bounds
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
        label=cbar_label,
    )

    cbar.ax.set_yticklabels(MONTHS_Y_AXIS_LABEL)
    cbar.ax.invert_yaxis()

    cbar.ax.tick_params(labelsize=9.0, length=0)


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

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cdutil
import matplotlib
import numpy as np
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter

from e3sm_diags.derivations.default_regions import regions_specs
from e3sm_diags.logger import custom_logger
from e3sm_diags.parameter.streamflow_parameter import StreamflowParameter
from e3sm_diags.plot.utils import _save_plot

matplotlib.use("Agg")
import matplotlib.colors as colors  # isort:skip  # noqa: E402
import matplotlib.pyplot as plt  # isort:skip  # noqa: E402

logger = custom_logger(__name__)

plotTitle = {"fontsize": 11.5}
plotSideTitle = {"fontsize": 9.5}

# Border padding relative to subplot axes for saving individual panels
# (left, bottom, width, height) in page coordinates
BORDER_PADDING = (-0.14, -0.06, 0.04, 0.08)

PROJECTION = ccrs.PlateCarree(central_longitude=0)
PROJECTION_FUNC = ccrs.PlateCarree


# Position and sizes of subplot axes in page coordinates (0 to 1)
# (left, bottom, width, height) in page coordinates
ANNUAL_MAP_PANEL_CFG = [
    (0.1691, 0.6810, 0.6465, 0.2258),
    (0.1691, 0.3961, 0.6465, 0.2258),
    (0.1691, 0.1112, 0.6465, 0.2258),
]


def plot_annual_map(export: np.ndarray, bias, parameter: StreamflowParameter):
    fig = plt.figure(figsize=parameter.figsize, dpi=parameter.dpi)
    fig.suptitle(parameter.main_title_annual_map, x=0.5, y=0.97, fontsize=15)

    _plot_panel_annual_map(0, fig, export, bias, parameter)
    _plot_panel_annual_map(1, fig, export, bias, parameter)
    _plot_panel_annual_map(2, fig, export, bias, parameter)

    # Set the output file name before saving the plot.
    parameter.output_file = parameter.output_file_annual_map
    _save_plot(fig, parameter, border_padding=BORDER_PADDING)

    plt.close()


def _plot_panel_annual_map(panel_index, fig, export, bias_array, parameter):
    if panel_index == 0:
        panel_type = "test"
    elif panel_index == 1:
        panel_type = "ref"
    elif panel_index == 2:
        panel_type = "bias"
    else:
        raise Exception("Invalid panel_index={}".format(panel_index))

    # Plot of streamflow gauges. Color -> peak month, marker size -> seasonality index.

    # Position and sizes of subplot axes in page coordinates (0 to 1)
    ax = fig.add_axes(ANNUAL_MAP_PANEL_CFG[panel_index], projection=PROJECTION)
    region_str = parameter.regions[0]
    region = regions_specs[region_str]
    if "domain" in region.keys():  # type: ignore
        # Get domain to plot
        domain = region["domain"]  # type: ignore
    else:
        # Assume global domain
        domain = cdutil.region.domain(latitude=(-90.0, 90, "ccb"))
    kargs = domain.components()[0].kargs
    # lon_west, lon_east, lat_south, lat_north = (0, 360, -90, 90)
    lon_west, lon_east, lat_south, lat_north = (-180, 180, -90, 90)
    if "longitude" in kargs:
        lon_west, lon_east, _ = kargs["longitude"]
    if "latitude" in kargs:
        lat_south, lat_north, _ = kargs["latitude"]
    lon_covered = lon_east - lon_west
    lon_step = _determine_tick_step(lon_covered)
    xticks = np.arange(lon_west, lon_east, lon_step)
    # Subtract 0.50 to get 0 W to show up on the right side of the plot.
    # If less than 0.50 is subtracted, then 0 W will overlap 0 E on the left side of the plot.
    # If a number is added, then the value won't show up at all.
    xticks = np.append(xticks, lon_east - 0.50)
    lat_covered = lat_north - lat_south
    lat_step = _determine_tick_step(lat_covered)
    yticks = np.arange(lat_south, lat_north, lat_step)
    yticks = np.append(yticks, lat_north)
    ax.set_extent([lon_west, lon_east, lat_south, lat_north], crs=PROJECTION)

    # Stream gauges
    color_list, value_min, value_max, norm = setup_annual_map(
        parameter, panel_type, bias_array
    )
    plot_gauges_annual_map(
        panel_type,
        export,
        bias_array,
        value_min,
        value_max,
        color_list,
        ax,
    )

    # Full world would be aspect 360/(2*180) = 1
    ax.set_aspect((lon_east - lon_west) / (2 * (lat_north - lat_south)))
    ax.coastlines(lw=0.3)
    ax.add_feature(cfeature.RIVERS)
    if panel_type == "test":
        title = parameter.test_title
    elif panel_type == "ref":
        title = parameter.reference_title
    elif panel_type == "bias":
        title = "Relative Bias"
    else:
        raise Exception("Invalid panel_type={}".format(panel_type))
    title = (None, title, None)
    if title[0] is not None:
        ax.set_title(title[0], loc="left", fontdict=plotSideTitle)
    if title[1] is not None:
        ax.set_title(title[1], fontdict=plotTitle)
    if title[2] is not None:
        ax.set_title(title[2], loc="right", fontdict=plotSideTitle)
    ax.set_xticks(xticks, crs=PROJECTION_FUNC())
    ax.set_yticks(yticks, crs=PROJECTION_FUNC())
    lon_formatter = LongitudeFormatter(zero_direction_label=True, number_format=".0f")
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.tick_params(labelsize=8.0, direction="out", width=1)
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")

    # Color bar
    # Position and sizes of subplot axes in page coordinates (0 to 1)
    # (left, bottom, width, height) in page coordinates
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
    else:
        raise Exception("Invalid panel_type={}".format(panel_type))
    cbar = fig.colorbar(
        matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm),
        cax=cbax,
        label=cbar_label,
        extend="both",
    )

    if panel_type in ["test", "ref"]:
        pass
    elif panel_type == "bias":
        step_size = (value_max - value_min) // 5
        ticks = np.arange(int(value_min), int(value_max) + step_size, step_size)
        cbar.ax.tick_params(labelsize=9.0, length=0)
        cbar.ax.set_yticklabels(ticks)
    else:
        raise Exception("Invalid panel_type={}".format(panel_type))


def setup_annual_map(parameter, panel_type, bias_array):
    # Continuous colormap
    colormap = plt.get_cmap("jet_r")
    color_list = list(map(lambda index: colormap(index)[:3], range(colormap.N)))
    if panel_type in ["test", "ref"]:
        value_min, value_max = 1, 1e4
        # https://matplotlib.org/3.2.1/tutorials/colors/colormapnorms.html
        norm = matplotlib.colors.LogNorm(vmin=value_min, vmax=value_max)
    elif panel_type == "bias":
        if parameter.print_statements:
            value_min = np.floor(np.min(bias_array))
            value_max = np.ceil(np.max(bias_array))
            logger.info(
                "Bias of mean annual discharge {} min={}, max={}".format(
                    panel_type, value_min, value_max
                )
            )

        value_min = -100
        value_max = 100
        norm = matplotlib.colors.Normalize()
    else:
        raise Exception("Invalid panel_type={}".format(panel_type))
    return color_list, value_min, value_max, norm


def plot_gauges_annual_map(
    panel_type, export, bias_array, value_min, value_max, color_list, ax
):
    # `export` is the array of gauges. Each gauge has multiple fields -- e.g., lat is index 7
    for gauge, i in zip(export, range(len(export))):
        if panel_type == "test":
            # Test mean annual discharge
            value = gauge[1]
        elif panel_type == "ref":
            # Ref mean annual discharge
            value = gauge[0]
        elif panel_type == "bias":
            # Bias
            value = bias_array[i]
        else:
            raise Exception("Invalid panel_type={}".format(panel_type))
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
        else:
            raise Exception("Invalid panel_type={}".format(panel_type))
        lat = gauge[7]
        lon = gauge[8]

        color = color_list[int(normalized_value * (len(color_list) - 1))]
        # https://scitools.org.uk/iris/docs/v1.9.2/examples/General/projections_and_annotations.html
        # Place a single marker point for each gauge.
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

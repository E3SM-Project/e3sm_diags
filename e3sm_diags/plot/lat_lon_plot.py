from typing import TYPE_CHECKING

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib
import xarray as xr
import xcdat as xc

from e3sm_diags.derivations.default_regions_xr import REGION_SPECS
from e3sm_diags.logger import _setup_child_logger
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.plot.utils import (
    DEFAULT_PANEL_CFG,
    _add_colorbar,
    _add_contour_plot,
    _add_grid_res_info,
    _add_min_mean_max_text,
    _add_rmse_corr_text,
    _configure_titles,
    _configure_x_and_y_axes,
    _get_c_levels_and_norm,
    _get_x_ticks,
    _get_y_ticks,
    _make_lon_cyclic,
    _save_plot,
)

if TYPE_CHECKING:
    from e3sm_diags.driver.lat_lon_driver import MetricsDict


matplotlib.use("Agg")
import matplotlib.pyplot as plt  # isort:skip  # noqa: E402

logger = _setup_child_logger(__name__)


def plot(
    parameter: CoreParameter,
    da_test: xr.DataArray,
    da_ref: xr.DataArray | None,
    da_diff: xr.DataArray | None,
    metrics_dict: MetricsDict,
):
    """Plot the variable's metrics generated for the lat_lon set.

    Parameters
    ----------
    parameter : CoreParameter
        The CoreParameter object containing plot configurations.
    da_test : xr.DataArray
        The test data.
    da_ref : xr.DataArray | None
        The optional reference data.
    da_diff : xr.DataArray | None
        The difference between `da_test` and `da_ref` (both are gridded to
        the lower resolution of the two beforehand).
    metrics_dict : Metrics
        The metrics.
    """
    fig = plt.figure(figsize=parameter.figsize, dpi=parameter.dpi)
    fig.suptitle(parameter.main_title, x=0.5, y=0.96, fontsize=18)

    # The variable units.
    units = metrics_dict["unit"]

    # Add the first subplot for test data.
    min1 = metrics_dict["test"]["min"]  # type: ignore
    mean1 = metrics_dict["test"]["mean"]  # type: ignore
    max1 = metrics_dict["test"]["max"]  # type: ignore

    _add_colormap(
        0,
        da_test,
        fig,
        parameter,
        parameter.test_colormap,
        parameter.contour_levels,
        title=(parameter.test_name_yrs, parameter.test_title, units),  # type: ignore
        metrics=(max1, mean1, min1),  # type: ignore
    )

    # Add the second and third subplots for ref data and the differences,
    # respectively.
    if da_ref is not None and da_diff is not None:
        min2 = metrics_dict["ref"]["min"]  # type: ignore
        mean2 = metrics_dict["ref"]["mean"]  # type: ignore
        max2 = metrics_dict["ref"]["max"]  # type: ignore

        _add_colormap(
            1,
            da_ref,
            fig,
            parameter,
            parameter.reference_colormap,
            parameter.contour_levels,
            title=(parameter.ref_name_yrs, parameter.reference_title, units),  # type: ignore
            metrics=(max2, mean2, min2),  # type: ignore
        )

        min3 = metrics_dict["diff"]["min"]  # type: ignore
        mean3 = metrics_dict["diff"]["mean"]  # type: ignore
        max3 = metrics_dict["diff"]["max"]  # type: ignore
        r = metrics_dict["misc"]["rmse"]  # type: ignore
        c = metrics_dict["misc"]["corr"]  # type: ignore

        _add_colormap(
            2,
            da_diff,
            fig,
            parameter,
            parameter.diff_colormap,
            parameter.diff_levels,
            title=(None, parameter.diff_title, units),  # type: ignore
            metrics=(max3, mean3, min3, r, c),  # type: ignore
        )

    _save_plot(fig, parameter)

    plt.close(fig)


def _add_colormap(
    subplot_num: int,
    var: xr.DataArray,
    fig: plt.Figure,
    parameter: CoreParameter,
    color_map: str,
    contour_levels: list[float],
    title: tuple[str | None, str, str],
    metrics: tuple[float, ...],
):
    """Adds a colormap containing the variable data and metrics to the figure.

    This function is used by:
      - `lat_lon_plot.py`
      - `aerosol_aeronet_plot.py` (when refactored).

    Parameters
    ----------
    subplot_num : int
        The subplot number.
    var : xr.DataArray
        The variable to plot.
    fig : plt.Figure
        The figure object to add the subplot to.
    parameter : CoreParameter
        The CoreParameter object containing plot configurations.
    color_map : str
        The colormap styling to use (e.g., "cet_rainbow.rgb").
    contour_levels : list[float]
        The map contour levels.
    title : tuple[str | None, str, str]
        A tuple of strings to form the title of the colormap, in the format
        (<optional> years, title, units).
    metrics : tuple[float, ...]
        A tuple of metrics for this subplot.
    """
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
    x_ticks = _get_x_ticks(lon_west, lon_east, is_global_domain, is_lon_full)

    lat_south, lat_north = lat_slice
    y_ticks = _get_y_ticks(lat_south, lat_north)

    # Get the cartopy projection based on region info.
    # --------------------------------------------------------------------------
    projection = ccrs.PlateCarree()
    if is_global_domain or is_lon_full:
        projection = ccrs.PlateCarree(central_longitude=180)

    # Get the figure Axes object using the projection above.
    # --------------------------------------------------------------------------
    ax = fig.add_axes(DEFAULT_PANEL_CFG[subplot_num], projection=projection)
    ax.set_extent([lon_west, lon_east, lat_south, lat_north], crs=projection)
    contour_plot = _add_contour_plot(
        ax, var, lon, lat, color_map, ccrs.PlateCarree(), norm, c_levels
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

    # Configure the titles, x and y axes, and colorbar.
    # --------------------------------------------------------------------------
    _configure_titles(ax, title)
    _configure_x_and_y_axes(
        ax, x_ticks, y_ticks, ccrs.PlateCarree(), parameter.current_set
    )
    _add_colorbar(fig, subplot_num, DEFAULT_PANEL_CFG, contour_plot, c_levels)

    # Add metrics text to the figure.
    # --------------------------------------------------------------------------
    _add_min_mean_max_text(fig, subplot_num, DEFAULT_PANEL_CFG, metrics)

    if len(metrics) == 5:
        _add_rmse_corr_text(fig, subplot_num, DEFAULT_PANEL_CFG, metrics)

    _add_grid_res_info(fig, subplot_num, region_key, lat, lon, DEFAULT_PANEL_CFG)

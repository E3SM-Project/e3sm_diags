from typing import List, Literal, Tuple, Union

import cartopy.crs as ccrs
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import xcdat as xc

from e3sm_diags.driver.utils.type_annotations import MetricsDict
from e3sm_diags.logger import custom_logger
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.plot.utils import (
    _add_colorbar,
    _add_contour_plot,
    _add_min_mean_max_text,
    _add_rmse_corr_text,
    _configure_titles,
    _get_c_levels_and_norm,
    _make_lon_cyclic,
    _save_plot,
)

matplotlib.use("Agg")
import matplotlib.path as mpath  # isort:skip  # noqa: E402

logger = custom_logger(__name__)

PLOT_SECONDARY_TITLE = 9.0

# Position and sizes of subplot axes in page coordinates (0 to 1)
PANEL_CFG = [
    (0.27, 0.65, 0.3235, 0.25),
    (0.27, 0.35, 0.3235, 0.25),
    (0.27, 0.05, 0.3235, 0.25),
]

# Border padding relative to subplot axes for saving individual panels
# (left, bottom, right, top) in page coordinates
BORDER_PADDING = (-0.02, -0.01, 0.14, 0.04)


def plot(
    parameter: CoreParameter,
    da_test: xr.DataArray,
    da_ref: xr.DataArray,
    da_diff: xr.DataArray,
    metrics_dict: MetricsDict,
):
    """Plot the variable's metrics generated for the lat_lon set.

    Parameters
    ----------
    parameter : CoreParameter
        The CoreParameter object containing plot configurations.
    da_test : xr.DataArray
        The test data.
    da_ref : xr.DataArray
        The reference data.
    da_diff : xr.DataArray
        The difference between `da_test` and `da_ref` (both are gridded to
        the lower resolution of the two beforehand).
    metrics_dict : Metrics
        The metrics.
    """
    fig = plt.figure(figsize=parameter.figsize, dpi=parameter.dpi)

    units = metrics_dict["units"]

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

    # Figure title
    fig.suptitle(parameter.main_title, x=0.5, y=0.97, fontsize=18)

    _save_plot(fig, parameter, PANEL_CFG, BORDER_PADDING)

    plt.close()


def _add_colormap(
    subplot_num: int,
    var: xr.DataArray,
    fig: plt.Figure,
    parameter: CoreParameter,
    color_map: str,
    contour_levels: List[float],
    title: Tuple[Union[str, None], str, str],
    metrics: Tuple[float, ...],
):
    """Adds a colormap containing the variable data and metrics to the figure.

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
        The colormap stylin˘g to use (e.g., "cet_rainbow.rgb").
    contour_levels : List[float]
        The map contour levels.
    title : Tuple[Union[str, None], str, str]
        A tuple of strings to form the title of the colormap, in the format
        (<optional> years, title, units).
    metrics : Tuple[float, ...]
        A tuple of metrics for this subplot.
    """
    var = _make_lon_cyclic(var)
    lat = xc.get_dim_coords(var, axis="Y")
    lon = xc.get_dim_coords(var, axis="X")

    var = var.squeeze()

    pole, proj = _get_pole_and_projection(parameter)

    # Configure contour levels and boundary norm.
    # --------------------------------------------------------------------------
    c_levels, norm = _get_c_levels_and_norm(contour_levels)

    # Configure the figure Axes object using the projection above.
    # --------------------------------------------------------------------------
    ax = fig.add_axes(PANEL_CFG[subplot_num], projection=proj)
    ax.set_global()

    ax.gridlines()
    if pole == "N":
        ax.set_extent([-180, 180, 50, 90], crs=ccrs.PlateCarree())
    elif pole == "S":
        ax.set_extent([-180, 180, -55, -90], crs=ccrs.PlateCarree())

    theta = np.linspace(0, 2 * np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    ax.set_boundary(circle, transform=ax.transAxes)

    # Configure contour plot.
    # --------------------------------------------------------------------------
    contour_plot = _add_contour_plot(
        ax, var, lon, lat, color_map, ccrs.PlateCarree(), norm, c_levels
    )

    ax.set_aspect("auto")
    ax.coastlines(lw=0.3)

    # Configure the titles, and colorbar.
    # --------------------------------------------------------------------------
    _configure_titles(ax, title, secondary_fontsize=PLOT_SECONDARY_TITLE)
    _add_colorbar(
        fig,
        subplot_num,
        PANEL_CFG,
        contour_plot,
        c_levels,
        rect=(0.35, 0.0354, 0.0326, 0.1792),
    )

    _add_min_mean_max_text(
        fig,
        subplot_num,
        PANEL_CFG,
        metrics,
        left_text_pos=(0.35, 0.225),
        right_text_pos=(0.45, 0.225),
    )

    if len(metrics) == 5:
        _add_rmse_corr_text(
            fig,
            subplot_num,
            PANEL_CFG,
            metrics,
            left_text_pos=(0.35, 0.0),
            right_text_pos=(0.45, 0.0),
        )


def _get_pole_and_projection(
    parameter: CoreParameter,
) -> Tuple[Literal["N", "S"], Union[ccrs.NorthPolarStereo, ccrs.SouthPolarStereo]]:
    var_region = parameter.var_region

    if var_region.find("N") != -1:
        pole = "N"
        proj = ccrs.NorthPolarStereo(central_longitude=0)
    elif var_region.find("S") != -1:
        pole = "S"
        proj = ccrs.SouthPolarStereo(central_longitude=0)
    else:
        raise RuntimeError(
            f"The variable region ('{var_region}') does not contain 'N' or 'S' for "
            "polar set plotting."
        )

    return pole, proj  # type: ignore

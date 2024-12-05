from typing import List, Tuple, Union

import matplotlib
import xarray as xr
import xcdat as xc

from e3sm_diags.driver.utils.type_annotations import MetricsDict
from e3sm_diags.logger import custom_logger
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.plot.utils import (
    DEFAULT_PANEL_CFG,
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
import matplotlib.pyplot as plt  # isort:skip  # noqa: E402

logger = custom_logger(__name__)


def plot(
    parameter: CoreParameter,
    da_test: xr.DataArray,
    da_ref: xr.DataArray,
    da_diff: xr.DataArray,
    metrics_dict: MetricsDict,
):
    # Create figure, projection
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

    fig.suptitle(parameter.main_title, x=0.5, y=0.96, fontsize=18)

    _save_plot(fig, parameter)

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
    var = _make_lon_cyclic(var)
    lon = xc.get_dim_coords(var, axis="X")
    plev = xc.get_dim_coords(var, axis="Z")

    var = var.squeeze()

    # Configure contour levels and boundary norm.
    # --------------------------------------------------------------------------
    c_levels, norm = _get_c_levels_and_norm(contour_levels)

    # Configure contour plot.
    # --------------------------------------------------------------------------
    ax = fig.add_axes(DEFAULT_PANEL_CFG[subplot_num])
    contour_plot = _add_contour_plot(
        ax, var, lon, plev, color_map, None, norm, c_levels
    )

    # Configure the aspect ratio and plot titles.
    # --------------------------------------------------------------------------
    ax.set_aspect("auto")

    _configure_titles(ax, title)

    # Configure x and y axis.
    # --------------------------------------------------------------------------
    ax.set_xticks([0, 60, 120, 180, 240, 300, 359.99])
    ax.set_xticklabels(["0", "60E", "120E", "180", "120W", "60W", "0"])

    ax.tick_params(labelsize=8.0, direction="out", width=1)
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")

    if parameter.plot_log_plevs:
        ax.set_yscale("log")

    if parameter.plot_plevs:
        plev_ticks = parameter.plevs
        plt.yticks(plev_ticks, plev_ticks)

    ax.invert_yaxis()
    ax.set_ylabel("Pressure (mb)")

    # Add and configure the color bar.
    # --------------------------------------------------------------------------
    _add_colorbar(fig, subplot_num, DEFAULT_PANEL_CFG, contour_plot, c_levels)

    # Add metrics text to the figure.
    # --------------------------------------------------------------------------
    _add_min_mean_max_text(fig, subplot_num, DEFAULT_PANEL_CFG, metrics)

    if len(metrics) == 5:
        _add_rmse_corr_text(fig, subplot_num, DEFAULT_PANEL_CFG, metrics)

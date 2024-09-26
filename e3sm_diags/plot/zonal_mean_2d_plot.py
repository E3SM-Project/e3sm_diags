from typing import List, Optional, Tuple

import matplotlib
import numpy as np
import xarray as xr
import xcdat as xc

from e3sm_diags.driver.utils.type_annotations import MetricsDict
from e3sm_diags.logger import custom_logger
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.parameter.zonal_mean_2d_parameter import DEFAULT_PLEVS
from e3sm_diags.plot.utils import (
    DEFAULT_PANEL_CFG,
    _add_colorbar,
    _add_contour_plot,
    _add_min_mean_max_text,
    _add_rmse_corr_text,
    _configure_titles,
    _configure_x_and_y_axes,
    _get_c_levels_and_norm,
    _save_plot,
)

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # isort:skip  # noqa: E402

logger = custom_logger(__name__)


# Configs for x axis ticks and x axis limits.
X_TICKS = np.array([-90, -60, -30, 0, 30, 60, 90])
X_LIM = -90, 90


def plot(
    parameter: CoreParameter,
    da_test: xr.DataArray,
    da_ref: xr.DataArray,
    da_diff: xr.DataArray,
    metrics_dict: MetricsDict,
):
    """Plot the variable's metrics generated by the zonal_mean_2d set.

    Parameters
    ----------
    parameter : CoreParameter
        The CoreParameter object containing plot configurations.
    da_test : xr.DataArray
        The test data.
    da_ref : xr.DataArray
        The reference data.
    da_diff : xr.DataArray
        The difference between `da_test` and `da_ref` (both are regridded to
        the lower resolution of the two beforehand).
    metrics_dict : Metrics
        The metrics.
    """
    fig = plt.figure(figsize=parameter.figsize, dpi=parameter.dpi)
    fig.suptitle(parameter.main_title, x=0.5, y=0.96, fontsize=18)

    # The variable units.
    units = metrics_dict["units"]

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
        title=(None, parameter.diff_title, da_diff.attrs["units"]),
        metrics=(max3, mean3, min3, r, c),  # type: ignore
    )

    _save_plot(fig, parameter)

    plt.close()


def _add_colormap(
    subplot_num: int,
    var: xr.DataArray,
    fig: plt.Figure,
    parameter: CoreParameter,
    color_map: str,
    contour_levels: List[float],
    title: Tuple[Optional[str], str, str],
    metrics: Tuple[float, ...],
):
    lat = xc.get_dim_coords(var, axis="Y")
    plev = xc.get_dim_coords(var, axis="Z")
    var = var.squeeze()

    # Configure contour levels
    # --------------------------------------------------------------------------
    c_levels, norm = _get_c_levels_and_norm(contour_levels)

    # Add the contour plot
    # --------------------------------------------------------------------------
    ax = fig.add_axes(DEFAULT_PANEL_CFG[subplot_num], projection=None)

    contour_plot = _add_contour_plot(
        ax, parameter, var, lat, plev, color_map, None, norm, c_levels
    )

    # Configure the aspect ratio and plot titles.
    # --------------------------------------------------------------------------
    ax.set_aspect("auto")
    _configure_titles(ax, title)

    # Configure x and y axis.
    # --------------------------------------------------------------------------
    _configure_x_and_y_axes(ax, X_TICKS, None, None, parameter.current_set)
    ax.set_xlim(X_LIM)

    if parameter.plot_log_plevs:
        ax.set_yscale("log")

    if parameter.plot_plevs:
        plev_ticks = parameter.plevs
        plt.yticks(plev_ticks, plev_ticks)

    # For default plevs, specify the pressure axis and show the 50 mb tick
    # at the top.
    if (
        not parameter.plot_log_plevs
        and not parameter.plot_plevs
        and parameter.plevs == DEFAULT_PLEVS
    ):
        plev_ticks = parameter.plevs
        new_ticks = [plev_ticks[0]] + plev_ticks[1::2]
        new_ticks = [int(x) for x in new_ticks]
        plt.yticks(new_ticks, new_ticks)

    plt.ylabel("pressure (mb)")
    ax.invert_yaxis()

    # Add and configure the color bar.
    # --------------------------------------------------------------------------
    _add_colorbar(fig, subplot_num, DEFAULT_PANEL_CFG, contour_plot, c_levels)

    # Add metrics text.
    # --------------------------------------------------------------------------
    # Min, Mean, Max
    _add_min_mean_max_text(fig, subplot_num, DEFAULT_PANEL_CFG, metrics)

    if len(metrics) == 5:
        _add_rmse_corr_text(fig, subplot_num, DEFAULT_PANEL_CFG, metrics)
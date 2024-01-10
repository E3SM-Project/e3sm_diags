from typing import List, Tuple

import matplotlib
import xarray as xr
import xcdat as xc

from e3sm_diags.driver.utils.type_annotations import MetricsDict
from e3sm_diags.logger import custom_logger
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.parameter.zonal_mean_2d_parameter import DEFAULT_PLEVS
from e3sm_diags.plot import get_colormap
from e3sm_diags.plot.utils import (
    PANEL,
    _add_min_mean_max_text,
    _add_rmse_corr_text,
    _configure_cbar,
    _configure_titles,
    _save_plot,
)

matplotlib.use("Agg")
import matplotlib.colors as colors  # isort:skip  # noqa: E402
import matplotlib.pyplot as plt  # isort:skip  # noqa: E402

logger = custom_logger(__name__)


# Configs for x axis ticks and x axis limits.
X_TICKS = [-90, -60, -30, 0, 30, 60, 90]
X_LIM = -90, 90


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
    da_ref : xr.DataArray | None
        The optional reference data.
    ds_diff : xr.DataArray | None
        The difference between ``ds_test_regrid`` and ``ds_ref_regrid``.
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

    plt.close()


def _add_colormap(
    subplot_num: int,
    var: xr.DataArray,
    fig: plt.Figure,
    parameter: CoreParameter,
    color_map: str,
    contour_levels: List[float],
    title: Tuple[str | None, str, str],
    metrics: Tuple[float, ...],
):
    lat = xc.get_dim_coords(var, axis="Y")
    plev = xc.get_dim_coords(var, axis="Z")
    var = var.squeeze()

    # Configure contour levels
    # --------------------------------------------------------------------------
    c_levels = None
    norm = None

    if len(contour_levels) > 0:
        c_levels = [-1.0e8] + contour_levels + [1.0e8]
        norm = colors.BoundaryNorm(boundaries=c_levels, ncolors=256)

    # Add the contour plot
    # --------------------------------------------------------------------------
    ax = fig.add_axes(PANEL[subplot_num], projection=None)
    color_map = get_colormap(color_map, parameter)

    p1 = ax.contourf(
        lat,
        plev,
        var,
        norm=norm,
        levels=c_levels,
        cmap=color_map,
        extend="both",
    )

    # Configure the aspect ratio.
    # --------------------------------------------------------------------------
    ax.set_aspect("auto")

    # Configure the titles.
    # --------------------------------------------------------------------------
    ax = _configure_titles(ax, title)

    # Configure x and y axis.
    # --------------------------------------------------------------------------
    ax.set_xticks(X_TICKS)
    ax.set_xlim(X_LIM)

    ax.tick_params(labelsize=8.0, direction="out", width=1)

    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")

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
    cbax = fig.add_axes(
        (PANEL[subplot_num][0] + 0.6635, PANEL[subplot_num][1] + 0.0215, 0.0326, 0.1792)
    )

    cbar = fig.colorbar(p1, cax=cbax)
    cbar = _configure_cbar(cbar, c_levels)

    # Add metrics text.
    # --------------------------------------------------------------------------
    # Min, Mean, Max
    fig = _add_min_mean_max_text(subplot_num, fig, metrics)

    if len(metrics) == 5:
        fig = _add_rmse_corr_text(subplot_num, fig, metrics)

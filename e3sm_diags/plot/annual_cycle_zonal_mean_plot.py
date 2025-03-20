from typing import List, Optional, Tuple

import matplotlib
import numpy as np
import xarray as xr
import xcdat as xc
from cartopy.mpl.ticker import LatitudeFormatter

from e3sm_diags.logger import _setup_child_logger
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.plot.utils import (
    DEFAULT_PANEL_CFG,
    _add_colorbar,
    _add_contour_plot,
    _configure_titles,
    _get_c_levels_and_norm,
    _save_plot,
)

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # isort:skip  # noqa: E402

logger = _setup_child_logger(__name__)


# Configs for x axis ticks and x axis limits.
X_TICKS = ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"]
Y_TICKS = np.array([-90, -60, -30, 0, 30, 60, 90])
Y_LIM = -90, 90


def plot(
    parameter: CoreParameter,
    da_test: xr.DataArray,
    da_ref: xr.DataArray,
    da_diff: xr.DataArray,
):
    """Plot the variable's metrics generated by the annual_cycle_zonal_mean set.

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
    """
    fig = plt.figure(figsize=parameter.figsize, dpi=parameter.dpi)
    fig.suptitle(parameter.main_title, x=0.5, y=0.96, fontsize=18)

    units = da_test.units

    _add_colormap(
        0,
        da_test,
        fig,
        parameter,
        parameter.test_colormap,
        parameter.contour_levels,
        title=(parameter.test_name_yrs, parameter.test_title, units),
    )

    _add_colormap(
        1,
        da_ref,
        fig,
        parameter,
        parameter.reference_colormap,
        parameter.contour_levels,
        title=(parameter.ref_name_yrs, parameter.reference_title, units),
    )

    _add_colormap(
        2,
        da_diff,
        fig,
        parameter,
        parameter.diff_colormap,
        parameter.diff_levels,
        title=(None, parameter.diff_title, da_diff.attrs["units"]),
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
):
    lat = xc.get_dim_coords(var, axis="Y")
    time = xc.get_dim_coords(var, axis="T")

    var = var.squeeze()

    # Configure contour levels
    # --------------------------------------------------------------------------
    c_levels, norm = _get_c_levels_and_norm(contour_levels)

    # Add the contour plot
    # --------------------------------------------------------------------------
    ax = fig.add_axes(DEFAULT_PANEL_CFG[subplot_num], projection=None)
    var = var.transpose(lat.name, time.name)
    contour_plot = _add_contour_plot(
        ax, var, time, lat, color_map, None, norm, c_levels
    )

    # Configure the aspect ratio and plot titles.
    # --------------------------------------------------------------------------
    ax.set_aspect("auto")
    _configure_titles(ax, title)

    # Configure x and y axis.
    # --------------------------------------------------------------------------
    plt.xticks(time, X_TICKS)
    lat_formatter = LatitudeFormatter()
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.tick_params(labelsize=8.0, direction="out", width=1)
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")

    # Add and configure the color bar.
    # --------------------------------------------------------------------------
    _add_colorbar(fig, subplot_num, DEFAULT_PANEL_CFG, contour_plot, c_levels)

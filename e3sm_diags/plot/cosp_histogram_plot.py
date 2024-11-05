from typing import List, Tuple, Union

import matplotlib
import numpy as np
import xarray as xr

from e3sm_diags.logger import custom_logger
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.plot.utils import _get_colormap, _save_plot

matplotlib.use("Agg")
import matplotlib.colors as colors  # isort:skip  # noqa: E402
import matplotlib.pyplot as plt  # isort:skip  # noqa: E402

logger = custom_logger(__name__)

MAIN_TITLE_FONT = 11.5
SECONDARY_TITLE_FONT = 9.5


# Position and sizes of subplot axes in page coordinates (0 to 1)
PANEL_CFG = [
    (0.1691, 0.6810, 0.6465, 0.2150),
    (0.1691, 0.3961, 0.6465, 0.2150),
    (0.1691, 0.1112, 0.6465, 0.2150),
]

# Border padding relative to subplot axes for saving individual panels
# (left, bottom, right, top) in page coordinates
BORDER_PADDING = (-0.10, -0.05, 0.13, 0.033)


def plot(
    parameter: CoreParameter,
    da_test: xr.DataArray,
    da_ref: xr.DataArray,
    da_diff: xr.DataArray,
):
    fig = plt.figure(figsize=parameter.figsize, dpi=parameter.dpi)

    _add_colormap(
        0,
        da_test,
        fig,
        parameter.contour_levels,
        "rainbow",
        title=(parameter.test_name_yrs, parameter.test_title, da_test.units),
    )
    _add_colormap(
        1,
        da_ref,
        fig,
        parameter.contour_levels,
        "rainbow",
        title=(parameter.ref_name_yrs, parameter.reference_title, da_test.units),
    )
    _add_colormap(
        2,
        da_diff,
        fig,
        parameter.diff_levels,
        parameter.diff_colormap,
        title=(parameter.diff_name, parameter.diff_title, da_test.units),
    )

    # Figure title
    fig.suptitle(parameter.main_title, x=0.5, y=0.96, fontsize=18)

    _save_plot(fig, parameter, PANEL_CFG, BORDER_PADDING)

    plt.close()


def _add_colormap(
    subplot_num: int,
    var: xr.DataArray,
    fig: plt.Figure,
    contour_levels: List[float],
    color_map: str,
    title: Tuple[Union[str, None], str, str],
):
    # Contour levels
    levels = None
    norm = None
    if len(contour_levels) > 0:
        levels = [-1.0e8] + contour_levels + [1.0e8]
        norm = colors.BoundaryNorm(boundaries=levels, ncolors=256)

    # Contour plot
    ax = fig.add_axes(PANEL_CFG[subplot_num])

    color_map = _get_colormap(color_map)
    p1 = plt.pcolormesh(var, cmap=color_map, norm=norm)
    # Calculate 3 x 3 grids for cloud fraction for nine cloud class
    # Place cloud fraction of each cloud class in plot:
    cld_3x3 = np.zeros((3, 3))
    for j in range(3):
        for i in range(3):
            if "MISR" in str(var.name):
                if j == 0:
                    cld_3x3[j, i] = var[0:6, 2 * i : 2 * i + 2].sum()
                    ax.text(
                        i * 2 + 1,
                        3,
                        "%.1f" % cld_3x3[j, i],
                        horizontalalignment="center",
                        verticalalignment="center",
                        fontsize=25,
                    )
                elif j == 1:
                    cld_3x3[j, i] = var[6:9, 2 * i : 2 * i + 2].sum()
                    ax.text(
                        i * 2 + 1,
                        7.5,
                        "%.1f" % cld_3x3[j, i],
                        horizontalalignment="center",
                        verticalalignment="center",
                        fontsize=25,
                    )
                elif j == 2:
                    cld_3x3[j, i] = var[9:, 2 * i : 2 * i + 2].sum()
                    ax.text(
                        i * 2 + 1,
                        12,
                        "%.1f" % cld_3x3[j, i],
                        horizontalalignment="center",
                        verticalalignment="center",
                        fontsize=25,
                    )

            else:
                if j == 2:
                    cld_3x3[j, i] = var[4:7, 2 * i : 2 * i + 2].sum()
                    ax.text(
                        i * 2 + 1,
                        j * 2 + 1.5,
                        "%.1f" % cld_3x3[j, i],
                        horizontalalignment="center",
                        verticalalignment="center",
                        fontsize=25,
                    )
                else:
                    cld_3x3[j, i] = var[2 * j : 2 * j + 2, 2 * i : 2 * i + 2].sum()
                    ax.text(
                        i * 2 + 1,
                        j * 2 + 1,
                        "%.1f" % cld_3x3[j, i],
                        horizontalalignment="center",
                        verticalalignment="center",
                        fontsize=25,
                    )

    cld_3x3.sum()
    # Place vertical/horizonal line to separate cloud class
    plt.axvline(x=2, linewidth=2, color="k")
    plt.axvline(x=4, linewidth=2, color="k")

    if "MISR" in str(var.name):
        plt.axhline(y=6, linewidth=2, color="k")
        plt.axhline(y=9, linewidth=2, color="k")
        yticks = [0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 7, 9, 11, 13, 15, 17, 23]
        ylabels = ["%.1f" % i for i in yticks]
        yticks_position = np.linspace(0, 15, 16)
        plt.yticks(yticks_position, ylabels)
        ax.set_ylabel("Cloud Top Height (km)")
    else:
        plt.axhline(y=2, linewidth=2, color="k")
        plt.axhline(y=4, linewidth=2, color="k")
        yticks = [1000.0, 800.0, 680.0, 560.0, 440.0, 310.0, 180.0, 0.0]
        ylabels = ["%.1f" % i for i in yticks]
        ax.set_yticklabels(ylabels)
        ax.set_ylabel("Cloud Top Pressure (mb)")

    xticks = [0.3, 1.3, 3.6, 9.4, 23, 60, 379]
    xlabels = ["%.1f" % i for i in xticks]
    ax.set_xticklabels(xlabels)
    ax.set_xlabel("Cloud Optical Thickness")

    if title[0] is not None:
        ax.set_title(title[0], loc="left", fontdict={"fontsize": SECONDARY_TITLE_FONT})
    if title[1] is not None:
        ax.set_title(title[1], fontdict={"fontsize": MAIN_TITLE_FONT})

    ax.set_title("%", loc="right", fontdict={"fontsize": SECONDARY_TITLE_FONT})

    # Color bar
    cbax = fig.add_axes(
        (
            PANEL_CFG[subplot_num][0] + 0.6635,
            PANEL_CFG[subplot_num][1] + 0.0215,
            0.0326,
            0.1792,
        )
    )
    cbar = fig.colorbar(p1, cax=cbax, extend="both")

    if levels is None:
        cbar.ax.tick_params(labelsize=9.0, length=0)

    else:
        cbar.set_ticks(levels[1:-1])
        labels = ["%4.1f" % level for level in levels[1:-1]]
        cbar.ax.set_yticklabels(labels, ha="right")
        cbar.ax.tick_params(labelsize=9.0, pad=25, length=0)

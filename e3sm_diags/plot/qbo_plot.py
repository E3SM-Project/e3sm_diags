from __future__ import annotations

from typing import List, Literal, TypedDict

import matplotlib
import numpy as np
import xcdat as xc

from e3sm_diags.logger import custom_logger
from e3sm_diags.parameter.qbo_parameter import QboParameter
from e3sm_diags.plot.utils import _save_plot

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # isort:skip  # noqa: E402

logger = custom_logger(__name__)

PANEL_CFG = [
    (0.075, 0.70, 0.6, 0.225),
    (0.075, 0.425, 0.6, 0.225),
    (0.725, 0.425, 0.2, 0.5),
    (0.075, 0.075, 0.85, 0.275),
]

LABEL_SIZE = 14
CMAP = plt.cm.RdBu_r


class XAxis(TypedDict):
    axis_range: list[int]
    axis_scale: Literal["linear"]
    label: str
    data: np.ndarray
    data_label: str | None
    data2: np.ndarray | None
    data2_label: str | None


class YAxis(TypedDict):
    axis_range: list[int]
    axis_scale: Literal["linear", "log"]
    label: str
    data: np.ndarray
    data_label: str | None
    data2: np.ndarray | None
    data2_label: str | None


class ZAxis(TypedDict):
    data: np.ndarray


def plot(parameter: QboParameter, test_dict, ref_dict):
    fig = plt.figure(figsize=(14, 14))

    test_z_axis = xc.get_dim_coords(test_dict["qbo"], axis="Z")
    ref_z_axis = xc.get_dim_coords(ref_dict["qbo"], axis="Z")

    months = np.minimum(ref_dict["qbo"].shape[0], test_dict["qbo"].shape[0])
    x_test, y_test = np.meshgrid(np.arange(0, months), test_z_axis)
    x_ref, y_ref = np.meshgrid(np.arange(0, months), ref_z_axis)

    color_levels0 = np.arange(-50, 51, 100.0 / 20.0)

    # Panel 0 (Top Left)
    x: XAxis = dict(
        axis_range=[0, months],
        axis_scale="linear",
        label=" ",
        data=x_test,
        data_label=None,
        data2=None,
        data2_label=None,
    )
    y: YAxis = dict(
        axis_range=[100, 1],
        axis_scale="log",
        label="hPa",
        data=y_test,
        data_label=None,
        data2=None,
        data2_label=None,
    )
    z: ZAxis = dict(data=test_dict["qbo"].T[:, :months])
    title = "{} U [{}] 5S-5N ({})".format(test_dict["name"], "m/s", parameter.test_yrs)
    _add_color_map(
        0,
        fig,
        "contourf",
        title,
        x,
        y,
        z=z,
        plot_colors=CMAP,
        color_levels=color_levels0,
        color_ticks=[-50, -25, -5, 5, 25, 50],
    )

    # Panel 1 (Middle Left)
    x = dict(
        axis_range=[0, months],
        axis_scale="linear",
        label="month",
        data=x_ref,
        data_label=None,
        data2=None,
        data2_label=None,
    )
    y = dict(
        axis_range=[100, 1],
        axis_scale="log",
        label="hPa",
        data=y_ref,
        data_label=None,
        data2=None,
        data2_label=None,
    )
    z = dict(data=ref_dict["qbo"].T[:, :months])
    title = "{} U [{}] 5S-5N ({})".format(ref_dict["name"], "m/s", parameter.ref_yrs)
    _add_color_map(
        1,
        fig,
        "contourf",
        title,
        x,
        y,
        z=z,
        plot_colors=CMAP,
        color_levels=color_levels0,
        color_ticks=[-50, -25, -5, 5, 25, 50],
    )

    # Panel 2 (Top/Middle Right)
    x = dict(
        axis_range=[0, 30],
        axis_scale="linear",
        label="Amplitude (m/s)",
        data=test_dict["amplitude"][:],
        data_label=test_dict["name"],
        data2=ref_dict["amplitude"][:],
        data2_label=ref_dict["name"],
    )
    y = dict(
        axis_range=[100, 1],
        axis_scale="log",
        label="Pressure (hPa)",
        data=test_z_axis[:],
        data_label=None,
        data2=ref_z_axis[:],
        data2_label=None,
    )
    title = "QBO Amplitude \n (period = 20-40 months)"
    _add_color_map(2, fig, "line", title, x, y)

    # Panel 3 (Bottom)
    x = dict(
        axis_range=[0, 50],
        axis_scale="linear",
        label="Period (months)",
        data=test_dict["period_new"],
        data_label=test_dict["name"],
        data2=ref_dict["period_new"],
        data2_label=ref_dict["name"],
    )
    y = dict(
        axis_range=[-1, 25],
        axis_scale="linear",
        label="Amplitude (m/s)",
        data=test_dict["amplitude_new"],
        data_label=None,
        data2=ref_dict["amplitude_new"],
        data2_label=None,
    )
    title = "QBO Spectral Density (Eq. 18-22 hPa zonal winds)"
    _add_color_map(3, fig, "line", title, x, y)

    plt.tight_layout()

    # Figure title
    fig.suptitle(parameter.main_title, x=0.5, y=0.97, fontsize=15)

    # Save figure
    _save_plot(fig, parameter, PANEL_CFG)

    plt.close()


def _add_color_map(
    subplot_num: int,
    fig: plt.Figure,
    plot_type: Literal["contourf", "line"],
    title: str,
    x: XAxis,
    y: YAxis,
    z: ZAxis | None = None,
    plot_colors: plt.cm.ColormapRegistry | None = None,
    color_levels: np.ndarray | None = None,
    color_ticks: List[int] | None = None,
):
    # x,y,z should be of the form:
    # dict(axis_range=None, axis_scale=None, data=None, data_label=None, data2=None, data2_label=None, label=None)

    # Create new figure axis using dimensions from panel (hard coded)
    ax = fig.add_axes(PANEL_CFG[subplot_num])
    # Plot either a contourf or line plot
    if plot_type == "contourf":
        if z is None:
            raise RuntimeError(f"Must set `z` arg to use plot_type={plot_type}.")

        p1 = ax.contourf(
            x["data"], y["data"], z["data"], color_levels, cmap=plot_colors
        )
        cbar = plt.colorbar(p1, ticks=color_ticks)
        cbar.ax.tick_params(labelsize=LABEL_SIZE)

    if plot_type == "line":
        (p1,) = ax.plot(x["data"], y["data"], "-ok")
        (p2,) = ax.plot(x["data2"], y["data2"], "--or")

        plt.grid("on")
        ax.legend(
            (p1, p2),
            (x["data_label"], x["data2_label"]),
            loc="upper right",
            fontsize=LABEL_SIZE,
        )

    ax.set_title(title, size=LABEL_SIZE, weight="demi")
    ax.set_xlabel(x["label"], size=LABEL_SIZE)
    ax.set_ylabel(y["label"], size=LABEL_SIZE)

    plt.yscale(y["axis_scale"])
    plt.ylim([y["axis_range"][0], y["axis_range"][1]])
    plt.yticks(size=LABEL_SIZE)
    plt.xscale(x["axis_scale"])
    plt.xlim([x["axis_range"][0], x["axis_range"][1]])
    plt.xticks(size=LABEL_SIZE)

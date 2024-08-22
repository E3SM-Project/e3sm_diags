import matplotlib
import numpy as np
import scipy.stats

from e3sm_diags.logger import custom_logger
from e3sm_diags.parameter.streamflow_parameter import StreamflowParameter
from e3sm_diags.plot.utils import _save_plot

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # isort:skip  # noqa: E402

logger = custom_logger(__name__)

# Position and sizes of subplot axes in page coordinates (0 to 1)
# (left, bottom, width, height) in page coordinates
ANNUAL_SCATTER_PANEL_CFG = [(0.0900, 0.2000, 0.7200, 0.6000)]


def plot_annual_scatter(xs, ys, zs, parameter: StreamflowParameter):
    fig = plt.figure(figsize=parameter.figsize, dpi=parameter.dpi)

    ax = fig.add_axes(ANNUAL_SCATTER_PANEL_CFG[0])
    cmap = plt.get_cmap("jet")
    ax.scatter(xs, ys, label="Scatterplot", marker="o", s=10, c=zs, cmap=cmap)
    r, _ = scipy.stats.pearsonr(xs, ys)
    r2 = r * r
    r2_str = "{0:.2f}".format(r2)
    bounds = [0.01, 100000]
    ax.plot(bounds, bounds, color="red", linestyle="-")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(
        "{} streamflow ($m^3$/$s$)".format(parameter.reference_title),
        fontsize=12,
    )
    ax.set_ylabel("{} streamflow ($m^3$/$s$)".format(parameter.test_title), fontsize=12)
    ax.set_xlim(bounds[0], bounds[1])
    ax.set_ylim(bounds[0], bounds[1])
    ax.tick_params(axis="both", labelsize=12)

    # Color bar
    # Position and sizes of subplot axes in page coordinates (0 to 1)
    # (left, bottom, width, height) in page coordinates
    cbax = fig.add_axes(
        (
            ANNUAL_SCATTER_PANEL_CFG[0][0] + 0.7535,
            ANNUAL_SCATTER_PANEL_CFG[0][1] + 0.0515,
            0.0326,
            0.1792 * 2,
        )
    )
    cbar_label = "Drainage area bias (%)"
    cbar = fig.colorbar(matplotlib.cm.ScalarMappable(cmap=cmap), cax=cbax)
    cbar.ax.set_ylabel(cbar_label, fontsize=12)

    zs_max = np.ceil(np.max(zs))
    zs_min = np.floor(np.min(zs))
    step_size = (zs_max - zs_min) // 5
    try:
        ticks = np.arange(zs_min, zs_max + step_size, step_size)
        cbar.ax.set_yticklabels(ticks)
    except ValueError:
        # `zs` has invalid values (likely from no area_upstream being found).
        # Just use default colorbar.
        pass
    cbar.ax.tick_params(labelsize=12.0, length=0)

    # Figure title
    if parameter.main_title_annual_scatter == "":
        main_title_annual_scatter = "Annual mean streamflow\n{} vs {}".format(
            parameter.test_title, parameter.reference_title
        )
    else:
        main_title_annual_scatter = parameter.main_title_annual_scatter
    ax.set_title(main_title_annual_scatter, loc="center", y=1.05, fontsize=15)

    legend_title = "$R^2$={}, (n={})".format(r2_str, xs.shape[0])
    ax.legend(handles=[], title=legend_title, loc="upper left", prop={"size": 12})

    # Set the output file name before saving the plot.
    parameter.output_file = parameter.output_file_annual_scatter
    _save_plot(fig, parameter)

    plt.close()

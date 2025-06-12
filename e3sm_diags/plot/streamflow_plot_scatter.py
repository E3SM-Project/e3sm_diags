import matplotlib
import numpy as np
import scipy.stats

from e3sm_diags.logger import _setup_child_logger
from e3sm_diags.parameter.streamflow_parameter import StreamflowParameter
from e3sm_diags.plot.utils import _save_plot

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # isort:skip  # noqa: E402

logger = _setup_child_logger(__name__)

# Position and sizes of subplot axes in page coordinates (0 to 1)
# (left, bottom, width, height) in page coordinates
ANNUAL_SCATTER_PANEL_CFG = [(0.0900, 0.2000, 0.7200, 0.6000)]


def plot_annual_scatter(parameter: StreamflowParameter, export_data: np.ndarray):
    """Plot the streamflow annual scatter.

    Parameters
    ----------
    parameter : StreamflowParameter
        The streamflow parameter.
    export_data : np.ndarray
        The export data.
    """
    fig = plt.figure(figsize=parameter.figsize, dpi=parameter.dpi)

    ann_mean_ref = export_data[:, 0]
    ann_mean_test = export_data[:, 1]
    pct_drainage_area_bias = export_data[:, 2]

    # Get the figure Axes object and configure axes.
    # --------------------------------------------------------------------------
    ax = fig.add_axes(ANNUAL_SCATTER_PANEL_CFG[0])

    bounds = [0.01, 100000]
    ax.plot(bounds, bounds, color="red", linestyle="-")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(
        f"{parameter.reference_title} streamflow ($m^3$/$s$)",
        fontsize=12,
    )
    ax.set_ylabel(f"{parameter.test_title} streamflow ($m^3$/$s$)", fontsize=12)
    ax.set_xlim(bounds[0], bounds[1])
    ax.set_ylim(bounds[0], bounds[1])
    ax.tick_params(axis="both", labelsize=12)

    # Configure the title.
    # --------------------------------------------------------------------------
    if parameter.main_title_annual_scatter == "":
        main_title_annual_scatter = (
            f"Annual mean streamflow\n{parameter.test_title} vs "
            f"{parameter.reference_title}"
        )
    else:
        main_title_annual_scatter = parameter.main_title_annual_scatter

    ax.set_title(main_title_annual_scatter, loc="center", y=1.05, fontsize=15)

    # Configure the legend.
    # --------------------------------------------------------------------------
    r, _ = scipy.stats.pearsonr(ann_mean_ref, ann_mean_test)
    r2 = r * r
    r2_str = "{0:.2f}".format(r2)

    legend_title = f"$R^2$={r2_str}, (n={ann_mean_ref.shape[0]})"
    ax.legend(handles=[], title=legend_title, loc="upper left", prop={"size": 12})

    # Configure the color map.
    # --------------------------------------------------------------------------
    cmap = plt.get_cmap("jet")
    ax.scatter(
        ann_mean_ref,
        ann_mean_test,
        label="Scatterplot",
        marker="o",
        s=10,
        c=pct_drainage_area_bias,
        cmap=cmap,
    )

    # Configure the colorbar.
    # --------------------------------------------------------------------------
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
    cbar.ax.tick_params(labelsize=12.0, length=0)

    pct_drainage_area_max = np.ceil(np.max(pct_drainage_area_bias))
    pct_drainage_area_min = np.floor(np.min(pct_drainage_area_bias))
    step_size = (pct_drainage_area_max - pct_drainage_area_min) // 5

    try:
        ticks = np.arange(
            pct_drainage_area_min, pct_drainage_area_max + step_size, step_size
        )
        cbar.ax.set_yticklabels(ticks)
    except ValueError:
        # `pct_drainage_area_bias` has invalid values (likely from no area_upstream being found).
        # Just use default colorbar.
        pass

    # NOTE: Need to set the output filename to the name of the specific
    # streamflow plot before saving the plot, otherwise the filename will
    # be blank.
    parameter.output_file = parameter.output_file_annual_scatter
    _save_plot(fig, parameter)

    plt.close(fig)

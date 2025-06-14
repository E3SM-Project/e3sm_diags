import matplotlib
import numpy as np
import xarray as xr

from e3sm_diags.driver.utils.type_annotations import MetricsDict
from e3sm_diags.logger import _setup_child_logger
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.plot.lat_lon_plot import _add_colormap
from e3sm_diags.plot.utils import _save_main_plot, _save_single_subplot

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # isort:skip  # noqa: E402

logger = _setup_child_logger(__name__)

# Plot scatter plot
# Position and sizes of subplot axes in page coordinates (0 to 1)
# (left, bottom, width, height) in page coordinates
PANEL_CFG = [
    (0.09, 0.40, 0.72, 0.30),
    (0.19, 0.2, 0.62, 0.30),
]
# Border padding relative to subplot axes for saving individual panels
# (left, bottom, right, top) in page coordinates.
BORDER_PADDING_COLORMAP = (-0.06, 0.25, 0.13, 0.25)
BORDER_PADDING_SCATTER = (-0.08, -0.04, 0.15, 0.04)


def plot(
    parameter: CoreParameter,
    da_test: xr.DataArray,
    test_site_arr: np.ndarray,
    ref_site_arr: np.ndarray,
    metrics_dict: MetricsDict,
):
    """Plot the test variable's metrics generated for the aerosol_aeronet set.

    Parameters
    ----------
    parameter : CoreParameter
        The CoreParameter object containing plot configurations.
    da_test : xr.DataArray
        The test data.
    test_site : np.ndarray
        The array containing values for the test site.
    ref_site : np.ndarray
        The array containing values for the ref site.
    metrics_dict : MetricsDict
        The metrics.
    """
    fig = plt.figure(figsize=parameter.figsize, dpi=parameter.dpi)
    fig.suptitle(parameter.var_id, x=0.5, y=0.97)

    # Add the colormap subplot for test data.
    min = metrics_dict["min"]
    mean = metrics_dict["mean"]
    max = metrics_dict["max"]

    _add_colormap(
        0,
        da_test,
        fig,
        parameter,
        parameter.test_colormap,
        parameter.contour_levels,
        title=(parameter.test_name_yrs, None, None),  # type: ignore
        metrics=(max, mean, min),  # type: ignore
    )

    # Add the scatter plot.
    ax = fig.add_axes(PANEL_CFG[1])
    ax.set_title(f"{parameter.var_id} from AERONET sites")

    # Define 1:1 line, and x, y axis limits.
    if parameter.var_id == "AODVIS":
        x1 = np.arange(0.01, 3.0, 0.1)
        y1 = np.arange(0.01, 3.0, 0.1)
        plt.xlim(0.03, 1)
        plt.ylim(0.03, 1)
    else:
        x1 = np.arange(0.0001, 1.0, 0.01)
        y1 = np.arange(0.0001, 1.0, 0.01)
        plt.xlim(0.001, 0.3)
        plt.ylim(0.001, 0.3)

    plt.loglog(x1, y1, "-k", linewidth=0.5)
    plt.loglog(x1, y1 * 0.5, "--k", linewidth=0.5)
    plt.loglog(x1 * 0.5, y1, "--k", linewidth=0.5)

    corr = np.corrcoef(ref_site_arr, test_site_arr)
    xmean = np.mean(ref_site_arr)
    ymean = np.mean(test_site_arr)
    ax.text(
        0.3,
        0.9,
        f"Mean (test): {ymean:.3f} \n Mean (ref): {xmean:.3f}\n Corr: {corr[0, 1]:.2f}",
        horizontalalignment="right",
        verticalalignment="top",
        transform=ax.transAxes,
    )

    # Configure axis ticks.
    plt.tick_params(axis="both", which="major")
    plt.tick_params(axis="both", which="minor")

    # Configure axis labels
    plt.xlabel(f"ref: {parameter.ref_name_yrs}")
    plt.ylabel(f"test: {parameter.test_name_yrs}")

    plt.loglog(ref_site_arr, test_site_arr, "kx", markersize=3.0, mfc="none")

    _save_plot_aerosol_aeronet(fig, parameter)

    plt.close(fig)


def _save_plot_aerosol_aeronet(fig, parameter):
    """Save aerosol_aeronet plots using the _save_single_subplot helper function.

    This function handles the special case where different border padding is needed
    for each panel by calling _save_single_subplot twice with panel-specific
    configurations (BORDER_PADDING_COLORMAP for panel 0, BORDER_PADDING_SCATTER for panel 1).
    """
    # Save the main plot
    _save_main_plot(parameter)

    # Save subplots with different border padding by calling general function
    # for each panel individually
    if parameter.output_format_subplot:
        # Save colormap panel (panel 0) with its specific border padding
        _save_single_subplot(fig, parameter, 0, PANEL_CFG[0], BORDER_PADDING_COLORMAP)

        # Save scatter panel (panel 1) with its specific border padding
        _save_single_subplot(fig, parameter, 1, PANEL_CFG[1], BORDER_PADDING_SCATTER)

from typing import TYPE_CHECKING

import matplotlib
import numpy as np
import xarray as xr

from e3sm_diags.logger import _setup_child_logger
from e3sm_diags.parameter.area_mean_time_series_parameter import (
    AreaMeanTimeSeriesParameter,
)
from e3sm_diags.plot.utils import _save_plot

if TYPE_CHECKING:
    from e3sm_diags.driver.area_mean_time_series_driver import RefsTestMetrics

matplotlib.use("agg")
import matplotlib.pyplot as plt  # isort:skip  # noqa: E402

logger = _setup_child_logger(__name__)

# Data is a list based on the len of the regions parameter. Each element is a
# tuple (test data, ref data and metrics) for that region.
LINE_COLOR = ["r", "b", "g", "m", "c", "y"]

# Position and sizes of subplot axes in page coordinates (0 to 1). The
# dimensions [left, bottom, width, height] of the new axes. All quantities are
# in fractions of figure width and height.
PANEL_CFG = [
    (0.1, 0.70, 0.25, 0.225),
    (0.4, 0.70, 0.25, 0.225),
    (0.7, 0.70, 0.25, 0.225),
    (0.1, 0.38, 0.25, 0.225),
    (0.4, 0.38, 0.25, 0.225),
    (0.7, 0.38, 0.25, 0.225),
    (0.1, 0.06, 0.25, 0.225),
    (0.4, 0.06, 0.25, 0.225),
    (0.7, 0.06, 0.25, 0.225),
]

# Border padding relative to subplot axes for saving individual panels.
# (left, bottom, right, top) in page coordinates.
BORDER_PADDING = (-0.047, -0.06, 0.006, 0.03)


def plot(
    var: str,
    parameter: AreaMeanTimeSeriesParameter,
    metrics_dict: dict[str, RefsTestMetrics],
):
    # Create the figure.
    fig = plt.figure(figsize=(17.0, 10.0), dpi=parameter.dpi)

    test_name_yrs = parameter.test_name_yrs
    test_name = test_name_yrs.split("(")[0].replace(" ", "")
    if test_name == "":
        test_name = "test data"

    for idx_region, ds_region in enumerate(metrics_dict.values()):
        # Add the figure axis object by region.
        ax = fig.add_axes(PANEL_CFG[idx_region])

        # Plot the test data.
        test_var = ds_region.test[var]
        test_label = _get_mean_and_std_label(test_var, test_name)
        ax.plot(test_var, "k", linewidth=2, label=test_label)

        # Plot the list of reference data (if they exist).
        refs = ds_region.refs
        for idx_ref, ref_var in enumerate(refs):
            ref_label = _get_mean_and_std_label(ref_var, ref_var.attrs["ref_name"])
            ax.plot(ref_var, LINE_COLOR[idx_ref], linewidth=2, label=ref_label)

        # Perform truncation division to accomodate long time records to get
        # the step sizes for the X axis ticks.
        start_time = int(parameter.start_yr)  # type: ignore
        end_time = int(parameter.end_yr)  # type: ignore
        num_year = end_time - start_time + 1

        if num_year > 19:
            stepsize = num_year // 10
        else:
            stepsize = 1

        # Configure X and Y axes.
        x = np.arange(num_year)
        ax.set_xticks(x[::stepsize])
        x_ticks_labels = np.arange(start_time, end_time + 1, stepsize)
        ax.set_xticklabels(x_ticks_labels, rotation=45, fontsize=8)
        ax.set_xlabel("Year")

        ax.set_ylabel(var + " (" + test_var.units + ")")

        # Configure legend and title.
        ax.legend(loc="best", prop={"size": 7})
        ax.set_title(parameter.regions[idx_region], fontsize=10)

    # Figure title.
    fig.suptitle(
        "Annual mean " + var + " time series over regions " + parameter.test_name_yrs,
        x=0.3,
        y=0.98,
        fontsize=15,
    )

    parameter.output_file = var
    _save_plot(fig, parameter, PANEL_CFG, BORDER_PADDING)

    plt.close(fig)


def _get_mean_and_std_label(var: xr.DataArray, name: str) -> str:
    """Get the plot label using the mean and std deviation of the yearly mean.

    Parameters
    ----------
    var : xr.DataArray
        The test or reference variable.
    name : str
        The test or reference label name.

    Returns
    -------
    str
        The plot label.
    """
    mean = np.mean(var)
    std = np.std(var)
    label = name + f"(mean: {mean:.2f}, std: {std:.3f})"

    return label

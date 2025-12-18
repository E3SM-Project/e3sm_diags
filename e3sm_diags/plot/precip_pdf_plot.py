"""
Plotting functions for precipitation PDF diagnostics.

Creates frequency and amount PDF plots for precipitation data.

Based on original work by Chris Terai.
Modified to integrate into E3SM Diags.
"""
from __future__ import annotations

import os

import matplotlib
import xarray as xr

from e3sm_diags.driver.utils.io import _get_output_dir
from e3sm_diags.logger import _setup_child_logger
from e3sm_diags.parameter.core_parameter import CoreParameter

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # isort:skip  # noqa: E402

logger = _setup_child_logger(__name__)

# Plot title and side title configurations
PLOT_TITLE = {"fontsize": 13}
PLOT_SIDE_TITLE = {"fontsize": 11}


def plot(
    parameter: CoreParameter,
    test_pdf: xr.Dataset,
    ref_pdf: xr.Dataset,
    region: str,
):
    """Plot precipitation PDFs for test and reference data.

    Parameters
    ----------
    parameter : CoreParameter
        The CoreParameter object containing plot configurations
    test_pdf : xr.Dataset
        Test data PDF
    ref_pdf : xr.Dataset
        Reference data PDF
    region : str
        Region name for plot title
    """
    # Create figure with two subplots (frequency and amount PDFs)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12))

    bin_centers = test_pdf["bin_centers"].values

    # Plot 1: Frequency PDF (df/dlog(P))
    ax1.plot(
        bin_centers,
        test_pdf["FREQPDF"].values,
        color="tab:blue",
        linewidth=2,
        label=parameter.test_name_yrs,
    )
    ax1.plot(
        bin_centers,
        ref_pdf["FREQPDF"].values,
        color="tab:brown",
        linewidth=2,
        label=parameter.ref_name_yrs,
    )

    ax1.set_ylabel("df/dlog(P)", fontsize=14)
    ax1.set_xlabel("P (mm/day)", fontsize=14)
    ax1.set_xscale("log")
    ax1.tick_params(axis="both", labelsize=12)
    ax1.legend(loc="upper right", fontsize=12)
    ax1.grid(True, alpha=0.3)
    ax1.set_title(f"{region} Precipitation Frequency PDF (daily)", fontsize=13)

    # Plot 2: Amount PDF (dA/dlog(P))
    # This shows the contribution to total precipitation amount from each bin
    ax2.plot(
        bin_centers,
        test_pdf["AMNTPDF"].values,
        color="tab:blue",
        linewidth=2,
        label=parameter.test_name_yrs,
    )
    ax2.plot(
        bin_centers,
        ref_pdf["AMNTPDF"].values,
        color="tab:brown",
        linewidth=2,
        label=parameter.ref_name_yrs,
    )

    ax2.set_ylabel("dA/dlog(P)", fontsize=14)
    ax2.set_xlabel("P (mm/day)", fontsize=14)
    ax2.set_xscale("log")
    ax2.tick_params(axis="both", labelsize=12)
    ax2.legend(loc="upper right", fontsize=12)
    ax2.grid(True, alpha=0.3)
    ax2.set_title(f"{region} Precipitation Amount PDF (daily)", fontsize=13)

    plt.tight_layout()

    # Save the plot
    _save_plot(fig, parameter)

    plt.close()


def _save_plot(fig: plt.figure, parameter: CoreParameter):
    """Save the plot using the figure object and parameter configs.

    Parameters
    ----------
    fig : plt.figure
        The plot figure
    parameter : CoreParameter
        The CoreParameter with file configurations
    """
    for f in parameter.output_format:
        f = f.lower().split(".")[-1]
        fnm = os.path.join(
            _get_output_dir(parameter),
            parameter.output_file + "." + f,
        )
        plt.savefig(fnm, dpi=150, bbox_inches="tight")
        logger.info(f"Plot saved in: {fnm}")

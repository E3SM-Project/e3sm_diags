from __future__ import print_function

import os

import matplotlib
import numpy as np

from e3sm_diags.driver.utils.general import get_output_dir
from e3sm_diags.logger import custom_logger

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # isort:skip  # noqa: E402

logger = custom_logger(__name__)

plotTitle = {"fontsize": 11.5}
plotSideTitle = {"fontsize": 9.5}


def plot(test_site, ref_site, parameter):
    # Plot scatter plot
    # Set global matplotlib parameters
    plt.clf()
    plt.rcParams["axes.linewidth"] = 0.5

    # Create figure, projection
    panel = [(0.0900, 0.2000, 0.7200, 0.6000)]
    fig = plt.figure(figsize=(2.5, 2), dpi=300)

    # Figure title
    title = f"Comparison of {parameter.var_id} from AERONET sites"
    fig.suptitle(title, x=0.5, y=0.97, fontsize=8)
    # ax = plt.axes()
    ax = fig.add_axes(panel[0])

    # define 1:1 line

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

    # add text
    corr = np.corrcoef(ref_site, test_site)
    xmean = np.mean(ref_site)
    ymean = np.mean(test_site)
    ax.text(0.02, 0.98, "{} {:.3f}".format("mean(Model): ", ymean), fontsize=5)
    ax.text(0.02, 0.9, "{} {:.3f}".format("mean(Reference): ", xmean), fontsize=5)
    ax.text(0.02, 0.82, "{} {:.2f}".format("corr: ", corr[0, 1]), fontsize=5)

    # axis ticks
    plt.tick_params(axis="both", which="major", labelsize=5, length=1.6, width=0.5)
    plt.tick_params(axis="both", which="minor", labelsize=5, length=1.0, width=0.5)

    # axis labels
    plt.xlabel(parameter.ref_name_yrs, fontsize=5)
    plt.ylabel(parameter.test_name_yrs, fontsize=5)

    plt.loglog(ref_site, test_site, "kx", markersize=3.0, mfc="none", label="Z03")

    # legend
    plt.legend(frameon=False, prop={"size": 5})

    for f in parameter.output_format:
        f = f.lower().split(".")[-1]
        fnm = os.path.join(
            get_output_dir(parameter.current_set, parameter),
            f"{parameter.output_file}" + "." + f,
        )
        plt.savefig(fnm, bbox_inches="tight", dpi=300)
        logger.info(f"Plot saved in: {fnm}")

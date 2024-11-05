from __future__ import annotations

import math
import os
import warnings
from typing import TYPE_CHECKING, List, TypedDict

import matplotlib
import numpy as np
from matplotlib.gridspec import GridSpec

from e3sm_diags.driver.utils.diurnal_cycle_xr import _fft_all_grid
from e3sm_diags.driver.utils.io import _get_output_dir
from e3sm_diags.logger import custom_logger
from e3sm_diags.parameter.arm_diags_parameter import ARMDiagsParameter

matplotlib.use("agg")

import matplotlib.pyplot as plt  # isort:skip  # noqa: E402

logger = custom_logger(__name__)

if TYPE_CHECKING:
    from e3sm_diags.driver.arm_diags_driver import RefsTestMetrics


MONTHS = [
    "Jan",
    "Feb",
    "Mar",
    "Apr",
    "May",
    "Jun",
    "Jul",
    "Aug",
    "Sep",
    "Oct",
    "Nov",
    "Dec",
]


# Region information for convection onset statistics.
class Stats(TypedDict):
    cwv_max: int
    cwv_min: int
    bin_width: float
    sitename: str


class RegionStats(TypedDict):
    twpc1: Stats
    twpc2: Stats
    twpc3: Stats
    sgpc1: Stats


REGION_INFO: RegionStats = {
    "twpc1": {
        "cwv_max": 69,
        "cwv_min": 28,
        "bin_width": 1.5,
        "sitename": "Manus Island",
    },
    "twpc2": {"cwv_max": 70, "cwv_min": 28, "bin_width": 2.0, "sitename": "Nauru"},
    "twpc3": {"cwv_max": 85, "cwv_min": 28, "bin_width": 2.0, "sitename": "Darwin"},
    "sgpc1": {"cwv_max": 75, "cwv_min": 20, "bin_width": 2.0, "sitename": "SGP"},
}

# Precipitation threshold for convection onset, default 0.5 (in mm/hr).
PRECIP_THRESHOLD = 0.5


def _plot_diurnal_cycle(parameter: ARMDiagsParameter, vars_to_data: RefsTestMetrics):
    """Plot the diurnal cycle of Total Precipitation Rate.

    Parameters:
    -----------
    parameter : ARMDiagsParameter
        The ARMDiagsParameter object containing the parameters for plotting.
    vars_to_data : RefsTestMetrics
        The ordered dictionary containing the variables and their corresponding
        data.
    """
    test = vars_to_data.test[0]
    ref = vars_to_data.refs[0][0]
    lst = vars_to_data.misc[0]
    t_conv = lst[0][0]

    output_file_name = parameter.output_file + "-" + "diurnal-cycle"

    fig = plt.figure()
    ax = fig.add_axes([0.15, 0.1, 0.8, 0.8])

    for index in range(2):
        if index == 0:
            data = test
            line_c = "k"
            data_name = parameter.test_name_yrs
        else:
            data = ref
            line_c = "r"
            data_name = parameter.ref_name

        time_freq = len(data)
        res = int(24 / time_freq)
        c, maxvalue, tmax = _fft_all_grid(data, np.array([0]))

        # Configure x and y axes.
        # ----------------------------------------------------------------------
        x_axis1 = np.linspace(0, 48 - res, time_freq * 2)
        ax.plot(x_axis1, np.concatenate((data, data)), "." + line_c, label=data_name)

        x_axis2 = np.linspace(0, 48 - res, time_freq * 2 * 3)
        w = 2.0 * np.pi / 24
        yax = (c + maxvalue[0] * np.sin(w * x_axis2 + np.pi / 2 - tmax[0] * w))[0]
        ax.plot(x_axis2, yax, line_c, label="First harmonic")

        plt.xlim([24 - t_conv, 47 - t_conv + 1])
        plt.ylim([0, 5.5])
        plt.xlabel("local solar time [hr]")
        plt.ylabel("Total Precipitation Rate" + " (" + parameter.var_units + ")")  # type: ignore

        x_axis3 = np.arange(24 - t_conv, 47 - t_conv, 3)
        x_ticks_hrs = ["0h", "3h", "6h", "9h", "12h", "15h", "18h", "21h"]
        plt.xticks(x_axis3, x_ticks_hrs)

        # Configure legend and title.
        # ----------------------------------------------------------------------
        plt.legend(loc="upper right")
        plt.title(output_file_name.replace("-", " "))

    _save_plots(parameter, output_file_name, parameter.output_format)
    plt.close()


def _plot_diurnal_cycle_zt(parameter: ARMDiagsParameter, vars_to_data: RefsTestMetrics):
    """Plot the diurnal cycle of cloud fraction for each month.

    Parameters:
    -----------
    parameter : ARMDiagsParameter
        The ARMDiagsParameter object containing the parameters for plotting.
    vars_to_data : RefsTestMetrics
        The ordered dictionary containing the variables and their corresponding
        data.
    """
    ref = vars_to_data.refs[0]
    test = vars_to_data.test
    lst = vars_to_data.misc

    for index in range(2):
        fig, axs = plt.subplots(
            4,
            3,
            figsize=(15, 12),
            facecolor="w",
            edgecolor="k",
            sharex=True,
            sharey=True,
        )
        fig.subplots_adjust(hspace=0.4, wspace=0.1)

        axs = axs.ravel()
        t_conv = lst[0][0][0]

        for imon in range(12):
            if index == 0:
                title = parameter.ref_name
                data = ref
                data_name = "ref"
                if "armdiags" in title:
                    data = data[:, :, ::-1]

            else:
                title = parameter.test_name_yrs
                data = test
                data_name = "test"

            time_freq = data.shape[1]

            # Configure x and y axes.
            # ------------------------------------------------------------------
            yy = np.linspace(0, 48, time_freq * 2)
            xx = np.linspace(100, 1000, 37)
            x, y = np.meshgrid(xx, yy)
            data_con = np.concatenate((data[imon, :, :], data[imon, :, :]), axis=0)
            im = axs[imon].pcolormesh(
                y, x, data_con[:, :], vmin=0, vmax=30, cmap="jet", shading="auto"
            )
            axs[imon].set_title(MONTHS[imon])
            plt.xlim([24 - t_conv, 47 - t_conv])
            xax = np.arange(24 - t_conv, 47 - t_conv, 3)

            my_xticks = ["0", "3", "6", "9", "12", "15", "18", "21"]
            plt.xticks(xax, my_xticks)
            axs[imon].xaxis.set_tick_params(labelbottom=True)

            axs[imon].set_xlabel("Local time (hr)")

        # Configure y label.
        # ----------------------------------------------------------------------)
        for ax in axs[::3]:
            ax.set_ylabel("Pressure (mb)")
        axs[0].invert_yaxis()

        # Configure titles.
        # ----------------------------------------------------------------------)
        site = parameter.output_file.split("-")[-1]
        suptitle = "Cloud Fraction Monthly Diurnal Cycle " + site + "\n" + title

        plt.suptitle(suptitle, fontsize=20)
        fig.subplots_adjust(right=0.8)

        # Configure colorbar.
        # ----------------------------------------------------------------------)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
        fig.colorbar(im, cax=cbar_ax)
        plt.title("cl (%)")

        output_file_name = parameter.output_file + "-" + data_name
        _save_plots(parameter, output_file_name, parameter.output_format)

        plt.close()


def _plot_convection_onset_statistics(
    parameter: ARMDiagsParameter,
    region: str,
    test_pr: np.ndarray,
    test_prw: np.ndarray,
    ref_pr: np.ndarray,
    ref_prw: np.ndarray,
):
    """Plot the convection onset statistics.

    Parameters
    ----------
    parameter : ARMDiagsParameter
        The ARMDiagsParameter object containing the parameters for the plot.
    region : str
        The region for which the plot is generated.
    test_pr : np.ndarray
        The test precipitation data.
    test_prw : np.ndarray
        The test precipitable water data.
    ref_pr : np.ndarray
        The reference precipitation data.
    ref_prw : np.ndarray
        The reference precipitable water data.

    Notes
    -----
    - Original code: Kathleen Schiro, python version 22 Dec 2016, University of
        California Dept. of Atmospheric and Oceanic Sciences
    - Modifications: Baird Langenbrunner, Yi-Hung Kuo
    - Modifications: Jill Zhang, Cheng Tao
    - Scientific supervision: Prof. J David Neelin

    For related publications and research information see the Neelin group
    webpage  http://www.atmos.ucla.edu/~csi/.
    """
    region_info: Stats | None = REGION_INFO.get(region)  # type: ignore
    if region_info is None:
        raise ValueError(f"Invalid region: {region}")

    cwv_max = region_info["cwv_max"]
    cwv_min = region_info["cwv_min"]
    bin_width = region_info["bin_width"]
    site_name = region_info["sitename"]

    fig, axes = plt.subplots(1, 3, figsize=(12, 3))
    fig.subplots_adjust(wspace=0.3)
    title = ""

    for index in range(2):
        if index == 0:
            precip = test_pr
            cwv = test_prw
            data_name = "Test: " + parameter.test_name_yrs
            line_color = ["black", "grey"]
            time_interval = parameter.time_interval
        else:
            precip = ref_pr
            cwv = ref_prw
            data_name = "Ref: " + parameter.ref_name
            line_color = ["blue", "steelblue"]
            time_interval = 1

        number_of_bins = int(np.ceil((cwv_max - cwv_min) / bin_width))
        bin_center = np.arange(
            (cwv_min + (bin_width / 2)),
            (cwv_max - (bin_width / 2)) + bin_width,
            bin_width,
        )
        if len(bin_center) != number_of_bins:
            bin_center = np.arange(
                (cwv_min + (bin_width / 2)), (cwv_max - (bin_width / 2)), bin_width
            )

        # Define variables for binning
        bin_index = np.zeros([number_of_bins, cwv.size])
        precip_binned = np.empty([number_of_bins, cwv.size]) * np.nan
        precip_counts = np.zeros([number_of_bins, cwv.size])

        # FIXME: Why are we ignoring warnings here?
        warnings.filterwarnings("ignore")
        # Bin the data by CWV value as specified above
        for i in range(0, number_of_bins):
            tmp1 = np.where(cwv > cwv_min + (i * bin_width))
            bin_index[i, tmp1] = 1
            tmp2 = np.where(cwv > cwv_min + (i * bin_width) + bin_width)
            bin_index[i, tmp2] = 0

        for i in range(0, number_of_bins):
            tmp1 = np.where(bin_index[i, :] == 1)
            precip_binned[i, tmp1] = precip[tmp1]
            tmp2 = np.where(bin_index[i, :] != 1)
            precip_binned[i, tmp2] = np.nan

        for i in range(0, number_of_bins):
            tmp1 = np.where(precip_binned[i, :] >= PRECIP_THRESHOLD)
            precip_counts[i, tmp1] = 1
            for j in range(0, cwv.size):
                if np.isnan(precip_binned[i, j]):
                    precip_counts[i, j] = np.nan

        # Create binned arrays
        hist_cwv = np.empty([number_of_bins, 1]) * np.nan
        hist_precip_points = np.empty([number_of_bins, 1]) * np.nan

        pr_binned_mean = np.empty([number_of_bins, 1]) * np.nan
        pr_binned_std = np.empty([number_of_bins, 1]) * np.nan
        pr_probability = np.empty([number_of_bins, 1]) * np.nan

        errorbar_precip_points = np.empty([number_of_bins, 1]) * np.nan
        errorbar_precip = np.empty([number_of_bins, 1]) * np.nan
        errorbar_precip_binom = np.empty([number_of_bins, 2]) * np.nan

        # Fill binned arrays
        hist_cwv = bin_index.sum(axis=1)
        hist_cwv[hist_cwv <= 1] = 0

        hist_precip_points = np.nansum(precip_counts, axis=1)
        hist_precip_points[hist_precip_points <= 1] = 0

        pr_binned_mean = np.nanmean(precip_binned, axis=1)
        pr_binned_std = np.nanstd(precip_binned, axis=1)

        r = np.empty([1, number_of_bins]) * np.nan
        r = np.sum(~np.isnan(precip_counts), axis=1)
        pr_probability = np.nansum(precip_counts, axis=1) / r

        freq_cwv = (hist_cwv / bin_width) / np.nansum(hist_cwv)
        freq_precipitating_points = hist_precip_points / bin_width / np.nansum(hist_cwv)

        for i in range(0, number_of_bins):
            errorbar_precip[i] = pr_binned_std[i] / math.sqrt(hist_cwv[i])
            errorbar_precip_points[i] = (
                math.sqrt(hist_precip_points[i])
                / np.nansum(hist_cwv / bin_width)
                / bin_width
            )
            z = 0.675

            phat = hist_precip_points[i] / hist_cwv[i]

            errorbar_precip_binom[i, 0] = z * math.sqrt(phat * (1 - phat) / hist_cwv[i])
            errorbar_precip_binom[i, 1] = z * math.sqrt(phat * (1 - phat) / hist_cwv[i])

        # General plot configurations.
        # ----------------------------------------------------------------------
        axes_fontsize = 12  # size of font in all plots
        marker_size = 40  # size of markers in scatter plots
        xtick_pad = 10  # padding between x tick labels and actual plot
        bin_width = (np.max(bin_center) - np.min(bin_center)) / number_of_bins

        # Figure 1.
        # --------------------------------------------------------------------------
        ax1 = axes[0]
        xulim = 5 * np.ceil(np.max(np.round(bin_center + bin_width / 2)) / 5)
        xllim = 5 * np.floor(np.min(np.round(bin_center - bin_width / 2)) / 5)

        ax1.tick_params(labelsize=axes_fontsize)
        ax1.tick_params(axis="x", pad=10)
        ax1.errorbar(
            bin_center,
            pr_binned_mean,
            yerr=errorbar_precip.squeeze(),
            ls="none",
            color="black",
        )

        ax1.scatter(
            bin_center,
            pr_binned_mean,
            edgecolor="none",
            facecolor=line_color[0],
            s=marker_size,
            clip_on=True,
            zorder=3,
            label=data_name.split(":")[0],
        )
        ax1.set_xlim(xllim - 10, cwv_max)
        ax1.set_ylim(0, 3)
        ax1.set_ylabel("Precip (mm/hr)", fontsize=axes_fontsize)
        ax1.set_xlabel("CWV (mm)", fontsize=axes_fontsize)
        ax1.set_axisbelow(True)
        legend_handles, legend_labels = ax1.get_legend_handles_labels()
        ax1.legend(legend_handles, legend_labels, loc="upper left", frameon=False)

        # Figure 2 (probability pickup)
        # ----------------------------------------------------------------------
        ax2 = axes[1]
        xulim = 5 * np.ceil(np.max(np.round(bin_center + bin_width / 2)) / 5)
        xllim = 5 * np.floor(np.min(np.round(bin_center - bin_width / 2)) / 5)
        # ax2.set_xlim(xllim-10,xulim+15)
        ax2.tick_params(labelsize=axes_fontsize)
        ax2.errorbar(
            bin_center,
            pr_probability,
            yerr=errorbar_precip_binom.T,
            fmt="none",
            color="black",
        )
        ax2.tick_params(axis="x", pad=xtick_pad)
        ax2.scatter(
            bin_center,
            pr_probability,
            marker="d",
            s=marker_size,
            edgecolor="none",
            facecolor=line_color[0],
            zorder=3,
            label=data_name.split(":")[0],
        )
        ax2.set_xlim(xllim - 10, cwv_max)
        ax2.set_ylim(0, 1)
        ax2.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
        ax2.set_ylabel("Probability of Precip.", fontsize=axes_fontsize)
        ax2.set_xlabel("CWV (mm)", fontsize=axes_fontsize)
        ax2.set_axisbelow(True)

        legend_handles, legend_labels = ax2.get_legend_handles_labels()
        ax2.legend(legend_handles, legend_labels, loc="upper left", frameon=False)
        title = (
            title
            + data_name
            + ":  "
            + str(time_interval)
            + " hrly("
            + line_color[0]
            + ")\n"
        )

        # Figure 3 (non-normalized PDF).
        # ----------------------------------------------------------------------
        ax3 = axes[2]
        ax3.set_yscale("log")

        xulim = 5 * np.ceil(np.max(np.round(bin_center + bin_width / 2)) / 5)
        xllim = 5 * np.floor(np.min(np.round(bin_center - bin_width / 2)) / 5)
        # ax3.set_xlim(xllim-10,xulim+15)
        ax3.set_xlim(xllim - 10, cwv_max)
        ax3.set_xticks(
            np.arange(np.ceil(xllim / 10) * 10 - 10, np.ceil(xulim / 10) * 10 + 15, 15)
        )
        # low_lim = -6.0
        low_lim = -4.0
        ax3.set_ylim(10**low_lim, 100)
        ax3.set_yticks(10 ** np.arange(low_lim, 2, dtype="float64"))
        ax3.tick_params(labelsize=axes_fontsize)
        ax3.tick_params(axis="x", pad=xtick_pad)
        freq_precipitating_points[freq_precipitating_points == 0] = np.nan
        freq_cwv[freq_cwv == 0] = np.nan

        ax3.scatter(
            bin_center,
            freq_cwv,
            color=line_color[0],
            label=data_name.split(":")[0] + ": all",
        )
        ax3.scatter(
            bin_center,
            freq_precipitating_points,
            edgecolor="none",
            facecolor=line_color[1],
            s=marker_size,
            zorder=3,
            label=data_name.split(":")[0] + ": precip $>$ 0.5 mm/hr ",
        )
        ax3.set_ylabel("PDF", fontsize=axes_fontsize)
        ax3.set_xlabel("CWV (mm)", fontsize=axes_fontsize)
        ax3.set_axisbelow(True)

        # create legend
        legend_handles, legend_labels = ax3.get_legend_handles_labels()
        ax3.legend(
            legend_handles,
            legend_labels,
            loc="upper left",
            bbox_to_anchor=(0.1, 0.95),
            fontsize=9,
            scatterpoints=1,
            handlelength=0,
            labelspacing=0,
            borderpad=0,
            borderaxespad=0,
            frameon=False,
        )

    plt.suptitle(
        "Convection Onset Metrics" + " at " + site_name, y=1.15, fontweight="bold"
    )
    plt.title(title, ha="left", x=-2, y=0.98)

    # Save the figure.
    output_file_name = parameter.output_file
    _save_plots(
        parameter,
        output_file_name,
        parameter.output_format,
        transparent=True,
        bbox_inches="tight",
    )

    plt.close()


def _get_seasonal_mean(data: np.ndarray) -> np.ndarray:
    """Calculate annual mean and seasonal mean of input data (mean of 12 month)

    TODO: Use climo_xr to get weighted seasonal mean

    Parameters:
    -----------
    data : np.ndarray
        Input data array.

    Returns:
    --------
    np.ndarray
        Array containing the annual mean and seasonal means.
    """
    ac = data
    ac_dec = np.concatenate((ac, ac))[11:23]
    season = np.nanmean(ac_dec.reshape(-1, 3), axis=1)
    ann = np.nanmean(season)

    return np.hstack([ann, season.data])


def _plot_annual_cycle(
    parameter: ARMDiagsParameter, var: str, vars_to_data: RefsTestMetrics
):
    """Plot the annual cycle of a variable.

    Parameters:
    -----------
    parameter : ARMDiagsParameter
        The parameter object containing information about the plot.
    var : str
        The variable to plot.
    vars_to_data : RefsTestMetrics
        A dictionary mapping variable names to data.
    """
    line_color = ["r", "b", "g", "m"]
    fig = plt.figure()

    ax1 = fig.add_axes([0.15, 0.1, 0.8, 0.8])  # Create axes
    xax = np.arange(1, 13, 1)

    refs = vars_to_data.refs
    test = vars_to_data.test
    ax1.plot(xax, test, "k", linewidth=2, label="Test: " + parameter.test_name_yrs)

    test_season = _get_seasonal_mean(test)

    for i_ref, ref in enumerate(refs):
        ref_season = _get_seasonal_mean(ref)
        ax1.plot(
            xax,
            ref,
            line_color[i_ref],
            linewidth=2,
            label="Ref: " + parameter.ref_name,
        )

    my_xticks = ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"]

    plt.xticks(xax, my_xticks)
    plt.xlim(1, 12)
    ymin, ymax = plt.gca().get_ylim()
    plt.ylim(0.8 * ymin, 1.2 * ymax)

    plt.xlabel("Month")
    plt.legend(loc="best", prop={"size": 10})

    if var == "PRECT":
        plt.ylabel("Total Precipitation Rate" + " (" + parameter.var_units + ")")  # type: ignore
    else:
        plt.ylabel(parameter.var_name + " (" + parameter.var_units + ")")  # type: ignore

    # Add a table at the bottom of the axes
    bias = test_season - ref_season
    cell_text = np.round(np.vstack((test_season, ref_season, bias)), 2)
    collabel = ("ANN", "DJF", "MAM", "JJA", "SON")
    rows = ("Test", "Ref", "Bias")
    plt.table(
        cellText=cell_text,
        rowLabels=rows,
        colLabels=collabel,
        alpha=0.8,
        bbox=[0.15, 0.0, 0.8, 0.2],
    )

    output_file_name = parameter.output_file
    plt.title(output_file_name.replace("-", " "))
    _save_plots(parameter, output_file_name, parameter.output_format)

    plt.close()


def _plot_aerosol_activation(
    parameter: ARMDiagsParameter,
    region: str,
    variable: str,
    test_a_num: np.ndarray,
    test_ccn_num: np.ndarray,
    ref_a_num: np.ndarray,
    ref_ccn_num: np.ndarray,
):
    """Plot the aerosol activation.

    Parameters
    ----------
    parameter : ARMDiagsParameter
        The ARMDiagsParameter object containing the parameters for the plot.
    region : str
        The region for which the plot is generated.
    variable : str
        The variable for which the plot is generated.
    test_a_num : np.ndarray
        The test aerosol number concentration data.
    test_ccn_num : np.ndarray
        The test CCN number concentration data.
    ref_a_num : np.ndarray
        The reference aerosol number concentration data.
    ref_ccn_num : np.ndarray
        The reference CCN number concentration data.

    Notes
    -----
    Program for generate aerosol-to-ccn activate metric
    - Original code: Xiaojian Zheng, Cheng Tao
    - Modifications: Jill Zhang

    For related publications and research information see
    https://github.com/ARM-DOE/arm-gcm-diagnostics/blob/master/docs/ARM_DIAGS_v3_TechReport.pdf #
    """
    _subplot_aerosol_ccn(
        parameter,
        region,
        variable,
        test_a_num,
        test_ccn_num,
    )
    _subplot_aerosol_ccn(
        parameter, region, variable, ref_a_num, ref_ccn_num, test=False
    )
    return


def _subplot_aerosol_ccn(
    parameter: ARMDiagsParameter,
    region: str,
    variable: str,
    a_num: np.ndarray,
    ccn_num: np.ndarray,
    test: bool = True,
):
    """
    Plot the aerosol activation.

    Parameters
    ----------
    parameter : ARMDiagsParameter
        The ARMDiagsParameter object containing the parameters for the plot.
    region : str
        The region for which the plot is generated.
    variable : str
        The variable for which the plot is generated.
    a_num : np.ndarray
        The aerosol number concentration data.
    ccn_num : np.ndarray
        The CCN number concentration data.
    test : bool, optional
        Whether the data is from the test model or the reference model, by
        default True.
    """
    # Bulk aerosol vs. ccn
    if region == "sgpc1":
        ccn_num_pedge = np.arange(0, 6200, 100)
        a_num_pedge = np.arange(0, 6200, 100)
        pvmax = 6000
    elif region == "enac1":
        ccn_num_pedge = np.arange(0, 1020, 20)
        a_num_pedge = np.arange(0, 1020, 20)
        pvmax = 1000
    else:
        msg = "Aerosol activation at Site: {} is not supported yet".format(region)
        raise RuntimeError(msg)

    a_num = np.array(a_num)
    ccn_num = np.array(ccn_num)

    ratio_all = ccn_num / a_num
    ratio_mean = np.nanmean(ratio_all)
    ratio_std = np.nanstd(ratio_all)

    output_str = "test"
    if test is False:
        output_str = "ref"

    if parameter.ref_name == "armdiags" and test is False:
        data_name = "OBS"
    else:
        data_name = "Test Model"

    fig = plt.figure(figsize=(12, 10))

    fsize = 30
    xysize = 30
    gspec = GridSpec(ncols=1, nrows=1, figure=fig)
    ax1 = fig.add_subplot(gspec[0])
    ax1.set_title(
        f"{region.upper()} Bulk Aerosol Activation ({data_name})", fontsize=fsize
    )
    h2d02, xeg02, yeg02, im02 = plt.hist2d(
        a_num, ccn_num, bins=[a_num_pedge, ccn_num_pedge], cmap="turbo", density=True
    )
    ax1.plot([0, pvmax], [0, pvmax], "r", lw=3)
    ax1.text(
        0.02,
        0.9,
        "Ratio = " + "%.2f" % ratio_mean + r"$\pm$" + "%.2f" % ratio_std,
        color="r",
        ha="left",
        va="center",
        transform=ax1.transAxes,
        fontsize=xysize,
    )
    ax1.set_xlabel("Aerosol Num. Conc. (# $cm^{-3}$)", fontsize=xysize)
    ax1.set_ylabel(
        f"CCN Num. Conc. @0.{variable[-1]}%SS (# $cm^{-3}$)", fontsize=xysize
    )
    ax1.tick_params(
        labelsize=xysize, length=10, width=2, direction="out", which="major"
    )
    ax1.tick_params(length=7, width=3, direction="out", which="minor")
    cb1 = plt.colorbar()
    cb1.ax.tick_params(labelsize=13)
    cb1.set_label("Probability Density", fontsize=15)
    for axis in ["top", "bottom", "left", "right"]:
        ax1.spines[axis].set_linewidth(2)
    plt.subplots_adjust(left=0.16, right=1.01, bottom=0.11, top=0.94, hspace=0.15)

    output_file_name = f"{parameter.output_file}-{output_str}"
    _save_plots(
        parameter,
        output_file_name,
        parameter.output_format,
        transparent=True,
        bbox_inches="tight",
    )

    plt.close()


def _save_plots(
    parameter: ARMDiagsParameter,
    output_file_name: str,
    output_format: List[str],
    transparent: bool = False,
    bbox_inches: str | None = None,
):
    """
    Save the generated plots in the specified output formats.

    Parameters
    ----------
    parameter : ARMDiagsParameter
        The ARMDiagsParameter object containing the parameters for the plot.
    output_file_name : str
        The name of the output file.
    output_format : List[str]
        The list of output formats to save the plots in.
    transparent : bool, optional
        Whether the plot should be transparent, by default False.
    bbox_inches : str, optional
        The bounding box in inches, by default None.
    """
    for f in output_format:
        f = f.lower().split(".")[-1]
        fnm = os.path.join(_get_output_dir(parameter), output_file_name + "." + f)
        plt.savefig(fnm, transparent=transparent, bbox_inches=bbox_inches)

        fnm = os.path.join(
            _get_output_dir(parameter),
            output_file_name + "." + f,
        )
        logger.info(f"Plot saved in: {fnm}")

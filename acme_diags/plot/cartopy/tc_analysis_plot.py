import os

import cartopy.crs as ccrs
import matplotlib
import numpy as np

from acme_diags.driver.utils.general import get_output_dir

matplotlib.use("agg")
import matplotlib.pyplot as plt  # isort:skip  # noqa: E402


def plot(test, ref, parameter, basin_dict):
    test_metrics = test["metrics"]
    ref_metrics = ref["metrics"]
    fig, axes = plt.subplots(2, 3, figsize=(12, 6), sharex=True, sharey=True)
    fig.subplots_adjust(hspace=0.4, wspace=0.15)
    axes = axes.ravel()

    ace_ref = []
    ace_test = []
    num_ref = []
    num_test = []
    for ifig, (basin, basin_info) in enumerate(basin_dict.items()):
        ace_ref.append(ref_metrics[basin][0])
        ace_test.append(test_metrics[basin][0])
        num_ref.append(ref_metrics[basin][3])
        num_test.append(test_metrics[basin][3])
        ax = axes[ifig]
        ax.plot(np.arange(1, 7), ref_metrics[basin][1], "k", linewidth=2, label="Ref")
        ax.plot(
            np.arange(1, 7), test_metrics[basin][1], "--k", linewidth=2, label="Test"
        )
        ax.legend()
        ax.set_title(basin_info[0])
        ax.set_ylabel("Fraction")
        plt.xticks([1, 2, 3, 4, 5, 6], ["Cat0", "Cat1", "Cat2", "Cat3", "Cat4", "Cat5"])
        ax.xaxis.set_tick_params(labelbottom=True)

    plt.ylim(-0.05, 0.9)
    plt.title("TC intensity")
    output_file_name = "tc-intensity"
    for f in parameter.output_format:
        f = f.lower().split(".")[-1]
        fnm = os.path.join(
            get_output_dir(parameter.current_set, parameter),
            output_file_name + "." + f,
        )
        plt.savefig(fnm)
        print("Plot saved in: " + fnm)
        plt.close()

    # TC frequency of each basins
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(111)

    N = 6
    ind = np.arange(N)  # the x locations for the groups
    width = 0.27

    # ref_vals = [pdf_at_obs, pdf_wp_obs, pdf_ep_obs, pdf_ni_obs, pdf_si_obs, pdf_sp_obs]
    ref_vals = num_ref / np.sum(num_ref)
    rects1 = ax.bar(ind + width / 2, ref_vals, width, color="darkgrey")
    test_vals = num_test / np.sum(num_test)
    rects2 = ax.bar(ind - width / 2, test_vals, width, color="black")
    print(
        # "total number of storms",
        # num_storms,
        "total number based on 6 basins",
        np.sum(num_test),
        num_test,
    )

    ax.set_xticks(ind)
    ax.set_xticklabels(
        (
            "North Atlantic",
            "Northwest Pacific",
            "Eastern Pacific",
            "North Indian",
            "South Indian",
            "South Pacific",
        )
    )

    # TODO num_year
    num_year = 10
    ax.legend(
        (rects2[0], rects1[0]),
        (
            "Observations (~{})".format(int(np.sum(num_ref) / num_year)),
            "E3SM-HR (~{})".format(int(np.sum(num_test) / num_year)),
        ),
    )
    ax.set_ylabel("Fraction")
    ax.set_title("Relative frequency of TCs for each ocean basins")

    output_file_name = "tc-frequency"
    for f in parameter.output_format:
        f = f.lower().split(".")[-1]
        fnm = os.path.join(
            get_output_dir(parameter.current_set, parameter),
            output_file_name + "." + f,
        )
        plt.savefig(fnm)
        print("Plot saved in: " + fnm)
        plt.close()

    # plt.ylim(-.05,.9)
    fig1 = plt.figure(figsize=(12, 6))
    ax = fig1.add_subplot(111)

    N = 6
    ind = np.arange(N)  # the x locations for the groups
    width = 0.27

    ref_vals = ace_ref / np.sum(ace_ref)
    rects2 = ax.bar(ind - width / 2, ref_vals, width, color="black")
    test_vals = ace_test / np.sum(ace_test)
    rects1 = ax.bar(ind + width / 2, test_vals, width, color="darkgrey")

    ax.set_xticks(ind)
    ax.set_xticklabels(
        (
            "North Atlantic",
            "Northwest Pacific",
            "Eastern Pacific",
            "North Indian",
            "South Indian",
            "South Pacific",
        )
    )

    ax.legend((rects2[0], rects1[0]), ("ref", "test"))
    ax.set_ylabel("Fraction")
    ax.set_title(
        "Distribution of accumulated cyclone energy (ACE) among various ocean basins"
    )
    output_file_name = "ace-distribution"
    for f in parameter.output_format:
        f = f.lower().split(".")[-1]
        fnm = os.path.join(
            get_output_dir(parameter.current_set, parameter),
            output_file_name + "." + f,
        )
        plt.savefig(fnm)
        print("Plot saved in: " + fnm)
        plt.close()

    fig, axes = plt.subplots(2, 3, figsize=(12, 6), sharex=True, sharey=True)
    fig.subplots_adjust(hspace=0.4, wspace=0.15)
    axes = axes.ravel()

    for ifig, (basin, basin_info) in enumerate(basin_dict.items()):
        ax = axes[ifig]
        ax.plot(np.arange(1, 13), ref_metrics[basin][2], "k", linewidth=2, label="Obs")
        ax.plot(
            np.arange(1, 13),
            test_metrics[basin][2],
            "--k",
            linewidth=2,
            label="E3SM-HR",
        )
        ax.legend()
        ax.set_title(basin_info[0])
        ax.set_ylabel("Fraction")
        plt.xticks(
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
            ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"],
        )
        ax.xaxis.set_tick_params(labelbottom=True)
    output_file_name = "tc-frequency-annual-cycle"
    for f in parameter.output_format:
        f = f.lower().split(".")[-1]
        fnm = os.path.join(
            get_output_dir(parameter.current_set, parameter),
            output_file_name + "." + f,
        )
        plt.savefig(fnm)
        print("Plot saved in: " + fnm)
        plt.close()

    ##########################################################
    # Plotting tracks density
    # file_name = "cyclones_all_obs.nc"
    # fin = cdms2.open(os.path.join(test_data_path, file_name))
    # var = fin("density", squeeze=1) / 38

    clevs = np.arange(0, 1.1, 0.1)

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(10, 10),
        subplot_kw={"projection": ccrs.PlateCarree(central_longitude=180)},
    )

    for ax in axes:
        if ax == axes[0]:
            var = test["cyclone_density"]  # /num_year
        else:
            var = ref["cyclone_density"]  # /num_year
        ax.coastlines()
        p1 = ax.contourf(
            var.getLongitude(),
            var.getLatitude(),
            var,
            transform=ccrs.PlateCarree(),
            extend="both",
            levels=clevs,
            cmap="jet",
        )

        cbax = fig.add_axes([0.92, 0.1, 0.01, 0.7])
        cbar = fig.colorbar(p1, cax=cbax)
        cbar.ax.tick_params(labelsize=9.0, length=0)

    output_file_name = "tc-density-map"
    for f in parameter.output_format:
        f = f.lower().split(".")[-1]
        fnm = os.path.join(
            get_output_dir(parameter.current_set, parameter),
            output_file_name + "." + f,
        )
        plt.savefig(fnm)
        print("Plot saved in: " + fnm)
        plt.close()

    ##########################################################
    # African Easterly Waves (AEWs) counts
    # file_name = "aew_era5_stitch_hist_5e-6.nc"
    # fin = cdms2.open(os.path.join(reference_data_path, file_name))
    # var = fin("density", squeeze=1)

    clevs = np.arange(0, 15.1, 1)
    fig, axes = plt.subplots(
        2,
        1,
        figsize=(10, 10),
        subplot_kw={"projection": ccrs.PlateCarree(central_longitude=180)},
    )
    for ax in axes:
        if ax == axes[0]:
            var = test["aew_density"]  # /num_year
        else:
            var = ref["aew_density"]  # /num_year
        ax.coastlines()
        p1 = ax.contourf(
            var.getLongitude(),
            var.getLatitude(),
            var,
            transform=ccrs.PlateCarree(),
            extend="both",
            levels=clevs,
            cmap="jet",
        )

        cbax = fig.add_axes([0.92, 0.1, 0.01, 0.7])
        cbar = fig.colorbar(p1, cax=cbax)
        cbar.ax.tick_params(labelsize=9.0, length=0)

    output_file_name = "aew-density-map"
    for f in parameter.output_format:
        f = f.lower().split(".")[-1]
        fnm = os.path.join(
            get_output_dir(parameter.current_set, parameter),
            output_file_name + "." + f,
        )
        plt.savefig(fnm)
        print("Plot saved in: " + fnm)
        plt.close()

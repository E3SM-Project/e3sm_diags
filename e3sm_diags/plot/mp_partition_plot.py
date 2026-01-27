import matplotlib

from e3sm_diags.logger import _setup_child_logger
from e3sm_diags.plot.utils import _save_main_plot

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # isort:skip  # noqa: E402

logger = _setup_child_logger(__name__)

plotTitle = {"fontsize": 11.5}
plotSideTitle = {"fontsize": 9.5}


def plot(metrics_dict, parameter):
    fig, ax = plt.subplots()

    ax.plot(
        metrics_dict["test"]["T"],
        metrics_dict["test"]["LCF"],
        color="black",
        label=parameter.test_name_yrs,
        zorder=1,
    )

    if parameter.run_type == "model-vs-model":
        ax.plot(
            metrics_dict["ref"]["T"],
            metrics_dict["ref"]["LCF"],
            color="blue",
            label=parameter.ref_name_yrs,
            zorder=1,
        )

    # Plot multiple sets of reference datasets as background
    # Results from Hu et al. 2010
    ax.plot(
        metrics_dict["obs"]["T"],
        metrics_dict["obs"]["LCF"],
        color="red",
        label="Hu et al. 2010",
        zorder=2,
    )
    # Results from McCoy et al. 2016
    ax.hlines(
        y=0.5,
        xmin=254,
        xmax=258,
        linewidth=5,
        color="r",
        label="254K < T5050 CALIPSO < 258K \n McCoy et al. 2016",
    )
    # Results from E3SM.v2.LR.historical
    # ax.plot(
    #     metrics_dict["E3SM.v2.LR.historical"]["T"],
    #     metrics_dict["E3SM.v2.LR.historical"]["LCF"],
    #     color="green",
    #     label="E3SM v2.LR.historical(1985-2014)",
    #     zorder=1,
    # )
    # Results from CMIP model results from McCoy et al. 2015
    for idx, imod in enumerate(list(metrics_dict["cmip5"].keys())):
        cmip5 = metrics_dict["cmip5"]
        if idx == 0:
            ax.plot(
                cmip5[imod]["T"],
                cmip5[imod]["LCF"],
                linewidth=1,
                color="lightgrey",
                label="CMIP5",
                zorder=-1,
            )
        else:
            ax.plot(
                cmip5[imod]["T"],
                cmip5[imod]["LCF"],
                linewidth=1,
                color="lightgrey",
                zorder=-1,
            )  # , label=cmip5['MODEL'][imod][0][0])

    ax.hlines(y=0.5, xmin=220, xmax=280, linewidth=1, linestyles="--", color="grey")
    ax.set_ylabel("Liquid Condensate Fraction")
    ax.set_xlim(220, 280)
    ax.set_ylim(0, 1.2)
    ax.set_xlabel("Temperature (K)")
    ax.legend(loc="upper left")
    ax.set_title("Mixed-phase Partition LCF [30S - 70S]", pad=20)

    # Add variable information as subtitle below title
    if hasattr(parameter, "mp_variables_used"):
        var_text = parameter.mp_variables_used

        # Add reference variables for model-vs-model mode
        if parameter.run_type == "model-vs-model" and hasattr(
            parameter, "mp_ref_variables_used"
        ):
            # Check if variables are the same
            if parameter.mp_variables_used == parameter.mp_ref_variables_used:
                var_text = parameter.mp_variables_used
            else:
                var_text = f"Test {parameter.mp_variables_used}; Ref {parameter.mp_ref_variables_used}"

        plt.figtext(0.5, 0.90, var_text, ha="center", fontsize=8)

    _save_main_plot(parameter)

    plt.close(fig)

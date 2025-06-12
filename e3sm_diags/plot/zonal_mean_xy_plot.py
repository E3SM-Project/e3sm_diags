import matplotlib
import xarray as xr
import xcdat as xc

from e3sm_diags.logger import _setup_child_logger
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.plot.utils import _save_plot

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # isort:skip  # noqa: E402

logger = _setup_child_logger(__name__)

# Plot title and side title configurations.
PLOT_TITLE = {"fontsize": 12.5}
PLOT_SIDE_TITLE = {"fontsize": 11.5}

# Position and sizes of subplot axes in page coordinates (0 to 1)
PANEL_CONFIGS = [
    (0.1500, 0.5500, 0.7500, 0.3000),
    (0.1500, 0.1300, 0.7500, 0.3000),
]

# Border padding relative to subplot axes for saving individual PANELs
# (left, bottom, right, top) in page coordinates
BORDER_PADDING = (-0.14, -0.06, 0.04, 0.08)


def plot(
    parameter: CoreParameter,
    da_test: xr.DataArray,
    da_ref: xr.DataArray,
    da_diff: xr.DataArray,
):
    """Plot the variable's metrics generated for the zonal_mean_xy set.

    Parameters
    ----------
    parameter : CoreParameter
        The CoreParameter object containing plot configurations.
    da_test : xr.DataArray
        The test data.
    da_ref : xr.DataArray
        The reference data.
    da_diff : xr.DataArray
        The difference data.
    """
    # Create figure
    fig = plt.figure(figsize=parameter.figsize, dpi=parameter.dpi)
    lat_test = xc.get_dim_coords(da_test, axis="Y")
    lat_ref = xc.get_dim_coords(da_ref, axis="Y")
    lat_diff = xc.get_dim_coords(da_diff, axis="Y")
    long_name = da_test.attrs["long_name"]

    # Top PANEL
    ax1 = fig.add_axes(PANEL_CONFIGS[0])
    ax1.plot(lat_test, da_test.values, "k", linewidth=2)
    ax1.plot(
        lat_ref,
        da_ref.values,
        "r",
        linewidth=2,
    )
    ax1.set_xticks([-90, -60, -30, 0, 30, 60, 90])
    ax1.set_xlim(-90, 90)
    ax1.tick_params(labelsize=11.0, direction="out", width=1)
    ax1.xaxis.set_ticks_position("bottom")
    ax1.set_ylabel(long_name + " (" + da_test.units + ")")

    test_title = "Test" if parameter.test_title == "" else parameter.test_title
    test_title += " : {}".format(parameter.test_name_yrs)
    ref_title = (
        "Reference" if parameter.reference_title == "" else parameter.reference_title
    )
    ref_title += " : {}".format(parameter.ref_name_yrs)
    fig.text(
        PANEL_CONFIGS[0][0],
        PANEL_CONFIGS[0][1] + PANEL_CONFIGS[0][3] + 0.03,
        test_title,
        ha="left",
        fontdict=PLOT_SIDE_TITLE,
        color="black",
    )
    fig.text(
        PANEL_CONFIGS[0][0],
        PANEL_CONFIGS[0][1] + PANEL_CONFIGS[0][3] + 0.01,
        ref_title,
        ha="left",
        fontdict=PLOT_SIDE_TITLE,
        color="red",
    )

    # Bottom PANEL
    ax2 = fig.add_axes(PANEL_CONFIGS[1])
    ax2.plot(lat_diff, da_diff.values, "k", linewidth=2)
    ax2.axhline(y=0, color="0.5")
    ax2.set_title(parameter.diff_title, fontdict=PLOT_SIDE_TITLE, loc="center")
    ax2.set_xticks([-90, -60, -30, 0, 30, 60, 90])
    ax2.set_xlim(-90, 90)
    ax2.tick_params(labelsize=11.0, direction="out", width=1)
    ax2.xaxis.set_ticks_position("bottom")
    ax2.set_ylabel(long_name + " (" + da_test.units + ")")

    # Figure title
    fig.suptitle(parameter.main_title, x=0.5, y=0.95, fontsize=18)

    _save_plot(fig, parameter, PANEL_CONFIGS, BORDER_PADDING)

    plt.close(fig)

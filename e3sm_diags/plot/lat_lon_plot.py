from __future__ import annotations

from typing import TYPE_CHECKING

import matplotlib
import xarray as xr

from e3sm_diags.logger import custom_logger
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.plot.utils import _add_colormap, _save_plot

if TYPE_CHECKING:
    from e3sm_diags.driver.lat_lon_driver import MetricsDict


matplotlib.use("Agg")
import matplotlib.pyplot as plt  # isort:skip  # noqa: E402

logger = custom_logger(__name__)


def plot(
    parameter: CoreParameter,
    da_test: xr.DataArray,
    da_ref: xr.DataArray | None,
    da_diff: xr.DataArray | None,
    metrics_dict: MetricsDict,
):
    """Plot the variable's metrics generated for the lat_lon set.

    Parameters
    ----------
    parameter : CoreParameter
        The CoreParameter object containing plot configurations.
    da_test : xr.DataArray
        The test data.
    da_ref : xr.DataArray | None
        The optional reference data.
    da_diff : xr.DataArray | None
        The difference between `da_test` and `da_ref` (both are gridded to
        the lower resolution of the two beforehand).
    metrics_dict : Metrics
        The metrics.
    """
    fig = plt.figure(figsize=parameter.figsize, dpi=parameter.dpi)
    fig.suptitle(parameter.main_title, x=0.5, y=0.96, fontsize=18)

    # The variable units.
    units = metrics_dict["unit"]

    # Add the first subplot for test data.
    min1 = metrics_dict["test"]["min"]  # type: ignore
    mean1 = metrics_dict["test"]["mean"]  # type: ignore
    max1 = metrics_dict["test"]["max"]  # type: ignore

    _add_colormap(
        0,
        da_test,
        fig,
        parameter,
        parameter.test_colormap,
        parameter.contour_levels,
        title=(parameter.test_name_yrs, parameter.test_title, units),  # type: ignore
        metrics=(max1, mean1, min1),  # type: ignore
    )

    # Add the second and third subplots for ref data and the differences,
    # respectively.
    if da_ref is not None and da_diff is not None:
        min2 = metrics_dict["ref"]["min"]  # type: ignore
        mean2 = metrics_dict["ref"]["mean"]  # type: ignore
        max2 = metrics_dict["ref"]["max"]  # type: ignore

        _add_colormap(
            1,
            da_ref,
            fig,
            parameter,
            parameter.reference_colormap,
            parameter.contour_levels,
            title=(parameter.ref_name_yrs, parameter.reference_title, units),  # type: ignore
            metrics=(max2, mean2, min2),  # type: ignore
        )

        min3 = metrics_dict["diff"]["min"]  # type: ignore
        mean3 = metrics_dict["diff"]["mean"]  # type: ignore
        max3 = metrics_dict["diff"]["max"]  # type: ignore
        r = metrics_dict["misc"]["rmse"]  # type: ignore
        c = metrics_dict["misc"]["corr"]  # type: ignore

        _add_colormap(
            2,
            da_diff,
            fig,
            parameter,
            parameter.diff_colormap,
            parameter.diff_levels,
            title=(None, parameter.diff_title, units),  # type: ignore
            metrics=(max3, mean3, min3, r, c),  # type: ignore
        )

    _save_plot(fig, parameter)

    plt.close()

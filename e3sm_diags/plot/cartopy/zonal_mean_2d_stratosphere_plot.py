import xarray as xr

from e3sm_diags.driver.utils.type_annotations import MetricsDict
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.plot.cartopy.zonal_mean_2d_plot import plot as base_plot


def plot(
    parameter: CoreParameter,
    da_test: xr.DataArray,
    da_ref: xr.DataArray,
    da_diff: xr.DataArray,
    metrics_dict: MetricsDict,
):
    return base_plot(parameter, da_test, da_ref, da_diff, metrics_dict)

import copy

from e3sm_diags.driver.zonal_mean_2d_driver import run_diag as base_run_diag
from e3sm_diags.parameter.zonal_mean_2d_parameter import ZonalMean2dParameter
from e3sm_diags.parameter.zonal_mean_2d_stratosphere_parameter import (
    DEFAULT_PLEVS,
    ZonalMean2dStratosphereParameter,
)

DEFAULT_PLEVS = copy.deepcopy(DEFAULT_PLEVS)


def run_diag(
    parameter: ZonalMean2dStratosphereParameter,
) -> ZonalMean2dParameter:
    return base_run_diag(parameter, default_plevs=DEFAULT_PLEVS)

from typing import TYPE_CHECKING

from e3sm_diags.driver.lat_lon_driver import run_diag as base_run_diag

if TYPE_CHECKING:
    from e3sm_diags.parameter.core_parameter import CoreParameter
    from e3sm_diags.parameter.lat_lon_land_parameter import LatLonLandParameter


def run_diag(parameter: LatLonLandParameter) -> CoreParameter:
    return base_run_diag(parameter)

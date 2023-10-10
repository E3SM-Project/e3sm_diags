from __future__ import annotations

from typing import TYPE_CHECKING

from e3sm_diags.driver.lat_lon_driver import run_diag as base_run_diag

if TYPE_CHECKING:
    from e3sm_diags.parameter.core_parameter import CoreParameter
    from e3sm_diags.parameter.lat_lon_river_parameter import LatLonRiverParameter


def run_diag(parameter: LatLonRiverParameter) -> CoreParameter:
    return base_run_diag(parameter)

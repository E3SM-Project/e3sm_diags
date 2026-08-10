from __future__ import annotations

from e3sm_diags.plot.cartopy_utils import SPHERICAL_GLOBE, plate_carree


def test_plate_carree_uses_spherical_globe() -> None:
    projection = plate_carree(central_longitude=180.0)

    assert projection.globe is SPHERICAL_GLOBE
    assert projection.globe.ellipse is None
    assert projection.globe.semimajor_axis == 6_378_137.0
    assert projection.globe.semiminor_axis == 6_378_137.0
    assert projection.proj4_params["lon_0"] == 180.0

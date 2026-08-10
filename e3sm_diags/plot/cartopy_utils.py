from __future__ import annotations

import cartopy.crs as ccrs

SPHERICAL_GLOBE = ccrs.Globe(
    # PROJ 9.8 can reject Cartopy's implicit ellipsoid-to-sphere conversion
    # while projecting global PlateCarree geometries. Use an explicit sphere
    # for the affected PlateCarree axes and transforms.
    ellipse=None,
    semimajor_axis=6_378_137.0,
    semiminor_axis=6_378_137.0,
)


def plate_carree(central_longitude: float = 0.0) -> ccrs.PlateCarree:
    """Create a Plate Carrée CRS backed by a spherical globe.

    This is a workaround to Cartopy's compatibility issue with PROJ 9.8.
    Refer to https://github.com/SciTools/cartopy/pull/2653

    Parameters
    ----------
    central_longitude : float, optional
        The central longitude of the projection, by default 0.0.

    Returns
    -------
    cartopy.crs.PlateCarree
        A Plate Carrée coordinate reference system that uses the WGS84
        semi-major axis as a spherical radius.
    """
    return ccrs.PlateCarree(
        central_longitude=central_longitude,
        globe=SPHERICAL_GLOBE,
    )

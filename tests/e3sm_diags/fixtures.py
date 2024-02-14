from typing import Literal

import cftime
import numpy as np
import xarray as xr

time_decoded = xr.DataArray(
    data=np.array(
        [
            cftime.DatetimeGregorian(2000, 1, 16, 12, 0, 0, 0, has_year_zero=False),
            cftime.DatetimeGregorian(2000, 2, 15, 12, 0, 0, 0, has_year_zero=False),
            cftime.DatetimeGregorian(2000, 3, 16, 12, 0, 0, 0, has_year_zero=False),
            cftime.DatetimeGregorian(2000, 4, 16, 0, 0, 0, 0, has_year_zero=False),
            cftime.DatetimeGregorian(2000, 5, 16, 12, 0, 0, 0, has_year_zero=False),
        ],
        dtype=object,
    ),
    dims=["time"],
    attrs={
        "axis": "T",
        "long_name": "time",
        "standard_name": "time",
    },
)
lat = xr.DataArray(
    data=np.array([-90, -88.75, 88.75, 90]),
    dims=["lat"],
    attrs={"units": "degrees_north", "axis": "Y", "standard_name": "latitude"},
)

lon = xr.DataArray(
    data=np.array([0, 1.875, 356.25, 358.125]),
    dims=["lon"],
    attrs={"units": "degrees_east", "axis": "X", "standard_name": "longitude"},
)

lev = xr.DataArray(
    data=[800, 600, 400, 200],
    dims=["lev"],
    attrs={"units": "mb", "positive": "down", "axis": "Z"},
)


def generate_lev_dataset(
    long_name: Literal["hybrid", "pressure", "isobaric"], pressure_vars: bool = True
) -> xr.Dataset:
    """Generate a dataset with a Z axis ("lev").

    Parameters
    ----------
    long_name : {"hybrid", "pressure", "isobaric"}
        The long name attribute for the Z axis coordinates.
    pressure_vars : bool, optional
        Whether or not to include variables ps, hyam, or hybm, by default True.

    Returns
    -------
    xr.Dataset
    """
    ds = xr.Dataset(
        data_vars={
            "so": xr.DataArray(
                name="so",
                data=np.ones((5, 4, 4, 4)),
                coords={"time": time_decoded, "lev": lev, "lat": lat, "lon": lon},
                attrs={"units": "ppt"},
            ),
        },
        coords={
            "lat": lat.copy(),
            "lon": lon.copy(),
            "time": time_decoded.copy(),
            "lev": lev.copy(),
        },
    )

    ds["time"].encoding["calendar"] = "standard"

    ds = ds.bounds.add_missing_bounds(axes=["X", "Y", "Z", "T"])

    ds["lev"].attrs["axis"] = "Z"
    ds["lev"].attrs["bounds"] = "lev_bnds"
    ds["lev"].attrs["long_name"] = long_name

    if pressure_vars:
        ds["ps"] = xr.DataArray(
            name="ps",
            data=np.ones((5, 4, 4)),
            coords={"time": ds.time, "lat": ds.lat, "lon": ds.lon},
            attrs={"long_name": "surface_pressure", "units": "Pa"},
        )

        if long_name == "hybrid":
            ds["hyam"] = xr.DataArray(
                name="hyam",
                data=np.ones((4)),
                coords={"lev": ds.lev},
                attrs={"long_name": "hybrid A coefficient at layer midpoints"},
            )
            ds["hybm"] = xr.DataArray(
                name="hybm",
                data=np.ones((4)),
                coords={"lev": ds.lev},
                attrs={"long_name": "hybrid B coefficient at layer midpoints"},
            )

    return ds

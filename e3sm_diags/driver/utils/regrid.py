from typing import List

import cf_xarray as cfxr  # noqa: F401
import xarray as xr
import xcdat as xc


def has_z_axis_coords(data_var: xr.DataArray) -> bool:
    """Checks whether the data variable has a Z axis.

    Parameters
    ----------
    data_var : xr.DataArray
        The data variable.

    Returns
    -------
    bool
        True if data variable has Z axis, else False.
    """
    try:
        get_z_axis_coords(data_var)
        return True
    except KeyError:
        return False


def get_z_axis_coords(data_var: xr.DataArray) -> xr.DataArray:
    """Gets the Z axis coordinates.

    Returns True if:
    - Data variable has a "Z" axis in the cf-xarray mapping dict
    - A coordinate has a matching "positive" attribute ("up" or "down")
    - A coordinate has a matching "name"
    - # TODO: conditional for valid pressure units with "Pa"

    Parameters
    ----------
    data_var : xr.DataArray
        The data variable.

    Returns
    -------
    xr.DataArray
        The Z axis coordinates.

    Notes
    -----
    Based on
    - https://cdms.readthedocs.io/en/latest/_modules/cdms2/avariable.html#AbstractVariable.getLevel
    - https://cdms.readthedocs.io/en/latest/_modules/cdms2/axis.html#AbstractAxis.isLevel
    """
    try:
        z_coords = xc.get_dim_coords(data_var, axis="Z")
        return z_coords
    except KeyError:
        pass

    for coord in data_var.coords.values():
        if coord.name in ["lev", "plev", "depth"]:
            return coord

    raise KeyError(
        f"No Z axis coordinate were found in the '{data_var.name}' "
        "Make sure the variable has Z axis coordinates"
    )


def select_region(
    region: str,
    var: xr.DataArray,
    land_frac: xr.DataArray,
    ocean_frac: xr.DataArray,
    regrid_tool: str,
    regrid_method: str,
) -> xr.DataArray:
    """Select desired regions for the variable.

    Parameters
    ----------
    region : str.
        TODO: The region, options include "global"....
    var : xr.DataArray
        The variable.
    land_frac : xr.DataArray
        The land mask.
    ocean_frac : xr.DataArray
        The ocean mask.
    regrid_tool : str
        The regridding tool to use. Options include "xesmf" and "regrid2".
    regrid_method: str
        TODO: The regridding method to use. Options include: "conservation",
        "bilinear", ....

    Returns
    -------
    xr.DataArray
        The variable subsetted by region(s).
    """

    pass


def convert_to_pressure_levels(
    dataset: xr.Dataset,
    data_var: xr.DataArray,
    plevs: List[float],
) -> xr.DataArray:
    # 1. Get the z coordinates (mv_plv)
    # 2. If "hybrid" in the `long_name`, convert hybrid to pressure
    # 3. If "pressure" or "isobaric" in `long_name`, convert pressure to plevs
    # 4. Raise ValueError if vertical level is neither hybrid nor pressure
    # 5. Return pressure levels
    z_levels = get_z_axis_coords(data_var)
    z_long_name = z_levels.attrs.get("long_name")

    if z_long_name is None:
        raise KeyError(
            f"The vertical level ({z_levels.name}) for '{data_var.name}' does "
            "not have a 'long_name' attribute to determine whether it hybrid "
            "or pressure."
        )

    z_long_name = z_long_name.lower()
    if "hybrid" in z_long_name:
        p_levels = _hybrid_to_plevs(dataset, data_var, plevs)
    elif "pressure" in z_long_name or "isobaric" in z_long_name:
        p_levels = _pressure_to_plevs(data_var, plevs)
    else:
        raise ValueError(
            f"The vertical level ({z_levels.name}) for '{data_var.name}' is "
            "not hybrid or pressure. Its long name must either be 'hybrid', "
            "'pressure', or 'isobaric'."
        )

    return p_levels


def _hybrid_to_plevs(
    dataset: xr.Dataset,
    var: xr.DataArray,
    plevs: List[float],
) -> xr.DataArray:
    ds = dataset.copy()
    ps = ds.get("ps")
    hyam = ds.get("hyam")
    hybm = ds.get("hybm")

    if ps is None or hyam is None or hybm is None:
        raise KeyError(
            f"The dataset for {var.name} does not contain 'ps', 'hyam' "
            "and/or 'hybm' to reconstruct pressure from hybrid."
        )

    # Convert units from 'Pa' to mb
    p0 = 1000.0  # noqa: F841
    ps = ps / 100.0

    # Reconstruct pressure from hybrid
    # Make sure z is positive down
    pass


def _pressure_to_plevs(var: xr.DataArray, plevs: List[float]) -> xr.DataArray:
    pass

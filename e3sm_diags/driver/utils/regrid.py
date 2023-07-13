from __future__ import annotations

from typing import List

import xarray as xr
import xcdat as xc


def has_z_axis(data_var: xr.DataArray) -> bool:
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
        get_z_axis(data_var)
        return True
    except KeyError:
        return False


def get_z_axis(data_var: xr.DataArray) -> xr.DataArray:
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


def regrid_z_axis_to_plevs(
    dataset: xr.Dataset,
    var_key: str,
    plevs: List[int] | List[float],
) -> xr.Dataset:
    """Regrid a variable's Z axis to the desired pressure levels.

    The Z axis (e.g., 'lev') must either by hybrid or pressure, which is
    determined by the "long_name" attribute. Valid "long_name" values include
    "hybrid", "isobaric", and "pressure".

    Parameters
    ----------
    dataset : xr.Dataset
        The dataset with the variable on a Z axis.
    var_key : str
        The variable key.
    plevs : List[int] | List[float]
        A 1-D array of floats or integers representing output pressure levels
        in mb units. This parameter is usually set by ``CoreParameter.plevs``
        attribute. For example, ``plevs=[850.0, 200.0]``.

    Returns
    -------
    xr.Dataset
        The dataset with the variables's Z axis regridded to the desired
        pressure levels.

    Raises
    ------
    KeyError
        If the Z axis has no "long_name" attribute to determine whether it is
        hybrid or pressure.
    ValueError
        If the Z axis "long_name" attribute is not "hybrid", "isobaric",
        or "pressure".
    """
    ds = dataset.copy()
    dv = ds[var_key]

    z_axis = get_z_axis(dv)
    z_long_name = z_axis.attrs.get("long_name")

    if z_long_name is None:
        raise KeyError(
            f"The vertical level ({z_axis.name}) for '{dv.name}' does "
            "not have a 'long_name' attribute to determine whether it is hybrid "
            "or pressure."
        )

    z_long_name = z_long_name.lower()

    if "pressure" in z_long_name or "isobaric" in z_long_name:
        ds_plevs = _pressure_to_plevs(ds, dv, plevs)
    elif "hybrid" in z_long_name:
        ds_plevs = _hybrid_to_plevs(ds, dv.name, plevs)
    else:
        raise ValueError(
            f"The vertical level ({z_axis.name}) for '{dv.name}' is "
            "not hybrid or pressure. Its long name must either be 'hybrid', "
            "'pressure', or 'isobaric'."
        )

    return ds_plevs


def _pressure_to_plevs(
    dataset: xr.Dataset,
    var: xr.DataArray,
    plevs: List[int] | List[float],
) -> xr.DataArray:
    """Regrid pressure data to desired pressure level(s) in mb units.

    Parameters
    ----------
    dataset : xr.Dataset
        The dataset with a variable using pressure data.
    var_key : xr.DataArray.
        The variable.
    plevs : List[int] | List[float]
        A 1-D array of floats or integers representing output pressure levels
        in mb units. This parameter is usually set by ``CoreParameter.plevs``
        attribute. For example, ``plevs=[850.0, 200.0]``.

    Returns
    -------
    xr.DataArray
        The variable with a Z axis regridded to pressure levels (mb units).
    """
    ds = dataset.copy()

    # Create the output pressure grid to regrid to using the `plevs` array.
    pressure_grid = xc.create_grid(z=xc.create_axis("lev", plevs))

    # Get the pressure data ("ps").
    ps = ds.get("ps")
    if ps is None:
        raise KeyError(
            f"The dataset for '{var}' does not contain 'ps' (to convert to the "
            "desired pressure level(s)."
        )

    # Convert ps from 'Pa' to 'mb' if it is not already in 'mb'.
    ps = _convert_units_from_pa_to_mb(ps)

    # Perform the vertical regridding using log linear method.
    result = ds.regridder.vertical(
        var.name,
        output_grid=pressure_grid,
        tool="xgcm",
        method="log",
        target_data=ps,
    )

    return result


def _hybrid_to_plevs(
    dataset: xr.Dataset,
    var_key: str,
    plevs: List[int] | List[float],
) -> xr.DataArray:
    """Regrid the variable's hybrid-sigma levels to the desired pressure levels.

    Steps:

        1. Convert hybrid-sigma levels to pressure data
        2. Regrid the pressure data to the pressure levels (plevs)
        3. Return the dataset with the Z axis on pressure levels

    Parameters
    ----------
    dataset : xr.Dataset
        The dataset with the variable using hybrid-sigma levels.
    var_key : var_key.
        The variable key.
    plevs : List[int] | List[float]
        A 1-D array of floats or integers representing output pressure levels
        in mb units. This parameter is usually set by ``CoreParameter.plevs``
        attribute. For example, ``plevs=[850.0, 200.0]``.

    Returns
    -------
    xr.DataArray
        The variable with a Z axis regridded to pressure levels (mb units).
    """
    ds = dataset.copy()

    # Convert hybrid-sigma levels to pressure data.
    # Formula: p(k) = hyam(k) * p0 + hybm(k) * ps
    pressure_data = _hybrid_to_pressure(ds, var_key)

    # Create the output pressure grid to regrid to using the `plevs` array.
    pressure_grid = xc.create_grid(z=xc.create_axis("lev", plevs))

    # TODO: Do we need to make sure z is positive down and why?

    # 4. Regrid the variable using the pressure grid, pressure data,
    # and the log linear method.
    result = ds.regridder.vertical(
        var_key,
        output_grid=pressure_grid,
        tool="xgcm",
        method="log",
        target_data=pressure_data,
    )

    return result


def _hybrid_to_pressure(dataset: xr.Dataset, var_key: str) -> xr.DataArray:
    """Convert hybrid-sigma levels to pressure data (mb).

    Formula: p(k) = hyam(k) * p0 + hybm(k) * ps
      * "hyam" - 1-D array equal to hybrid A coefficients.
      * "p0" - Scalar numeric value equal to surface reference pressure with
          the same units as "ps" (mb).
      * "hybm" - 1-D array equal to hybrid B coefficients.
      * "ps" - 2-D array equal to surface pressure data in Pa and converted to
          mb.

    Parameters
    ----------
    dataset : xr.Dataset
        The dataset containing the variable and hybrid levels.
    var_key : str
        The variable key.

    Returns
    -------
    xr.DataArray
        The variable's pressure data.

    Raises
    ------
    KeyError
        If the dataset does not contain pressure data ('ps') or any of the
        hybrid levels ('hyam' and 'hymb').
    """
    ds = dataset.copy()
    p0 = 1000.0  # in mb, which is 100000 in Pa (mb * 100)
    ps = ds.get("ps", ds.get("PS", None))
    hyam = ds.get("hyam")
    hybm = ds.get("hybm")

    if ps is None or hyam is None or hybm is None:
        raise KeyError(
            f"The dataset for '{var_key}' does not contain 'ps', 'hyam' "
            "and/or 'hybm' to reconstruct them to pressure data."
        )

    ps = _convert_units_from_pa_to_mb(ps)

    # Formula: p(k) = hyam(k) * p0 + hybm(k) * ps
    pressure_data = hyam * p0 + hybm * ps
    pressure_data.attrs["units"] = "mb"

    return pressure_data


def _convert_units_from_pa_to_mb(ps: xr.DataArray) -> xr.DataArray:
    """Convert ps from Pa (Pascal pressure) to mb (millibars).

    The more common unit on weather maps is mb.

    Parameters
    ----------
    ps : xr.DataArray
        2-D array equal to surface pressure data.

    Returns
    -------
    xr.DataArray
        ps in 'mb' units (if not already).

    Raises
    ------
    ValueError
        If 'ps' variable has no 'units' attribute.
    ValueError
        If 'ps' variable has units not in 'mb' or 'Pa'.
    """
    ps_units = ps.attrs.get("units")

    if ps_units is None:
        raise ValueError(
            "'{ps.name}' has no 'units' attribute to determine if data is in 'mb' or "
            "'Pa' units."
        )

    if ps_units == "Pa":
        ps = ps / 100.0
        ps.attrs["units"] = "mb"
    elif ps_units == "mb":
        pass
    else:
        raise ValueError(
            f"'{ps.name}' should be in 'mb' or 'Pa' (which gets converted to 'mb'), "
            f"not {ps_units}."
        )

    return ps

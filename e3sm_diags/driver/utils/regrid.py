from typing import List, Optional

import numpy as np  # noqa: F401
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


def convert_z_axis_to_pressure_levels(
    dataset: xr.Dataset,
    var_key: str,
    plevs: List[float],
) -> xr.DataArray:
    """Converts a variable's Z axis (vertical) to pressure levels.

    The Z axis must either by hybrid or pressure, which is determined by
    the "long_name" attribute ("hybrid", "isobaric" and "pressure").

    Parameters
    ----------
    dataset : xr.Dataset
        The dataset.
    var_key : str
        The variable key.
    plevs : List[float]
        A 1-D array of floats representing output pressure levels. This
        parameter is usually set by ``CoreParameter.plevs`` attribute.

    Returns
    -------
    xr.DataArray
        The variable with a Z axis using pressure levels.

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
    data_var = ds[var_key]

    z_levels = get_z_axis_coords(data_var)
    z_long_name = z_levels.attrs.get("long_name")

    if z_long_name is None:
        raise KeyError(
            f"The vertical level ({z_levels.name}) for '{data_var.name}' does "
            "not have a 'long_name' attribute to determine whether it is hybrid "
            "or pressure."
        )

    z_long_name = z_long_name.lower()
    if "pressure" in z_long_name or "isobaric" in z_long_name:
        p_levels = _pressure_to_plevs(ds, data_var, plevs)
    elif "hybrid" in z_long_name:
        p_levels = _hybrid_to_plevs(dataset, data_var, plevs)
    else:
        raise ValueError(
            f"The vertical level ({z_levels.name}) for '{data_var.name}' is "
            "not hybrid or pressure. Its long name must either be 'hybrid', "
            "'pressure', or 'isobaric'."
        )

    return p_levels


def _pressure_to_plevs(
    ds: xr.Dataset, var_key: str, plevs: List[float]
) -> xr.DataArray:
    """Convert pressure coordinates to desired pressure level(s)

    Parameters
    ----------
    ds : xr.Dataset
        The dataset.
    var_key : str
        The variable.
    plevs : List[float]
        A 1-D array of floats representing output pressure levels. This
        parameter is usually set by ``CoreParameter.plevs`` attribute.

    Returns
    -------
    xr.DataArray
        The variable with a Z axis using pressure levels (mb).
    """
    ds_cp = ds.copy()
    dv = ds_cp[var_key].copy()

    # Get the pressure data ("ps").
    ps = ds_cp.get("ps")

    if ps is None:
        raise KeyError(
            f"The dataset for '{var_key}' does not contain 'ps' (to convert to the "
            "desired pressure level(s)."
        )

    if dv.attrs["units"] == "Pa":
        ps = _convert_units_from_pa_to_mb(ps)

    # 3. Create the pressure grid, which is the output grid for interpolation.
    pressure_grid = xc.create_grid(z=xc.create_axis("lev", plevs))

    result = ds.regridder.vertical(var_key, pressure_grid, method="log", target_data=ps)

    return result


def _hybrid_to_plevs(
    dataset: xr.Dataset,
    var_key: str,
    plevs: List[float],
    lev_dim: Optional[str] = None,
) -> xr.DataArray:
    """Convert the variable's hybrid-sigma levels to the desired pressure levels.

    Parameters
    ----------
    dataset : xr.Dataset
        The dataset containing a Z axis.
    var_key : str
        The variable key.
    plevs : List[float]
        A 1-D array of floats representing output pressure levels. This
        parameter is usually set by ``CoreParameter.plevs`` attribute.
    lev_dim : Optional[str], optional
        The key for the level (Z) dimension, by default None. This function will
        try to map to the Z dimension, but if it cannot then this parameter
        must be set.

    Returns
    -------
    xr.DataArray
        The variable with a Z axis using pressure levels (mb).
    """
    ds = dataset.copy()

    # 2. Convert hybrid-sigma levels to pressure coordinates.
    pressure_coords = _hybrid_to_pressure(ds, var_key)

    # 3. Create the pressure grid, which is the output grid for interpolation.
    pressure_grid = xc.create_grid(z=xc.create_axis("lev", plevs))

    # TODO: Do we need to make sure z is positive down and why?

    # 4. Regrid the variable using the pressure grid, pressure coordinates,
    # and the log linear method.
    result = ds.regridder.vertical(
        var_key, pressure_grid, method="log", target_data=pressure_coords
    )

    return result


def _hybrid_to_pressure(ds: xr.Dataset, var_key: str) -> xr.DataArray:
    """Convert hybrid-sigma levels to pressure coordinates (mb).

    Formula: p(k) = hya(k) * p0 + hyb(k) * ps
      * "hya" - 1-D array equal to hybrid A coefficients
      * "p0" - Scalar numeric value equal to surface reference pressure with
          the same units as ps.
      * "hyb" - 1-D array equal to hybrid B coefficients
      * "ps" - 2-D array equal to surface pressure data in Pa or hPA (mb)

    Parameters
    ----------
    ds : xr.Dataset
        The dataset containing the variable and hybrid levels.
    var_key : str
        The variable key.

    Returns
    -------
    xr.DataArray
        The variable's pressure coordinates.

    Raises
    ------
    KeyError
        If the dataset does not contain any of the hybrid levels.
    """
    ds_cp = ds.copy()

    # Get the pressure data ("ps") and hybrid levels ("hya" and "hyb").
    ps = ds_cp.get("ps")
    hya = ds_cp.get("a")
    hyb = ds_cp.get("b")

    if ps is None or hya is None or hyb is None:
        raise KeyError(
            f"The dataset for '{var_key}' does not contain 'ps', 'hya' "
            "and/or 'hyb' to reconstruct them to pressure coordinates."
        )

    # Convert the hybrid levels to pressure coordinates.
    p0 = 1000.0
    ps = _convert_units_from_pa_to_mb(ps)
    pressure = hya * p0 + hyb * ps
    pressure.attrs["units "] = "mb"

    return pressure


def _convert_units_from_pa_to_mb(da: xr.DataArray) -> xr.DataArray:
    """Convert Pa (Pascal pressure units) converted to mb (millibars).

    The more common unit on weather maps is mb.

    Parameters
    ----------
    da : xr.DataArray
        The variable with Pa units.

    Returns
    -------
    xr.DataArray
        The variable with mb units.
    """
    da = da / 100.0
    da.attrs["units"] = "mb"

    return da

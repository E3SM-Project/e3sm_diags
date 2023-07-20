from __future__ import annotations

from typing import List, Literal, Optional, Tuple

import xarray as xr
import xcdat as xc

# Valid hybrid-sigma levels keys that can be found in datasets.
HYBRID_SIGMA_KEYS = {
    "p0": ("p0", "P0"),
    "ps": ("ps", "PS"),
    "hyam": ("hyam", "hya", "a"),
    "hybm": ("hybm", "hyb", "b"),
}


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
    return None


def regrid_to_lower_res(
    ds_a: xr.Dataset,
    ds_b: xr.DataArray,
    var_key: str,
    tool: Literal["esmf", "xesmf", "regrid2"],
    method: str,
) -> Tuple[xr.DataArray, xr.DataArray]:
    """Horizontally regrid two DataArray using the lower resolution of the two.

    A variable has a lower resolution if it has less latitude coordinates,
    and vice versa. If both resolutions are the same, no regridding will happen.

    Parameters
    ----------
    ds_a : xr.DataArray
        The first Dataset containing ``var_key``.
    ds_b : xr.DataArray
        The second Dataset containing ``var_key``.
    var_key : str
        The key of the variable in both datasets to regrid.
    tool : Literal["esmf", "xesmf", "regrid2"]
        The regridding tool to use. Note, "esmf" is accepted for backwards
        compatibility with e3sm_diags and is simply updated to "xesmf".
    method : str
        The regridding method to use. Refer to [1]_ for more information on
        these options.

        esmf/xesmf options:
          - "bilinear"
          - "conservative"
          - "conservative_normed"
          - "patch"
          - "nearest_s2d"
          - "nearest_d2s"

        regrid2 options:
          - "conservative"

    Returns
    -------
    Tuple[xr.DataArray, xr.DataArray]
        A tuple of both DataArrays regridded to the lower resolution of the two.

    References
    ----------
    .. [1] https://xcdat.readthedocs.io/en/stable/generated/xarray.Dataset.regridder.horizontal.html
    """
    # TODO: Accept "esmf" as `tool` value for now because `CoreParameter`
    # defines `self.regrid_tool="esmf"` by default and
    # `e3sm_diags.driver.utils.general.regrid_to_lower_res()` accepts "esmf".
    # Once this function is deprecated, we can remove "esmf" as an option here
    # and update `CoreParameter.regrid_tool` to "xesmf"`.
    if tool == "esmf":
        tool = "xesmf"

    lat_a = xc.get_dim_coords(ds_a[var_key], axis="Y")
    lat_b = xc.get_dim_coords(ds_b[var_key], axis="Y")

    is_a_lower_res = len(lat_a) < len(lat_b)
    is_b_lower_res = len(lat_b) < len(lat_a)

    if is_a_lower_res:
        output_grid = ds_a.regridder.grid
        ds_b_regrid = ds_b.regridder.horizontal(
            var_key, output_grid, tool=tool, method=method
        )

        return ds_a, ds_b_regrid
    elif is_b_lower_res:
        output_grid = ds_b.regridder.grid
        ds_a_regrid = ds_a.regridder.horizontal(
            var_key, output_grid, tool=tool, method=method
        )

        return ds_a_regrid, ds_b

    return ds_a, ds_b


def regrid_z_axis_to_plevs(
    dataset: xr.Dataset,
    var_key: str,
    plevs: List[int] | List[float],
) -> xr.Dataset:
    """Regrid a variable's Z axis to the desired pressure levels (mb units).

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
        pressure levels (mb units).

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

    # Hybrid must be the first conditional statement because the long_name attr
    # can be "hybrid sigma pressure coordinate" which includes "hybrid" and
    # "pressure".
    if "hybrid" in z_long_name:
        ds_plevs = _hybrid_to_plevs(ds, var_key, plevs)
    elif "pressure" in z_long_name or "isobaric" in z_long_name:
        ds_plevs = _pressure_to_plevs(ds, var_key, plevs)
    else:
        raise ValueError(
            f"The vertical level ({z_axis.name}) for '{dv.name}' is "
            "not hybrid or pressure. Its long name must either include 'hybrid', "
            "'pressure', or 'isobaric'."
        )

    return ds_plevs


def _hybrid_to_plevs(
    dataset: xr.Dataset,
    var_key: str,
    plevs: List[int] | List[float],
) -> xr.DataArray:
    """Regrid a variable's hybrid-sigma levels to the desired pressure levels.

    Steps:
        1. Create the output pressure grid using ``plevs``.
        2. Convert hybrid-sigma levels to pressure coordinates.
        3. Regrid the pressure coordinates to the output pressure grid (plevs).

    Parameters
    ----------
    dataset : xr.Dataset
        The dataset with the variable using hybrid-sigma levels.
    var_key : var_key.
        The variable key.
    plevs : List[int] | List[float]
        A 1-D array of floats or integers representing output pressure levels
        in mb units. For example, ``plevs=[850.0, 200.0]``. This parameter is
        usually set by the ``CoreParameter.plevs`` attribute.

    Returns
    -------
    xr.DataArray
        The variable with a Z axis regridded to pressure levels (mb units).
    """
    # TODO: Do we need to convert the Z axis to mb units if it is in PA?
    ds = dataset.copy()

    pressure_grid = xc.create_grid(z=xc.create_axis("lev", plevs))

    pressure_coords = _hybrid_to_pressure(ds, var_key)

    result = ds.regridder.vertical(
        var_key,
        output_grid=pressure_grid,
        tool="xgcm",
        method="log",
        target_data=pressure_coords,
    )

    return result


def _hybrid_to_pressure(dataset: xr.Dataset, var_key: str) -> xr.DataArray:
    """Regrid hybrid-sigma levels to pressure coordinates (mb).

    Formula: p(k) = hyam(k) * p0 + hybm(k) * ps
      * p: pressure data (mb).
      * hyam: 1-D array equal to hybrid A coefficients.
      * p0: Scalar numeric value equal to surface reference pressure with
          the same units as "ps" (mb).
      * hybm: 1-D array equal to hybrid B coefficients.
      * ps: 2-D array equal to surface pressure data (mb, converted from Pa).

    Parameters
    ----------
    dataset : xr.Dataset
        The dataset containing the variable and hybrid levels.
    var_key : str
        The variable key.

    Returns
    -------
    xr.DataArray
        The variable with a Z axis on pressure coordinates.

    Raises
    ------
    KeyError
        If the dataset does not contain pressure data (ps) or any of the
        hybrid levels (hyam, hymb).

    Notes
    -----
    This function is equivalent to `geocat.comp.interp_hybrid_to_pressure()`
    and `cdutil.vertical.reconstructPressureFromHybrid()`.
    """
    ds = dataset.copy()

    # p0 is statically set to mb (1000) instead of retrieved from the dataset
    # because the pressure data should be in mb.
    p0 = 1000
    ps = _get_hybrid_sigma_level(ds, "ps")
    hyam = _get_hybrid_sigma_level(ds, "hyam")
    hybm = _get_hybrid_sigma_level(ds, "hybm")

    if ps is None or hyam is None or hybm is None:
        raise KeyError(
            f"The dataset for '{var_key}' does not contain hybrid-sigma level 'ps', "
            "'hyam' and/or 'hybm' to use for reconstructing to pressure data."
        )

    ps = _convert_units_to_mb(ps)

    pressure_data = hyam * p0 + hybm * ps
    pressure_data.attrs["units"] = "mb"

    return pressure_data


def _get_hybrid_sigma_level(
    ds: xr.Dataset, name: Literal["ps", "p0", "hyam", "hybm"]
) -> Optional[xr.DataArray]:
    """Get the hybrid-sigma level xr.DataArray from the xr.Dataset.

    This functions retrieves the valid keys for the specified hybrid-sigma
    level and loops over them. A dictionary look-up is performed and the first
    match is passed. If there are no matches, None is returned.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset.
    name : {"ps", "p0", "hyam", "hybm"}
        The name of the hybrid-sigma level to get.

    Returns
    -------
    Optional[xr.DataArray]
        The hybrid-sigma level xr.DataArray if found or None.
    """
    keys = HYBRID_SIGMA_KEYS[name]

    for key in keys:
        da = ds.get(key)

        if da is not None:
            return da

    return None


def _pressure_to_plevs(
    dataset: xr.Dataset,
    var_key: str,
    plevs: List[int] | List[float],
) -> xr.DataArray:
    """Regrids pressure coordinates to the desired pressure level(s).

    Parameters
    ----------
    dataset : xr.Dataset
        The dataset with a variable using pressure data.
    var_key : str
        The variable key.
    plevs : List[int] | List[float]
        A 1-D array of floats or integers representing output pressure levels
        in mb units. This parameter is usually set by ``CoreParameter.plevs``
        attribute. For example, ``plevs=[850.0, 200.0]``.

    Returns
    -------
    xr.DataArray
        The variable with a Z axis on pressure levels (mb).
    """
    ds = dataset.copy()

    # Create the output pressure grid to regrid to using the `plevs` array.
    pressure_grid = xc.create_grid(z=xc.create_axis("lev", plevs))

    # Convert pressure coordinates to mb if it is not already in mb.
    lev_key = xc.get_dim_keys(ds[var_key], axis="Z")
    ds[lev_key] = _convert_units_to_mb(ds[lev_key])

    result = ds.regridder.vertical(
        var_key,
        output_grid=pressure_grid,
        tool="xgcm",
        method="log",
    )

    return result


def _convert_units_to_mb(da: xr.DataArray) -> xr.DataArray:
    """Convert DataArray to mb (millibars) if not in mb.

    Unit conversion formulas:
      * mb = Pa / 100
      * Pa = (mb * 100)

    The more common unit on weather maps is mb.

    Parameters
    ----------
    da : xr.DataArray
        An xr.DataArray, usually for "lev" or "ps".

    Returns
    -------
    xr.DataArray
        The DataArray in mb units.

    Raises
    ------
    ValueError
        If ``da`` DataArray has no 'units' attribute.
    ValueError
        If ``da`` DataArray has units not in mb or Pa.
    """
    units = da.attrs.get("units")

    if units is None:
        raise ValueError(
            "'{ps.name}' has no 'units' attribute to determine if data is in 'mb' or "
            "'Pa' units."
        )

    if units == "mb":
        pass
    elif units == "Pa":
        with xr.set_options(keep_attrs=True):
            da = da / 100.0

        da.attrs["units"] = "mb"
    else:
        raise ValueError(
            f"'{da.name}' should be in 'mb' or 'Pa' (which gets converted to 'mb'), "
            f"not '{units}'."
        )

    return da

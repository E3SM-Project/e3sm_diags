"""
This module defines general utilities for deriving variables, including unit
conversion functions, renaming variables, etc.
"""

import cf_units
import xarray as xr


def rename(new_name: str):
    """Given the new name, just return it."""
    # FIXME: This function does nothing.
    # Related issue: https://github.com/E3SM-Project/e3sm_diags/issues/796
    return new_name


def aplusb(var1: xr.DataArray, var2: xr.DataArray, target_units=None):
    """Returns var1 + var2. If both of their units are not the same,
    it tries to convert both of their units to target_units"""

    if target_units is not None:
        var1 = convert_units(var1, target_units)
        var2 = convert_units(var2, target_units)

    return var1 + var2


def convert_units(var: xr.DataArray, target_units: str) -> xr.DataArray:  # noqa: C901
    """Convert the units of a variable to the target units.

    Parameters
    ----------
    var : xr.DataArray
        The input data array with units to be converted.
    target_units : str
        The target units to convert the data array to.

    Returns
    -------
    xr.DataArray
        A new data array with the units converted to the target units.

    Notes
    -----
    The function handles various specific non-CF compliant unit conversions
    based on the variable name and current units. Otherwise it will use the
    ``cf_units`` library to attempt to perform the conversion if the units are
    recognized by UDUNITS-2 (must be CF-compliant).
    """
    var_new = var.copy()

    with xr.set_options(keep_attrs=True):
        if var_new.attrs.get("units") is None:
            if var_new.name == "SST":
                var_new.attrs["units"] = target_units
            elif var_new.name == "ICEFRAC":
                var_new.attrs["units"] = target_units
                var_new = 100.0 * var_new
            elif var_new.name == "AODVIS":
                var_new.attrs["units"] = target_units
            elif var_new.name == "AODDUST":
                var_new.attrs["units"] = target_units
        elif var_new.name == "FAREA_BURNED":
            var_new = var_new * 1e9
            var_new.attrs["units"] = target_units
        elif var_new.attrs["units"] == "gC/m^2":
            var_new = var_new / 1000.0
            var_new.attrs["units"] = target_units
        elif var_new.name == "FLOODPLAIN_VOLUME" and target_units == "km3":
            var_new = var_new / 1.0e9
            var_new.attrs["units"] = target_units
        elif var_new.name == "AOD_550_ann":
            var_new.attrs["units"] = target_units
        elif var_new.name == "AOD_550":
            var_new.attrs["units"] = target_units
        elif var_new.attrs["units"] == "C" and target_units == "DegC":
            var_new.attrs["units"] = target_units
        elif var_new.attrs["units"] == "N/m2" and target_units == "N/m^2":
            var_new.attrs["units"] = target_units
        elif (
            var_new.name == "AODVIS"
            or var_new.name == "AOD_550_ann"
            or var_new.name == "TOTEXTTAU"
        ):
            var_new.attrs["units"] = target_units
        elif var_new.attrs["units"] == "fraction":
            var_new = 100.0 * var_new
            var_new.attrs["units"] = target_units
        elif var_new.attrs["units"] == "mb":
            var_new.attrs["units"] = target_units
        elif var_new.attrs["units"] == "gpm":  # geopotential meter
            var_new = var_new / 9.8 / 100  # convert to hecto meter
            var_new.attrs["units"] = target_units
        elif var_new.attrs["units"] == "Pa/s":
            var_new = var_new / 100.0 * 24 * 3600
            var_new.attrs["units"] = target_units
        elif var_new.attrs["units"] == "mb/day":
            var_new.attrs["units"] = target_units
        elif var_new.name == "prw" and var_new.attrs["units"] == "cm":
            var_new = var_new * 10.0  # convert from 'cm' to 'kg/m2' or 'mm'
            var_new.attrs["units"] = target_units
        elif var_new.attrs["units"] in ["gC/m^2/s"] and target_units == "g*/m^2/day":
            var_new = var_new * 24 * 3600
            var_new.attrs["units"] = var_new.attrs["units"][0:7] + "day"
        elif (
            var_new.attrs["units"] in ["gN/m^2/s", "gP/m^2/s"]
            and target_units == "mg*/m^2/day"
        ):
            var_new = var_new * 24 * 3600 * 1000.0
            var_new.attrs["units"] = "m" + var_new.attrs["units"][0:7] + "day"
        elif var_new.attrs["units"] in ["gN/m^2/day", "gP/m^2/day", "gC/m^2/day"]:
            pass
        elif var_new.attrs["units"] == "gC/m^2" and target_units == "kgC/m^2":
            var_new = var_new / 1000.0
            var_new.attrs["units"] = target_units
        elif var_new.attrs["units"] == "gC/m^2/s" and target_units == "kgC/m^2":
            var_new = var_new / 1000.0
            var_new.attrs["units"] = target_units
        elif var_new.attrs["units"] == "gC/m^2/s" and target_units == "kgC/m^2/s":
            var_new = var_new / 1000.0
            var_new.attrs["units"] = target_units
        else:
            original_udunit = cf_units.Unit(var_new.attrs["units"])
            target_udunit = cf_units.Unit(target_units)

            var_new.values = original_udunit.convert(var_new.values, target_udunit)
            var_new.attrs["units"] = target_units

    return var_new


def _apply_land_sea_mask(
    var: xr.DataArray, var_mask: xr.DataArray, lower_limit: float
) -> xr.DataArray:
    """Apply a land or sea mask on the variable.

    Parameters
    ----------
    var : xr.DataArray
        The variable.
    var_mask : xr.DataArray
        The variable mask ("LANDFRAC" or "OCNFRAC").
    lower_limit : float
        Update the mask variable with a lower limit. All values below the
        lower limit will be masked.

    Returns
    -------
    xr.DataArray
        The masked variable.
    """
    cond = var_mask > lower_limit
    masked_var = var.where(cond=cond, drop=False)

    return masked_var

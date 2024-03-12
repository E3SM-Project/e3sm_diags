"""
This module defines general utilities for deriving variables, including unit
conversion functions, renaming variables, etc.
"""
import xarray as xr
from genutil import udunits


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


def convert_units(var: xr.DataArray, target_units: str):  # noqa: C901
    if var.attrs.get("units") is None:
        if var.name == "SST":
            var.attrs["units"] = target_units
        elif var.name == "ICEFRAC":
            var.attrs["units"] = target_units
            var = 100.0 * var
        elif var.name == "AODVIS":
            var.attrs["units"] = target_units
        elif var.name == "AODDUST":
            var.attrs["units"] = target_units
    elif var.name == "FAREA_BURNED":
        var = var * 1e9
        var.attrs["units"] = target_units
    elif var.attrs["units"] == "gC/m^2":
        var = var / 1000.0
        var.attrs["units"] = target_units
    elif var.name == "FLOODPLAIN_VOLUME" and target_units == "km3":
        var = var / 1.0e9
        var.attrs["units"] = target_units
    elif var.name == "AOD_550_ann":
        var.attrs["units"] = target_units
    elif var.name == "AOD_550":
        var.attrs["units"] = target_units
    elif var.attrs["units"] == "C" and target_units == "DegC":
        var.attrs["units"] = target_units
    elif var.attrs["units"] == "N/m2" and target_units == "N/m^2":
        var.attrs["units"] = target_units
    elif var.name == "AODVIS" or var.name == "AOD_550_ann" or var.name == "TOTEXTTAU":
        var.attrs["units"] = target_units
    elif var.attrs["units"] == "fraction":
        var = 100.0 * var
        var.attrs["units"] = target_units
    elif var.attrs["units"] == "mb":
        var.attrs["units"] = target_units
    elif var.attrs["units"] == "gpm":  # geopotential meter
        var = var / 9.8 / 100  # convert to hecto meter
        var.attrs["units"] = target_units
    elif var.attrs["units"] == "Pa/s":
        var = var / 100.0 * 24 * 3600
        var.attrs["units"] = target_units
    elif var.attrs["units"] == "mb/day":
        var = var
        var.attrs["units"] = target_units
    elif var.name == "prw" and var.attrs["units"] == "cm":
        var = var * 10.0  # convert from 'cm' to 'kg/m2' or 'mm'
        var.attrs["units"] = target_units
    elif var.attrs["units"] in ["gC/m^2/s"] and target_units == "g*/m^2/day":
        var = var * 24 * 3600
        var.attrs["units"] = var.attrs["units"][0:7] + "day"
    elif (
        var.attrs["units"] in ["gN/m^2/s", "gP/m^2/s"] and target_units == "mg*/m^2/day"
    ):
        var = var * 24 * 3600 * 1000.0
        var.attrs["units"] = "m" + var.attrs["units"][0:7] + "day"
    elif var.attrs["units"] in ["gN/m^2/day", "gP/m^2/day", "gC/m^2/day"]:
        pass
    else:
        # FIXME: Replace genutil.udunits module.
        temp = udunits(1.0, var.attrs["units"])
        coeff, offset = temp.how(target_units)

        # Keep all of the attributes except the units.
        with xr.set_options(keep_attrs=True):
            var = coeff * var + offset

        var.attrs["units"] = target_units

    return var


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

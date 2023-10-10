"""
This module defines general utilities for deriving variables, including unit
conversion functions, renaming variables, etc.
"""
from typing import TYPE_CHECKING, Optional, Tuple

import MV2
import numpy as np
import xarray as xr
from genutil import udunits

if TYPE_CHECKING:
    from cdms2.axis import FileAxis
    from cdms2.fvariable import FileVariable


def rename(new_name: str):
    """Given the new name, just return it."""
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


def adjust_prs_val_units(
    prs: "FileAxis", prs_val: float, prs_val0: Optional[float]
) -> float:
    """Adjust the prs_val units based on the prs.id"""
    # FIXME: Refactor this function to operate on xr.Dataset/xr.DataArray.
    # COSP v2 cosp_pr in units Pa instead of hPa as in v1
    # COSP v2 cosp_htmisr in units m instead of km as in v1
    adjust_ids = {"cosp_prs": 100, "cosp_htmisr": 1000}

    if prs_val0:
        prs_val = prs_val0
        if prs.id in adjust_ids.keys() and max(prs.getData()) > 1000:
            prs_val = prs_val * adjust_ids[prs.id]

    return prs_val


def determine_cloud_level(
    prs_low: float,
    prs_high: float,
    low_bnds: Tuple[int, int],
    high_bnds: Tuple[int, int],
) -> str:
    """Determines the cloud type based on prs values and the specified boundaries"""
    # Threshold for cloud top height: high cloud (<440hPa or > 7km), midlevel cloud (440-680hPa, 3-7 km) and low clouds (>680hPa, < 3km)
    if prs_low in low_bnds and prs_high in high_bnds:
        return "middle cloud fraction"
    elif prs_low in low_bnds:
        return "high cloud fraction"
    elif prs_high in high_bnds:
        return "low cloud fraction"
    else:
        return "total cloud fraction"


def cosp_bin_sum(
    cld: "FileVariable",
    prs_low0: Optional[float],
    prs_high0: Optional[float],
    tau_low0: Optional[float],
    tau_high0: Optional[float],
):
    # FIXME: Refactor this function to operate on xr.Dataset/xr.DataArray.
    """sum of cosp bins to calculate cloud fraction in specified cloud top pressure / height and
    cloud thickness bins, input variable has dimension (cosp_prs,cosp_tau,lat,lon)/(cosp_ht,cosp_tau,lat,lon)
    """
    prs: FileAxis = cld.getAxis(0)
    tau: FileAxis = cld.getAxis(1)

    prs_low: float = adjust_prs_val_units(prs, prs[0], prs_low0)
    prs_high: float = adjust_prs_val_units(prs, prs[-1], prs_high0)

    if prs_low0 is None and prs_high0 is None:
        prs_lim = "total cloud fraction"

    tau_high, tau_low, tau_lim = determine_tau(tau, tau_low0, tau_high0)

    if cld.id == "FISCCP1_COSP":  # ISCCP model
        cld_bin = cld(cosp_prs=(prs_low, prs_high), cosp_tau=(tau_low, tau_high))
        simulator = "ISCCP"
    if cld.id == "CLISCCP":  # ISCCP obs
        cld_bin = cld(isccp_prs=(prs_low, prs_high), isccp_tau=(tau_low, tau_high))

    if cld.id == "CLMODIS":  # MODIS
        prs_lim = determine_cloud_level(prs_low, prs_high, (440, 44000), (680, 68000))
        simulator = "MODIS"

        if prs.id == "cosp_prs":  # Model
            cld_bin = cld(
                cosp_prs=(prs_low, prs_high), cosp_tau_modis=(tau_low, tau_high)
            )
        elif prs.id == "modis_prs":  # Obs
            cld_bin = cld(modis_prs=(prs_low, prs_high), modis_tau=(tau_low, tau_high))

    if cld.id == "CLD_MISR":  # MISR model
        if max(prs) > 1000:  # COSP v2 cosp_htmisr in units m instead of km as in v1
            cld = cld[
                1:, :, :, :
            ]  # COSP v2 cosp_htmisr[0] equals to 0 instead of -99 as in v1, therefore cld needs to be masked manually
        cld_bin = cld(cosp_htmisr=(prs_low, prs_high), cosp_tau=(tau_low, tau_high))
        prs_lim = determine_cloud_level(prs_low, prs_high, (7, 7000), (3, 3000))
        simulator = "MISR"
    if cld.id == "CLMISR":  # MISR obs
        cld_bin = cld(misr_cth=(prs_low, prs_high), misr_tau=(tau_low, tau_high))

    cld_bin_sum = MV2.sum(MV2.sum(cld_bin, axis=1), axis=0)

    try:
        cld_bin_sum.long_name = simulator + ": " + prs_lim + " with " + tau_lim
        # cld_bin_sum.long_name = "{}: {} with {}".format(simulator, prs_lim, tau_lim)
    except BaseException:
        pass
    return cld_bin_sum


def determine_tau(
    tau: "FileAxis", tau_low0: Optional[float], tau_high0: Optional[float]
):
    # FIXME: Refactor this function to operate on xr.Dataset/xr.DataArray.
    tau_low = tau[0]
    tau_high = tau[-1]

    if tau_low0 is None and tau_high0:
        tau_high = tau_high0
        tau_lim = "tau <" + str(tau_high0)
    elif tau_high0 is None and tau_low0:
        tau_low = tau_low0
        tau_lim = "tau >" + str(tau_low0)
    elif tau_low0 is None and tau_high0 is None:
        tau_lim = str(tau_low) + "< tau < " + str(tau_high)
    else:
        tau_low = tau_low0
        tau_high = tau_high0
        tau_lim = str(tau_low) + "< tau < " + str(tau_high)

    return tau_high, tau_low, tau_lim


def cosp_histogram_standardize(cld: "FileVariable"):
    # TODO: Refactor this function to operate on xr.Dataset/xr.DataArray.
    """standarize cloud top pressure and cloud thickness bins to dimensions that
    suitable for plotting, input variable has dimention (cosp_prs,cosp_tau)"""
    prs = cld.getAxis(0)
    tau = cld.getAxis(1)

    prs[0]
    prs_high = prs[-1]
    tau[0]
    tau_high = tau[-1]

    prs_bounds = getattr(prs, "bounds")
    if prs_bounds is None:
        cloud_prs_bounds = np.array(
            [
                [1000.0, 800.0],
                [800.0, 680.0],
                [680.0, 560.0],
                [560.0, 440.0],
                [440.0, 310.0],
                [310.0, 180.0],
                [180.0, 0.0],
            ]
        )  # length 7
        prs.setBounds(np.array(cloud_prs_bounds, dtype=np.float32))

    tau_bounds = getattr(tau, "bounds")
    if tau_bounds is None:
        cloud_tau_bounds = np.array(
            [
                [0.3, 1.3],
                [1.3, 3.6],
                [3.6, 9.4],
                [9.4, 23],
                [23, 60],
                [60, 379],
            ]
        )  # length 6
        tau.setBounds(np.array(cloud_tau_bounds, dtype=np.float32))

    if cld.id == "FISCCP1_COSP":  # ISCCP model
        cld_hist = cld(cosp_tau=(0.3, tau_high))
    if cld.id == "CLISCCP":  # ISCCP obs
        cld_hist = cld(isccp_tau=(0.3, tau_high))

    if cld.id == "CLMODIS":  # MODIS
        try:
            cld_hist = cld(cosp_tau_modis=(0.3, tau_high))  # MODIS model
        except BaseException:
            cld_hist = cld(modis_tau=(0.3, tau_high))  # MODIS obs

    if cld.id == "CLD_MISR":  # MISR model
        if max(prs) > 1000:  # COSP v2 cosp_htmisr in units m instead of km as in v1
            cld = cld[
                1:, :, :, :
            ]  # COSP v2 cosp_htmisr[0] equals to 0 instead of -99 as in v1, therefore cld needs to be masked manually
            prs_high = 1000.0 * prs_high
        cld_hist = cld(cosp_tau=(0.3, tau_high), cosp_htmisr=(0, prs_high))
    if cld.id == "CLMISR":  # MISR obs
        cld_hist = cld(misr_tau=(0.3, tau_high), misr_cth=(0, prs_high))

    return cld_hist

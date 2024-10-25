"""This module defines formula functions used for deriving variables.

The function arguments usually accept variables represented by `xr.DataArray`.
NOTE: If a function involves arithmetic between two or more `xr.DataArray`,
the arithmetic should be wrapped with `with xr.set_options(keep_attrs=True)`
to keep attributes on the resultant `xr.DataArray`.
"""
from typing import List

import numpy as np
import xarray as xr

from e3sm_diags.derivations.utils import convert_units

AVOGADRO_CONST = 6.022e23
AIR_DENS = 1.225  # standard air density 1.225kg/m3


def sum_vars(vars: List[xr.DataArray]) -> xr.DataArray:
    """Sum DataArrays using Python's `.sum()` and preserve attrs.

    Pythons sum iterates over the iterable (the list of DataArrays) and
    adds all elements, which is different from NumPy which performs a sum
    reduction over an axis/axes. This function ensures the DataArray attributes
    are perserved by invoking the `.sum()` call within the context of
    `xr.set_options()`.

    Parameters
    ----------
    vars : List[xr.DataArray]
        A list of variables.

    Returns
    -------
    xr.DataArray
        The sum of the variables
    """
    with xr.set_options(keep_attrs=True):
        result: xr.DataArray = sum(vars)  # type: ignore

    return result


def qflxconvert_units(var: xr.DataArray):
    with xr.set_options(keep_attrs=True):
        if (
            var.attrs["units"] == "kg/m2/s"
            or var.attrs["units"] == "kg m-2 s-1"
            or var.attrs["units"] == "mm/s"
        ):
            # need to find a solution for units not included in udunits
            # var = convert_units( var, 'kg/m2/s' )
            var = var * 3600.0 * 24  # convert to mm/day
            var.attrs["units"] = "mm/day"
        elif var.attrs["units"] == "mm/hr":
            var = var * 24.0
            var.attrs["units"] = "mm/day"

    return var


def qsat(temp: xr.DataArray, surfp: xr.DataArray) -> xr.DataArray:
    # Function to calculate saturation specific humidity based on air
    # temperature and surface pressure, following:
    # https://confluence.ecmwf.int/pages/viewpage.action?pageId=171411214
    # Input: temperature (temp) with units K and surface pressure (surfp) with
    # units Pa.
    Rdry = 287.0597
    Rvap = 461.5250
    # Constants for Teten’s formula: for saturation over water:
    # a1 = 611.21 Pa, a3 = 17.502 and a4 = 32.19 K, at T0  = 273.16 K.
    a1 = 611.21
    a3 = 17.502
    a4 = 32.19
    T0 = 273.16

    # Calculation of saturation water vapour pressure (sat_wvp) from Teten's
    # formula/
    sat_wvp: xr.DataArray = a1 * np.exp(a3 * (temp - T0) / (temp - a4))  # type: ignore

    # Calculation of saturation specific humidity at 2m qsat  (equal to huss)
    # with units g/kg.
    qsat: xr.DataArray = (
        (Rdry / Rvap) * sat_wvp / (surfp - ((1 - Rdry / Rvap) * sat_wvp)) * 1000.0
    )

    # Reset axes, which were dropped during calculation
    qsat.attrs["units"] = "g/kg"
    qsat.attrs["id"] = "QREFHT"
    qsat.attrs["long_name"] = "Specific Humidity"

    return qsat


def w_convert_q(var: xr.DataArray):
    if var.attrs["units"] == "mol/mol":
        var = (
            var * 18.0 / 28.97 * 1000.0
        )  # convert from volume mixing ratio to mass mixing ratio in units g/kg
        var.attrs["units"] = "g/kg"
        var.attrs["long_name"] = "H2OLNZ (radiation)"
    return var


def molec_convert_units(vars: List[xr.DataArray], molar_weight: float) -> xr.DataArray:
    """Sum the list of variables and convert the molecular units.

    Parameters
    ----------
    vars : List[xr.DataArray]
        The list of variables.
    molar_weight : float
        The molar weight to use for the conversion.

    Returns
    -------
    xr.DataArray
        The final result.
    """
    result = sum_vars(vars)

    # Convert molec/cm2/s to kg/m2/s
    if result.attrs["units"] == "molec/cm2/s":
        result = result / AVOGADRO_CONST * molar_weight * 10.0
        result.attrs["units"] = "kg/m2/s"

    return result


def a_num_sum(a1: xr.DataArray, a2: xr.DataArray, a3: xr.DataArray):
    # Calculate: total aerosol number concentration (#/cm3)

    with xr.set_options(keep_attrs=True):
        var = (a1 + a2 + a3) * AIR_DENS / 1e6
    var.name = "a_num"
    var["units"] = "/cm3"
    var["long_name"] = "aerosol number concentration"
    return var


def so4_mass_sum(a1: xr.DataArray, a2: xr.DataArray):
    # Calculate: SO4 mass conc. (ng/m3) (< 1um)
    with xr.set_options(keep_attrs=True):
        var = (a1 + a2) * AIR_DENS * 1e9
    var.name = "so4_mass"
    var.units = "\u03bcg/m3"
    var.long_name = "SO4 mass conc."
    return var


def qflx_convert_to_lhflx(
    qflx: xr.DataArray,
    precc: xr.DataArray,
    precl: xr.DataArray,
    precsc: xr.DataArray,
    precsl: xr.DataArray,
):
    # A more precise formula to close atmospheric energy budget:
    # LHFLX is modified to account for the latent energy of frozen precipitation.
    # LHFLX = (Lv+Lf)*QFLX - Lf*1.e3*(PRECC+PRECL-PRECSC-PRECSL)
    # Constants, from AMWG diagnostics
    Lv = 2.501e6
    Lf = 3.337e5
    var = (Lv + Lf) * qflx - Lf * 1.0e3 * (precc + precl - precsc - precsl)
    var.attrs["units"] = "W/m2"
    var.attrs["long_name"] = "Surface latent heat flux"
    return var


def qflx_convert_to_lhflx_approxi(var: xr.DataArray):
    # QFLX units: kg/((m^2)*s)
    # Multiply by the latent heat of condensation/vaporization (in J/kg)
    # kg/((m^2)*s) * J/kg = J/((m^2)*s) = (W*s)/((m^2)*s) = W/(m^2)
    with xr.set_options(keep_attrs=True):
        new_var = var * 2.5e6

    new_var.name = "LHFLX"
    new_var.attrs["units"] = "W/m2"
    new_var.attrs["long_name"] = "Surface latent heat flux"

    return new_var


def pminuse_1(
    precc: xr.DataArray, precl: xr.DataArray, qflx: xr.DataArray
) -> xr.DataArray:
    var_prect = prect(precc, precl)
    var_qflx = qflxconvert_units(qflx)

    with xr.set_options(keep_attrs=True):
        var = var_prect - var_qflx

    var_final = pminuse_convert_units(var)

    return var_final


def pminuse_2(pr: xr.DataArray, evspsbl: xr.DataArray) -> xr.DataArray:
    with xr.set_options(keep_attrs=True):
        var = pr + evspsbl

    var_final = pminuse_convert_units(var)

    return var_final


def pminuse_3(pr: xr.DataArray, evspsbl: xr.DataArray) -> xr.DataArray:
    with xr.set_options(keep_attrs=True):
        var = pr - evspsbl

    var_final = pminuse_convert_units(var)

    return var_final


def pminuse_convert_units(var: xr.DataArray):
    units = var.attrs.get("units")

    matching_units = ["kg/m2/s", "kg m-2 s-1", "kg/s/m^2"]
    if units in matching_units:
        # need to find a solution for units not included in udunits
        # var = convert_units( var, 'kg/m2/s' )
        var = var * 3600.0 * 24  # convert to mm/day

    var.attrs["units"] = "mm/day"
    var.attrs["long_name"] = "precip. flux - evap. flux"
    return var


def prect(precc: xr.DataArray, precl: xr.DataArray):
    """Total precipitation flux = convective + large-scale"""
    with xr.set_options(keep_attrs=True):
        var = precc + precl

    var = convert_units(var, "mm/day")
    var.name = "PRECT"
    var.attrs["long_name"] = "Total precipitation rate (convective + large-scale)"
    return var


def prect_frac(precc: xr.DataArray, precl: xr.DataArray):
    """convective precipitation fraction = convective /(convective + large-scale)"""
    with xr.set_options(keep_attrs=True):
        var = precc / (precc + precl) * 100.0

    var.attrs["units"] = "%"
    var.attrs["long_name"] = "convective precipitation fraction"

    return var


def precst(precc: xr.DataArray, precl: xr.DataArray):
    """Total precipitation flux = convective + large-scale"""
    with xr.set_options(keep_attrs=True):
        var = precc + precl

    var = convert_units(var, "mm/day")
    var.name = "PRECST"
    var.attrs["long_name"] = "Total snowfall flux (convective + large-scale)"
    return var


def tref_range(tmax: xr.DataArray, tmin: xr.DataArray):
    """TREF daily range = TREFMXAV - TREFMNAV"""
    var = tmax - tmin
    var.name = "TREF_range"
    var.attrs["units"] = "K"
    var.attrs["long_name"] = "Surface Temperature Daily Range"
    return var


def tauxy(taux: xr.DataArray, tauy: xr.DataArray):
    """tauxy = (taux^2 + tauy^2)sqrt"""
    with xr.set_options(keep_attrs=True):
        var = (taux**2 + tauy**2) ** 0.5

    var = convert_units(var, "N/m^2")
    var.name = "TAUXY"
    var.attrs["long_name"] = "Total surface wind stress"
    return var


def fp_uptake(a: xr.DataArray, b: xr.DataArray):
    """plant uptake of soil mineral N"""
    var = a / b
    var.name = "FP_UPTAKE"
    var.attrs["units"] = "dimensionless"
    var.attrs["long_name"] = "Plant uptake of soil mineral N"
    return var


def albedo(rsdt: xr.DataArray, rsut: xr.DataArray):
    """TOA (top-of-atmosphere) albedo, rsut / rsdt, unit is nondimension"""
    var = rsut / rsdt
    var.name = "ALBEDO"
    var.attrs["units"] = "dimensionless"
    var.attrs["long_name"] = "TOA albedo"
    return var


def albedoc(rsdt: xr.DataArray, rsutcs: xr.DataArray):
    """TOA (top-of-atmosphere) albedo clear-sky, rsutcs / rsdt, unit is nondimension"""
    var = rsutcs / rsdt
    var.name = "ALBEDOC"
    var.attrs["units"] = "dimensionless"
    var.attrs["long_name"] = "TOA albedo clear-sky"

    var = _replace_inf_with_nan(var)
    return var


def albedo_srf(rsds: xr.DataArray, rsus: xr.DataArray):
    """Surface albedo, rsus / rsds, unit is nondimension"""
    var = rsus / rsds
    var.name = "ALBEDOC_SRF"
    var.attrs["units"] = "dimensionless"
    var.attrs["long_name"] = "Surface albedo"
    return var


def rst(rsdt: xr.DataArray, rsut: xr.DataArray):
    """TOA (top-of-atmosphere) net shortwave flux"""
    with xr.set_options(keep_attrs=True):
        var = rsdt - rsut

    var.name = "FSNTOA"
    var.attrs["long_name"] = "TOA net shortwave flux"
    return var


def rstcs(rsdt: xr.DataArray, rsutcs: xr.DataArray):
    """TOA (top-of-atmosphere) net shortwave flux clear-sky"""
    with xr.set_options(keep_attrs=True):
        var = rsdt - rsutcs

    var.name = "FSNTOAC"
    var.attrs["long_name"] = "TOA net shortwave flux clear-sky"
    return var


def swcfsrf(fsns: xr.DataArray, fsnsc: xr.DataArray):
    """Surface shortwave cloud forcing"""
    with xr.set_options(keep_attrs=True):
        var = fsns - fsnsc

    var.name = "SCWFSRF"
    var.attrs["long_name"] = "Surface shortwave cloud forcing"
    return var


def lwcfsrf(flns: xr.DataArray, flnsc: xr.DataArray):
    """Surface longwave cloud forcing, for ACME model, upward is postitive for LW , for ceres, downward is postive for both LW and SW"""
    with xr.set_options(keep_attrs=True):
        var = -(flns - flnsc)

    var.name = "LCWFSRF"
    var.attrs["long_name"] = "Surface longwave cloud forcing"
    return var


def swcf(fsntoa: xr.DataArray, fsntoac: xr.DataArray):
    """TOA shortwave cloud forcing"""
    with xr.set_options(keep_attrs=True):
        var = fsntoa - fsntoac

    var.name = "SWCF"
    var.attrs["long_name"] = "TOA shortwave cloud forcing"
    return var


def lwcf(flntoa: xr.DataArray, flntoac: xr.DataArray):
    """TOA longwave cloud forcing"""
    with xr.set_options(keep_attrs=True):
        var = flntoa - flntoac

    var.name = "LWCF"
    var.attrs["long_name"] = "TOA longwave cloud forcing"
    return var


def netcf2(swcf: xr.DataArray, lwcf: xr.DataArray):
    """TOA net cloud forcing"""
    with xr.set_options(keep_attrs=True):
        var = swcf + lwcf

    var.name = "NETCF"
    var.attrs["long_name"] = "TOA net cloud forcing"
    return var


def netcf4(
    fsntoa: xr.DataArray,
    fsntoac: xr.DataArray,
    flntoa: xr.DataArray,
    flntoac: xr.DataArray,
):
    """TOA net cloud forcing"""
    with xr.set_options(keep_attrs=True):
        var = fsntoa - fsntoac + flntoa - flntoac

    var.name = "NETCF"
    var.attrs["long_name"] = "TOA net cloud forcing"
    return var


def netcf2srf(swcf: xr.DataArray, lwcf: xr.DataArray):
    """Surface net cloud forcing"""
    with xr.set_options(keep_attrs=True):
        var = swcf + lwcf

    var.name = "NETCF_SRF"
    var.attrs["long_name"] = "Surface net cloud forcing"
    return var


def netcf4srf(
    fsntoa: xr.DataArray,
    fsntoac: xr.DataArray,
    flntoa: xr.DataArray,
    flntoac: xr.DataArray,
):
    """Surface net cloud forcing"""
    with xr.set_options(keep_attrs=True):
        var = fsntoa - fsntoac + flntoa - flntoac

    var.name = "NETCF4SRF"
    var.attrs["long_name"] = "Surface net cloud forcing"
    return var


def fldsc(ts: xr.DataArray, flnsc: xr.DataArray):
    """Clearsky Surf LW downwelling flux"""
    with xr.set_options(keep_attrs=True):
        var = 5.67e-8 * ts**4 - flnsc

    var.name = "FLDSC"
    var.attrs["units"] = "W/m2"
    var.attrs["long_name"] = "Clearsky Surf LW downwelling flux"
    return var


def restom(fsnt: xr.DataArray, flnt: xr.DataArray):
    """TOM(top of model) Radiative flux"""
    with xr.set_options(keep_attrs=True):
        var = fsnt - flnt

    var.name = "RESTOM"
    var.attrs["long_name"] = "TOM(top of model) Radiative flux"
    return var


def restom3(swdn: xr.DataArray, swup: xr.DataArray, lwup: xr.DataArray):
    """TOM(top of model) Radiative flux"""
    with xr.set_options(keep_attrs=True):
        var = swdn - swup - lwup

    var.long_name = "TOM(top of model) Radiative flux"

    return var


def restoa(fsnt: xr.DataArray, flnt: xr.DataArray):
    """TOA(top of atmosphere) Radiative flux"""
    with xr.set_options(keep_attrs=True):
        var = fsnt - flnt

    var.name = "RESTOA"
    var.attrs["long_name"] = "TOA(top of atmosphere) Radiative flux"
    return var


def flus(flds: xr.DataArray, flns: xr.DataArray):
    """Surface Upwelling LW Radiative flux"""
    with xr.set_options(keep_attrs=True):
        var = flns + flds

    var.name = "FLUS"
    var.attrs["long_name"] = "Upwelling longwave flux at surface"
    return var


def fsus(fsds: xr.DataArray, fsns: xr.DataArray):
    """Surface Up-welling SW Radiative flux"""
    with xr.set_options(keep_attrs=True):
        var = fsds - fsns

    var.name = "FSUS"
    var.attrs["long_name"] = "Upwelling shortwave flux at surface"
    return var


def netsw(rsds: xr.DataArray, rsus: xr.DataArray):
    """Surface SW Radiative flux"""
    with xr.set_options(keep_attrs=True):
        var = rsds - rsus

    var.name = "FSNS"
    var.attrs["long_name"] = "Surface SW Radiative flux"
    return var


def netlw(rlds: xr.DataArray, rlus: xr.DataArray):
    """Surface LW Radiative flux"""
    with xr.set_options(keep_attrs=True):
        var = -(rlds - rlus)

    var.name = "NET_FLUX_SRF"
    var.attrs["long_name"] = "Surface LW Radiative flux"
    return var


def netflux4(
    fsns: xr.DataArray, flns: xr.DataArray, lhflx: xr.DataArray, shflx: xr.DataArray
):
    """Surface Net flux"""
    with xr.set_options(keep_attrs=True):
        var = fsns - flns - lhflx - shflx

    var.name = "NET_FLUX_SRF"
    var.attrs["long_name"] = "Surface Net flux"
    return var


def netflux6(
    rsds: xr.DataArray,
    rsus: xr.DataArray,
    rlds: xr.DataArray,
    rlus: xr.DataArray,
    hfls: xr.DataArray,
    hfss: xr.DataArray,
):
    """Surface Net flux"""
    with xr.set_options(keep_attrs=True):
        var = rsds - rsus + (rlds - rlus) - hfls - hfss

    var.name = "NET_FLUX_SRF"
    var.attrs["long_name"] = "Surface Net flux"
    return var


def _replace_inf_with_nan(var: xr.DataArray) -> xr.DataArray:
    """Replaces `np.inf` with `np.nan`.

    This function is useful for division arithmetic where divide by zero might
    occur. For example, in `albedoc()`, there is reference file where `rsdt`
    contains 0s. This  function divides `rsutcs / rsdt`. `rsdt` contains 0s,
    which results in divide by zeros.

    - CDAT/cdms2 replaces these values with the floats from `rsutcs`, but they
    are masked so they will be outputted as `np.nan`.
    - Xarray and NumPy replaces these values with `np.inf`. We replace `np.inf`
    with `np.nan` to maintain the behavior of the CDAT-based code.
    - Related ref file: '/global/cfs/cdirs/e3sm/diagnostics/observations/Atm
    /climatology/ceres_ebaf_toa_v4.1/ceres_ebaf_toa_v4.1_JJA_200106_201808_climo.nc'

    Parameters
    ----------
    var : xr.DataArray
        The variable containing `np.inf`.

    Returns
    -------
    xr.DataArray
        The variable with `np.inf` replaced with `np.nan`.
    """
    var_new = xr.where(var != np.inf, var, np.nan, keep_attrs=True)

    return var_new


def aero_burden_fxn(var: xr.DataArray) -> xr.DataArray:
    """Scale the aerosol burden by 1e6.

    Parameters
    ----------
    var : xr.DataArray
        The input burden in kg/m2.

    Returns
    -------
    xr.DataArray
        The output burden in 1e-6 kg/m2.
    """
    with xr.set_options(keep_attrs=True):
        burden = var * 1e6

    burden.attrs["units"] = "1e-6 kg/m2"

    return burden


def aero_mass_fxn(var: xr.DataArray) -> xr.DataArray:
    """Scale the given mass by 1e12.

    Parameters
    ----------
    var : xr.DataArray
        The input mass in kg/kg.

    Returns
    -------
    xr.DataArray
        The aerosol mass concentration in 1e-12 kg/kg units.
    """
    with xr.set_options(keep_attrs=True):
        mass = var * 1e12

    mass.attrs["units"] = "1e-12 kg/kg"

    return mass


def incldtop_cdnc(cdnc: xr.DataArray, lcc: xr.DataArray) -> xr.DataArray:
    """Return the in-cloud cloud droplet number concentration at cloud top.

    Parameters
    ----------
    cdnc : xr.DataArray
        Cloud droplet number concentration in 1/m3.
    lcc : xr.DataArray
        Liquid cloud fraction.

    Returns
    -------
    xr.DataArray
        In-cloud cdnc at cloud top in 1/cm3.
    """
    with xr.set_options(keep_attrs=True):
        var = cdnc * 1e-6 / lcc

    var.attrs["units"] = "1/cm3"
    var.attrs["long_name"] = "In-cloud-top CDNC"

    return var


def cldtop_cdnc(cdnc: xr.DataArray) -> xr.DataArray:
    """Return the in-grid cloud droplet number concentration at cloud top.

    Parameters
    ----------
    cdnc : xr.DataArray
        Cloud droplet number concentration in 1/m3.

    Returns
    -------
    xr.DataArray
        In-grid cdnc at cloud top in 1/cm3.
    """
    with xr.set_options(keep_attrs=True):
        var = cdnc * 1e-6

    var.attrs["units"] = "1/cm3"
    var.attrs["long_name"] = "In-grid cloud-top CDNC"

    return var


def incldtop_icnc(icnc: xr.DataArray, icc: xr.DataArray) -> xr.DataArray:
    """Return the in-cloud ice crystal number concentration at cloud top.

    Parameters
    ----------
    icnc : xr.DataArray
        Ice crystal number concentration in 1/m3.
    icc : xr.DataArray
        ice cloud fraction.

    Returns
    -------
    xr.DataArray
        In-cloud cdnc at cloud top in 1/cm3.
    """
    with xr.set_options(keep_attrs=True):
        var = icnc * 1e-6 / icc

    var.attrs["units"] = "1/cm3"
    var.attrs["long_name"] = "In-cloud-top ICNC"

    return var


def cldtop_icnc(icnc: xr.DataArray) -> xr.DataArray:
    """Return the in-grid ice crystal number concentration at cloud top.

    Parameters
    ----------
    icnc : xr.DataArray
        Cloud crystal number concentration in 1/m3.

    Returns
    -------
    xr.DataArray
        In-grid icnc at cloud top in 1/cm3.
    """
    with xr.set_options(keep_attrs=True):
        var = icnc * 1e-6

    var.attrs["units"] = "1/cm3"
    var.attrs["long_name"] = "In-grid cloud-top ICNC"

    return var


def incld_lwp(lwp: xr.DataArray, lcc: xr.DataArray) -> xr.DataArray:
    """Return the in-cloud liquid water path (LWP).

    Parameters
    ----------
    lwp : xr.DataArray
        Liquid water path in kg/m2.
    lcc : xr.DataArray
        Liquid cloud fraction.

    Returns
    -------
    xr.DataArray
        In-cloud liquid water path in g/cm3.
    """
    with xr.set_options(keep_attrs=True):
        var = 1e3 * lwp / lcc

    var.attrs["units"] = "g/cm3"
    var.attrs["long_name"] = "In-cloud LWP"

    return var


def cld_lwp(lwp: xr.DataArray) -> xr.DataArray:
    """Return the grid-mean-cloud LWP in g/cm3.

    Parameters
    ----------
    lwp : xr.DataArray
        Liquid Water Path (LWP) value.

    Returns
    -------
    xr.DataArray
        Grid-mean-cloud LWP in g/cm3.
    """
    with xr.set_options(keep_attrs=True):
        var = 1e3 * lwp

    var.attrs["units"] = "g/cm3"
    var.attrs["long_name"] = "In-grid LWP"

    return var


def incld_iwp(iwp: xr.DataArray, icc: xr.DataArray) -> xr.DataArray:
    """Return the in-cloud ice water path (IWP).

    Parameters
    ----------
    iwp : xr.DataArray
        Ice water path in kg/m2.
    icc : xr.DataArray
        Ice cloud fraction.

    Returns
    -------
    xr.DataArray
        In-cloud IWP in g/cm3.
    """
    with xr.set_options(keep_attrs=True):
        var = 1e3 * iwp / icc

    var.attrs["units"] = "g/cm3"
    var.attrs["long_name"] = "In-cloud IWP"

    return var


def cld_iwp(iwp: xr.DataArray) -> xr.DataArray:
    """Return the in-grid ice water path (IWP).

    Parameters
    ----------
    iwp : xr.DataArray
        Ice water path in kg/m2.

    Returns
    -------
    xr.DataArray
        In-grid IWP in g/cm3.
    """
    with xr.set_options(keep_attrs=True):
        var = 1e3 * iwp

    var.attrs["units"] = "g/cm3"
    var.attrs["long_name"] = "In-grid IWP"

    return var


def erf_tot(fsnt: xr.DataArray, flnt: xr.DataArray) -> xr.DataArray:
    """
    Calculate the total effective radiative forcing (ERFtot).

    Parameters
    ----------
    fsnt : xr.DataArray
        The incoming sw radiation at the top of the atmosphere.
    flnt : xr.DataArray
        The outgoing lw radiation at the top of the atmosphere.

    Returns
    -------
    xr.DataArray
        The ERFtot which represents the total erf.

    See Ghan 2013 for derivation of ERF decomposition: https://doi.org/10.5194/acp-13-9971-2013
    """
    with xr.set_options(keep_attrs=True):
        var = fsnt - flnt

    var.attrs["units"] = "W/m2"
    var.attrs["long_name"] = "ERFtot: total effect"
    return var


def erf_ari(
    fsnt: xr.DataArray, flnt: xr.DataArray, fsnt_d1: xr.DataArray, flnt_d1: xr.DataArray
) -> xr.DataArray:
    """
    Calculate aerosol--radiation interactions (ARI) part of effective radiative forcing (ERF).

    Parameters
    ----------
    fsnt : xr.DataArray
        Net solar flux at the top of the atmosphere.
    flnt : xr.DataArray
        Net longwave flux at the top of the atmosphere.
    fsnt_d1 : xr.DataArray
        fsnt without aerosols.
    flnt_d1 : xr.DataArray
        flnt without aerosols.

    Returns
    -------
    xr.DataArray
        ERFari (aka, direct effect) in W/m2.

    See Ghan 2013 for derivation of ERF decomposition: https://doi.org/10.5194/acp-13-9971-2013
    """
    with xr.set_options(keep_attrs=True):
        var = (fsnt - flnt) - (fsnt_d1 - flnt_d1)

    var.attrs["units"] = "W/m2"
    var.attrs["long_name"] = "ERFari: direct effect"

    return var


def erf_aci(
    fsnt_d1: xr.DataArray,
    flnt_d1: xr.DataArray,
    fsntc_d1: xr.DataArray,
    flntc_d1: xr.DataArray,
) -> xr.DataArray:
    """
    Calculate aerosol--cloud interactions (ACI) part of effectie radiative forcing (ERF)

    Parameters
    ----------
    fsnt_d1 : xr.DataArray
        Downward shortwave radiation toa without aerosols.
    flnt_d1 : xr.DataArray
        Upward longwave radiation toa without aerosols.
    fsntc_d1 : xr.DataArray
        fsnt_d1 without clouds.
    flntc_d1 : xr.DataArray
        flnt_d1 without clouds.

    Returns
    -------
    xr.DataArray
        ERFaci (aka, indirect effect) in W/m2.

    Notes
    -----
    See Ghan 2013 for derivation of ERF decomposition: https://doi.org/10.5194/acp-13-9971-2013
    """
    with xr.set_options(keep_attrs=True):
        var = (fsnt_d1 - flnt_d1) - (fsntc_d1 - flntc_d1)

    var.attrs["units"] = "W/m2"
    var.attrs["long_name"] = "ERFaci: indirect effect"

    return var


def erf_res(fsntc_d1: xr.DataArray, flntc_d1: xr.DataArray) -> xr.DataArray:
    """
    Calculate the residual effect (RES) part of effective radiative forcin g.

    Parameters
    ----------
    fsntc_d1 : xr.DataArray
        Downward solar radiation at the top of the atmosphere with neither
        clouds nor aerosols.
    flntc_d1 : xr.DataArray
        Upward longwave radiation at the top of the atmosphere with neither
        clouds nor aerosols.

    Returns
    -------
    xr.DataArray
        ERFres (aka, surface effect) in W/m2.

    Notes
    -----
    See Ghan 2013 for derivation of ERF decomposition: https://doi.org/10.5194/acp-13-9971-2013
    """
    with xr.set_options(keep_attrs=True):
        var = fsntc_d1 - flntc_d1

    var.attrs["units"] = "W/m2"
    var.attrs["long_name"] = "ERFres: residual effect"

    return var

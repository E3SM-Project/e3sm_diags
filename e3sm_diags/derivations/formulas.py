"""This module defines formula functions used for deriving variables.

The function arguments usually accept variables represented by `xr.DataArray`.
NOTE: If a function involves arithmetic between two or more `xr.DataArray`,
the arithmetic should be wrapped with `with xr.set_options(keep_attrs=True)`
to keep attributes on the resultant `xr.DataArray`.
"""
import xarray as xr

from e3sm_diags.derivations.utils import convert_units

AVOGADRO_CONST = 6.022e23


def qflxconvert_units(var: xr.DataArray):
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


def w_convert_q(var: xr.DataArray):
    if var.attrs["units"] == "mol/mol":
        var = (
            var * 18.0 / 28.97 * 1000.0
        )  # convert from volume mixing ratio to mass mixing ratio in units g/kg
        var.attrs["units"] = "g/kg"
        var.attrs["long_name"] = "H2OLNZ (radiation)"
    return var


def molec_convert_units(var: xr.DataArray, molar_weight: float):
    # Convert molec/cm2/s to kg/m2/s
    if var.attrs["units"] == "molec/cm2/s":
        var = var / AVOGADRO_CONST * molar_weight * 10.0
        var.attrs["units"] == "kg/m2/s"
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
    return new_var


def pminuse_convert_units(var: xr.DataArray):
    if (
        var.attrs["units"] == "kg/m2/s"
        or var.attrs["units"] == "kg m-2 s-1"
        or var.attrs["units"] == "kg/s/m^2"
    ):
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

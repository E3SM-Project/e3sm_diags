"""This module defines formula functions used for deriving variables."""
import xarray as xr

from e3sm_diags.derivations.utils import convert_units

AVOGADRO_CONST = 6.022e23


def qflxconvert_units(var):
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


def w_convert_q(var):
    if var.attrs["units"] == "mol/mol":
        var = (
            var * 18.0 / 28.97 * 1000.0
        )  # convert from volume mixing ratio to mass mixing ratio in units g/kg
        var.attrs["units"] = "g/kg"
        var.attrs["long_name"] = "H2OLNZ (radiation)"
    return var


def molec_convert_units(var, molar_weight):
    # Convert molec/cm2/s to kg/m2/s
    if var.attrs["units"] == "molec/cm2/s":
        var = var / AVOGADRO_CONST * molar_weight * 10.0
        var.attrs["units"] == "kg/m2/s"
    return var


def qflx_convert_to_lhflx(qflx, precc, precl, precsc, precsl):
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


def qflx_convert_to_lhflx_approxi(var):
    # QFLX units: kg/((m^2)*s)
    # Multiply by the latent heat of condensation/vaporization (in J/kg)
    # kg/((m^2)*s) * J/kg = J/((m^2)*s) = (W*s)/((m^2)*s) = W/(m^2)
    with xr.set_options(keep_attrs=True):
        new_var = var * 2.5e6

    new_var.name = "LHFLX"
    return new_var


def pminuse_convert_units(var):
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


def prect(precc, precl):
    """Total precipitation flux = convective + large-scale"""
    with xr.set_options(keep_attrs=True):
        var = precc + precl

    var = convert_units(var, "mm/day")
    var.name = "PRECT"
    var.attrs["long_name"] = "Total precipitation rate (convective + large-scale)"
    return var


def precst(precc, precl):
    """Total precipitation flux = convective + large-scale"""
    with xr.set_options(keep_attrs=True):
        var = precc + precl

    var = convert_units(var, "mm/day")
    var.name = "PRECST"
    var.attrs["long_name"] = "Total snowfall flux (convective + large-scale)"
    return var


def tref_range(tmax, tmin):
    """TREF daily range = TREFMXAV - TREFMNAV"""
    var = tmax - tmin
    var.name = "TREF_range"
    var.attrs["units"] = "K"
    var.attrs["long_name"] = "Surface Temperature Daily Range"
    return var


def tauxy(taux, tauy):
    """tauxy = (taux^2 + tauy^2)sqrt"""
    with xr.set_options(keep_attrs=True):
        var = (taux**2 + tauy**2) ** 0.5

    var = convert_units(var, "N/m^2")
    var.name = "TAUXY"
    var.attrs["long_name"] = "Total surface wind stress"
    return var


def fp_uptake(a, b):
    """plant uptake of soil mineral N"""
    var = a / b
    var.name = "FP_UPTAKE"
    var.attrs["units"] = "dimensionless"
    var.attrs["long_name"] = "Plant uptake of soil mineral N"
    return var


def albedo(rsdt, rsut):
    """TOA (top-of-atmosphere) albedo, rsut / rsdt, unit is nondimension"""
    var = rsut / rsdt
    var.name = "ALBEDO"
    var.attrs["units"] = "dimensionless"
    var.attrs["long_name"] = "TOA albedo"
    return var


def albedoc(rsdt, rsutcs):
    """TOA (top-of-atmosphere) albedo clear-sky, rsutcs / rsdt, unit is nondimension"""
    var = rsutcs / rsdt
    var.name = "ALBEDOC"
    var.attrs["units"] = "dimensionless"
    var.attrs["long_name"] = "TOA albedo clear-sky"
    return var


def albedo_srf(rsds, rsus):
    """Surface albedo, rsus / rsds, unit is nondimension"""
    var = rsus / rsds
    var.name = "ALBEDOC_SRF"
    var.attrs["units"] = "dimensionless"
    var.attrs["long_name"] = "Surface albedo"
    return var


def rst(rsdt, rsut):
    """TOA (top-of-atmosphere) net shortwave flux"""
    with xr.set_options(keep_attrs=True):
        var = rsdt - rsut

    var.name = "FSNTOA"
    var.attrs["long_name"] = "TOA net shortwave flux"
    return var


def rstcs(rsdt, rsutcs):
    """TOA (top-of-atmosphere) net shortwave flux clear-sky"""
    with xr.set_options(keep_attrs=True):
        var = rsdt - rsutcs

    var.name = "FSNTOAC"
    var.attrs["long_name"] = "TOA net shortwave flux clear-sky"
    return var


def swcfsrf(fsns, fsnsc):
    """Surface shortwave cloud forcing"""
    with xr.set_options(keep_attrs=True):
        var = fsns - fsnsc

    var.name = "SCWFSRF"
    var.attrs["long_name"] = "Surface shortwave cloud forcing"
    return var


def lwcfsrf(flns, flnsc):
    """Surface longwave cloud forcing, for ACME model, upward is postitive for LW , for ceres, downward is postive for both LW and SW"""
    with xr.set_options(keep_attrs=True):
        var = -(flns - flnsc)

    var.name = "LCWFSRF"
    var.attrs["long_name"] = "Surface longwave cloud forcing"
    return var


def swcf(fsntoa, fsntoac):
    """TOA shortwave cloud forcing"""
    with xr.set_options(keep_attrs=True):
        var = fsntoa - fsntoac

    var.name = "SWCF"
    var.attrs["long_name"] = "TOA shortwave cloud forcing"
    return var


def lwcf(flntoa, flntoac):
    """TOA longwave cloud forcing"""
    with xr.set_options(keep_attrs=True):
        var = flntoa - flntoac

    var.name = "LWCF"
    var.attrs["long_name"] = "TOA longwave cloud forcing"
    return var


def netcf2(swcf, lwcf):
    """TOA net cloud forcing"""
    with xr.set_options(keep_attrs=True):
        var = swcf + lwcf

    var.name = "NETCF"
    var.attrs["long_name"] = "TOA net cloud forcing"
    return var


def netcf4(fsntoa, fsntoac, flntoa, flntoac):
    """TOA net cloud forcing"""
    with xr.set_options(keep_attrs=True):
        var = fsntoa - fsntoac + flntoa - flntoac

    var.name = "NETCF"
    var.attrs["long_name"] = "TOA net cloud forcing"
    return var


def netcf2srf(swcf, lwcf):
    """Surface net cloud forcing"""
    with xr.set_options(keep_attrs=True):
        var = swcf + lwcf

    var.name = "NETCF_SRF"
    var.attrs["long_name"] = "Surface net cloud forcing"
    return var


def netcf4srf(fsntoa, fsntoac, flntoa, flntoac):
    """Surface net cloud forcing"""
    with xr.set_options(keep_attrs=True):
        var = fsntoa - fsntoac + flntoa - flntoac

    var.name = "NETCF4SRF"
    var.attrs["long_name"] = "Surface net cloud forcing"
    return var


def fldsc(ts, flnsc):
    """Clearsky Surf LW downwelling flux"""
    with xr.set_options(keep_attrs=True):
        var = 5.67e-8 * ts**4 - flnsc

    var.name = "FLDSC"
    var.attrs["units"] = "W/m2"
    var.attrs["long_name"] = "Clearsky Surf LW downwelling flux"
    return var


def restom(fsnt, flnt):
    """TOM(top of model) Radiative flux"""
    with xr.set_options(keep_attrs=True):
        var = fsnt - flnt

    var.name = "RESTOM"
    var.attrs["long_name"] = "TOM(top of model) Radiative flux"
    return var


def restoa(fsnt, flnt):
    """TOA(top of atmosphere) Radiative flux"""
    with xr.set_options(keep_attrs=True):
        var = fsnt - flnt

    var.name = "RESTOA"
    var.attrs["long_name"] = "TOA(top of atmosphere) Radiative flux"
    return var


def flus(flds, flns):
    """Surface Upwelling LW Radiative flux"""
    with xr.set_options(keep_attrs=True):
        var = flns + flds

    var.name = "FLUS"
    var.attrs["long_name"] = "Upwelling longwave flux at surface"
    return var


def fsus(fsds, fsns):
    """Surface Up-welling SW Radiative flux"""
    with xr.set_options(keep_attrs=True):
        var = fsds - fsns

    var.name = "FSUS"
    var.attrs["long_name"] = "Upwelling shortwave flux at surface"
    return var


def netsw(rsds, rsus):
    """Surface SW Radiative flux"""
    with xr.set_options(keep_attrs=True):
        var = rsds - rsus

    var.name = "FSNS"
    var.attrs["long_name"] = "Surface SW Radiative flux"
    return var


def netlw(rlds, rlus):
    """Surface LW Radiative flux"""
    with xr.set_options(keep_attrs=True):
        var = -(rlds - rlus)

    var.name = "NET_FLUX_SRF"
    var.attrs["long_name"] = "Surface LW Radiative flux"
    return var


def netflux4(fsns, flns, lhflx, shflx):
    """Surface Net flux"""
    with xr.set_options(keep_attrs=True):
        var = fsns - flns - lhflx - shflx

    var.name = "NET_FLUX_SRF"
    var.attrs["long_name"] = "Surface Net flux"
    return var


def netflux6(rsds, rsus, rlds, rlus, hfls, hfss):
    """Surface Net flux"""
    with xr.set_options(keep_attrs=True):
        var = rsds - rsus + (rlds - rlus) - hfls - hfss

    var.name = "NET_FLUX_SRF"
    var.attrs["long_name"] = "Surface Net flux"
    return var

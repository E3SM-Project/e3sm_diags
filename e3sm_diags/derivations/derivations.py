"""This module stores definitions for derived variables.

`DERIVED_VARIABLES` is a dictionary that stores definitions for derived
variables. The driver uses the Dataset class to search for available variable
keys and attempts to map them to a formula function to calculate a derived
variable.

For example to derive 'PRECT':
  1. In `DERIVED_VARIABLE` there is an entry for 'PRECT'.
  2. The netCDF file does not have a 'PRECT' variable, but has the 'PRECC'
     and 'PRECT' variables.
  3. 'PRECC' and 'PRECL' are used to derive `PRECT` by passing the
     data for these variables to the formula function 'prect()'.
"""
from collections import OrderedDict
from typing import Callable, Dict, Tuple, Union

from e3sm_diags.derivations.formulas import (
    a_num_sum,
    aero_burden_fxn,
    aero_mass_fxn,
    albedo,
    albedo_srf,
    albedoc,
    cld_iwp,
    cld_lwp,
    cldtop_cdnc,
    cldtop_icnc,
    erf_aci,
    erf_ari,
    erf_res,
    erf_tot,
    fldsc,
    flus,
    fp_uptake,
    fsus,
    incld_iwp,
    incld_lwp,
    incldtop_cdnc,
    incldtop_icnc,
    lwcf,
    lwcfsrf,
    molec_convert_units,
    netcf2,
    netcf2srf,
    netcf4,
    netcf4srf,
    netflux4,
    netflux6,
    netlw,
    netsw,
    pminuse_1,
    pminuse_2,
    pminuse_3,
    pminuse_convert_units,
    precst,
    prect,
    prect_frac,
    qflx_convert_to_lhflx,
    qflx_convert_to_lhflx_approxi,
    qflxconvert_units,
    qsat,
    restoa,
    restom,
    restom3,
    rst,
    rstcs,
    so4_mass_sum,
    sum_vars,
    swcf,
    swcfsrf,
    tauxy,
    tref_range,
    w_convert_q,
)
from e3sm_diags.derivations.formulas_cosp import (
    cosp_bin_sum,
    cosp_histogram_standardize,
)
from e3sm_diags.derivations.utils import (
    _apply_land_sea_mask,
    aplusb,
    convert_units,
    rename,
)

# A type annotation ordered dictionary that maps a tuple of source variable(s)
# to a derivation function.
DerivedVariableMap = Union[
    OrderedDict[Tuple[str, ...], Callable], Dict[Tuple[str, ...], Callable]
]

# A type annotation for a dictionary mapping the key of a derived variable
# to an ordered dictionary that maps a tuple of source variable(s) to a
# derivation function.
DerivedVariablesMap = Dict[str, DerivedVariableMap]

# A list of functions that need the target variable key as an argument to use
# for further operations. This variable is used in the Dataset class variable
# derivation method.
FUNC_NEEDS_TARGET_VAR = [cosp_bin_sum, cosp_histogram_standardize]


# TODO: Replace OrderedDict with normal dictionary and remove lambda calls
# that aren't necessary (e.g., `rename()`).
DERIVED_VARIABLES: DerivedVariablesMap = {
    "PRECT": {
        ("PRECT",): lambda pr: convert_units(rename(pr), target_units="mm/day"),
        ("pr",): lambda pr: qflxconvert_units(rename(pr)),
        ("PRECC", "PRECL"): lambda precc, precl: prect(precc, precl),
        ("sat_gauge_precip",): rename,
        ("PrecipLiqSurfMassFlux", "PrecipIceSurfMassFlux"): prect,  # EAMxx
    },
    "PRECST": {
        ("prsn",): lambda prsn: qflxconvert_units(rename(prsn)),
        ("PRECSC", "PRECSL"): precst,
        ("VapWaterPath",): lambda prw: convert_units(
            rename(prw), target_units="kg/m2"
        ),  # EAMxx
    },
    "PRECC": {
        ("PRECC",): lambda pr: convert_units(pr, target_units="mm/day"),
        ("prc",): rename,
    },
    "PRECL": {
        ("PRECL",): lambda pr: convert_units(rename(pr), target_units="mm/day"),
    },
    "PRECC_Frac": {
        ("PRECC", "PRECL"): lambda precc, precl: prect_frac(precc, precl),
    },
    # Sea Surface Temperature: Degrees C
    # Temperature of the water, not the air. Ignore land.
    "SST": OrderedDict(
        [
            # lambda sst: convert_units(rename(sst),target_units="degC")),
            (("sst",), rename),
            (
                ("TS", "OCNFRAC"),
                lambda ts, ocnfrac: _apply_land_sea_mask(
                    convert_units(ts, target_units="degC"),
                    ocnfrac,
                    lower_limit=0.9,
                ),
            ),
            (("SST",), lambda sst: convert_units(sst, target_units="degC")),
        ]
    ),
    "TMQ": OrderedDict(
        [
            (("PREH2O",), rename),
            (
                ("prw",),
                lambda prw: convert_units(rename(prw), target_units="kg/m2"),
            ),
        ]
    ),
    "SOLIN": {("rsdt",): rename, ("SW_flux_dn_at_model_top",): rename},  # EAMxx
    "ALBEDO": {
        ("ALBEDO",): rename,
        ("SOLIN", "FSNTOA"): lambda solin, fsntoa: albedo(solin, solin - fsntoa),
        ("rsdt", "rsut"): albedo,
        (
            "SW_flux_up_at_model_top",
            "SW_clrsky_flux_up_at_model_top",
        ): swcf,  # EAMxx
    },
    "ALBEDOC": OrderedDict(
        [
            (("ALBEDOC",), rename),
            (
                ("SOLIN", "FSNTOAC"),
                lambda solin, fsntoac: albedoc(solin, solin - fsntoac),
            ),
            (("rsdt", "rsutcs"), lambda rsdt, rsutcs: albedoc(rsdt, rsutcs)),
        ]
    ),
    "ALBEDO_SRF": OrderedDict(
        [
            (("ALBEDO_SRF",), rename),
            (("rsds", "rsus"), lambda rsds, rsus: albedo_srf(rsds, rsus)),
            (
                ("FSDS", "FSNS"),
                lambda fsds, fsns: albedo_srf(fsds, fsds - fsns),
            ),
        ]
    ),
    # Pay attention to the positive direction of SW and LW fluxes
    "SWCF": OrderedDict(
        [
            (("SWCF",), rename),
            (
                ("toa_net_sw_all_mon", "toa_net_sw_clr_mon"),
                lambda net_all, net_clr: swcf(net_all, net_clr),
            ),
            (
                ("toa_net_sw_all_mon", "toa_net_sw_clr_t_mon"),
                lambda net_all, net_clr: swcf(net_all, net_clr),
            ),
            (("toa_cre_sw_mon",), rename),
            (
                ("FSNTOA", "FSNTOAC"),
                lambda fsntoa, fsntoac: swcf(fsntoa, fsntoac),
            ),
            (("rsut", "rsutcs"), lambda rsutcs, rsut: swcf(rsut, rsutcs)),
        ]
    ),
    "SWCFSRF": OrderedDict(
        [
            (("SWCFSRF",), rename),
            (
                ("sfc_net_sw_all_mon", "sfc_net_sw_clr_mon"),
                lambda net_all, net_clr: swcfsrf(net_all, net_clr),
            ),
            (
                ("sfc_net_sw_all_mon", "sfc_net_sw_clr_t_mon"),
                lambda net_all, net_clr: swcfsrf(net_all, net_clr),
            ),
            (("sfc_cre_net_sw_mon",), rename),
            (("FSNS", "FSNSC"), lambda fsns, fsnsc: swcfsrf(fsns, fsnsc)),
        ]
    ),
    "LWCF": OrderedDict(
        [
            (("LWCF",), rename),
            (
                ("toa_net_lw_all_mon", "toa_net_lw_clr_mon"),
                lambda net_all, net_clr: lwcf(net_clr, net_all),
            ),
            (
                ("toa_net_lw_all_mon", "toa_net_lw_clr_t_mon"),
                lambda net_all, net_clr: lwcf(net_clr, net_all),
            ),
            (("toa_cre_lw_mon",), rename),
            (
                ("FLNTOA", "FLNTOAC"),
                lambda flntoa, flntoac: lwcf(flntoa, flntoac),
            ),
            (("rlut", "rlutcs"), lambda rlutcs, rlut: lwcf(rlut, rlutcs)),
        ]
    ),
    "LWCFSRF": OrderedDict(
        [
            (("LWCFSRF",), rename),
            (
                ("sfc_net_lw_all_mon", "sfc_net_lw_clr_mon"),
                lambda net_all, net_clr: lwcfsrf(net_clr, net_all),
            ),
            (
                ("sfc_net_lw_all_mon", "sfc_net_lw_clr_t_mon"),
                lambda net_all, net_clr: lwcfsrf(net_clr, net_all),
            ),
            (("sfc_cre_net_lw_mon",), rename),
            (("FLNS", "FLNSC"), lambda flns, flnsc: lwcfsrf(flns, flnsc)),
        ]
    ),
    "NETCF": OrderedDict(
        [
            (
                (
                    "toa_net_sw_all_mon",
                    "toa_net_sw_clr_mon",
                    "toa_net_lw_all_mon",
                    "toa_net_lw_clr_mon",
                ),
                lambda sw_all, sw_clr, lw_all, lw_clr: netcf4(
                    sw_all, sw_clr, lw_all, lw_clr
                ),
            ),
            (
                (
                    "toa_net_sw_all_mon",
                    "toa_net_sw_clr_t_mon",
                    "toa_net_lw_all_mon",
                    "toa_net_lw_clr_t_mon",
                ),
                lambda sw_all, sw_clr, lw_all, lw_clr: netcf4(
                    sw_all, sw_clr, lw_all, lw_clr
                ),
            ),
            (
                ("toa_cre_sw_mon", "toa_cre_lw_mon"),
                lambda swcf, lwcf: netcf2(swcf, lwcf),
            ),
            (("SWCF", "LWCF"), lambda swcf, lwcf: netcf2(swcf, lwcf)),
            (
                ("FSNTOA", "FSNTOAC", "FLNTOA", "FLNTOAC"),
                lambda fsntoa, fsntoac, flntoa, flntoac: netcf4(
                    fsntoa, fsntoac, flntoa, flntoac
                ),
            ),
        ]
    ),
    "NETCF_SRF": OrderedDict(
        [
            (
                (
                    "sfc_net_sw_all_mon",
                    "sfc_net_sw_clr_mon",
                    "sfc_net_lw_all_mon",
                    "sfc_net_lw_clr_mon",
                ),
                lambda sw_all, sw_clr, lw_all, lw_clr: netcf4srf(
                    sw_all, sw_clr, lw_all, lw_clr
                ),
            ),
            (
                (
                    "sfc_net_sw_all_mon",
                    "sfc_net_sw_clr_t_mon",
                    "sfc_net_lw_all_mon",
                    "sfc_net_lw_clr_t_mon",
                ),
                lambda sw_all, sw_clr, lw_all, lw_clr: netcf4srf(
                    sw_all, sw_clr, lw_all, lw_clr
                ),
            ),
            (
                ("sfc_cre_sw_mon", "sfc_cre_lw_mon"),
                lambda swcf, lwcf: netcf2srf(swcf, lwcf),
            ),
            (
                ("FSNS", "FSNSC", "FLNSC", "FLNS"),
                lambda fsns, fsnsc, flnsc, flns: netcf4srf(fsns, fsnsc, flnsc, flns),
            ),
        ]
    ),
    "FLNS": OrderedDict(
        [
            (
                ("sfc_net_lw_all_mon",),
                lambda sfc_net_lw_all_mon: -sfc_net_lw_all_mon,
            ),
            (("rlds", "rlus"), lambda rlds, rlus: netlw(rlds, rlus)),
        ]
    ),
    "FLNSC": OrderedDict(
        [
            (
                ("sfc_net_lw_clr_mon",),
                lambda sfc_net_lw_clr_mon: -sfc_net_lw_clr_mon,
            ),
            (
                ("sfc_net_lw_clr_t_mon",),
                lambda sfc_net_lw_clr_mon: -sfc_net_lw_clr_mon,
            ),
        ]
    ),
    "FLDS": OrderedDict([(("rlds",), rename)]),
    "FLUS": OrderedDict(
        [
            (("rlus",), rename),
            (("FLDS", "FLNS"), lambda FLDS, FLNS: flus(FLDS, FLNS)),
        ]
    ),
    "FLDSC": OrderedDict(
        [
            (("rldscs",), rename),
            (("TS", "FLNSC"), lambda ts, flnsc: fldsc(ts, flnsc)),
        ]
    ),
    "FSNS": OrderedDict(
        [
            (("sfc_net_sw_all_mon",), rename),
            (("rsds", "rsus"), lambda rsds, rsus: netsw(rsds, rsus)),
        ]
    ),
    "FSNSC": OrderedDict(
        [
            (("sfc_net_sw_clr_mon",), rename),
            (("sfc_net_sw_clr_t_mon",), rename),
        ]
    ),
    "FSDS": OrderedDict([(("rsds",), rename)]),
    "FSUS": OrderedDict(
        [
            (("rsus",), rename),
            (("FSDS", "FSNS"), lambda FSDS, FSNS: fsus(FSDS, FSNS)),
        ]
    ),
    "FSUSC": OrderedDict([(("rsuscs",), rename)]),
    "FSDSC": OrderedDict([(("rsdscs",), rename), (("rsdsc",), rename)]),
    # Net surface heat flux: W/(m^2)
    "NET_FLUX_SRF": OrderedDict(
        [
            # A more precise formula to close atmospheric surface budget, than the second entry.
            (
                ("FSNS", "FLNS", "QFLX", "PRECC", "PRECL", "PRECSC", "PRECSL", "SHFLX"),
                lambda fsns, flns, qflx, precc, precl, precsc, precsl, shflx: netflux4(
                    fsns,
                    flns,
                    qflx_convert_to_lhflx(qflx, precc, precl, precsc, precsl),
                    shflx,
                ),
            ),
            (
                ("FSNS", "FLNS", "LHFLX", "SHFLX"),
                lambda fsns, flns, lhflx, shflx: netflux4(fsns, flns, lhflx, shflx),
            ),
            (
                ("FSNS", "FLNS", "QFLX", "SHFLX"),
                lambda fsns, flns, qflx, shflx: netflux4(
                    fsns, flns, qflx_convert_to_lhflx_approxi(qflx), shflx
                ),
            ),
            (
                ("rsds", "rsus", "rlds", "rlus", "hfls", "hfss"),
                lambda rsds, rsus, rlds, rlus, hfls, hfss: netflux6(
                    rsds, rsus, rlds, rlus, hfls, hfss
                ),
            ),
        ]
    ),
    "FLUT": {("rlut",): rename, ("LW_flux_up_at_model_top",): rename},
    "FSUTOA": {("rsut",): rename},
    "FSUTOAC": {("rsutcs",): rename},
    "FLNT": {("FLNT",): rename},
    "FLUTC": {("rlutcs",): rename, ("LW_clrsky_flux_up_at_model_top",): rename},
    "FSNTOA": {
        ("FSNTOA",): rename,
        ("rsdt", "rsut"): rst,
        ("SW_flux_dn_at_model_top", "SW_flux_up_at_model_top"): rst,  # EAMxx
    },
    "FSNTOAC": {
        # Note: CERES_EBAF data in amwg obs sets misspells "units" as "lunits"
        ("FSNTOAC",): rename,
        ("rsdt", "rsutcs"): rstcs,
        ("SW_flux_dn_at_model_top", "SW_clrsky_flux_up_at_model_top"): rstcs,  # EAMxx
    },
    "RESTOM": {
        ("RESTOA",): rename,
        ("toa_net_all_mon",): rename,
        ("FSNT", "FLNT"): restom,
        (
            "SW_flux_dn_at_model_top",
            "SW_flux_up_at_model_top",
            "LW_flux_up_at_model_top",
        ): restom3,  # EAMxx
        ("rtmt",): rename,
    },
    "RESTOA": {
        ("RESTOM",): rename,
        ("toa_net_all_mon",): rename,
        ("FSNT", "FLNT"): restoa,
        ("rtmt",): rename,
    },
    "PRECT_LAND": OrderedDict(
        [
            (("PRECIP_LAND",), rename),
            # 0.5 just to match amwg
            (
                ("PRECC", "PRECL", "LANDFRAC"),
                lambda precc, precl, landfrac: _apply_land_sea_mask(
                    prect(precc, precl), landfrac, lower_limit=0.5
                ),
            ),
        ]
    ),
    "Z3": OrderedDict(
        [
            (
                ("zg",),
                lambda zg: convert_units(rename(zg), target_units="hectometer"),
            ),
            (("Z3",), lambda z3: convert_units(z3, target_units="hectometer")),
        ]
    ),
    "PSL": {
        ("PSL",): lambda psl: convert_units(psl, target_units="mbar"),
        ("psl",): lambda psl: convert_units(psl, target_units="mbar"),
        ("SeaLevelPressure",): lambda psl: convert_units(
            psl, target_units="mbar"
        ),  # EAMxx
    },
    "T": {
        ("ta",): rename,
        ("T",): lambda t: convert_units(t, target_units="K"),
        ("T_2m",): lambda t: convert_units(t, target_units="DegC"),  # EAMxx
    },
    "U": OrderedDict(
        [
            (("ua",), rename),
            (("U",), lambda u: convert_units(u, target_units="m/s")),
        ]
    ),
    "V": OrderedDict(
        [
            (("va",), rename),
            (("V",), lambda u: convert_units(u, target_units="m/s")),
        ]
    ),
    "TREFHT": OrderedDict(
        [
            (("TREFHT",), lambda t: convert_units(t, target_units="DegC")),
            (
                ("TREFHT_LAND",),
                lambda t: convert_units(t, target_units="DegC"),
            ),
            (("tas",), lambda t: convert_units(t, target_units="DegC")),
        ]
    ),
    # Surface water flux: kg/((m^2)*s)
    "QFLX": {
        ("evspsbl",): rename,
        ("QFLX",): qflxconvert_units,
        ("surf_evap",): qflxconvert_units,  # EAMxx
    },
    # Surface latent heat flux: W/(m^2)
    "LHFLX": {
        ("hfls",): rename,
        ("QFLX",): qflx_convert_to_lhflx_approxi,
        ("surface_upward_latent_heat_flux",): rename,  # EAMxx "s^-3 kg"
    },
    "SHFLX": {
        ("hfss",): rename,
        ("surf_sens_flux",): rename,  # EAMxx
    },
    "TGCLDLWP_OCN": OrderedDict(
        [
            (
                ("TGCLDLWP_OCEAN",),
                lambda x: convert_units(x, target_units="g/m^2"),
            ),
            (
                ("TGCLDLWP", "OCNFRAC"),
                lambda tgcldlwp, ocnfrac: _apply_land_sea_mask(
                    convert_units(tgcldlwp, target_units="g/m^2"),
                    ocnfrac,
                    lower_limit=0.65,
                ),
            ),
        ]
    ),
    "PRECT_OCN": OrderedDict(
        [
            (
                ("PRECT_OCEAN",),
                lambda x: convert_units(x, target_units="mm/day"),
            ),
            (
                ("PRECC", "PRECL", "OCNFRAC"),
                lambda a, b, ocnfrac: _apply_land_sea_mask(
                    aplusb(a, b, target_units="mm/day"),
                    ocnfrac,
                    lower_limit=0.65,
                ),
            ),
        ]
    ),
    "PREH2O_OCN": OrderedDict(
        [
            (("PREH2O_OCEAN",), lambda x: convert_units(x, target_units="mm")),
            (
                ("TMQ", "OCNFRAC"),
                lambda preh2o, ocnfrac: _apply_land_sea_mask(
                    preh2o, ocnfrac, lower_limit=0.65
                ),
            ),
        ]
    ),
    "CLDHGH": OrderedDict(
        [(("CLDHGH",), lambda cldhgh: convert_units(cldhgh, target_units="%"))]
    ),
    "CLDLOW": OrderedDict(
        [(("CLDLOW",), lambda cldlow: convert_units(cldlow, target_units="%"))]
    ),
    "CLDMED": OrderedDict(
        [(("CLDMED",), lambda cldmed: convert_units(cldmed, target_units="%"))]
    ),
    "CLDTOT": OrderedDict(
        [
            (("clt",), rename),
            (
                ("CLDTOT",),
                lambda cldtot: convert_units(cldtot, target_units="%"),
            ),
        ]
    ),
    "CLOUD": OrderedDict(
        [
            (("cl",), rename),
            (
                ("CLOUD",),
                lambda cldtot: convert_units(cldtot, target_units="%"),
            ),
        ]
    ),
    # below for COSP output
    # CLIPSO
    "CLDHGH_CAL": OrderedDict(
        [
            (
                ("CLDHGH_CAL",),
                lambda cldhgh: convert_units(cldhgh, target_units="%"),
            )
        ]
    ),
    "CLDLOW_CAL": OrderedDict(
        [
            (
                ("CLDLOW_CAL",),
                lambda cldlow: convert_units(cldlow, target_units="%"),
            )
        ]
    ),
    "CLDMED_CAL": OrderedDict(
        [
            (
                ("CLDMED_CAL",),
                lambda cldmed: convert_units(cldmed, target_units="%"),
            )
        ]
    ),
    "CLDTOT_CAL": OrderedDict(
        [
            (
                ("CLDTOT_CAL",),
                lambda cldtot: convert_units(cldtot, target_units="%"),
            )
        ]
    ),
    # ISCCP
    "CLDTOT_TAU1.3_ISCCP": {
        ("FISCCP1_COSP",): cosp_bin_sum,
        ("CLISCCP",): cosp_bin_sum,
    },
    "CLDTOT_TAU1.3_9.4_ISCCP": {
        ("FISCCP1_COSP",): cosp_bin_sum,
        ("CLISCCP",): cosp_bin_sum,
    },
    "CLDTOT_TAU9.4_ISCCP": {
        ("FISCCP1_COSP",): cosp_bin_sum,
        ("CLISCCP",): cosp_bin_sum,
    },
    # MODIS
    "CLDTOT_TAU1.3_MODIS": {
        ("CLMODIS",): cosp_bin_sum,
    },
    "CLDTOT_TAU1.3_9.4_MODIS": {
        ("CLMODIS",): cosp_bin_sum,
    },
    "CLDTOT_TAU9.4_MODIS": {
        ("CLMODIS",): cosp_bin_sum,
    },
    "CLDHGH_TAU1.3_MODIS": {
        ("CLMODIS",): cosp_bin_sum,
    },
    "CLDHGH_TAU1.3_9.4_MODIS": {
        ("CLMODIS",): cosp_bin_sum,
    },
    "CLDHGH_TAU9.4_MODIS": {
        ("CLMODIS",): cosp_bin_sum,
    },
    # MISR
    "CLDTOT_TAU1.3_MISR": {
        ("CLD_MISR",): cosp_bin_sum,
        ("CLMISR",): cosp_bin_sum,
    },
    "CLDTOT_TAU1.3_9.4_MISR": {
        ("CLD_MISR",): cosp_bin_sum,
        ("CLMISR",): cosp_bin_sum,
    },
    "CLDTOT_TAU9.4_MISR": {
        ("CLD_MISR",): cosp_bin_sum,
        ("CLMISR",): cosp_bin_sum,
    },
    "CLDLOW_TAU1.3_MISR": {
        ("CLD_MISR",): cosp_bin_sum,
        ("CLMISR",): cosp_bin_sum,
    },
    "CLDLOW_TAU1.3_9.4_MISR": {
        ("CLD_MISR",): cosp_bin_sum,
        ("CLMISR",): cosp_bin_sum,
    },
    "CLDLOW_TAU9.4_MISR": {
        ("CLD_MISR",): cosp_bin_sum,
        ("CLMISR",): cosp_bin_sum,
    },
    # COSP cloud fraction joint histogram
    "COSP_HISTOGRAM_MISR": {
        ("CLD_MISR",): cosp_histogram_standardize,
        ("CLMISR",): cosp_histogram_standardize,
    },
    "COSP_HISTOGRAM_MODIS": {
        ("CLMODIS",): cosp_histogram_standardize,
    },
    "COSP_HISTOGRAM_ISCCP": {
        ("FISCCP1_COSP",): cosp_histogram_standardize,
        ("CLISCCP",): cosp_histogram_standardize,
    },
    "ICEFRAC": {
        ("ICEFRAC",): lambda icefrac: convert_units(icefrac, target_units="%"),
    },
    "RELHUM": {
        ("hur",): lambda hur: convert_units(hur, target_units="%"),
        ("RELHUM",): lambda relhum: convert_units(relhum, target_units="%"),
    },
    "OMEGA": {
        ("wap",): lambda wap: convert_units(wap, target_units="mbar/day"),
        ("OMEGA",): lambda omega: convert_units(omega, target_units="mbar/day"),
    },
    "Q": {
        ("hus",): lambda q: convert_units(rename(q), target_units="g/kg"),
        ("Q",): lambda q: convert_units(rename(q), target_units="g/kg"),
        ("SHUM",): lambda shum: convert_units(shum, target_units="g/kg"),
    },
    "H2OLNZ": {
        ("hus",): lambda q: convert_units(rename(q), target_units="g/kg"),
        ("H2OLNZ",): lambda h2o: w_convert_q(h2o),
    },
    "TAUXY": {
        ("TAUX", "TAUY"): tauxy,
        ("tauu", "tauv"): tauxy,
    },
    "AODVIS": {
        ("od550aer",): rename,
        ("AODVIS",): lambda aod: convert_units(
            rename(aod), target_units="dimensionless"
        ),
        ("AOD_550",): lambda aod: convert_units(
            rename(aod), target_units="dimensionless"
        ),
        ("TOTEXTTAU",): lambda aod: convert_units(
            rename(aod), target_units="dimensionless"
        ),
        ("AOD_550_ann",): lambda aod: convert_units(
            rename(aod), target_units="dimensionless"
        ),
    },
    "AODABS": {("abs550aer",): rename},
    "AODDUST": {
        ("AODDUST",): lambda aod: convert_units(
            rename(aod), target_units="dimensionless"
        ),
    },
    # Surface temperature: Degrees C
    # (Temperature of the surface (land/water) itself, not the air)
    "TS": {
        ("ts",): rename,
        ("surf_radiative_T",): rename,  # EAMxx
    },
    "PS": {("ps",): rename},
    "U10": {("sfcWind",): rename},
    "QREFHT": {
        ("QREFHT",): lambda q: convert_units(q, target_units="g/kg"),
        ("huss",): lambda q: convert_units(q, target_units="g/kg"),
        ("d2m", "sp"): qsat,
    },
    "TAUX": {
        ("tauu",): lambda tauu: -tauu,
        ("surf_mom_flux_U",): lambda tauu: -tauu,  # EAMxx
    },
    "TAUY": {
        ("tauv",): lambda tauv: -tauv,
        ("surf_mom_flux_V",): lambda tauv: -tauv,  # EAMxx
    },
    "CLDICE": {("cli",): rename},
    "TGCLDIWP": {("clivi",): rename},
    "CLDLIQ": {("clw",): rename},
    "TGCLDCWP": {("clwvi",): rename},
    "O3": {("o3",): rename},
    "PminusE": {
        ("PminusE",): pminuse_convert_units,
        ("PRECC", "PRECL", "QFLX"): pminuse_1,
        ("F_prec", "F_evap"): pminuse_2,
        ("pr", "evspsbl"): pminuse_3,
    },
    "TREFMNAV": {
        ("TREFMNAV",): lambda t: convert_units(t, target_units="DegC"),
        ("tasmin",): lambda t: convert_units(t, target_units="DegC"),
    },
    "TREFMXAV": {
        ("TREFMXAV",): lambda t: convert_units(t, target_units="DegC"),
        ("tasmax",): lambda t: convert_units(t, target_units="DegC"),
    },
    "TREF_range": {
        ("TREFMXAV", "TREFMNAV"): lambda tmax, tmin: tref_range(tmax, tmin),
        ("tasmax", "tasmin"): lambda tmax, tmin: tref_range(tmax, tmin),
    },
    "TCO": {("TCO",): rename},
    "SCO": {("SCO",): rename},
    "bc_DDF": {
        ("bc_DDF",): rename,
        ("bc_a?DDF", "bc_c?DDF"): sum_vars,
    },
    "bc_SFWET": {
        ("bc_SFWET",): rename,
        ("bc_a?SFWET", "bc_c?SFWET"): sum_vars,
    },
    "SFbc": {
        ("SFbc",): rename,
        ("SFbc_a?",): sum_vars,
    },
    "bc_CLXF": {
        ("bc_CLXF",): rename,
        ("bc_a?_CLXF",): lambda x: molec_convert_units(x, 12),
    },
    "Mass_bc": {("Mass_bc",): rename},
    "dst_DDF": {
        ("dst_DDF",): rename,
        ("dst_a?DDF", "dst_c?DDF"): sum_vars,
    },
    "dst_SFWET": {
        ("dst_SFWET",): rename,
        ("dst_a?SFWET", "dst_c?SFWET"): sum_vars,
    },
    "SFdst": {("SFdst",): rename, ("SFdst_a?",): sum_vars},
    "Mass_dst": {("Mass_dst",): rename},
    "mom_DDF": {
        ("mom_DDF",): rename,
        ("mom_a?DDF", "mom_c?DDF"): sum_vars,
    },
    "mom_SFWET": {
        ("mom_SFWET",): rename,
        ("mom_a?SFWET", "mom_c?SFWET"): sum_vars,
    },
    "SFmom": {("SFmom",): rename, ("SFmom_a?",): sum_vars},
    "Mass_mom": {("Mass_mom",): rename},
    "ncl_DDF": {
        ("ncl_DDF",): rename,
        ("ncl_a?DDF", "ncl_c?DDF"): sum_vars,
    },
    "ncl_SFWET": {
        ("ncl_SFWET",): rename,
        ("ncl_a?SFWET", "ncl_c?SFWET"): sum_vars,
    },
    "SFncl": {("SFncl",): rename, ("SFncl_a?",): sum_vars},
    "Mass_ncl": {("Mass_ncl",): rename},
    "so4_DDF": {
        ("so4_DDF",): rename,
        ("so4_a?DDF", "so4_c?DDF"): sum_vars,
    },
    "so4_SFWET": {
        ("so4_SFWET",): rename,
        ("so4_a?SFWET", "so4_c?SFWET"): sum_vars,
    },
    "so4_CLXF": {
        ("so4_CLXF",): rename,
        ("so4_a?_CLXF",): lambda x: molec_convert_units(x, 115.0),
    },
    "SFso4": {
        ("SFso4",): rename,
        ("SFso4_a?",): sum_vars,
    },
    "Mass_so4": {("Mass_so4",): rename},
    "soa_DDF": {
        ("soa_DDF",): rename,
        ("soa_a?DDF", "soa_c?DDF"): sum_vars,
    },
    "soa_SFWET": {
        ("soa_SFWET",): rename,
        ("soa_a?SFWET", "soa_c?SFWET"): sum_vars,
    },
    "SFsoa": {("SFsoa",): rename, ("SFsoa_a?",): sum_vars},
    "Mass_soa": {("Mass_soa",): rename},
    "pom_DDF": {
        ("pom_DDF",): rename,
        ("pom_a?DDF", "pom_c?DDF"): sum_vars,
    },
    "pom_SFWET": {
        ("pom_SFWET",): rename,
        ("pom_a?SFWET", "pom_c?SFWET"): sum_vars,
    },
    "SFpom": {
        ("SFpom",): rename,
        ("SFpom_a?",): sum_vars,
    },
    "pom_CLXF": {
        ("pom_CLXF",): rename,
        ("pom_a?_CLXF",): lambda x: molec_convert_units(x, 12.0),
    },
    "Mass_pom": {("Mass_pom",): rename},
    # total aerosol number concentration (#/CC)
    "a_num": {
        ("cpc",): rename,
        # Aerosol concentration from Aitken, Accumu., and Coarse mode
        (
            "num_a1",
            "num_a2",
            "num_a3",
        ): lambda a1, a2, a3: a_num_sum(a1, a2, a3),
    },
    # total so4 mass concentration (ng/m3)
    "so4_mass": {
        ("sulfate",): rename,
        # Aerosol concentration from Aitken, Accumu., and Coarse mode
        (
            "so4_a1",
            "so4_a2",
        ): lambda a1, a2: so4_mass_sum(a1, a2),
    },
    # CCN 0.1%SS concentration (1/CC)
    "ccn01": {
        ("ccn01",): rename,
        ("CCN3",): rename,
    },
    # CCN 0.2%SS concentration (1/CC)
    "ccn02": {
        ("ccn02",): rename,
        ("CCN4",): rename,
    },
    # CCN 0.5%SS concentration (1/CC)
    "ccn05": {
        ("ccn05",): rename,
        ("CCN5",): rename,
    },
    # Land variables
    "SOILWATER_10CM": {("mrsos",): rename},
    "SOILWATER_SUM": {("mrso",): rename},
    "SOILICE_SUM": {("mrfso",): rename},
    "QRUNOFF": OrderedDict(
        [
            (("QRUNOFF",), lambda qrunoff: qflxconvert_units(qrunoff)),
            (("mrro",), lambda qrunoff: qflxconvert_units(qrunoff)),
        ]
    ),
    "QINTR": OrderedDict([(("prveg",), rename)]),
    "QVEGE": OrderedDict(
        [
            (("QVEGE",), lambda qevge: qflxconvert_units(rename(qevge))),
            (("evspsblveg",), lambda qevge: qflxconvert_units(rename(qevge))),
        ]
    ),
    "QVEGT": OrderedDict(
        [
            (("QVEGT",), lambda qevgt: qflxconvert_units(rename(qevgt))),
        ]
    ),
    "QSOIL": OrderedDict(
        [
            (("QSOIL",), lambda qsoil: qflxconvert_units(rename(qsoil))),
            (("evspsblsoi",), lambda qsoil: qflxconvert_units(rename(qsoil))),
        ]
    ),
    "QDRAI": OrderedDict(
        [
            (("QDRAI",), lambda q: qflxconvert_units(rename(q))),
        ]
    ),
    "QINFL": OrderedDict(
        [
            (("QINFL",), lambda q: qflxconvert_units(rename(q))),
        ]
    ),
    "QIRRIG_GRND": OrderedDict(
        [
            (("QIRRIG_GRND",), lambda q: qflxconvert_units(rename(q))),
        ]
    ),
    "QIRRIG_ORIG": OrderedDict(
        [
            (("QIRRIG_ORIG",), lambda q: qflxconvert_units(rename(q))),
        ]
    ),
    "QIRRIG_REAL": OrderedDict(
        [
            (("QIRRIG_REAL",), lambda q: qflxconvert_units(rename(q))),
        ]
    ),
    "QIRRIG_SURF": OrderedDict(
        [
            (("QIRRIG_SURF",), lambda q: qflxconvert_units(rename(q))),
        ]
    ),
    "QIRRIG_WM": OrderedDict(
        [
            (("QIRRIG_WM",), lambda q: qflxconvert_units(rename(q))),
        ]
    ),
    "QOVER": OrderedDict(
        [
            (("QOVER",), lambda q: qflxconvert_units(rename(q))),
            (("mrros",), lambda q: qflxconvert_units(rename(q))),
        ]
    ),
    "QRGWL": OrderedDict(
        [
            (("QRGWL",), lambda q: qflxconvert_units(rename(q))),
        ]
    ),
    "RAIN": OrderedDict(
        [
            (("RAIN",), lambda rain: qflxconvert_units(rename(rain))),
        ]
    ),
    "SNOW": OrderedDict(
        [
            (("SNOW",), lambda snow: qflxconvert_units(rename(snow))),
        ]
    ),
    "TRAN": OrderedDict([(("tran",), rename)]),
    "TSOI": OrderedDict([(("tsl",), rename)]),
    "LAI": OrderedDict([(("lai",), rename)]),
    # Additional land variables requested by BGC evaluation
    "FAREA_BURNED": OrderedDict(
        [
            (
                ("FAREA_BURNED",),
                lambda v: convert_units(v, target_units="proportionx10^9"),
            )
        ]
    ),
    "FLOODPLAIN_VOLUME": OrderedDict(
        [(("FLOODPLAIN_VOLUME",), lambda v: convert_units(v, target_units="km3"))]
    ),
    "TLAI": OrderedDict([(("TLAI",), rename)]),
    "EFLX_LH_TOT": OrderedDict([(("EFLX_LH_TOT",), rename)]),
    "GPP": OrderedDict(
        [(("GPP",), lambda v: convert_units(v, target_units="g*/m^2/day"))]
    ),
    "HR": OrderedDict(
        [(("HR",), lambda v: convert_units(v, target_units="g*/m^2/day"))]
    ),
    "NBP": OrderedDict(
        [(("NBP",), lambda v: convert_units(v, target_units="g*/m^2/day"))]
    ),
    "NPP": OrderedDict(
        [(("NPP",), lambda v: convert_units(v, target_units="g*/m^2/day"))]
    ),
    "TOTVEGC": OrderedDict(
        [(("TOTVEGC",), lambda v: convert_units(v, target_units="kgC/m^2"))]
    ),
    "TOTSOMC": OrderedDict(
        [(("TOTSOMC",), lambda v: convert_units(v, target_units="kgC/m^2"))]
    ),
    "TOTSOMN": OrderedDict([(("TOTSOMN",), rename)]),
    "TOTSOMP": OrderedDict([(("TOTSOMP",), rename)]),
    "FPG": OrderedDict([(("FPG",), rename)]),
    "FPG_P": OrderedDict([(("FPG_P",), rename)]),
    "TBOT": OrderedDict([(("TBOT",), rename)]),
    "CPOOL": OrderedDict(
        [(("CPOOL",), lambda v: convert_units(v, target_units="kgC/m^2"))]
    ),
    "LEAFC": OrderedDict(
        [(("LEAFC",), lambda v: convert_units(v, target_units="kgC/m^2"))]
    ),
    "SR": OrderedDict([(("SR",), lambda v: convert_units(v, target_units="kgC/m^2"))]),
    "RH2M": OrderedDict([(("RH2M",), rename)]),
    "DENIT": OrderedDict(
        [(("DENIT",), lambda v: convert_units(v, target_units="mg*/m^2/day"))]
    ),
    "GROSS_NMIN": OrderedDict(
        [(("GROSS_NMIN",), lambda v: convert_units(v, target_units="mg*/m^2/day"))]
    ),
    "GROSS_PMIN": OrderedDict(
        [(("GROSS_PMIN",), lambda v: convert_units(v, target_units="mg*/m^2/day"))]
    ),
    "NDEP_TO_SMINN": OrderedDict(
        [
            (
                ("NDEP_TO_SMINN",),
                lambda v: convert_units(v, target_units="mg*/m^2/day"),
            )
        ]
    ),
    "NFIX_TO_SMINN": OrderedDict(
        [
            (
                ("NFIX_TO_SMINN",),
                lambda v: convert_units(v, target_units="mg*/m^2/day"),
            )
        ]
    ),
    "PLANT_NDEMAND_COL": OrderedDict(
        [
            (
                ("PLANT_NDEMAND_COL",),
                lambda v: convert_units(v, target_units="mg*/m^2/day"),
            )
        ]
    ),
    "PLANT_PDEMAND_COL": OrderedDict(
        [
            (
                ("PLANT_PDEMAND_COL",),
                lambda v: convert_units(v, target_units="mg*/m^2/day"),
            )
        ]
    ),
    "SMINN_TO_PLANT": OrderedDict(
        [
            (
                ("SMINN_TO_PLANT",),
                lambda v: convert_units(v, target_units="mg*/m^2/day"),
            )
        ]
    ),
    "SMINP_TO_PLANT": OrderedDict(
        [
            (
                ("SMINP_TO_PLANT",),
                lambda v: convert_units(v, target_units="mg*/m^2/day"),
            )
        ]
    ),
    "SMIN_NO3_LEACHED": OrderedDict(
        [
            (
                ("SMIN_NO3_LEACHED",),
                lambda v: convert_units(v, target_units="mg*/m^2/day"),
            )
        ]
    ),
    "FP_UPTAKE": OrderedDict(
        [
            (("FP_UPTAKE",), rename),
            (
                ("SMINN_TO_PLANT", "PLANT_NDEMAND_COL"),
                lambda a, b: fp_uptake(a, b),
            ),
        ]
    ),
    # Ocean variables
    "tauuo": {("tauuo",): rename},
    "tos": {("tos",): rename},
    "thetaoga": {("thetaoga",): rename},
    "hfsifrazil": {("hfsifrazil",): rename},
    "sos": {("sos",): rename},
    "soga": {("soga",): rename},
    "tosga": {("tosga",): rename},
    "wo": {("wo",): rename},
    "thetao": {("thetao",): rename},
    "masscello": {("masscello",): rename},
    "wfo": {("wfo",): rename},
    "tauvo": {("tauvo",): rename},
    "vo": {("vo",): rename},
    "hfds": {("hfds",): rename},
    "volo": {("volo",): rename},
    "uo": {("uo",): rename},
    "zos": {("zos",): rename},
    "tob": {("tob",): rename},
    "sosga": {("sosga",): rename},
    "sfdsi": {("sfdsi",): rename},
    "zhalfo": {("zhalfo",): rename},
    "masso": {("masso",): rename},
    "so": {("so",): rename},
    "sob": {("sob",): rename},
    "mlotst": {("mlotst",): rename},
    "fsitherm": {("fsitherm",): rename},
    "msftmz": {("msftmz",): rename},
    # sea ice variables
    "sitimefrac": {("sitimefrac",): rename},
    "siconc": {("siconc",): rename},
    "sisnmass": {("sisnmass",): rename},
    "sisnthick": {("sisnthick",): rename},
    "simass": {("simass",): rename},
    "sithick": {("sithick",): rename},
    "siu": {("siu",): rename},
    "sitemptop": {("sitemptop",): rename},
    "siv": {("siv",): rename},
}


# Names of 2D aerosol burdens, including cloud-borne aerosols
aero_burden_list = [
    "ABURDENDUST",
    "ABURDENSO4",
    "ABURDENSO4_STR",
    "ABURDENSO4_TRO",
    "ABURDENPOM",
    "ABURDENMOM",
    "ABURDENSOA",
    "ABURDENBC",
    "ABURDENSEASALT",
]

# Add burden vars to DERIVED_VARIABLES
for aero_burden_item in aero_burden_list:
    DERIVED_VARIABLES[f"_{aero_burden_item}"] = OrderedDict(
        [((aero_burden_item,), aero_burden_fxn)]
    )


# Names of 2D mass slices of aerosol species
# Also add 3D masses while at it (if available)
aero_mass_list = []
for aero_name in ["dst", "mom", "pom", "so4", "soa", "ncl", "bc"]:
    for aero_lev in ["_srf", "_200", "_330", "_500", "_850", ""]:
        # Note that the empty string (last entry) will get the 3D mass fields
        aero_mass_list.append(f"Mass_{aero_name}{aero_lev}")


# Add burden vars to DERIVED_VARIABLES
for aero_mass_item in aero_mass_list:
    DERIVED_VARIABLES[f"_{aero_mass_item}"] = OrderedDict(
        [((aero_mass_item,), aero_mass_fxn)]
    )

# Add all the output_aerocom_aie.F90 variables to aero_rename_list
# components/eam/src/physics/cam/output_aerocom_aie.F90
aero_aerocom_list = [
    "angstrm",
    "aerindex",
    "cdr",
    "cdnc",
    "cdnum",
    "icnum",
    "clt",
    "lcc",
    "lwp",
    "iwp",
    "icr",
    "icc",
    "cod",
    "ccn",
    "ttop",
    "htop",
    "ptop",
    "autoconv",
    "accretn",
    "icnc",
    "rh700",
    "rwp",
    "intccn",
    "colrv",
    "lwp2",
    "iwp2",
    "lwpbf",
    "iwpbf",
    "cdnumbf",
    "icnumbf",
    "aod400",
    "aod700",
    "colccn.1",
    "colccn.3",
    "ccn.1bl",
    "ccn.3bl",
]

# Add aerocom vars to DERIVED_VARIABLES
for aero_aerocom_item in aero_aerocom_list:
    DERIVED_VARIABLES[aero_aerocom_item] = OrderedDict([((aero_aerocom_item,), rename)])

# add cdnc, icnc, lwp, iwp to DERIVED_VARIABLES
DERIVED_VARIABLES.update(
    {
        "in_cloud_cdnc": {("cdnc", "lcc"): incldtop_cdnc},
        "in_grid_cdnc": {("cdnc",): cldtop_cdnc},
        "in_cloud_icnc": {("icnc", "icc"): incldtop_icnc},
        "in_grid_icnc": {("icnc",): cldtop_icnc},
        "in_cloud_lwp": {("lwp", "lcc"): incld_lwp},
        "in_grid_lwp": {("lwp",): cld_lwp},
        "in_cloud_iwp": {("iwp", "icc"): incld_iwp},
        "in_grid_iwp": {("iwp",): cld_iwp},
    }
)


DERIVED_VARIABLES.update(
    {
        "ERFtot": {("FSNT", "FLNT"): erf_tot},
        "ERFari": {("FSNT", "FLNT", "FSNT_d1", "FLNT_d1"): erf_ari},
        "ERFaci": {("FSNT_d1", "FLNT_d1", "FSNTC_d1", "FLNTC_d1"): erf_aci},
        "ERFres": {("FSNTC_d1", "FLNTC_d1"): erf_res},
    }
)

# Add more AOD terms
# Note that AODVIS and AODDUST are already added elsewhere
aero_aod_list = [
    "AODBC",
    "AODPOM",
    "AODMOM",
    "AODSO4",
    "AODSO4_STR",
    "AODSO4_TRO",
    "AODSS",
    "AODSOA",
]

# Add aod vars to DERIVED_VARIABLES
for aero_aod_item in aero_aod_list:
    DERIVED_VARIABLES[aero_aod_item] = {(aero_aod_item,): rename}

# Add 3D variables related to aerosols and chemistry
# Note that O3 is already added above
# Note that 3D mass vars are already added by the empty string above ""
# Note that it is possible to create on-the-fly slices from these variables with
# a function of the form:
# def aero_3d_slice(var, lev):
#     return var[lev, :, :]
aero_chem_list = ["DMS", "H2O2", "H2SO4", "NO3", "OH", "SO2"]

# Add aero/chem vars to DERIVED_VARIABLES
for aero_chem_item in aero_chem_list:
    DERIVED_VARIABLES[aero_chem_item] = {(aero_chem_item,): rename}

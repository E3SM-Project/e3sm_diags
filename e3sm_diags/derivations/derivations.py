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
from typing import Callable, Dict, Tuple

from e3sm_diags.derivations.formulas import (
    albedo,
    albedo_srf,
    albedoc,
    fldsc,
    flus,
    fp_uptake,
    fsus,
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
    pminuse_convert_units,
    precst,
    prect,
    qflx_convert_to_lhflx,
    qflx_convert_to_lhflx_approxi,
    qflxconvert_units,
    restoa,
    restom,
    rst,
    rstcs,
    swcf,
    swcfsrf,
    tauxy,
    tref_range,
    w_convert_q,
)
from e3sm_diags.derivations.utils import (
    _apply_land_sea_mask,
    aplusb,
    convert_units,
    cosp_bin_sum,
    cosp_histogram_standardize,
    rename,
)

# A type annotation ordered dictionary that maps a tuple of source variable(s)
# to a derivation function.
DerivedVariableMap = OrderedDict[Tuple[str, ...], Callable]

# A type annotation for a dictionary mapping the key of a derived variable
# to an ordered dictionary that maps a tuple of source variable(s) to a
# derivation function.
DerivedVariablesMap = Dict[str, DerivedVariableMap]


DERIVED_VARIABLES: DerivedVariablesMap = {
    "PRECT": OrderedDict(
        [
            (
                ("PRECT",),
                lambda pr: convert_units(rename(pr), target_units="mm/day"),
            ),
            (("pr",), lambda pr: qflxconvert_units(rename(pr))),
            (("PRECC", "PRECL"), lambda precc, precl: prect(precc, precl)),
        ]
    ),
    "PRECST": OrderedDict(
        [
            (("prsn",), lambda prsn: qflxconvert_units(rename(prsn))),
            (
                ("PRECSC", "PRECSL"),
                lambda precsc, precsl: precst(precsc, precsl),
            ),
        ]
    ),
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
    "SOLIN": OrderedDict([(("rsdt",), rename)]),
    "ALBEDO": OrderedDict(
        [
            (("ALBEDO",), rename),
            (
                ("SOLIN", "FSNTOA"),
                lambda solin, fsntoa: albedo(solin, solin - fsntoa),
            ),
            (("rsdt", "rsut"), lambda rsdt, rsut: albedo(rsdt, rsut)),
        ]
    ),
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
    "FLUT": OrderedDict([(("rlut",), rename)]),
    "FSUTOA": OrderedDict([(("rsut",), rename)]),
    "FSUTOAC": OrderedDict([(("rsutcs",), rename)]),
    "FLNT": OrderedDict([(("FLNT",), rename)]),
    "FLUTC": OrderedDict([(("rlutcs",), rename)]),
    "FSNTOA": OrderedDict(
        [
            (("FSNTOA",), rename),
            (("rsdt", "rsut"), lambda rsdt, rsut: rst(rsdt, rsut)),
        ]
    ),
    "FSNTOAC": OrderedDict(
        [
            # Note: CERES_EBAF data in amwg obs sets misspells "units" as "lunits"
            (("FSNTOAC",), rename),
            (("rsdt", "rsutcs"), lambda rsdt, rsutcs: rstcs(rsdt, rsutcs)),
        ]
    ),
    "RESTOM": OrderedDict(
        [
            (("RESTOA",), rename),
            (("toa_net_all_mon",), rename),
            (("FSNT", "FLNT"), lambda fsnt, flnt: restom(fsnt, flnt)),
            (("rtmt",), rename),
        ]
    ),
    "RESTOA": OrderedDict(
        [
            (("RESTOM",), rename),
            (("toa_net_all_mon",), rename),
            (("FSNT", "FLNT"), lambda fsnt, flnt: restoa(fsnt, flnt)),
            (("rtmt",), rename),
        ]
    ),
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
    "PSL": OrderedDict(
        [
            (("PSL",), lambda psl: convert_units(psl, target_units="mbar")),
            (("psl",), lambda psl: convert_units(psl, target_units="mbar")),
        ]
    ),
    "T": OrderedDict(
        [
            (("ta",), rename),
            (("T",), lambda t: convert_units(t, target_units="K")),
        ]
    ),
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
    "QFLX": OrderedDict(
        [
            (("evspsbl",), rename),
            (("QFLX",), lambda qflx: qflxconvert_units(qflx)),
        ]
    ),
    # Surface latent heat flux: W/(m^2)
    "LHFLX": OrderedDict(
        [
            (("hfls",), rename),
            (("QFLX",), lambda qflx: qflx_convert_to_lhflx_approxi(qflx)),
        ]
    ),
    "SHFLX": OrderedDict([(("hfss",), rename)]),
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
    "CLDTOT_TAU1.3_ISCCP": OrderedDict(
        [
            (
                ("FISCCP1_COSP",),
                lambda cld: convert_units(
                    cosp_bin_sum(cld, None, None, 1.3, None), target_units="%"
                ),
            ),
            (
                ("CLISCCP",),
                lambda cld: convert_units(
                    cosp_bin_sum(cld, None, None, 1.3, None), target_units="%"
                ),
            ),
        ]
    ),
    "CLDTOT_TAU1.3_9.4_ISCCP": OrderedDict(
        [
            (
                ("FISCCP1_COSP",),
                lambda cld: convert_units(
                    cosp_bin_sum(cld, None, None, 1.3, 9.4), target_units="%"
                ),
            ),
            (
                ("CLISCCP",),
                lambda cld: convert_units(
                    cosp_bin_sum(cld, None, None, 1.3, 9.4), target_units="%"
                ),
            ),
        ]
    ),
    "CLDTOT_TAU9.4_ISCCP": OrderedDict(
        [
            (
                ("FISCCP1_COSP",),
                lambda cld: convert_units(
                    cosp_bin_sum(cld, None, None, 9.4, None), target_units="%"
                ),
            ),
            (
                ("CLISCCP",),
                lambda cld: convert_units(
                    cosp_bin_sum(cld, None, None, 9.4, None), target_units="%"
                ),
            ),
        ]
    ),
    # MODIS
    "CLDTOT_TAU1.3_MODIS": OrderedDict(
        [
            (
                ("CLMODIS",),
                lambda cld: convert_units(
                    cosp_bin_sum(cld, None, None, 1.3, None), target_units="%"
                ),
            ),
        ]
    ),
    "CLDTOT_TAU1.3_9.4_MODIS": OrderedDict(
        [
            (
                ("CLMODIS",),
                lambda cld: convert_units(
                    cosp_bin_sum(cld, None, None, 1.3, 9.4), target_units="%"
                ),
            ),
        ]
    ),
    "CLDTOT_TAU9.4_MODIS": OrderedDict(
        [
            (
                ("CLMODIS",),
                lambda cld: convert_units(
                    cosp_bin_sum(cld, None, None, 9.4, None), target_units="%"
                ),
            ),
        ]
    ),
    "CLDHGH_TAU1.3_MODIS": OrderedDict(
        [
            (
                ("CLMODIS",),
                lambda cld: convert_units(
                    cosp_bin_sum(cld, 440, 0, 1.3, None), target_units="%"
                ),
            ),
        ]
    ),
    "CLDHGH_TAU1.3_9.4_MODIS": OrderedDict(
        [
            (
                ("CLMODIS",),
                lambda cld: convert_units(
                    cosp_bin_sum(cld, 440, 0, 1.3, 9.4), target_units="%"
                ),
            ),
        ]
    ),
    "CLDHGH_TAU9.4_MODIS": OrderedDict(
        [
            (
                ("CLMODIS",),
                lambda cld: convert_units(
                    cosp_bin_sum(cld, 440, 0, 9.4, None), target_units="%"
                ),
            ),
        ]
    ),
    # MISR
    "CLDTOT_TAU1.3_MISR": OrderedDict(
        [
            (
                ("CLD_MISR",),
                lambda cld: convert_units(
                    cosp_bin_sum(cld, None, None, 1.3, None), target_units="%"
                ),
            ),
            (
                ("CLMISR",),
                lambda cld: convert_units(
                    cosp_bin_sum(cld, None, None, 1.3, None), target_units="%"
                ),
            ),
        ]
    ),
    "CLDTOT_TAU1.3_9.4_MISR": OrderedDict(
        [
            (
                ("CLD_MISR",),
                lambda cld: convert_units(
                    cosp_bin_sum(cld, None, None, 1.3, 9.4), target_units="%"
                ),
            ),
            (
                ("CLMISR",),
                lambda cld: convert_units(
                    cosp_bin_sum(cld, None, None, 1.3, 9.4), target_units="%"
                ),
            ),
        ]
    ),
    "CLDTOT_TAU9.4_MISR": OrderedDict(
        [
            (
                ("CLD_MISR",),
                lambda cld: convert_units(
                    cosp_bin_sum(cld, None, None, 9.4, None), target_units="%"
                ),
            ),
            (
                ("CLMISR",),
                lambda cld: convert_units(
                    cosp_bin_sum(cld, None, None, 9.4, None), target_units="%"
                ),
            ),
        ]
    ),
    "CLDLOW_TAU1.3_MISR": OrderedDict(
        [
            (
                ("CLD_MISR",),
                lambda cld: convert_units(
                    cosp_bin_sum(cld, 0, 3, 1.3, None), target_units="%"
                ),
            ),
            (
                ("CLMISR",),
                lambda cld: convert_units(
                    cosp_bin_sum(cld, 0, 3, 1.3, None), target_units="%"
                ),
            ),
        ]
    ),
    "CLDLOW_TAU1.3_9.4_MISR": OrderedDict(
        [
            (
                ("CLD_MISR",),
                lambda cld: convert_units(
                    cosp_bin_sum(cld, 0, 3, 1.3, 9.4), target_units="%"
                ),
            ),
            (
                ("CLMISR",),
                lambda cld: convert_units(
                    cosp_bin_sum(cld, 0, 3, 1.3, 9.4), target_units="%"
                ),
            ),
        ]
    ),
    "CLDLOW_TAU9.4_MISR": OrderedDict(
        [
            (
                ("CLD_MISR",),
                lambda cld: convert_units(
                    cosp_bin_sum(cld, 0, 3, 9.4, None), target_units="%"
                ),
            ),
            (
                ("CLMISR",),
                lambda cld: convert_units(
                    cosp_bin_sum(cld, 0, 3, 9.4, None), target_units="%"
                ),
            ),
        ]
    ),
    # COSP cloud fraction joint histogram
    "COSP_HISTOGRAM_MISR": OrderedDict(
        [
            (
                ("CLD_MISR",),
                lambda cld: cosp_histogram_standardize(rename(cld)),
            ),
            (("CLMISR",), lambda cld: cosp_histogram_standardize(rename(cld))),
        ]
    ),
    "COSP_HISTOGRAM_MODIS": OrderedDict(
        [
            (
                ("CLMODIS",),
                lambda cld: cosp_histogram_standardize(rename(cld)),
            ),
        ]
    ),
    "COSP_HISTOGRAM_ISCCP": OrderedDict(
        [
            (
                ("FISCCP1_COSP",),
                lambda cld: cosp_histogram_standardize(rename(cld)),
            ),
            (
                ("CLISCCP",),
                lambda cld: cosp_histogram_standardize(rename(cld)),
            ),
        ]
    ),
    "ICEFRAC": OrderedDict(
        [
            (
                ("ICEFRAC",),
                lambda icefrac: convert_units(icefrac, target_units="%"),
            )
        ]
    ),
    "RELHUM": OrderedDict(
        [
            (("hur",), lambda hur: convert_units(hur, target_units="%")),
            (
                ("RELHUM",),
                lambda relhum: convert_units(relhum, target_units="%"),
            )
            # (('RELHUM',), rename)
        ]
    ),
    "OMEGA": OrderedDict(
        [
            (
                ("wap",),
                lambda wap: convert_units(wap, target_units="mbar/day"),
            ),
            (
                ("OMEGA",),
                lambda omega: convert_units(omega, target_units="mbar/day"),
            ),
        ]
    ),
    "Q": OrderedDict(
        [
            (
                ("hus",),
                lambda q: convert_units(rename(q), target_units="g/kg"),
            ),
            (("Q",), lambda q: convert_units(rename(q), target_units="g/kg")),
            (("SHUM",), lambda shum: convert_units(shum, target_units="g/kg")),
        ]
    ),
    "H2OLNZ": OrderedDict(
        [
            (
                ("hus",),
                lambda q: convert_units(rename(q), target_units="g/kg"),
            ),
            (("H2OLNZ",), lambda h2o: w_convert_q(h2o)),
        ]
    ),
    "TAUXY": OrderedDict(
        [
            (("TAUX", "TAUY"), lambda taux, tauy: tauxy(taux, tauy)),
            (("tauu", "tauv"), lambda taux, tauy: tauxy(taux, tauy)),
        ]
    ),
    "AODVIS": OrderedDict(
        [
            (("od550aer",), rename),
            (
                ("AODVIS",),
                lambda aod: convert_units(rename(aod), target_units="dimensionless"),
            ),
            (
                ("AOD_550",),
                lambda aod: convert_units(rename(aod), target_units="dimensionless"),
            ),
            (
                ("TOTEXTTAU",),
                lambda aod: convert_units(rename(aod), target_units="dimensionless"),
            ),
            (
                ("AOD_550_ann",),
                lambda aod: convert_units(rename(aod), target_units="dimensionless"),
            ),
        ]
    ),
    "AODABS": OrderedDict([(("abs550aer",), rename)]),
    "AODDUST": OrderedDict(
        [
            (
                ("AODDUST",),
                lambda aod: convert_units(rename(aod), target_units="dimensionless"),
            )
        ]
    ),
    # Surface temperature: Degrees C
    # (Temperature of the surface (land/water) itself, not the air)
    "TS": OrderedDict([(("ts",), rename)]),
    "PS": OrderedDict([(("ps",), rename)]),
    "U10": OrderedDict([(("sfcWind",), rename)]),
    "QREFHT": OrderedDict([(("huss",), rename)]),
    "PRECC": OrderedDict([(("prc",), rename)]),
    "TAUX": OrderedDict([(("tauu",), lambda tauu: -tauu)]),
    "TAUY": OrderedDict([(("tauv",), lambda tauv: -tauv)]),
    "CLDICE": OrderedDict([(("cli",), rename)]),
    "TGCLDIWP": OrderedDict([(("clivi",), rename)]),
    "CLDLIQ": OrderedDict([(("clw",), rename)]),
    "TGCLDCWP": OrderedDict([(("clwvi",), rename)]),
    "O3": OrderedDict([(("o3",), rename)]),
    "PminusE": OrderedDict(
        [
            (("PminusE",), lambda pminuse: pminuse_convert_units(pminuse)),
            (
                (
                    "PRECC",
                    "PRECL",
                    "QFLX",
                ),
                lambda precc, precl, qflx: pminuse_convert_units(
                    prect(precc, precl) - pminuse_convert_units(qflx)
                ),
            ),
            (
                ("F_prec", "F_evap"),
                lambda pr, evspsbl: pminuse_convert_units(pr + evspsbl),
            ),
            (
                ("pr", "evspsbl"),
                lambda pr, evspsbl: pminuse_convert_units(pr - evspsbl),
            ),
        ]
    ),
    "TREFMNAV": OrderedDict(
        [
            (("TREFMNAV",), lambda t: convert_units(t, target_units="DegC")),
            (("tasmin",), lambda t: convert_units(t, target_units="DegC")),
        ]
    ),
    "TREFMXAV": OrderedDict(
        [
            (("TREFMXAV",), lambda t: convert_units(t, target_units="DegC")),
            (("tasmax",), lambda t: convert_units(t, target_units="DegC")),
        ]
    ),
    "TREF_range": OrderedDict(
        [
            (
                (
                    "TREFMXAV",
                    "TREFMNAV",
                ),
                lambda tmax, tmin: tref_range(tmax, tmin),
            ),
            (
                (
                    "tasmax",
                    "tasmin",
                ),
                lambda tmax, tmin: tref_range(tmax, tmin),
            ),
        ]
    ),
    "TCO": OrderedDict([(("TCO",), rename)]),
    "SCO": OrderedDict([(("SCO",), rename)]),
    "bc_DDF": OrderedDict(
        [
            (("bc_DDF",), rename),
            (
                (
                    "bc_a?DDF",
                    "bc_c?DDF",
                ),
                lambda *x: sum(x),
            ),
        ]
    ),
    "bc_SFWET": OrderedDict(
        [
            (("bc_SFWET",), rename),
            (
                (
                    "bc_a?SFWET",
                    "bc_c?SFWET",
                ),
                lambda *x: sum(x),
            ),
        ]
    ),
    "SFbc": OrderedDict(
        [
            (("SFbc",), rename),
            (("SFbc_a?",), lambda *x: sum(x)),
        ]
    ),
    "bc_CLXF": OrderedDict(
        [
            (("bc_CLXF",), rename),
            (("bc_a?_CLXF",), lambda *x: molec_convert_units(sum(x), 12.0)),
        ]
    ),
    "Mass_bc": OrderedDict(
        [
            (("Mass_bc",), rename),
        ]
    ),
    "dst_DDF": OrderedDict(
        [
            (("dst_DDF",), rename),
            (
                (
                    "dst_a?DDF",
                    "dst_c?DDF",
                ),
                lambda *x: sum(x),
            ),
        ]
    ),
    "dst_SFWET": OrderedDict(
        [
            (("dst_SFWET",), rename),
            (
                (
                    "dst_a?SFWET",
                    "dst_c?SFWET",
                ),
                lambda *x: sum(x),
            ),
        ]
    ),
    "SFdst": OrderedDict(
        [
            (("SFdst",), rename),
            (("SFdst_a?",), lambda *x: sum(x)),
        ]
    ),
    "Mass_dst": OrderedDict(
        [
            (("Mass_dst",), rename),
        ]
    ),
    "mom_DDF": OrderedDict(
        [
            (("mom_DDF",), rename),
            (
                (
                    "mom_a?DDF",
                    "mom_c?DDF",
                ),
                lambda *x: sum(x),
            ),
        ]
    ),
    "mom_SFWET": OrderedDict(
        [
            (("mom_SFWET",), rename),
            (
                (
                    "mom_a?SFWET",
                    "mom_c?SFWET",
                ),
                lambda *x: sum(x),
            ),
        ]
    ),
    "SFmom": OrderedDict(
        [
            (("SFmom",), rename),
            (("SFmom_a?",), lambda *x: sum(x)),
        ]
    ),
    "Mass_mom": OrderedDict(
        [
            (("Mass_mom",), rename),
        ]
    ),
    "ncl_DDF": OrderedDict(
        [
            (("ncl_DDF",), rename),
            (
                (
                    "ncl_a?DDF",
                    "ncl_c?DDF",
                ),
                lambda *x: sum(x),
            ),
        ]
    ),
    "ncl_SFWET": OrderedDict(
        [
            (("ncl_SFWET",), rename),
            (
                (
                    "ncl_a?SFWET",
                    "ncl_c?SFWET",
                ),
                lambda *x: sum(x),
            ),
        ]
    ),
    "SFncl": OrderedDict(
        [
            (("SFncl",), rename),
            (("SFncl_a?",), lambda *x: sum(x)),
        ]
    ),
    "Mass_ncl": OrderedDict(
        [
            (("Mass_ncl",), rename),
        ]
    ),
    "so4_DDF": OrderedDict(
        [
            (("so4_DDF",), rename),
            (
                (
                    "so4_a?DDF",
                    "so4_c?DDF",
                ),
                lambda *x: sum(x),
            ),
        ]
    ),
    "so4_SFWET": OrderedDict(
        [
            (("so4_SFWET",), rename),
            (
                (
                    "so4_a?SFWET",
                    "so4_c?SFWET",
                ),
                lambda *x: sum(x),
            ),
        ]
    ),
    "so4_CLXF": OrderedDict(
        [
            (("so4_CLXF",), rename),
            (
                ("so4_a?_CLXF",),
                lambda *x: molec_convert_units(sum(x), 115.0),
            ),
        ]
    ),
    "SFso4": OrderedDict(
        [
            (("SFso4",), rename),
            (("SFso4_a?",), lambda *x: sum(x)),
        ]
    ),
    "Mass_so4": OrderedDict(
        [
            (("Mass_so4",), rename),
        ]
    ),
    "soa_DDF": OrderedDict(
        [
            (("soa_DDF",), rename),
            (
                (
                    "soa_a?DDF",
                    "soa_c?DDF",
                ),
                lambda *x: sum(x),
            ),
        ]
    ),
    "soa_SFWET": OrderedDict(
        [
            (("soa_SFWET",), rename),
            (
                (
                    "soa_a?SFWET",
                    "soa_c?SFWET",
                ),
                lambda *x: sum(x),
            ),
        ]
    ),
    "SFsoa": OrderedDict(
        [
            (("SFsoa",), rename),
            (("SFsoa_a?",), lambda *x: sum(x)),
        ]
    ),
    "Mass_soa": OrderedDict(
        [
            (("Mass_soa",), rename),
        ]
    ),
    "pom_DDF": OrderedDict(
        [
            (("pom_DDF",), rename),
            (
                (
                    "pom_a?DDF",
                    "pom_c?DDF",
                ),
                lambda *x: sum(x),
            ),
        ]
    ),
    "pom_SFWET": OrderedDict(
        [
            (("pom_SFWET",), rename),
            (
                (
                    "pom_a?SFWET",
                    "pom_c?SFWET",
                ),
                lambda *x: sum(x),
            ),
        ]
    ),
    "SFpom": OrderedDict(
        [
            (("SFpom",), rename),
            (("SFpom_a?",), lambda *x: sum(x)),
        ]
    ),
    "pom_CLXF": OrderedDict(
        [
            (("pom_CLXF",), rename),
            (("pom_a?_CLXF",), lambda *x: molec_convert_units(sum(x), 12.0)),
        ]
    ),
    "Mass_pom": OrderedDict(
        [
            (("Mass_pom",), rename),
        ]
    ),
    # Land variables
    "SOILWATER_10CM": OrderedDict([(("mrsos",), rename)]),
    "SOILWATER_SUM": OrderedDict([(("mrso",), rename)]),
    "SOILICE_SUM": OrderedDict([(("mrfso",), rename)]),
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
    "tauuo": OrderedDict([(("tauuo",), rename)]),
    "tos": OrderedDict([(("tos",), rename)]),
    "thetaoga": OrderedDict([(("thetaoga",), rename)]),
    "hfsifrazil": OrderedDict([(("hfsifrazil",), rename)]),
    "sos": OrderedDict([(("sos",), rename)]),
    "soga": OrderedDict([(("soga",), rename)]),
    "tosga": OrderedDict([(("tosga",), rename)]),
    "wo": OrderedDict([(("wo",), rename)]),
    "thetao": OrderedDict([(("thetao",), rename)]),
    "masscello": OrderedDict([(("masscello",), rename)]),
    "wfo": OrderedDict([(("wfo",), rename)]),
    "tauvo": OrderedDict([(("tauvo",), rename)]),
    "vo": OrderedDict([(("vo",), rename)]),
    "hfds": OrderedDict([(("hfds",), rename)]),
    "volo": OrderedDict([(("volo",), rename)]),
    "uo": OrderedDict([(("uo",), rename)]),
    "zos": OrderedDict([(("zos",), rename)]),
    "tob": OrderedDict([(("tob",), rename)]),
    "sosga": OrderedDict([(("sosga",), rename)]),
    "sfdsi": OrderedDict([(("sfdsi",), rename)]),
    "zhalfo": OrderedDict([(("zhalfo",), rename)]),
    "masso": OrderedDict([(("masso",), rename)]),
    "so": OrderedDict([(("so",), rename)]),
    "sob": OrderedDict([(("sob",), rename)]),
    "mlotst": OrderedDict([(("mlotst",), rename)]),
    "fsitherm": OrderedDict([(("fsitherm",), rename)]),
    "msftmz": OrderedDict([(("msftmz",), rename)]),
    # sea ice variables
    "sitimefrac": OrderedDict([(("sitimefrac",), rename)]),
    "siconc": OrderedDict([(("siconc",), rename)]),
    "sisnmass": OrderedDict([(("sisnmass",), rename)]),
    "sisnthick": OrderedDict([(("sisnthick",), rename)]),
    "simass": OrderedDict([(("simass",), rename)]),
    "sithick": OrderedDict([(("sithick",), rename)]),
    "siu": OrderedDict([(("siu",), rename)]),
    "sitemptop": OrderedDict([(("sitemptop",), rename)]),
    "siv": OrderedDict([(("siv",), rename)]),
}

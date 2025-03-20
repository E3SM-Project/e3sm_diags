# NOTE: This module uses the deprecated e3sm_diags.driver.utils.dataset.Dataset
# class, which was replaced by e3sm_diags.driver.utils.dataset_xr.Dataset.
from __future__ import annotations

import glob
import json
import os
from typing import TYPE_CHECKING  # , Optional

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from matplotlib.colors import BoundaryNorm, ListedColormap
from tropical_subseasonal_plot import plot
from tropical_subseasonal_viewer import create_viewer
from zwf import zwf_functions as wf

import e3sm_diags
from e3sm_diags.driver import utils
from e3sm_diags.logger import _setup_child_logger
from e3sm_diags.plot.mp_partition_plot import plot

logger = _setup_child_logger(__name__)

# Script to compute and plot spectral powers of a subseasonal tropical field in
#   zonal wavenumber-frequency space.  Both the plot files and files containing the
#   associated numerical data shown in the plots are created.

# Authors: Jim Benedict and Brian Medeiros
# Modified by Jill Zhang to integrate into E3SM Diags.


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]
    """Return index of [array] closest in value to [value]
    Example:
    array = [ 0.21069679  0.61290182  0.63425412  0.84635244  0.91599191  0.00213826
              0.17104965  0.56874386  0.57319379  0.28719469]
    print(find_nearest(array, value=0.5))
    # 0.568743859261

    """


def wf_analysis(x, **kwargs):
    """Return zonal wavenumber-frequency power spectra of x.  The returned spectra are:
    spec_sym:    Raw (non-normalized) power spectrum of the component of x that is symmetric about the equator.
    spec_asym:   Raw (non-normalized) power spectrum of the component of x that is antisymmetric about the equator.
    nspec_sym:   Normalized (by a smoothed red-noise background spectrum) power spectrum of the component of x that is symmetric about the equator.
    nspec_asym:  Normalized (by a smoothed red-noise background spectrum) power spectrum of the component of x that is antisymmetric about the equator.

    The NCL version of 'wkSpaceTime' smooths the symmetric and antisymmetric components
    along the frequency dimension using a 1-2-1 filter once.

    """
    # Get the "raw" spectral power
    # OPTIONAL kwargs:
    # segsize, noverlap, spd, latitude_bounds (tuple: (south, north)), dosymmetries, rmvLowFrq

    z2 = wf.spacetime_power(x, **kwargs)
    z2avg = z2.mean(dim="component")
    z2.loc[{"frequency": 0}] = np.nan  # get rid of spurious power at \nu = 0 (mean)

    # Following NCL's wkSpaceTime, apply one pass of a 1-2-1 filter along the frequency
    #   domain to the raw (non-normalized) spectra/um.
    #   Do not use 0 frequency when smoothing here.
    #   Use weights that sum to 1 to ensure that smoothing is conservative.
    z2s = wf.smoothFrq121(z2, 1)

    # The background is supposed to be derived from both symmetric & antisymmetric
    # Inputs to the background spectrum calculation should be z2avg
    background = wf.smoothBackground_wavefreq(z2avg)
    # separate components
    spec_sym = z2s[0, ...]
    spec_asy = z2s[1, ...]
    # normalize:  Following NCL's wkSpaceTime, use lightly smoothed version of spectra/um
    #             as numerator
    nspec_sym = spec_sym / background
    nspec_asy = spec_asy / background

    spec = xr.merge(
        [
            spec_sym.rename("spec_raw_sym"),
            spec_asy.rename("spec_raw_asy"),
            nspec_sym.rename("spec_norm_sym"),
            nspec_asy.rename("spec_norm_asy"),
            background.rename("spec_background"),
        ],
        compat="override",
    )
    spec_all = spec.drop("component")
    spec_all["spec_raw_sym"].attrs = {"component": "symmetric", "type": "raw"}
    spec_all["spec_raw_asy"].attrs = {"component": "antisymmetric", "type": "raw"}
    spec_all["spec_norm_sym"].attrs = {"component": "symmetric", "type": "normalized"}
    spec_all["spec_norm_asy"].attrs = {
        "component": "antisymmetric",
        "type": "normalized",
    }
    spec_all["spec_background"].attrs = {"component": "", "type": "background"}

    return spec_all


def calculate_spectrum(path, variable):
    var = xr.open_mfdataset(glob.glob(f"{test_data_path}/{variable}_*.nc")).sel(
        lat=slice(-15, 15)
    )[variable]

    # TODO: subset time

    # Unit conversion
    if var.name == "PRECT":
        if var.attrs["units"] == "m/s" or var.attrs["units"] == "m s{-1}":
            print(
                "\nBEFORE unit conversion: Max/min of data: "
                + str(var.values.max())
                + "   "
                + str(var.values.min())
            )
            var.values = (
                var.values * 1000.0 * 86400.0
            )  # convert m/s to mm/d, do not alter metadata (yet)
            var.attrs["units"] = "mm/d"  # adjust metadata to reflect change in units
            print(
                "\nAFTER unit conversion: Max/min of data: "
                + str(var.values.max())
                + "   "
                + str(var.values.min())
            )
    if var.name == "precipAvg":
        if var.attrs["units"] == "mm/hr":
            print(
                "\nBEFORE unit conversion: Max/min of data: "
                + str(var.values.max())
                + "   "
                + str(var.values.min())
            )
            var.values = (
                data.values * 24.0
            )  # convert mm/hr to mm/d, do not alter metadata (yet)
            var.attrs["units"] = "mm/d"  # adjust metadata to reflect change in units
            print(
                "\nAFTER unit conversion: Max/min of data: "
                + str(var.values.max())
                + "   "
                + str(var.values.min())
            )

    # Wavenumber Frequency Analysis
    spec_all = wf_analysis(var, **opt)
    # spec_all.to_netcdf(outDataDir + "/full_spec.nc")
    return spec_all


#
# Options ... right now these only go into wk.spacetime_power()
#
do_zooming = False  # Set to True to also make plots to zoom into MJO spectral region,
#   in addition to the default (larger) spectral region
latBound = (-15, 15)  # latitude bounds for analysis
spd = 1  # SAMPLES PER DAY
nDayWin = 96  # Wheeler-Kiladis [WK] temporal window length (days)
nDaySkip = -60  # time (days) between temporal windows [segments]
# negative means there will be overlapping temporal segments
twoMonthOverlap = -1 * nDaySkip

vari = "PRECT"
srcID = "model"
outDataDir = "/global/cfs/cdirs/e3sm/www/chengzhu/tests/tropical_diags"
outDataDir = "/Users/zhang40/Documents/repos/e3sm_diags/auxiliary_tools/tropical_subseasonal_diags/data"

opt = {
    "segsize": nDayWin,
    "noverlap": twoMonthOverlap,
    "spd": spd,
    "latitude_bounds": latBound,
    "dosymmetries": True,
    "rmvLowFrq": True,
}


# datapath = '/global/cfs/cdirs/e3sm/forsyth/E3SMv2/v2.LR.historical_0201/post/atm/180x360_aave/ts/daily/5yr'
datapath = "/Users/zhang40/Documents/e3sm_diags_data/e3sm_diags_test_data/E3SM_v2_daily"

from e3sm_diags.parameter.core_parameter import CoreParameter

parameter = CoreParameter()

test_data_path = (
    "/Users/zhang40/Documents/e3sm_diags_data/e3sm_diags_test_data/E3SM_v2_daily"
)
parameter.test_data_path = test_data_path
parameter.test_timeseries_input = True
parameter.test_start_yr = "2000"
parameter.test_end_yr = "2014"
parameter.ref_data_path = test_data_path
parameter.ref_timeseries_input = True
parameter.ref_start_yr = "2000"
parameter.ref_end_yr = "2014"
parameter.variables = ["PRECT"]
season = "ANN"

test_data = utils.dataset.Dataset(parameter, test=True)
parameter.test_name_yrs = utils.general.get_name_and_yrs(parameter, test_data, season)

ref_data = utils.dataset.Dataset(parameter, ref=True)
parameter.ref_name_yrs = utils.general.get_name_and_yrs(parameter, ref_data, season)

for variable in parameter.variables:
    # test = calculate_spectrum(parameter.test_data_path, variable)
    # test.to_netcdf("data/full_spec_test.nc")
    # ref = calculate_spectrum(parameter.ref_data_path, variable)
    # ref.to_netcdf("data/full_spec_ref.nc")
    # Below uses intermediate saved files for development
    test = xr.open_dataset(
        "/Users/zhang40/Documents/repos/e3sm_diags/auxiliary_tools/tropical_subseasonal_diags/data/full_spec_ref.nc"
    ).load()
    ref = xr.open_dataset(
        "/Users/zhang40/Documents/repos/e3sm_diags/auxiliary_tools/tropical_subseasonal_diags/data/full_spec_ref.nc"
    ).load()
    parameter.var_id = variable

    for diff_name in ["raw_sym", "raw_asy", "norm_sym", "norm_asy", "background"]:
        # Compute percentage difference
        diff = (
            100
            * (test[f"spec_{diff_name}"] - ref[f"spec_{diff_name}"])
            / ref[f"spec_{diff_name}"]
        )
        diff.name = f"spec_{diff_name}"
        diff.attrs.update(test[f"spec_{diff_name}"].attrs)
        parameter.spec_type = diff_name
        plot(parameter, test[f"spec_{diff_name}"], ref[f"spec_{diff_name}"], diff)
        if "norm" in diff_name:
            parameter.spec_type = f"{diff_name}_zoom"
            plot(
                parameter,
                test[f"spec_{diff_name}"],
                ref[f"spec_{diff_name}"],
                diff,
                do_zoom=True,
            )


display_name, url = create_viewer(".", parameter)
print("Viewer Created: ", url)

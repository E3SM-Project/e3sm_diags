from __future__ import annotations

import glob
from typing import TYPE_CHECKING  # , Optional

import numpy as np
import xarray as xr

from e3sm_diags.driver import utils
from e3sm_diags.driver.utils import zwf_functions as wf
from e3sm_diags.logger import custom_logger
from e3sm_diags.plot.cartopy.tropical_subseasonal_plot import plot

if TYPE_CHECKING:
    from e3sm_diags.parameter.tropical_subseasonal_parameter import (
        TropicalSubseasonalParameter,
    )


logger = custom_logger(__name__)

#   Script to compute and plot spectral powers of a subseasonal tropical field in
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


def calculate_spectrum(path, variable, start_year, end_year):
    latBound = (-15, 15)  # latitude bounds for analysis
    spd = 1  # SAMPLES PER DAY
    nDayWin = 96  # Wheeler-Kiladis [WK] temporal window length (days)
    nDaySkip = -60  # time (days) between temporal windows [segments]
    # negative means there will be overlapping temporal segments
    twoMonthOverlap = -1 * nDaySkip

    opt = {
        "segsize": nDayWin,
        "noverlap": twoMonthOverlap,
        "spd": spd,
        "latitude_bounds": latBound,
        "dosymmetries": True,
        "rmvLowFrq": True,
    }

    var = xr.open_mfdataset(glob.glob(f"{path}/{variable}_*.nc")).sel(
        lat=slice(-15, 15), time=slice(f"{start_year}-01-01", f"{end_year}-12-31")
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

    # Wavenumber Frequency Analysis
    spec_all = wf_analysis(var, **opt)
    # spec_all.to_netcdf(outDataDir + "/full_spec.nc")
    return spec_all


def run_diag(parameter: TropicalSubseasonalParameter) -> TropicalSubseasonalParameter:
    """Runs the wavenumber frequency analysis for tropical variability.

    :param parameter: Parameters for the run
    :type parameter: CoreParameter
    :raises ValueError: Invalid run type
    :return: Parameters for the run
    :rtype: CoreParameter
    """
    run_type = parameter.run_type
    # variables = parameter.variables
    season = "ANN"

    test_data = utils.dataset.Dataset(parameter, test=True)
    parameter.test_name_yrs = utils.general.get_name_and_yrs(
        parameter, test_data, season
    )

    ref_data = utils.dataset.Dataset(parameter, ref=True)
    parameter.ref_name_yrs = utils.general.get_name_and_yrs(parameter, ref_data, season)

    for variable in parameter.variables:
        test = calculate_spectrum(
            parameter.test_data_path,
            variable,
            parameter.test_start_yr,
            parameter.test_end_yr,
        )
        test.to_netcdf(f"{parameter.results_dir}/full_spec_test.nc")
        if run_type == "model_vs_model":
            ref = calculate_spectrum(
                parameter.reference_data_path,
                variable,
                parameter.ref_start_yr,
                parameter.ref_end_yr,
            )
        elif run_type == "model_vs_obs":
            if parameter.ref_start_yr == "":
                parameter.ref_name_yrs = parameter.reference_name
                # read precalculated data.
            else:
                ref_data_path = f"{parameter.reference_data_path}/{parameter.ref_name}"
                # parameter.ref_name_yrs = f"{parameter.ref_name}({parameter.test_start_yr}-{parameter.test_start_yr})"
                ref = calculate_spectrum(
                    ref_data_path,
                    variable,
                    parameter.ref_start_yr,
                    parameter.ref_end_yr,
                )
                ref.to_netcdf(
                    f"{parameter.results_dir}/full_spec_ref_{parameter.ref_name}.nc"
                )
        # ref = calculate_spectrum(parameter.test_data_path, variable)
        # test = xr.open_dataset(f"{parameter.results_dir}/full_spec_test.nc").load()
        # ref.to_netcdf(f"{parameter.results_dir}/full_spec_ref.nc")
        # ref = xr.open_dataset(f"{parameter.results_dir}/full_spec_ref.nc").load()
        # TODO save to netcdf
        parameter.var_id = variable
        for diff_name in ["raw_sym", "raw_asy", "norm_sym", "norm_asy", "background"]:
            # Compute percentage difference
            diff = (
                100
                * (test[f"spec_{diff_name}"] - ref[f"spec_{diff_name}"])
                / ref[f"spec_{diff_name}"]
            )
            print("diff_name****888")
            print(diff)
            diff.name = f"spec_{diff_name}"
            diff.attrs.update(test[f"spec_{diff_name}"].attrs)
            parameter.spec_type = diff_name
            parameter.output_file = f"{parameter.var_id}_{parameter.spec_type}_15N-15S"
            parameter.diff_title = "percent difference"
            plot(parameter, test[f"spec_{diff_name}"], ref[f"spec_{diff_name}"], diff)
            if "norm" in diff_name:
                parameter.spec_type = f"{diff_name}_zoom"
                parameter.output_file = (
                    f"{parameter.var_id}_{parameter.spec_type}_15N-15S"
                )
                plot(
                    parameter,
                    test[f"spec_{diff_name}"],
                    ref[f"spec_{diff_name}"],
                    diff,
                    do_zoom=True,
                )

    return parameter

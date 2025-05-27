"""
Script to compute and plot spectral powers of a subseasonal tropical field in
zonal wavenumber-frequency space.  Both the plot files and files containing the
associated numerical data shown in the plots are created.

Authors: Jim Benedict and Brian Medeiros
Modified by Jill Zhang to integrate into E3SM Diags.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Tuple

import numpy as np
import xarray as xr

from e3sm_diags.driver.utils import zwf_functions as wf
from e3sm_diags.driver.utils.climo_xr import ClimoFreq
from e3sm_diags.driver.utils.dataset_xr import Dataset
from e3sm_diags.logger import _setup_child_logger
from e3sm_diags.plot.tropical_subseasonal_plot import plot

if TYPE_CHECKING:
    from e3sm_diags.parameter.tropical_subseasonal_parameter import (
        TropicalSubseasonalParameter,
    )


logger = _setup_child_logger(__name__)


def run_diag(parameter: TropicalSubseasonalParameter) -> TropicalSubseasonalParameter:
    """Runs the wavenumber frequency analysis for tropical variability.

    :param parameter: Parameters for the run
    :type parameter: CoreParameter
    :raises ValueError: Invalid run type
    :return: Parameters for the run
    :rtype: CoreParameter
    """
    run_type = parameter.run_type
    season: ClimoFreq = "ANN"

    test_data = Dataset(parameter, data_type="test")
    ref_data = Dataset(parameter, data_type="ref")

    for variable in parameter.variables:
        # Get test dataset
        test_ds = test_data.get_time_series_dataset(variable, single_point=True)
        test_spectrum, test_start, test_end = calculate_spectrum(test_ds, variable)

        # Update parameters with actual time range
        parameter.test_start_yr = test_start
        parameter.test_end_yr = test_end
        parameter.test_name_yrs = test_data.get_name_yrs_attr(season)

        # Get reference dataset
        if run_type == "model_vs_model":
            ref_ds = ref_data.get_time_series_dataset(variable, single_point=True)
            ref_spectrum, ref_start, ref_end = calculate_spectrum(ref_ds, variable)
        elif run_type == "model_vs_obs":
            # TODO use pre-calculated spectral power
            # if parameter.ref_start_yr == "":
            #    parameter.ref_name_yrs = parameter.reference_name
            #    # read precalculated data.
            # else:
            ref_ds = ref_data.get_time_series_dataset(variable, single_point=True)
            ref_spectrum, ref_start, ref_end = calculate_spectrum(ref_ds, variable)

        parameter.ref_start_yr = ref_start
        parameter.ref_end_yr = ref_end
        parameter.ref_name_yrs = ref_data.get_name_yrs_attr(season)
        parameter.var_id = variable

        for diff_name in ["raw_sym", "raw_asy", "norm_sym", "norm_asy", "background"]:
            diff = (
                100
                * (
                    test_spectrum[f"spec_{diff_name}"]
                    - ref_spectrum[f"spec_{diff_name}"]
                )
                / ref_spectrum[f"spec_{diff_name}"]
            )
            diff.name = f"spec_{diff_name}"
            diff.attrs.update(test_spectrum[f"spec_{diff_name}"].attrs)

            parameter.spec_type = diff_name
            parameter.output_file = f"{parameter.var_id}_{parameter.spec_type}_15N-15S"
            parameter.diff_title = "percent difference"
            plot(
                parameter,
                test_spectrum[f"spec_{diff_name}"],
                ref_spectrum[f"spec_{diff_name}"],
                diff,
            )

            if "norm" in diff_name:
                parameter.spec_type = f"{diff_name}_zoom"
                parameter.output_file = (
                    f"{parameter.var_id}_{parameter.spec_type}_15N-15S"
                )
                plot(
                    parameter,
                    test_spectrum[f"spec_{diff_name}"],
                    ref_spectrum[f"spec_{diff_name}"],
                    diff,
                    do_zoom=True,
                )

    return parameter


def calculate_spectrum(ds: xr.Dataset, variable: str) -> Tuple[xr.Dataset, str, str]:
    """Calculate wavenumber-frequency power spectra for a variable.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset containing the variable of interest
    variable : str
        Name of the variable to analyze

    Returns
    -------
    Tuple[xr.Dataset, str, str]
        Tuple containing:
        - Dataset with spectral power components
        - Start year (as string)
        - End year (as string)
    """
    # latitude bounds for analysis
    latBound = (-15, 15)
    # SAMPLES PER DAY
    spd = 1
    # Wheeler-Kiladis [WK] temporal window length (days)
    nDayWin = 96
    # time (days) between temporal windows [segments]
    nDaySkip = -60
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

    # Extract variable data from dataset and subset to tropical latitudes
    var = ds[variable].sel(lat=slice(-15, 15))

    # Get actual time range from the data
    actual_start = var.time.dt.year.values[0]
    actual_end = var.time.dt.year.values[-1]

    # Wavenumber Frequency Analysis
    spec_all = wf_analysis(var, **opt)
    return spec_all, str(actual_start), str(actual_end)


def wf_analysis(x, **kwargs):
    """Calculate zonal wavenumber-frequency power spectra of x.

    Input:
    ----------
    x : xr.DataArray
        The input variable with dimension (ntime, nlat, nlon)
    **kwargs : dict
        A dictionary providing options for spectrum calculation
        Example:
        # Latitude bounds
        latBound = (-15, 15)
        # SAMPLES PER DAY
        spd = 1
        # Wheeler-Kiladis [WK] temporal window length (days)
        nDayWin = 96
        # time (days) between temporal windows [segments]
        nDaySkip = -60
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

    Returns
    -------
    spectra: xr.DataArray
    The returned spectra are:
    spec_raw_sym:    Raw (non-normalized) power spectrum of the component of x that is symmetric about the equator.
    spec_raw_asy:   Raw (non-normalized) power spectrum of the component of x that is antisymmetric about the equator.
    spec_norm_sym:   Normalized (by a smoothed red-noise background spectrum) power spectrum of the component of x that is symmetric about the equator.
    spec_norm_asy:  Normalized (by a smoothed red-noise background spectrum) power spectrum of the component of x that is antisymmetric about the equator.

    The NCL version of 'wkSpaceTime' smooths the symmetric and antisymmetric components
    along the frequency dimension using a 1-2-1 filter once.

    """
    # Get the "raw" spectral power
    # OPTIONAL kwargs:
    # segsize, noverlap, spd, latitude_bounds (tuple: (south, north)), dosymmetries, rmvLowFrq

    # Interpolate missing values along longitude before spectral analysis
    if np.any(np.isnan(x)):
        logger.info(
            "Interpolating missing values along longitude before spectral analysis"
        )
        x = x.interpolate_na(dim="lon", method="linear", fill_value="extrapolate")

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

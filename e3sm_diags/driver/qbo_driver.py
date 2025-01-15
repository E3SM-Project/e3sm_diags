from __future__ import annotations

import json
import os
from typing import TYPE_CHECKING, Dict, Literal, Tuple, TypedDict

import numpy as np
import pywt
import scipy.fftpack
import xarray as xr
import xcdat as xc
from scipy.signal import detrend

from e3sm_diags.driver.utils.dataset_xr import Dataset
from e3sm_diags.driver.utils.io import _get_output_dir, _write_to_netcdf
from e3sm_diags.driver.utils.regrid import _subset_on_region
from e3sm_diags.logger import custom_logger
from e3sm_diags.metrics.metrics import spatial_avg
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.plot.qbo_plot import plot

logger = custom_logger(__name__)

if TYPE_CHECKING:
    from e3sm_diags.parameter.qbo_parameter import QboParameter

# The region will always be 5S5N
REGION = "5S5N"

# Target power spectral vertical level for the wavelet diagnostic.
POW_SPEC_LEV = 20.0


class MetricsDict(TypedDict):
    qbo: xr.DataArray
    psd_sum: np.ndarray
    amplitude: np.ndarray
    period_new: np.ndarray
    psd_x_new: np.ndarray
    amplitude_new: np.ndarray
    wave_period: np.ndarray
    wavelet: np.ndarray
    name: str


def run_diag(parameter: QboParameter) -> QboParameter:
    variables = parameter.variables

    test_ds = Dataset(parameter, data_type="test")
    ref_ds = Dataset(parameter, data_type="ref")

    for var_key in variables:
        logger.info(f"Variable={var_key}")

        ds_test = test_ds.get_time_series_dataset(var_key)
        ds_ref = ref_ds.get_time_series_dataset(var_key)

        ds_test_region = _subset_on_region(ds_test, var_key, REGION)
        ds_ref_region = _subset_on_region(ds_ref, var_key, REGION)

        # Convert plevs of test and ref for unified units and direction
        ds_test_region = _unify_plev(ds_test_region, var_key)
        ds_ref_region = _unify_plev(ds_ref_region, var_key)

        # Dictionaries to store information on the variable including the name,
        # the averaged variable, and metrics.
        test_dict: MetricsDict = {}  # type: ignore
        ref_dict: MetricsDict = {}  # type: ignore

        # Diagnostic 1: average over longitude & latitude to produce time-height
        # array of u field.
        test_dict["qbo"] = _spatial_avg(ds_test_region, var_key)
        ref_dict["qbo"] = _spatial_avg(ds_ref_region, var_key)

        # Diagnostic 2: calculate and plot the amplitude of wind variations with a 20-40 month period
        test_dict["psd_sum"], test_dict["amplitude"] = _get_20to40month_fft_amplitude(
            test_dict["qbo"]
        )
        ref_dict["psd_sum"], ref_dict["amplitude"] = _get_20to40month_fft_amplitude(
            ref_dict["qbo"]
        )

        # Diagnostic 3: calculate the Power Spectral Density.
        # Pre-process data to average over lat,lon,height
        x_test = _get_power_spectral_density(test_dict["qbo"])
        x_ref = _get_power_spectral_density(ref_dict["qbo"])

        # Calculate the PSD and interpolate to period_new. Specify periods to
        # plot.
        test_dict["period_new"] = ref_dict["period_new"] = np.concatenate(
            (np.arange(2.0, 33.0), np.arange(34.0, 100.0, 2.0)), axis=0
        )
        test_dict["psd_x_new"], test_dict["amplitude_new"] = _get_psd_from_deseason(
            x_test, test_dict["period_new"]
        )
        ref_dict["psd_x_new"], ref_dict["amplitude_new"] = _get_psd_from_deseason(
            x_ref, ref_dict["period_new"]
        )

        # Diagnostic 4: calculate the Wavelet
        test_dict["wave_period"], test_dict["wavelet"] = _calculate_wavelet(
            test_dict["qbo"]
        )
        ref_dict["wave_period"], ref_dict["wavelet"] = _calculate_wavelet(
            ref_dict["qbo"]
        )

        parameter.var_id = var_key
        parameter.output_file = "qbo_diags"
        parameter.main_title = (
            f"QBO index, amplitude, and power spectral density for {var_key}"
        )
        # Get the years of the data.
        parameter.viewer_descr[var_key] = parameter.main_title

        parameter.test_yrs = f"{test_ds.start_yr}-{test_ds.end_yr}"
        parameter.ref_yrs = f"{ref_ds.start_yr}-{ref_ds.end_yr}"

        # Write the averaged variables to netCDF. Save with data type as
        # `qbo_test` and `qbo_ref` to match CDAT codebase for regression
        # testing of the `.nc` files.
        _write_to_netcdf(parameter, test_dict["qbo"], var_key, "qbo_test")  # type: ignore
        _write_to_netcdf(parameter, ref_dict["qbo"], var_key, "qbo_ref")  # type: ignore

        # Write the metrics to .json files.
        test_dict["name"] = test_ds._get_test_name()
        ref_dict["name"] = ref_ds._get_ref_name()

        _save_metrics_to_json(parameter, test_dict, "test")  # type: ignore
        _save_metrics_to_json(parameter, ref_dict, "ref")  # type: ignore

        # plot the results.
        plot(parameter, test_dict, ref_dict)

    return parameter


def _save_metrics_to_json(
    parameter: CoreParameter,
    var_dict: Dict[str, str | np.ndarray],
    dict_type: Literal["test", "ref"],
):
    output_dir = _get_output_dir(parameter)
    filename = parameter.output_file + f"_{dict_type}.json"
    abs_path = os.path.join(output_dir, filename)

    # Convert all metrics from `np.ndarray` to a Python list for serialization
    # to `.json`.
    metrics_dict = {k: v for k, v in var_dict.items() if k != "qbo"}

    for key in metrics_dict.keys():
        if key != "name":
            metrics_dict[key] = metrics_dict[key].tolist()  # type: ignore

    with open(abs_path, "w") as outfile:
        json.dump(metrics_dict, outfile, default=str)

    logger.info("Metrics saved in: {}".format(abs_path))


def _unify_plev(ds_region: xr.Dataset, var_key: str) -> xr.Dataset:
    """Convert the Z-axis (plev) with units Pa to hPa.

    This function also orders the data by plev in ascending order (same as model
    data).

    Parameters
    ----------
    ds_region : xr.Dataset
        The dataset for the region.

    Returns
    -------
    xr.Dataset
        The dataset for the region with a converted Z-axis.
    """
    ds_region_new = ds_region.copy()
    # The dataset can have multiple Z axes (e.g., "lev", "ilev"), so get the
    # Z axis from the variable directly.
    z_axis = xc.get_dim_coords(ds_region[var_key], axis="Z")
    z_dim = z_axis.name

    if z_axis.attrs["units"] == "Pa":
        ds_region_new[z_dim] = z_axis / 100.0
        ds_region_new[z_dim].attrs["units"] = "hPa"

    ds_region_new = ds_region_new.sortby(z_dim, ascending=True)

    return ds_region_new


def _spatial_avg(ds: xr.Dataset, var_key: str) -> xr.DataArray:
    """Process the U variable for time and height by averaging of lat and lon.

    Diagnostic 1: average over longitude & latitude to produce time-height
    array of u field.

    Richter, J. H., Chen, C. C., Tang, Q., Xie, S., & Rasch, P. J. (2019).
    Improved Simulation of the QBO in E3SMv1. Journal of Advances in Modeling
    Earth Systems, 11(11), 3403-3418.

    U = "Monthly mean zonal mean zonal wind averaged between 5S and 5N as a
    function of pressure and time" (p. 3406)

    Source: https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2019MS001763

    Parameters
    ----------
    ds_region : xr.Dataset
        The dataset.
    var_key : str
        The key of the variable.

    Returns
    -------
    xr.DataArray
        The averaged variable.
    """
    var_avg = spatial_avg(ds, var_key, axis=["X", "Y"], as_list=False)

    return var_avg  # type: ignore


def _get_20to40month_fft_amplitude(
    var_avg: xr.DataArray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Calculates the amplitude of wind variations in the 20 - 40 month period.

    Parameters
    ----------
    var_avg : xr.DataArray
        The spatially averaged variable.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        The psd and amplitude arrays.
    """
    qboN_arr = np.squeeze(var_avg.values)

    levelN = xc.get_dim_coords(var_avg, axis="Z")
    psd_sumN = np.zeros(levelN.shape)
    amplitudeN = np.zeros(levelN.shape)

    for ilev in np.arange(len(levelN)):
        # `qboN[:, ilev]` returns the entire 0th dimension for ilev in the 1st
        # dimension of the array.
        y_input = deseason(np.squeeze(qboN_arr[:, ilev]))
        y = scipy.fftpack.fft(y_input)
        n = len(y)

        frequency = np.arange(n / 2) / n
        period = 1 / frequency
        values = y[0 : int(np.floor(n / 2))]
        fyy = values * np.conj(values)

        # Choose the range 20 - 40 months that captures most QBOs (in nature).
        psd_sumN[ilev] = 2 * np.nansum(fyy[(period <= 40) & (period >= 20)])
        amplitudeN[ilev] = np.sqrt(2 * psd_sumN[ilev]) * (frequency[1] - frequency[0])

    return psd_sumN, amplitudeN


def _get_power_spectral_density(var_avg: xr.DataArray):
    # Average over vertical levels and horizontal area (units: hPa)
    level_bottom = 22
    level_top = 18

    z_dim = xc.get_dim_keys(var_avg, axis="Z")
    # Average over vertical
    try:
        average = var_avg.sel({z_dim: slice(level_top, level_bottom)})
    except Exception as err:
        raise Exception(
            "No levels found between {}hPa and {}hPa".format(level_top, level_bottom)
        ) from err

    x0 = np.nanmean(average.values, axis=1)

    # x0 should now be 1D
    return x0


def ceil_log2(x):
    """
    Given a number, calculate the exponent for the next power of 2.

    Example:
        ceil_log2(16) = 4
        ceil_log2(17) = 5
    """
    return np.ceil(np.log2(x)).astype("int")


def _get_psd_from_deseason(xraw, period_new):
    x_deseasoned = deseason(xraw)

    # Sampling frequency: assumes frequency of sampling = 1 month
    sampling_frequency = 1
    # Calculate the period as a function of frequency
    period0 = 1 / sampling_frequency
    L0 = len(xraw)
    NFFT0 = 2 ** ceil_log2(L0)

    # Apply fft on x_deseasoned with n = NFFT
    x0 = scipy.fftpack.fft(x_deseasoned, n=NFFT0) / L0
    # Frequency (cycles/month). Frequency will be increasing.
    frequency0 = sampling_frequency * np.arange(0, (NFFT0 / 2 + 1)) / NFFT0
    # Period (months/cycle). Calculate as a function of frequency. Thus, period will be decreasing.
    period0 = 1 / frequency0

    # Calculate amplitude as a function of frequency
    amplitude0 = 2 * abs(x0[0 : int(NFFT0 / 2 + 1)])
    # Calculate power spectral density as a function of frequency
    psd_x0 = amplitude0**2 / L0

    # Total spectral power
    # In the next code block, we will perform an interpolation using the period
    # (interpolating values of amplitude0_flipped and psd_x0_flipped from period0_flipped to period_new).
    # For that interpolation, we want the period to be increasing.
    # Therefore, we will flip the following values:
    period0_flipped = period0[::-1]  # type: ignore
    amplitude0_flipped = amplitude0[::-1]
    psd_x0_flipped = psd_x0[::-1]

    amplitude_new0 = np.interp(
        period_new, period0_flipped[:-1], amplitude0_flipped[:-1]
    )
    psd_x_new0 = np.interp(period_new, period0_flipped[:-1], psd_x0_flipped[:-1])

    return psd_x_new0, amplitude_new0


def deseason(xraw):
    # Calculates the deseasonalized data
    months_per_year = 12
    # Create array to hold climatological values and deseasonalized data
    # Create months_per_year x 1 array of zeros
    xclim = np.zeros((months_per_year, 1))
    # Create array with same shape as xraw
    x_deseasoned = np.zeros(xraw.shape)
    # Iterate through all 12 months.
    for month in np.arange(months_per_year):
        # `xraw[month::12]` will return the data for this month every year (12 months)
        # (i.e., from month until the end of xraw, get every 12th month)
        # Get the mean of this month, using data from every year, ignoring NaNs
        xclim[month] = np.nanmean(xraw[month::months_per_year])
    num_years = int(np.floor(len(x_deseasoned) / months_per_year))
    # Iterate through all years in x_deseasoned (same number as in xraw)
    for year in np.arange(num_years):
        year_index = year * months_per_year
        # Iterate through all months of the year
        for month in np.arange(months_per_year):
            month_index = year_index + month
            # Subtract the month's mean over num_years from xraw's data for this month in this year
            # i.e., get the difference between this month's value and it's "usual" value
            x_deseasoned[month_index] = xraw[month_index] - xclim[month]
    return x_deseasoned


def _calculate_wavelet(var: xr.DataArray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate the wavelet spectrum for a given data array at a specified power
    spectral level.

    Parameters
    ----------
    data : xr.DataArray
        The variable data.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        The wavelet period and wavelet array.
    """
    # Find the closest value for power spectral level in the list
    test_lev = xc.get_dim_coords(var, axis="Z")
    test_lev_list = list(test_lev)
    closest_lev = min(test_lev_list, key=lambda x: abs(x - POW_SPEC_LEV))
    closest_index = test_lev_list.index(closest_lev)

    # Grab target vertical level
    data_avg = var.values[:, closest_index]

    # Convert to anomalies
    data_avg = data_avg - data_avg.mean()

    # Detrend the data
    detrended_data = detrend(data_avg)

    wave_period, wavelet = _get_psd_from_wavelet(detrended_data)

    # Get square root values of wavelet spectra
    wavelet = np.sqrt(wavelet)

    return wave_period, wavelet


def _get_psd_from_wavelet(data: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate the power spectral density (PSD) of the data using a complex
    Mortlet wavelet spectrum of degree 6.

    Parameters
    ----------
    data : np.ndarray
        The data to calculate the PSD for.

    Returns
    -------
    Tuple(np.ndarray, np.ndarray)
        The period and PSD arrays.
    """
    deg = 6
    period = np.arange(1, 55 + 1)
    widths = deg / (2 * np.pi / period)

    [cfs, freq] = pywt.cwt(data, scales=widths, wavelet="cmor1.5-1.0")
    psd = np.mean(np.square(np.abs(cfs)), axis=1)
    period = 1 / freq

    return (period, psd)

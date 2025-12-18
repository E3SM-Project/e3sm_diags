"""
Driver for precipitation PDF diagnostics.

Computes and plots precipitation probability density functions (PDFs)
for various regions (tropics, CONUS, etc.) from daily precipitation data.

Based on original work by Chris Terai (2015, 2020).
Modified to integrate into E3SM Diags.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import xarray as xr
import xcdat as xc

from e3sm_diags.derivations.default_regions_xr import REGION_SPECS
from e3sm_diags.driver.utils.climo_xr import ClimoFreq
from e3sm_diags.driver.utils.dataset_xr import Dataset
from e3sm_diags.logger import _setup_child_logger
from e3sm_diags.plot.precip_pdf_plot import plot

if TYPE_CHECKING:
    from e3sm_diags.parameter.precip_pdf_parameter import PrecipPDFParameter


logger = _setup_child_logger(__name__)


def run_diag(parameter: PrecipPDFParameter) -> PrecipPDFParameter:
    """Runs the precipitation PDF analysis.

    Parameters
    ----------
    parameter : PrecipPDFParameter
        Parameters for the run

    Returns
    -------
    PrecipPDFParameter
        Parameters for the run with results
    """
    run_type = parameter.run_type
    season: ClimoFreq = "ANN"

    test_data = Dataset(parameter, data_type="test")
    ref_data = Dataset(parameter, data_type="ref")

    for variable in parameter.variables:
        # Get test dataset
        test_ds = test_data.get_time_series_dataset(variable, single_point=True)
        test_pdf, test_start, test_end = calculate_precip_pdf(test_ds, variable)

        # Update parameters with actual time range
        parameter.test_start_yr = test_start
        parameter.test_end_yr = test_end
        parameter.test_name_yrs = test_data.get_name_yrs_attr(season)

        # Get reference dataset
        if run_type == "model_vs_model":
            ref_ds = ref_data.get_time_series_dataset(variable, single_point=True)
            ref_pdf, ref_start, ref_end = calculate_precip_pdf(ref_ds, variable)
        elif run_type == "model_vs_obs":
            ref_ds = ref_data.get_time_series_dataset(variable, single_point=True)
            ref_pdf, ref_start, ref_end = calculate_precip_pdf(ref_ds, variable)

        parameter.ref_start_yr = ref_start
        parameter.ref_end_yr = ref_end
        parameter.ref_name_yrs = ref_data.get_name_yrs_attr(season)
        parameter.var_id = variable

        # Calculate regional PDFs and plot
        for region in parameter.regions:
            parameter.output_file = f"{parameter.var_id}_PDF_{region}"

            # Extract regional PDFs
            test_regional = extract_regional_pdf(test_pdf, region)
            ref_regional = extract_regional_pdf(ref_pdf, region)

            # Call plot function
            plot(parameter, test_regional, ref_regional, region)

    return parameter


def calculate_precip_pdf(ds: xr.Dataset, variable: str) -> tuple[xr.Dataset, str, str]:
    """Calculate precipitation PDFs from daily data.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset containing the variable of interest (already time-sliced)
    variable : str
        Variable name (typically "PRECT")

    Returns
    -------
    tuple[xr.Dataset, str, str]
        Tuple containing:
        - Dataset with frequency and amount PDFs
        - Start year (as string)
        - End year (as string)
    """
    # Extract variable data from dataset
    var = ds[variable]

    # Get actual time range from the data
    actual_start = var.time.dt.year.values[0]
    actual_end = var.time.dt.year.values[-1]

    logger.info(
        f"Loaded {variable} data from {actual_start} to {actual_end}"
    )

    # Unit conversion to mm/day
    if var.name == "PRECT":
        if var.attrs["units"] == "m/s" or var.attrs["units"] == "m s{-1}":
            logger.info(
                "\nBEFORE unit conversion: Max/min of data: "
                + str(var.values.max())
                + "   "
                + str(var.values.min())
            )
            var.values = (
                var.values * 1000.0 * 86400.0
            )  # convert m/s to mm/d, do not alter metadata (yet)
            var.attrs["units"] = "mm/d"  # adjust metadata to reflect change in units
            logger.info(
                "\nAFTER unit conversion: Max/min of data: "
                + str(var.values.max())
                + "   "
                + str(var.values.min())
            )

    # Create precipitation bins (following Terai's WC_diag_amwg.py)
    # Exponentially spaced bins: 0.1 * 1.07^n for n=0 to 129
    threshold = 0.1  # mm/day - minimum precipitation threshold
    bin_increase = 1.07  # multiplicative factor
    bin_edges = threshold * bin_increase ** np.arange(0, 130)  # 130 edges = 129 bins

    # Bin centers using arithmetic mean (Terai's bincentertype='arithmetic')
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0
    num_bins = len(bin_centers)

    # Calculate PDFs for each grid point
    freq_pdf, amnt_pdf = calculate_gridded_pdf(var, bin_edges, bin_centers)

    # Create output dataset
    pdf_ds = xr.Dataset(
        {
            "FREQPDF": (["lat", "lon", "bin"], freq_pdf),
            "AMNTPDF": (["lat", "lon", "bin"], amnt_pdf),
            "bin_centers": (["bin"], bin_centers),
            "bin_edges": (["bin_edge"], bin_edges),
        },
        coords={
            "lat": var.lat,
            "lon": var.lon,
            "bin": np.arange(num_bins),
            "bin_edge": np.arange(num_bins + 1),
        },
    )

    pdf_ds["FREQPDF"].attrs["long_name"] = "Frequency PDF"
    pdf_ds["AMNTPDF"].attrs["long_name"] = "Amount PDF"
    pdf_ds["bin_centers"].attrs["long_name"] = "Precipitation bin centers"
    pdf_ds["bin_centers"].attrs["units"] = "mm/day"

    # Add coordinate bounds for area-weighted averaging
    pdf_ds = pdf_ds.bounds.add_missing_bounds(axes=["X", "Y"])

    return pdf_ds, str(actual_start), str(actual_end)


def calculate_gridded_pdf(
    precip: xr.DataArray, bin_edges: np.ndarray, bin_centers: np.ndarray
) -> tuple:
    """Calculate frequency and amount PDFs for each grid point.

    Following Chris Terai's WC_diag_amwg.py formulation:
    - Frequency PDF: df/dlog(P) = precip_fraction / binwidth
    - Amount PDF: dP/dlog(P) = precip_fraction / binwidth * binmean

    Parameters
    ----------
    precip : xr.DataArray
        Daily precipitation data (time, lat, lon)
    bin_edges : np.ndarray
        Edges of precipitation bins
    bin_centers : np.ndarray
        Centers of precipitation bins (bin means)

    Returns
    -------
    tuple
        (freq_pdf, amnt_pdf) arrays of shape (lat, lon, nbins)
    """
    nlat, nlon = len(precip.lat), len(precip.lon)
    nbins = len(bin_edges) - 1

    freq_pdf = np.zeros((nlat, nlon, nbins))
    amnt_pdf = np.zeros((nlat, nlon, nbins))

    precip_values = precip.values  # (time, lat, lon)

    # Calculate bin widths in log space
    bin_widths_log = np.diff(np.log10(bin_edges))

    for ilat in range(nlat):
        for ilon in range(nlon):
            precip_timeseries = precip_values[:, ilat, ilon]

            # Count total valid data points (following Terai's approach)
            # Only count data >= 0.0 (filters out missing/negative values)
            valid_mask = precip_timeseries >= 0.0
            total_count = np.sum(valid_mask)

            if total_count == 0:
                continue

            # Filter array for performance (operate on smaller array in loop)
            precip_timeseries = precip_timeseries[valid_mask]

            # Calculate PDFs for each bin
            for ibin in range(nbins):
                # Find data in this bin (following Terai's approach)
                # For the last bin, include all data >= last edge (open-ended)
                if ibin < nbins - 1:
                    mask = (precip_timeseries >= bin_edges[ibin]) & (
                        precip_timeseries < bin_edges[ibin + 1]
                    )
                else:
                    # Last bin is open-ended (all data >= last edge)
                    mask = precip_timeseries >= bin_edges[ibin]

                count_in_bin = np.sum(mask)

                if count_in_bin > 0:
                    # Precipitation fraction in this bin
                    precip_fraction = count_in_bin / total_count

                    # Frequency PDF: df/dlog(P)
                    freq_pdf[ilat, ilon, ibin] = precip_fraction / bin_widths_log[ibin]

                    # Amount PDF: dP/dlog(P) = freq_pdf * bin_center
                    # This weights by the precipitation rate
                    amnt_pdf[ilat, ilon, ibin] = (
                        precip_fraction / bin_widths_log[ibin] * bin_centers[ibin]
                    )

    return freq_pdf, amnt_pdf


def extract_regional_pdf(pdf_ds: xr.Dataset, region: str) -> xr.Dataset:
    """Extract and average PDF over a specified region.

    Parameters
    ----------
    pdf_ds : xr.Dataset
        Gridded PDF dataset
    region : str
        Region name (e.g., '15S15N', 'TROPICS', 'CONUS', etc.)
        Must be a valid region from REGION_SPECS

    Returns
    -------
    xr.Dataset
        Regionally averaged PDF with area weighting
    """
    # Use built-in region specifications
    if region not in REGION_SPECS:
        logger.warning(
            f"Unknown region: {region}. Available regions: {list(REGION_SPECS.keys())}. "
            f"Using global domain."
        )
        regional_pdf = pdf_ds
    else:
        specs = REGION_SPECS[region]
        lat_range = specs.get("lat")
        lon_range = specs.get("lon")

        regional_pdf = pdf_ds.copy()

        # Apply latitude subsetting
        if lat_range is not None:
            lat_min, lat_max = lat_range
            regional_pdf = regional_pdf.sel(lat=slice(lat_min, lat_max))

        # Apply longitude subsetting
        if lon_range is not None:
            lon_min, lon_max = lon_range
            regional_pdf = regional_pdf.sel(lon=slice(lon_min, lon_max))

    # Calculate area-weighted mean
    # Get spatial weights (area weights based on grid cell sizes)
    weights = regional_pdf.spatial.get_weights(axis=["X", "Y"])

    # Apply area-weighted averaging for each variable
    freq_pdf_mean = regional_pdf["FREQPDF"].weighted(weights).mean(dim=["lat", "lon"])
    amnt_pdf_mean = regional_pdf["AMNTPDF"].weighted(weights).mean(dim=["lat", "lon"])

    result = xr.Dataset(
        {
            "FREQPDF": freq_pdf_mean,
            "AMNTPDF": amnt_pdf_mean,
            "bin_centers": regional_pdf["bin_centers"],
        }
    )

    return result

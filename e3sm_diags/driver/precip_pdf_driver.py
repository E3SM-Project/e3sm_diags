"""
Driver for precipitation PDF diagnostics.

Computes and plots precipitation probability density functions (PDFs)
for various regions (tropics, CONUS, etc.) from daily precipitation data.

Based on original work by Chris Terai (2015, 2020).
Modified to integrate into E3SM Diags.
"""
from __future__ import annotations

import glob
from typing import TYPE_CHECKING

import numpy as np
import xarray as xr

from e3sm_diags.derivations.default_regions_xr import REGION_SPECS
from e3sm_diags.driver.utils.climo_xr import ClimoFreq
from e3sm_diags.driver.utils.dataset_xr import Dataset
from e3sm_diags.logger import custom_logger
from e3sm_diags.plot.precip_pdf_plot import plot

if TYPE_CHECKING:
    from e3sm_diags.parameter.precip_pdf_parameter import PrecipPDFParameter


logger = custom_logger(__name__)


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
        # Read and process test data
        test_pdf = calculate_precip_pdf(
            parameter.test_data_path,
            variable,
            parameter.test_start_yr,
            parameter.test_end_yr,
        )
        parameter.test_name_yrs = test_data.get_name_yrs_attr(season)

        # Read and process reference data
        if run_type == "model_vs_model":
            ref_pdf = calculate_precip_pdf(
                parameter.reference_data_path,
                variable,
                parameter.ref_start_yr,
                parameter.ref_end_yr,
            )
        elif run_type == "model_vs_obs":
            ref_data_path = f"{parameter.reference_data_path}/{parameter.ref_name}"
            ref_pdf = calculate_precip_pdf(
                ref_data_path,
                variable,
                parameter.ref_start_yr,
                parameter.ref_end_yr,
            )
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


def calculate_precip_pdf(
    path: str, variable: str, start_year: str, end_year: str
) -> xr.Dataset:
    """Calculate precipitation PDFs from daily data.

    Parameters
    ----------
    path : str
        Path to daily precipitation data
    variable : str
        Variable name (typically "PRECT")
    start_year : str
        Start year
    end_year : str
        End year

    Returns
    -------
    xr.Dataset
        Dataset containing frequency and amount PDFs
    """
    # Read daily precipitation data (following tropical_subseasonal pattern)
    try:
        var = xr.open_mfdataset(glob.glob(f"{path}/{variable}_*.nc")).sel(
            time=slice(f"{start_year}-01-01", f"{end_year}-12-31")
        )[variable]
        actual_start = var.time.dt.year.values[0]
        actual_end = var.time.dt.year.values[-1]

        logger.info(
            f"Loaded {variable} data from {actual_start} to {actual_end} from {path}"
        )
    except OSError:
        logger.error(
            f"No files to open for {variable} within {start_year} and {end_year} from {path}."
        )
        raise

    # Unit conversion to mm/day (following tropical_subseasonal pattern)
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

    # Create precipitation bins (following GPCP convention)
    # Logarithmically spaced bins from 0.1 to 600 mm/day
    num_bins = 129
    bin_edges = np.logspace(np.log10(0.1), np.log10(600.0), num_bins + 1)
    bin_centers = np.sqrt(bin_edges[:-1] * bin_edges[1:])  # Geometric mean

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

    return pdf_ds


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

            # Remove missing values
            precip_timeseries = precip_timeseries[~np.isnan(precip_timeseries)]

            if len(precip_timeseries) == 0:
                continue

            # Count total valid data points
            total_count = len(precip_timeseries)

            # Calculate PDFs for each bin
            for ibin in range(nbins):
                # Find data in this bin
                mask = (precip_timeseries >= bin_edges[ibin]) & (
                    precip_timeseries < bin_edges[ibin + 1]
                )
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
        Regionally averaged PDF
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
    # Simple average for now (TODO: add proper area weighting)
    freq_pdf_mean = regional_pdf["FREQPDF"].mean(dim=["lat", "lon"])
    amnt_pdf_mean = regional_pdf["AMNTPDF"].mean(dim=["lat", "lon"])

    result = xr.Dataset(
        {
            "FREQPDF": freq_pdf_mean,
            "AMNTPDF": amnt_pdf_mean,
            "bin_centers": regional_pdf["bin_centers"],
        }
    )

    return result

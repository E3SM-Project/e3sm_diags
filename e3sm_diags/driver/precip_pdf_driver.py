"""
Driver for precipitation PDF diagnostics.

Computes and plots precipitation probability density functions (PDFs)
for various regions (tropics, CONUS, etc.) from daily precipitation data.

Based on original work by Chris Terai (2015, 2020).
Modified to integrate into E3SM Diags.
"""

from __future__ import annotations

import os
from typing import TYPE_CHECKING

import numpy as np
import xarray as xr
import xcdat as xc

from e3sm_diags.derivations.default_regions_xr import REGION_SPECS
from e3sm_diags.driver.utils.dataset_xr import Dataset
from e3sm_diags.driver.utils.io import _get_output_dir
from e3sm_diags.logger import _setup_child_logger
from e3sm_diags.plot.precip_pdf_plot import plot

if TYPE_CHECKING:
    from e3sm_diags.parameter.precip_pdf_parameter import PrecipPDFParameter


logger = _setup_child_logger(__name__)


def subset_by_season(ds: xr.Dataset, variable: str, season: str) -> xr.Dataset:
    """Subset dataset by season using xarray's optimized groupby.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset containing the variable of interest
    variable : str
        Variable name (not used but kept for API consistency)
    season : str
        Season to subset: 'DJF', 'MAM', 'JJA', 'SON', or 'ANN' for all months

    Returns
    -------
    xr.Dataset
        Dataset subset to the specified season

    Notes
    -----
    Uses xarray's built-in groupby('time.season') which is optimized and faster
    than manual month filtering. Season definitions:
    - DJF: December, January, February
    - MAM: March, April, May
    - JJA: June, July, August
    - SON: September, October, November
    """
    # If annual, return full dataset
    if season == "ANN":
        return ds

    # Validate season
    valid_seasons = ["DJF", "MAM", "JJA", "SON", "ANN"]
    if season not in valid_seasons:
        raise ValueError(f"Invalid season: {season}. Must be one of {valid_seasons}")

    # Use xarray's optimized groupby to subset by season
    # Create a dictionary of season groups for fast access
    season_groups = {s: data for s, data in ds.groupby("time.season")}

    # Get the requested season's data
    ds_subset = season_groups[season]

    logger.info(f"Subset data to {season}: {len(ds_subset.time)} time steps")

    return ds_subset


def load_cached_pdf(
    parameter: PrecipPDFParameter,
    variable: str,
    data_type: str,
    season: str = "ANN",
) -> tuple[xr.Dataset | None, str | None, str | None]:
    """Try to load a cached global PDF from a previous run.

    Parameters
    ----------
    parameter : PrecipPDFParameter
        Parameter object containing path information
    variable : str
        Variable name (e.g., "PRECT")
    data_type : str
        Type of data: "test" or "ref"
    season : str
        Season identifier (e.g., "ANN", "DJF", "MAM", "JJA", "SON")

    Returns
    -------
    tuple[xr.Dataset | None, str | None, str | None]
        If cache exists: (pdf_dataset, start_year, end_year)
        If cache doesn't exist: (None, None, None)
    """
    # Get the dataset name and time range for constructing cache filename
    if data_type == "test":
        dataset_name = parameter.test_name
        start_yr = parameter.test_start_yr
        end_yr = parameter.test_end_yr
    else:
        dataset_name = parameter.ref_name  # type: ignore[assignment]
        start_yr = parameter.ref_start_yr
        end_yr = parameter.ref_end_yr

    # Construct the expected cache filename with season
    output_dir = _get_output_dir(parameter)
    cache_filename = f"{variable}_PDF_global_{data_type}_{dataset_name}_{start_yr}-{end_yr}_{season}.nc"
    cache_filepath = os.path.join(output_dir, cache_filename)

    # Try to load the cached file
    try:
        if os.path.exists(cache_filepath):
            logger.info(f"Found cached PDF: {cache_filepath}")
            pdf_ds = xr.open_dataset(cache_filepath)

            # Extract time range from attributes
            cached_start = pdf_ds.attrs.get("start_year", start_yr)
            cached_end = pdf_ds.attrs.get("end_year", end_yr)

            return pdf_ds, str(cached_start), str(cached_end)
        else:
            return None, None, None
    except Exception as e:
        logger.warning(f"Failed to load cached PDF from {cache_filepath}: {e}")
        return None, None, None


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

    # Normalize ref_name to a list
    if isinstance(parameter.ref_name, str):
        if parameter.ref_name == "":
            # Default to both GPCP and IMERG when not specified
            ref_names = ["GPCP", "IMERG"]
            logger.info("No ref_name specified, using default: GPCP and IMERG")
        else:
            ref_names = [parameter.ref_name]
    else:
        ref_names = parameter.ref_name

    logger.info(f"Processing reference datasets: {ref_names}")

    test_data = Dataset(parameter, data_type="test")

    # Determine which seasons to process
    if parameter.season_subset:
        seasons_to_process = ["ANN", "DJF", "MAM", "JJA", "SON"]
        logger.info("Processing all months (ANN) and all seasons (DJF, MAM, JJA, SON)")
    else:
        seasons_to_process = ["ANN"]
        logger.info("Processing all months (ANN) only")

    for variable in parameter.variables:
        for season in seasons_to_process:
            logger.info(f"\n{'=' * 60}")
            logger.info(f"Processing {variable} for season: {season}")
            logger.info(f"{'=' * 60}\n")

            # Try to load cached test PDF first (performance optimization)
            test_pdf, test_start, test_end = load_cached_pdf(
                parameter, variable, "test", season
            )

            if test_pdf is None:
                # Cache miss - calculate test PDF from raw data
                logger.info(
                    f"Calculating test PDF for {variable} {season} (no cache found)"
                )
                test_ds = test_data.get_time_series_dataset(variable, single_point=True)

                # Subset by season
                test_ds_season = subset_by_season(test_ds, variable, season)

                test_pdf, test_start, test_end = calculate_precip_pdf(
                    test_ds_season, variable
                )

                # Save for future use
                save_global_pdf_to_netcdf(
                    parameter, test_pdf, variable, "test", test_start, test_end, season
                )
            else:
                # Cache hit - using cached PDF
                logger.info(
                    f"Using cached test PDF for {variable} {season} ({test_start}-{test_end})"
                )

            # Update parameters with actual time range
            parameter.test_start_yr = test_start  # type: ignore[assignment]
            parameter.test_end_yr = test_end  # type: ignore[assignment]
            parameter.test_name_yrs = test_data.get_name_yrs_attr(
                "ANN"
            )  # Use ANN for name
            parameter.var_id = variable

            # Process all reference datasets
            ref_pdfs = []
            ref_info = []  # Store (name, start_yr, end_yr, name_yrs) for each ref dataset

            for ref_name in ref_names:
                # Temporarily set ref_name for Dataset to work correctly
                original_ref_name = parameter.ref_name
                parameter.ref_name = ref_name

                # Create ref_data for this specific reference dataset
                ref_data = Dataset(parameter, data_type="ref")

                # Try to load cached reference PDF first (performance optimization)
                ref_pdf, ref_start, ref_end = load_cached_pdf(
                    parameter, variable, "ref", season
                )

                if ref_pdf is None:
                    # Cache miss - calculate reference PDF from raw data
                    logger.info(
                        f"Calculating reference PDF for {variable} {season} from {ref_name} (no cache found)"
                    )
                    try:
                        if run_type == "model_vs_model":
                            ref_ds = ref_data.get_time_series_dataset(
                                variable, single_point=True
                            )
                            ref_ds_season = subset_by_season(ref_ds, variable, season)
                            ref_pdf, ref_start, ref_end = calculate_precip_pdf(
                                ref_ds_season, variable
                            )
                        elif run_type == "model_vs_obs":
                            ref_ds = ref_data.get_time_series_dataset(
                                variable, single_point=True
                            )
                            ref_ds_season = subset_by_season(ref_ds, variable, season)
                            ref_pdf, ref_start, ref_end = calculate_precip_pdf(
                                ref_ds_season, variable
                            )

                        # Save for future use
                        save_global_pdf_to_netcdf(
                            parameter,
                            ref_pdf,  # type: ignore[arg-type]
                            variable,
                            "ref",
                            ref_start,  # type: ignore[arg-type]
                            ref_end,  # type: ignore[arg-type]
                            season,
                        )
                    except Exception as e:
                        logger.warning(
                            f"Failed to load/calculate reference PDF for {ref_name}: {e}. Skipping."
                        )
                        parameter.ref_name = original_ref_name
                        continue
                else:
                    # Cache hit - using cached PDF
                    logger.info(
                        f"Using cached reference PDF for {variable} {season} from {ref_name} "
                        f"({ref_start}-{ref_end})"
                    )

                # Store the PDF and metadata
                ref_pdfs.append(ref_pdf)
                ref_name_yrs = ref_data.get_name_yrs_attr("ANN")
                ref_info.append((ref_name, ref_start, ref_end, ref_name_yrs))

                # Restore original ref_name
                parameter.ref_name = original_ref_name

            # Skip plotting if no reference data was successfully loaded
            if not ref_pdfs:
                logger.warning(
                    f"No reference data available for {variable} {season}. Skipping plotting."
                )
                continue

            # Log what reference datasets were successfully loaded
            logger.info(
                f"Successfully loaded {len(ref_pdfs)} reference dataset(s): "
                f"{[info[0] for info in ref_info]}"
            )

            # Calculate regional PDFs and plot
            for region in parameter.regions:
                # Create output filename with multiple ref names
                ref_names_str = "_".join([info[0] for info in ref_info])
                parameter.output_file = (
                    f"{parameter.var_id}_PDF_{region}_{ref_names_str}_{season}"
                )

                # Extract regional PDFs
                test_regional = extract_regional_pdf(test_pdf, region)
                # Note: ref_pdfs list only contains successfully loaded PDFs (no None values)
                ref_regionals = [
                    extract_regional_pdf(ref_pdf, region)  # type: ignore[arg-type]
                    for ref_pdf in ref_pdfs
                ]

                # Call plot function with season information and all reference datasets
                plot(parameter, test_regional, ref_regionals, ref_info, region, season)

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

    logger.info(f"Loaded {variable} data from {actual_start} to {actual_end}")

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
        lat_range = specs.get("lat")  # type: ignore[attr-defined]
        lon_range = specs.get("lon")  # type: ignore[attr-defined]

        regional_pdf = pdf_ds.copy()

        # Apply longitude subsetting (following regrid.py pattern)
        if lon_range is not None:
            regional_pdf = regional_pdf.sortby("lon")

            # Swap longitude axis if needed (0-360 to -180-180)
            # Check if region uses negative longitudes but data uses 0-360
            is_lon_axis_diff = lon_range[0] < 0 and regional_pdf["lon"].values[0] >= 0
            if is_lon_axis_diff:
                logger.info(
                    f"Converting longitude from 0-360 to -180 to 180 for region {region}"
                )
                regional_pdf = xc.swap_lon_axis(regional_pdf, to=(-180, 180))

            regional_pdf = regional_pdf.sel(lon=slice(*lon_range))

        # Apply latitude subsetting
        if lat_range is not None:
            regional_pdf = regional_pdf.sortby("lat")
            regional_pdf = regional_pdf.sel(lat=slice(*lat_range))

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


def save_global_pdf_to_netcdf(
    parameter: PrecipPDFParameter,
    pdf_ds: xr.Dataset,
    variable: str,
    data_type: str,
    start_yr: str,
    end_yr: str,
    season: str = "ANN",
) -> None:
    """Save global PDF dataset to a netCDF file for offline use.

    Parameters
    ----------
    parameter : PrecipPDFParameter
        Parameter object containing output configuration
    pdf_ds : xr.Dataset
        Global gridded PDF dataset containing FREQPDF and AMNTPDF
    variable : str
        Variable name (e.g., "PRECT")
    data_type : str
        Type of data: "test" or "ref"
    start_yr : str
        Start year of the data
    end_yr : str
        End year of the data
    season : str
        Season identifier (e.g., "ANN", "DJF", "MAM", "JJA", "SON")
    """
    if not parameter.save_netcdf:
        return

    # Get output directory
    output_dir = _get_output_dir(parameter)

    # Construct filename with dataset name and season to avoid overwrites
    # e.g., PRECT_PDF_global_test_E3SMv2_1996-2004_DJF.nc
    #       PRECT_PDF_global_ref_GPCP_1DD_Daily_1996-2010_ANN.nc
    if data_type == "test":
        dataset_name = parameter.test_name
    else:
        dataset_name = parameter.ref_name  # type: ignore[assignment]

    filename = f"{variable}_PDF_global_{data_type}_{dataset_name}_{start_yr}-{end_yr}_{season}.nc"
    filepath = os.path.join(output_dir, filename)

    # Add metadata to the dataset
    pdf_ds.attrs["variable"] = variable
    pdf_ds.attrs["data_type"] = data_type
    pdf_ds.attrs["dataset_name"] = dataset_name
    pdf_ds.attrs["start_year"] = start_yr
    pdf_ds.attrs["end_year"] = end_yr
    pdf_ds.attrs["season"] = season
    pdf_ds.attrs["description"] = (
        f"Global gridded precipitation PDFs for {variable} from {dataset_name} "
        f"({start_yr}-{end_yr}) - Season: {season}"
    )

    # Save to netCDF
    pdf_ds.to_netcdf(filepath)

    logger.info(f"Global PDF {data_type} dataset saved to: {filepath}")

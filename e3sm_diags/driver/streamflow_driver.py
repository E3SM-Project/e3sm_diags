from __future__ import annotations

import csv
from typing import TYPE_CHECKING, List, Tuple

import numpy as np
import scipy.io
import xarray as xr
import xcdat as xc

from e3sm_diags.driver.utils.dataset_xr import Dataset
from e3sm_diags.logger import _setup_child_logger
from e3sm_diags.plot.streamflow_plot_map import plot_annual_map
from e3sm_diags.plot.streamflow_plot_scatter import plot_annual_scatter
from e3sm_diags.plot.streamflow_plot_seasonality import plot_seasonality_map

logger = _setup_child_logger(__name__)

if TYPE_CHECKING:
    from e3sm_diags.parameter.streamflow_parameter import StreamflowParameter


# Target physical search distance (0.5 degrees for LR)
TARGET_SEARCH_DISTANCE = 0.5
# The max area error (percent) for all plots.
MAX_AREA_ERROR = 20


def _calculate_resolution(ds: xr.Dataset) -> float:
    """Calculate the grid resolution from a dataset.

    Determines the resolution by calculating the spacing between
    the first two latitude points.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset containing lat/lon coordinates

    Returns
    -------
    float
        The detected grid resolution in degrees
    """
    try:
        # Get the latitude coordinates
        lat_dim = xc.get_dim_coords(ds, axis="Y")
        lat_vals = lat_dim.values

        # Calculate spacing between first two points
        if len(lat_vals) >= 2:
            resolution = abs(lat_vals[1] - lat_vals[0])
            logger.info(f"Detected grid resolution: {resolution} degrees")
            return resolution
        else:
            logger.warning("Not enough lat values to determine resolution")
            return 0.5
    except Exception as e:
        logger.warning(f"Could not detect grid resolution from dataset: {e}")
        logger.warning("Using default resolution of 0.5 degrees")
        return 0.5


def run_diag(parameter: StreamflowParameter) -> StreamflowParameter:
    """Get metrics for the streamflow set.

    Parameters
    ----------
    parameter : StreamflowParameter
        The parameter for the diagnostic.

    Returns
    -------
    StreamflowParameter
        The parameter for the diagnostic with the result (completed or failed).
    """
    gauges, is_ref_mat_file = _get_gauges(parameter)

    for var_key in parameter.variables:
        logger.info(f"Variable: {var_key}")

        # Get test data and resolution information
        test_array, area_upstream, resolution, search_radius = (
            _get_test_data_and_area_upstream(parameter, var_key)
        )
        ref_array = _get_ref_data(parameter, var_key, is_ref_mat_file)

        export_data = _generate_export_data(
            parameter,
            gauges,
            test_array,
            ref_array,
            area_upstream,
            is_ref_mat_file,
            resolution=resolution,
            search_radius=search_radius,
        )

        if parameter.test_title == "":
            parameter.test_title = parameter.test_name_yrs
        if parameter.reference_title == "":
            parameter.reference_title = parameter.ref_name_yrs

        # Plot the original ref and test (not regridded versions).
        plot_seasonality_map(parameter, export_data)
        plot_annual_map(parameter, export_data)
        plot_annual_scatter(parameter, export_data)

    return parameter


def _get_gauges(parameter: StreamflowParameter) -> Tuple[np.ndarray, bool]:
    """Get the gauges.

    Assume `model_vs_model` is an `nc` file and `model_vs_obs` is an `mat` file.
    If `model_vs_obs`, the metadata file of GSIM that has observed gauge lat lon
    and drainage area. This file includes 25765 gauges, which is a subset of the
    entire dataset (30959 gauges). The removed gauges are associated with very
    small drainage area (<1km2), which is not meaningful to be included.

    Parameters
    ----------
    parameter : StreamflowParameter
        The parameter.

    Returns
    -------
    Tuple[np.ndarray, bool]
        A tuple containing the gauges array and a boolean representing whether
        the reference file is a mat file (True) or not (False).

    Raises
    ------
    RuntimeError
        Non-GSIM reference file specified without using parameter `.gauges_path`
        attribute.
    RuntimeError
        Parameter run type is not supported.
    """
    ref_path = parameter.reference_data_path.rstrip("/")

    if parameter.run_type == "model_vs_model":
        is_ref_mat_file = False

        if parameter.gauges_path is None:
            raise RuntimeError(
                "To use a non-GSIM reference, please specify streamflow_param.gauges_path. "
                f"This might be {ref_path}/GSIM/GSIM_catchment_characteristics_all_1km2.csv"
            )

        else:
            gauges_path = parameter.gauges_path
    elif parameter.run_type == "model_vs_obs":
        is_ref_mat_file = True
        gauges_path = f"{ref_path}/GSIM/GSIM_catchment_characteristics_all_1km2.csv"
    else:
        raise RuntimeError(f"parameter.run_type={parameter.run_type} not supported")

    # Set path to the gauge metadata
    with open(gauges_path) as gauges_file:
        gauges_list = list(csv.reader(gauges_file))

    # Remove headers
    gauges_list.pop(0)
    gauges = np.array(gauges_list)

    return gauges, is_ref_mat_file


def _get_test_data_and_area_upstream(
    parameter: StreamflowParameter, var_key: str
) -> Tuple[np.ndarray, np.ndarray, float, int]:
    """Set up the test data and calculate grid resolution.

    Parameters
    ----------
    parameter : StreamflowParameter
        The parameter.
    var_key : str
        The key of the variable.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray, float, int]
        The test data, area upstream, detected resolution, and calculated search radius.
    """
    test_ds = Dataset(parameter, data_type="test")
    parameter.test_name_yrs = test_ds.get_name_yrs_attr()

    ds_test = test_ds.get_time_series_dataset(var_key)

    # Detect grid resolution from the dataset
    resolution = _calculate_resolution(ds_test)

    # Calculate search radius based on detected resolution
    search_radius = round(TARGET_SEARCH_DISTANCE / resolution)
    logger.info(f"Using search radius of {search_radius} grid cells")

    test_array = _get_var_data(ds_test, var_key)

    areatotal2 = ds_test["areatotal2"].values
    area_upstream = np.transpose(areatotal2, (1, 0)).astype(np.float64)

    return test_array, area_upstream, resolution, search_radius


def _get_ref_data(
    parameter: StreamflowParameter, var_key: str, is_ref_mat_file: bool
) -> np.ndarray:
    """Set up the reference data.

    Parameters
    ----------
    parameter : StreamflowParameter
        The parameter.
    var_key : str
        The key of the variable.
    is_ref_mat_file : bool
        If the reference data is from a mat file (True) or not (False).

    Returns
    -------
    np.ndarray
        The reference data.
    """
    ref_ds = Dataset(parameter, data_type="ref")

    if not is_ref_mat_file:
        parameter.ref_name_yrs = ref_ds.get_name_yrs_attr()

        ds_ref = ref_ds.get_time_series_dataset(var_key)
        ref_array = _get_var_data(ds_ref, var_key)
    else:
        # Load the observed streamflow dataset (GSIM)
        # the data has been reorganized to a 1380 * 30961 matrix. 1380 is the month
        # number from 1901.1 to 2015.12. 30961 include two columns for year and month plus
        # streamflow at 30959 gauge locations reported by GSIM
        ref_path = parameter.reference_data_path.rstrip("/")
        ref_mat_file = f"{ref_path}/GSIM/GSIM_198601_199512.mat"
        parameter.ref_name_yrs = ref_ds.get_name_yrs_attr(default_name="GSIM")

        ref_mat = scipy.io.loadmat(ref_mat_file)
        ref_array = ref_mat["GSIM"].astype(np.float64)

    return ref_array


def _get_var_data(ds: xr.Dataset, var_key: str) -> np.ndarray:
    """Get the variable data then subset on latitude and transpose.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset object.
    var_key : str
        The key of the variable.

    Returns
    -------
    np.ndarray
        The variable data.
    """
    da_var = ds[var_key].copy()
    lat_dim = xc.get_dim_keys(da_var, axis="Y")

    da_var_reg = da_var.sel({lat_dim: slice(-90, 90)})
    var_transposed = np.transpose(da_var_reg.values, (2, 1, 0))
    test_array = var_transposed.astype(np.float64)

    return test_array


def _generate_export_data(
    parameter: StreamflowParameter,
    gauges: np.ndarray,
    test_array: np.ndarray,
    ref_array: np.ndarray,
    area_upstream: np.ndarray,
    is_ref_mat_file: bool,
    resolution: float = 0.5,
    search_radius: int = 1,
) -> np.ndarray:
    """Generate the export data.

    Parameters
    ----------
    parameter : StreamflowParameter
        The parameter.
    gauges : np.ndarray
        The gauges.
    test_array : np.ndarray
        The test data.
    ref_array : np.ndarray
        The reference data.
    area_upstream : np.ndarray
        The area upstream.
    is_ref_mat_file : bool
        If the reference data is a mat file or not.

    Returns
    -------
    np.ndarray
        The export data as a 2D array with the format:
          - col 0: year
          - col 1: month
          - [i, 0]: annual mean for reference data
          - [i, 1]: annual mean for test data
          - [i, 2]: percentage of drainage area bias
          - [i, 3]: seasonality index of reference data
          - [i, 4]: peak month flow of reference data
          - [i, 5]: seasonality index of test data
          - [i, 6]: peak month flow of test data
          - [i, 7:9]: Lat lon index of reference data.
    Notes
    -----
    TODO: This function should be refactored to make it readable and
    maintainable. The number of code comments suggest that the code is not
    understandable and needs to be explained line by line.
    """
    # Center the reference lat lon to the grid center and use it
    # to create the shape for the export data array.
    bins = np.floor(gauges[:, 7:9].astype(np.float64) / resolution)
    lat_lon = (bins + 0.5) * resolution

    export_data = np.zeros((lat_lon.shape[0], 9))

    for i in range(lat_lon.shape[0]):
        if parameter.max_num_gauges and i >= parameter.max_num_gauges:
            break

        lat_ref = lat_lon[i, 1]
        lon_ref = lat_lon[i, 0]

        # Estimated drainage area (km^2) from ref
        area_ref = gauges[i, 13].astype(np.float64)
        drainage_area_error, lat_lon_ref = _get_drainage_area_error(
            lon_ref,
            lat_ref,
            area_upstream,
            area_ref,
            resolution=resolution,
            search_radius=search_radius,
        )

        if is_ref_mat_file:
            origin_id = gauges[i, 1].astype(np.int64)
            # Column 0 -- year
            # Column 1 -- month
            # Column origin_id + 1 -- the ref streamflow from gauge with the corresponding origin_id
            extracted = ref_array[:, [0, 1, origin_id + 1]]
            monthly_mean = np.zeros((12, 1)) + np.nan
            # For GSIM, shape is (1380,)
            month_array = extracted[:, 1]
            for month in range(12):
                # Add 1 to month to account for the months being 1-indexed
                month_array_boolean = month_array == month + 1
                s = np.sum(month_array_boolean)
                if s > 0:
                    # `extracted[:,1]`: for all x, examine `extracted[x,1]`
                    # `extracted[:,1] == m`: Boolean array where 0 means the item in position [x,1] is NOT m,
                    # and 1 means it is m
                    # Example:
                    # a = [[1,2,3,4],
                    #      [5,6,7,8],
                    #      [1,2,3,4]]
                    # a[:,1]: [[2],
                    #          [6],
                    #          [2]]
                    # a[:,1] == 2: [[1], # False is 0, True is 1
                    #               [0],
                    #               [1]]
                    # a[a[:,1] == 2, 2]: [[3],
                    #                     [3]]
                    monthly_mean[month] = np.nanmean(extracted[month_array_boolean, 2])
            # This is ref annual mean streamflow
            annual_mean_ref = np.mean(monthly_mean)
        if is_ref_mat_file and np.isnan(annual_mean_ref):
            # All elements of row i will be nan
            export_data[i, :] = np.nan
        else:
            if is_ref_mat_file:
                # Reshape extracted[:,2] into a 12 x ? matrix; -1 means to
                # calculate the size of the missing dimension.
                # Note that `np.reshape(extracted[:, 2], (12,-1))` will not work.
                # We do need to go from (12n x 1) to (12 x n).
                # `reshape` alone would make the first row [January of year 1, February of year 1,...]
                # (i.e., 12 sequential rows with n entries)
                # We actually want the first row to be [January of year 1, January of year 2,...]
                # (i.e., n sequential columns with 12 entries)
                # So, we use `reshape` to slice into n segments of length 12 and then we `transpose`.
                mmat = np.transpose(np.reshape(extracted[:, 2], (-1, 12)))
                mmat_id = np.sum(mmat, axis=0).transpose()
                if np.sum(~np.isnan(mmat_id), axis=0) > 0:
                    # There's at least one year of record
                    monthly = mmat[:, ~np.isnan(mmat_id)]
                else:
                    monthly = monthly_mean
                seasonality_index_ref, peak_month_ref = _get_seasonality(monthly)
            else:
                ref_lon = int(
                    1 + (lat_lon_ref[1] - (-180 + resolution / 2)) / resolution
                )
                ref_lat = int(
                    1 + (lat_lon_ref[0] - (-90 + resolution / 2)) / resolution
                )
                ref = np.squeeze(ref_array[ref_lon - 1, ref_lat - 1, :])
                # Note that `np.reshape(ref, (12,-1))` will not work.
                # We do need to go from (12n x 1) to (12 x n).
                # `reshape` alone would make the first row [January of year 1, February of year 1,...]
                # (i.e., 12 sequential rows with n entries)
                # We actually want the first row to be [January of year 1, January of year 2,...]
                # (i.e., n sequential columns with 12 entries)
                # So, we use `reshape` to slice into n segments of length 12 and then we `transpose`.
                mmat = np.transpose(np.reshape(ref, (-1, 12)))
                monthly_mean_ref = np.nanmean(mmat, axis=1)
                annual_mean_ref = np.mean(monthly_mean_ref)
                if np.isnan(annual_mean_ref) == 1:
                    # The identified grid is in the ocean
                    monthly = np.ones((12, 1))
                else:
                    monthly = mmat

                seasonality_index_ref, peak_month_ref = _get_seasonality(monthly)

            test_lon = int(1 + (lat_lon_ref[1] - (-180 + resolution / 2)) / resolution)
            test_lat = int(1 + (lat_lon_ref[0] - (-90 + resolution / 2)) / resolution)
            # For edison: 600x1
            test = np.squeeze(test_array[test_lon - 1, test_lat - 1, :])
            # For edison: 12x50
            # Note that `np.reshape(test, (12,-1))` will not work.
            # We do need to go from (12n x 1) to (12 x n).
            # `reshape` alone would make the first row [January of year 1, February of year 1,...]
            # (i.e., 12 sequential rows with n entries)
            # We actually want the first row to be [January of year 1, January of year 2,...]
            # (i.e., n sequential columns with 12 entries)
            # So, we use `reshape` to slice into n segments of length 12 and then we `transpose`.
            mmat = np.transpose(np.reshape(test, (-1, 12)))
            monthly_mean_test = np.nanmean(mmat, axis=1)
            annual_mean_test = np.mean(monthly_mean_test)

            if np.isnan(annual_mean_test) == 1:
                # The identified grid is in the ocean
                monthly = np.ones((12, 1))
            else:
                monthly = mmat

            seasonality_index_test, peak_month_test = _get_seasonality(monthly)

            # TODO: The export data structure should be turned into a dict.
            export_data[i, 0] = annual_mean_ref
            export_data[i, 1] = annual_mean_test
            # Convert from fraction to percetange.
            export_data[i, 2] = drainage_area_error * 100
            export_data[i, 3] = seasonality_index_ref
            export_data[i, 4] = peak_month_ref
            export_data[i, 5] = seasonality_index_test
            export_data[i, 6] = peak_month_test
            export_data[i, 7:9] = lat_lon_ref

    export_data = _remove_gauges_with_nan_flow(export_data, area_upstream)

    return export_data


def _get_drainage_area_error(
    lon_ref: float,
    lat_ref: float,
    area_upstream: np.ndarray,
    area_ref: float,
    resolution: float = 0.5,
    search_radius: int = 1,
) -> Tuple[np.ndarray, List[float]]:
    """Get the drainage area error.

    Parameters
    ----------
    lon_ref : float
        The reference longitude.
    lat_ref : float
        The reference latitude.
    area_upstream : np.ndarray
        The area upstream.
    area_ref : float
        The reference area.
    resolution : float, optional
        The grid resolution in degrees, by default 0.5
    search_radius : int, optional
        The search radius in grid cells, by default 1

    Returns
    -------
    Tuple[np.ndarray, list[float, float]]
        The drainage area error and reference lat/lon.
    """
    k_bound = len(range(-search_radius, search_radius + 1))
    k_bound *= k_bound

    area_test = np.zeros((k_bound, 1))
    error_test = np.zeros((k_bound, 1))
    lat_lon_test = np.zeros((k_bound, 2))
    k = 0

    for i in range(-search_radius, search_radius + 1):
        for j in range(-search_radius, search_radius + 1):
            x = int(
                1 + ((lon_ref + j * resolution) - (-180 + resolution / 2)) / resolution
            )
            y = int(
                1 + ((lat_ref + i * resolution) - (-90 + resolution / 2)) / resolution
            )
            area_test[k] = area_upstream[x - 1, y - 1] / 1000000
            error_test[k] = np.abs(area_test[k] - area_ref) / area_ref
            lat_lon_test[k, 0] = lat_ref + i * resolution
            lat_lon_test[k, 1] = lon_ref + j * resolution
            k += 1

    # The id of the center grid in the searching area
    center_id = (k_bound - 1) / 2

    lat_lon_ref = [lat_ref, lon_ref]
    drainage_area_error = error_test[int(center_id)]

    return drainage_area_error, lat_lon_ref


def _get_seasonality(monthly: np.ndarray) -> Tuple[int, float]:
    """Get the seasonality.

    Parameters
    ----------
    monthly : np.ndarray
        The monthly data.

    Returns
    -------
    Tuple[int, float]
        A tuple including the seasonality index and the peak flow month.

    Raises
    ------
    RuntimeError
        If the monthly data does not include 12 months.
    """
    monthly = monthly.astype(np.float64)

    # See https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2018MS001603 Equations 1 and 2
    if monthly.shape[0] != 12:
        raise RuntimeError(f"monthly.shape={monthly.shape} does not include 12 months.")

    num_years = monthly.shape[1]
    p_k = np.zeros((12, 1))

    # The total streamflow for each year (sum of Q_ij in the denominator of Equation 1, for all j)
    # 1 x num_years
    total_streamflow = np.sum(monthly, axis=0)
    for month in range(12):
        # The streamflow for this month in each year (Q_ij in the numerator of Equation 1, for all j)
        # 1 x num_years
        streamflow_month_all_years = monthly[month, :]
        # Proportion that this month contributes to streamflow that year.
        # 1 x num_years
        # For all i, divide streamflow_month_all_years[i] by total_streamflow[i]
        streamflow_proportion = np.divide(streamflow_month_all_years, total_streamflow)
        # The sum is the sum over j in Equation 1.
        # Dividing the sum of proportions by num_years gives the *average* proportion of annual streamflow during
        # this month.
        # Multiplying by 12 makes it so that Pk_i (`p_k[month]`) will be 1 if all months have equal streamflow and
        # 12 if all streamflow occurs in one month.
        # These steps produce the 12/n factor in Equation 1.
        p_k[month] = np.nansum(streamflow_proportion) * 12 / num_years

    # From Equation 2
    seasonality_index = np.max(p_k)
    # `p_k == np.max(p_k)` produces a Boolean matrix, True if the value (i.e., streamflow) is the max value.
    # `np.where(p_k == np.max(p_k))` produces the indices (i.e., months) where the max value is reached.
    peak_month = np.where(p_k == np.max(p_k))[0]
    # If more than one month has peak streamflow, simply define the peak month as the first one of the peak months.
    # Month 0 is January, Month 1 is February, and so on.
    peak_month = peak_month[0]

    return seasonality_index, peak_month  # type: ignore


def _remove_gauges_with_nan_flow(
    export_data: np.ndarray, area_upstream: np.ndarray | None
) -> np.ndarray:
    """Remove gauges with NaN flow.

    Gauges will only plotted  if they have a non-nan value for both test and ref
    data.

    Parameters
    ----------
    export_data : np.ndarray
        The export data.
    area_upstream : np.ndarray | None
        The optional area upstream.

    Returns
    -------
    np.ndarray
        The export with gauges that have NaN flow removed.
    """
    export_data_new = np.array(export_data)
    export_data_new = export_data_new[~np.isnan(export_data_new[:, 0]), :]
    export_data_new = export_data_new[~np.isnan(export_data_new[:, 1]), :]

    if area_upstream is not None:
        export_data_new = export_data_new[export_data_new[:, 2] <= MAX_AREA_ERROR, :]

    return export_data_new

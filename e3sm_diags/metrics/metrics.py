"""This module stores functions to calculate metrics using Xarray objects."""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import numpy as np
import xarray as xr
import xcdat as xc
import xskillscore as xs

from e3sm_diags.logger import _setup_child_logger

if TYPE_CHECKING:
    import uxarray as ux

logger = _setup_child_logger(__name__)

Axis = Literal["X", "Y", "Z"]
DEFAULT_AXIS: list[Axis] = ["X", "Y"]


def spatial_avg(
    ds: xr.Dataset, var_key: str, axis: list[Axis] = DEFAULT_AXIS, as_list: bool = True
) -> list[float] | xr.DataArray:
    """Compute a variable's weighted spatial average.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset containing the variable.
    var_key : str
        The key of the variable in the dataset.
    axis : list[Axis]
        The list of axes to use for the computation, by default ["X", "Y"].
        Valid axes including "X", "Y", and "Z".
        The key of the variable.
    axis : list[str]
        A list of axis strings, by default ["X", "Y"].
    as_list : bool
        Return the spatial average as a list of floats, by default True.
        If False, return an xr.DataArray. Must be True to be serializable for
        writing out to a `.json` metrics file.

    Returns
    -------
    list[float] | xr.DataArray
        The spatial average of the variable based on the specified axis.

    Raises
    ------
    ValueError
        If the axis argument contains an invalid value.

    Notes
    -----
    Replaces `e3sm_diags.metrics.mean`.
    """
    dv = ds[var_key].copy()
    weights = _get_weights(ds, var_key, axis)
    dims = _get_dims(dv, axis)

    results = dv.weighted(weights).mean(dims, keep_attrs=True)

    if as_list:
        return results.data.tolist()

    return results


def std(ds: xr.Dataset, var_key: str, axis: list[Axis] = DEFAULT_AXIS) -> list[float]:
    """Compute the weighted standard deviation for a variable.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset containing the variable.
    var_key : str
        The key of the variable.
    axis : list[Axis]
        The list of axes to use for the computation, by default ["X", "Y"].
        Valid strings include "X", "Y", and "Z".

    Returns
    -------
    list[float]
        The standard deviation of the variable based on the specified axis.

    Raises
    ------
    ValueError
        If the axis argument contains an invalid value.

    Notes
    -----
    Replaces `e3sm_diags.metrics.std`.
    """
    dv = ds[var_key].copy()

    weights = _get_weights(ds, var_key, axis)
    dims = _get_dims(dv, axis)

    result = dv.weighted(weights).std(dim=dims, keep_attrs=True)

    return result.data.tolist()


def correlation(
    ds_a: xr.Dataset, ds_b: xr.Dataset, var_key: str, axis: list[Axis] = DEFAULT_AXIS
) -> list[float]:
    """Compute the correlation coefficient between two variables.

    This function uses the Pearson correlation coefficient. Refer to [1]_ for
    more information.

    Parameters
    ----------
    ds_a : xr.Dataset
        The first dataset.
    ds_b : xr.Dataset
        The second dataset.
    var_key: str
        The key of the variable.
    axis : list[Axis]
        The list of axes to use for the computation, by default ["X", "Y"].
        Valid axes including "X", "Y", and "Z".

    Returns
    -------
    list[float]
        The weighted correlation coefficient.

    References
    ----------

    .. [1] https://en.wikipedia.org/wiki/Pearson_correlation_coefficient

    Notes
    -----
    Replaces `e3sm_diags.metrics.corr`.
    """
    var_a = ds_a[var_key]
    var_b = ds_b[var_key]

    # Dimensions, bounds, and coordinates should be identical between datasets,
    # so use the first dataset and first variable to get dimensions and weights.
    dims = _get_dims(var_a, axis)
    weights = _get_weights(ds_a, var_key, axis)

    result = xs.pearson_r(var_a, var_b, dim=dims, weights=weights, skipna=True)
    results_list = result.data.tolist()

    return results_list


def rmse(
    ds_a: xr.Dataset, ds_b: xr.Dataset, var_key: str, axis: list[Axis] = DEFAULT_AXIS
) -> list[float]:
    """Calculates the root mean square error (RMSE) between two variables.

    Parameters
    ----------
    ds_a : xr.Dataset
        The first dataset.
    ds_b : xr.Dataset
        The second dataset.
    var_key: str
        The key of the variable.
    axis : list[Axis]
        The list of axes to use for the computation, by default ["X", "Y"].
        Valid axes including "X", "Y", and "Z".

    Returns
    -------
    list[float]
        The root mean square error.

    Notes
    -----
    Replaces `e3sm_diags.metrics.rmse`.
    """
    var_a = ds_a[var_key]
    var_b = ds_b[var_key]

    # Dimensions, bounds, and coordinates should be identical between datasets,
    # so use the first dataset and first variable to get dimensions and weights.
    dims = _get_dims(var_a, axis)
    weights = _get_weights(ds_a, var_key, axis)

    result = xs.rmse(var_a, var_b, dim=dims, weights=weights, skipna=True)
    results_list = result.data.tolist()

    return results_list


def _get_weights(ds: xr.Dataset, var_key: str, axis: list[Axis] = DEFAULT_AXIS):
    """Get weights for the specified axes of the variable.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset containing the variable.
    var_key : str
        The key of the variable.
    axis : list[Axis]
        The list of axes to use for the computation, by default ["X", "Y"].
        Valid axes including "X", "Y", and "Z".

    Returns
    -------
    xr.DataArray
        Weights for the specified axis.

    Notes
    -----
    xCDAT's ``ds.spatial.get_weights()`` method only supports rectilinear grids
    ("X", "Y") as of v0.6.1. This function computes weights for "X" and Y"
    using xCDAT, then calculates "Z" weights separately before combining all
    weights into a single matrix.
    """
    spatial_wts = None
    vertical_wts = None

    spatial_axis = [key for key in axis if key != "Z"]
    if len(spatial_axis) > 0:
        spatial_wts = ds.spatial.get_weights(spatial_axis, data_var=var_key)
        spatial_wts = spatial_wts.fillna(0)

    if "Z" in axis:
        vertical_wts = _get_z_weights(ds, var_key)

    if spatial_wts is not None and vertical_wts is not None:
        return spatial_wts * vertical_wts
    elif spatial_wts is not None:
        return spatial_wts
    elif vertical_wts is not None:
        return vertical_wts


def _get_z_weights(ds: xr.Dataset, var_key: str) -> xr.DataArray:
    """Get the Z axis weights using Z axis bounds.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset containing a Z axis.
    var_key : str
        The key of the variable to get weights for.

    Returns
    -------
    xr.DataArray
        Weights for the Z axis.

    Raises
    ------
    RuntimeError
        If the dataset has no Z axis bounds.

    Notes
    -----
    xCDAT spatial average get_weights() method does not support the Z axis as of
    v0.6.1. This function is temporarily used until get_weights() supports the
    Z axis. The logic is the same as xCDAT.

    * Related issue: https://github.com/xCDAT/xcdat/issues/596
    * Source: https://github.com/xCDAT/xcdat/blob/main/xcdat/spatial.py#L479-L495C16
    """
    try:
        z_bnds = ds.bounds.get_bounds("Z", var_key)
    except KeyError as err:
        raise RuntimeError(
            f"The dataset for {var_key} has no Z axis bounds to get weights for the Z "
            "axis."
        ) from err

    weights = np.abs(z_bnds[:, 1] - z_bnds[:, 0])
    weights = weights.fillna(0)

    return weights


def _get_dims(da: xr.DataArray, axis: list[Axis]) -> list[str]:
    """Get the dimensions for an axis in an xarray.DataArray.

    The dimensions are passed to the ``dim`` argument in xarray or xarray-based
    computational APIs, such as ``.std()``.

    Parameters
    ----------
    da : xr.DataArray
        The array.
    axis : list[Axis]
        The list of axes to get dimensions for. Valid strings
        include "X", "Y", and "Z".

    Returns
    -------
    list[str]
        A list of dimensions.
    """
    dims = []

    for a in axis:
        dim_key = xc.get_dim_keys(da, axis=a)
        dims.append(dim_key)

    return dims


def native_rmse(uxds_a: "ux.UxDataset", uxds_b: "ux.UxDataset", var_key: str) -> float:
    """Calculate RMSE for native grid datasets using uxarray and xskillscore.

    Parameters
    ----------
    uxds_a : ux.UxDataset
        The first uxarray dataset.
    uxds_b : ux.UxDataset
        The second uxarray dataset.
    var_key : str
        The key of the variable.

    Returns
    -------
    float
        The root mean square error.

    Raises
    ------
    RuntimeError
        If RMSE calculation fails.
    """
    try:
        import xskillscore as xs

        var_a = uxds_a[var_key]
        var_b = uxds_b[var_key]

        # Get spatial dimensions
        var_dims = list(var_a.dims)
        spatial_dims = [
            dim for dim in var_dims if "face" in dim or "node" in dim or "edge" in dim
        ]
        if not spatial_dims:
            spatial_dims = [dim for dim in var_dims if dim != "time"]

        # Get appropriate weights
        weights = None
        if var_a._face_centered():
            weights = var_a.uxgrid.face_areas
        elif var_a._edge_centered():
            weights = var_a.uxgrid.edge_node_distances

        return xs.rmse(
            var_a, var_b, dim=spatial_dims, weights=weights, skipna=True
        ).item()

    except Exception as e:
        raise RuntimeError(f"Failed to calculate native grid RMSE: {e}") from e


def native_correlation(
    uxds_a: "ux.UxDataset", uxds_b: "ux.UxDataset", var_key: str
) -> float:
    """Calculate Pearson correlation for native grid datasets using uxarray and xskillscore.

    Parameters
    ----------
    uxds_a : ux.UxDataset
        The first uxarray dataset.
    uxds_b : ux.UxDataset
        The second uxarray dataset.
    var_key : str
        The key of the variable.

    Returns
    -------
    float
        The Pearson correlation coefficient.

    Raises
    ------
    RuntimeError
        If correlation calculation fails.
    """
    try:
        import xskillscore as xs

        var_a = uxds_a[var_key]
        var_b = uxds_b[var_key]

        # Get spatial dimensions
        var_dims = list(var_a.dims)
        spatial_dims = [
            dim for dim in var_dims if "face" in dim or "node" in dim or "edge" in dim
        ]
        if not spatial_dims:
            spatial_dims = [dim for dim in var_dims if dim != "time"]

        # Get appropriate weights
        weights = None
        if var_a._face_centered():
            weights = var_a.uxgrid.face_areas
        elif var_a._edge_centered():
            weights = var_a.uxgrid.edge_node_distances

        return xs.pearson_r(
            var_a, var_b, dim=spatial_dims, weights=weights, skipna=True
        ).item()

    except Exception as e:
        raise RuntimeError(f"Failed to calculate native grid correlation: {e}") from e

import xarray as xr
import xcdat as xc

from e3sm_diags.logger import custom_logger

logger = custom_logger(__name__)


def correlation():
    # https://github.com/CDAT/genutil/blob/59517ab54e65c03098502f63434270f96020f3eb/Lib/statistics.py#L731-L793
    # https://github.com/CDAT/genutil/blob/59517ab54e65c03098502f63434270f96020f3eb/Lib/statistics.py#L268-L279
    pass


def spatial_avg(ds: xr.Dataset, var_key: str, axis=["X", "Y"]) -> xr.DataArray:
    """Compute a variable's weighted spatial average.

    This function is intended to replace ``e3sm_diags.metrics.mean()``.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset containing the variable.
    var_key : str
        The key of the varible.
    axis : list, optional
        The axis to compute spatial average on, by default ["X", "Y"]. Options
        include "X" and "Y".

    Returns
    -------
    xr.DataArray
        The spatial average of the variable based on the specified axis.

    Raises
    ------
    ValueError
        If the axis argument contains an invalid value.
    """
    for k in axis:
        if k not in ["X", "Y"]:
            raise ValueError(
                f"The `axis` argument has an unsupported value ('{k}'). "
                "Supported values include: ['X'], ['Y'], ['X', 'Y']."
            )

    ds_avg = ds.spatial.average(var_key, axis=axis, weights="generate")

    return ds_avg[var_key]


def std_xr(ds: xr.Dataset, var_key: str, axis=["X", "Y"]) -> xr.DataArray:
    """
    Compute a variable's weighted standard deviation on spatial axes.

    This function is intended to replace ``e3sm_diags.metrics.std()``.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset containing the variable.
    var_key : str
        The key of the variable.
    axis : list, optional
        The spatial axis to compute standard deviation on, by
        default ["X", "Y"]. Options include "X" and "Y".

    Returns
    -------
    xr.DataArray
        The standard deviation of the variable based on the specified axis.

    Raises
    ------
    ValueError
        If the axis argument contains an invalid value.
    """
    for k in axis:
        if k not in ["X", "Y"]:
            raise ValueError(
                f"The `axis` argument has an unsupported value ('{k}'). "
                "Supported values include: ['X'], ['Y'], ['X', 'Y']."
            )

    dims = []
    dv = ds[var_key].copy()

    # Get the weights for the data variable based on the specified axes.
    weights = ds.spatial.get_weights(axis, data_var=var_key)

    # Get the dimensions related to the axis.
    for a in axis:
        dim_key = xc.get_dim_keys(dv, axis=a)
        dims.append(dim_key)

    # Calculate weighted standard deviation.
    dv_std = dv.weighted(weights).std(dim=dims, keep_attrs=True)

    return dv_std


def rmse(ds: xr.DataArray, model: xr.DataArray, obs: xr.DataArray, axis=["X", "Y"]):
    # TODO: Look at xskillscore
    # - https://xskillscore.readthedocs.io/en/stable/api/xskillscore.rmse.html
    # https://github.com/xarray-contrib/xskillscore/blob/88474c98ad6078ebce9624b97b5afa4af0ba6e03/xskillscore/core/np_deterministic.py#L588-L628
    # https://github.com/scikit-learn/scikit-learn/blob/364c77e04/sklearn/metrics/_regression.py#L382

    # dims = xc.get_dim_keys(ds[model.name], axis=axis)
    # weights = ds.spatial.get_weights(model.name, axis=axis)
    # result = xs.rsme(model, obs, dim=dims, weight=weights)

    pass

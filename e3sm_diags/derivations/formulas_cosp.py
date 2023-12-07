from typing import Dict, List, Literal, Optional, TypedDict

import numpy as np
import xarray as xr

# Type annotations for CLOUD_HIST_MAP.
CloudAxis = Literal["prs", "tau"]
CloudHistMapAttrs = TypedDict(
    "CloudHistMapAttrs",
    {"keys": List[str], "default_bnds": np.ndarray, "min_mask": float},
)

# A dictionary storing attributes for "prs" and "tau" cloud axes.
# 1. 'keys' -- a list of valid variable keys found in the dataset
# 2. 'default_bnds' -- an array of floats for default bounds values if missing
# 3. 'min_mask' -- the minimum value to standardize the axis on for histogram
#    plotting.
CLOUD_HIST_MAP: Dict[CloudAxis, CloudHistMapAttrs] = {
    "prs": {
        "keys": ["cosp_htmisr", "misr_cth"],
        "default_bnds": np.array(
            [
                [1000.0, 800.0],
                [800.0, 680.0],
                [680.0, 560.0],
                [560.0, 440.0],
                [440.0, 310.0],
                [310.0, 180.0],
                [180.0, 0.0],
            ]
        ),
        "min_mask": 0.0,
    },
    "tau": {
        "keys": ["cosp_tau", "cosp_tau_modis", "modis_tau", "misr_tau", "isccp_tau"],
        "default_bnds": np.array(
            [
                [0.3, 1.3],
                [1.3, 3.6],
                [3.6, 9.4],
                [9.4, 23],
                [23, 60],
                [60, 379],
            ]
        ),
        "min_mask": 0.3,
    },
}


def cosp_histogram_standardize(ds: xr.Dataset, var: xr.DataArray) -> xr.Dataset:
    """Standardize cloud top pressure and cloud thickness bins.

    This standardization makes the cloud variable data suitable for plotting.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset with the cloud variable, cloud axes (prs, tau), and
        cloud axes bounds (optional). If bounds don't exist, they will be added
        using default values (refer to ``CLOUD_HIST_MAP``).
    var : xr.DataArray
        The cloud variable to be standardized.

    Returns
    -------
    xr.Dataset
        The dataset with the standardized cloud variable.
    """
    ds_new = ds.copy()
    prs = _get_cloud_axis(ds_new, var, "prs")
    tau = _get_cloud_axis(ds_new, var, "tau")

    # Mask on the min and max of prs and/or tau, then subset by dropping masked
    # values.
    var_new_std = _subset_cloud_var(var.copy(), prs, tau)

    # Align the dataset dimensions with the standardized variabe. `xr.align`
    # returns both objects and element 0 is the xr.Dataset that is needed.
    ds_final: xr.Dataset = xr.align(ds_new, var_new_std)[0]  # type: ignore
    ds_final[var_new_std.name] = var_new_std

    ds_final = _add_missing_cloud_bnds(ds_final, prs, "prs")
    ds_final = _add_missing_cloud_bnds(ds_final, tau, "tau")

    return ds_final


def _get_cloud_axis(ds: xr.Dataset, var: xr.DataArray, axis: CloudAxis) -> xr.DataArray:
    da_axis = None

    keys = CLOUD_HIST_MAP[axis]["keys"]
    for key in keys:
        if key in var.dims:
            da_axis = ds[key]

            break

    if da_axis is None:
        raise KeyError(
            f"The '{axis}' axis is not in the '{var.name}' to "
            "standardize the cosp histogram."
        )

    return da_axis


def _subset_cloud_var(
    var: xr.DataArray, prs: xr.DataArray, tau: xr.DataArray
) -> xr.DataArray:
    prs_min = CLOUD_HIST_MAP["prs"]["min_mask"]
    prs_max = prs[-1].item()

    tau_min = CLOUD_HIST_MAP["tau"]["min_mask"]
    tau_max = tau[-1].item()

    # MISR model and MISR obs
    if var.name in ["CLD_MISR", "CLMISR"]:
        # COSP v2 cosp_htmisr in units m instead of km as in v1 and COSP v2
        # cosp_htmisr[0] equals to 0 instead of -99 as in v1 so the cloud
        # varable needs to be masked manually by slicing out the first index
        # on the cosp_htmisr axis.
        if var.name == "CLD_MISR" and max(prs) > 1000:
            var = var.isel({prs.name: slice(1, None)})
            prs_max = 1000.0 * prs_max

        cond = (prs >= prs_min) & (prs <= prs_max) & (tau >= tau_min) & (tau <= tau_max)
    # ISCCP model, # ISCCP obs, and MODIS
    elif var.name in ["FISCCP1_COSP", "CLISCCP", "CLMODIS"]:
        cond = (tau >= tau_min) & (tau <= tau_max)

    var_sub = var.where(cond, drop=True)

    return var_sub


def _add_missing_cloud_bnds(
    ds: xr.Dataset, axis_var: xr.DataArray, axis: CloudAxis
) -> xr.Dataset:
    bnds_key = axis_var.attrs.get("bounds")

    if bnds_key is None or bnds_key not in ds.data_vars.keys():
        bnds_data = CLOUD_HIST_MAP[axis]["default_bnds"]
        default_bnds = xr.DataArray(
            name=f"{axis_var.name}_bnds",
            data=bnds_data,
            dims=list(axis_var.dims) + ["bnds"],
            coords=axis_var.coords,
        )
        ds[default_bnds.name] = default_bnds

    return ds


def cosp_bin_sum(
    cld: xr.DataArray,
    prs_low0: Optional[float],
    prs_high0: Optional[float],
    tau_low0: Optional[float],
    tau_high0: Optional[float],
):
    pass

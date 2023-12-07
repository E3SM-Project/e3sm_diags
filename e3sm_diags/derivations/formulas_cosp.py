from typing import Dict, List, Literal, Tuple, TypedDict

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


# TODO: Update these to reflect E3SM variable name in derived variables dicitonary
# Add sumulator name to the list of variables
CLD_BIN_SUM_RANGE = {
    # ISSCCP
    "CLDTOT_TAU1.3_ISCCP": {"prs": (None, None), "tau": (1.3, None)},
    "CLDTOT_TAU1.3_9.4_ISCCP": {"prs": (None, None), "tau": (1.3, 9.4)},
    "CLDTOT_TAU9.4_ISCCP": {"prs": (None, None), "tau": (9.4, None)},
    # MODIS
    "CLDTOT_TAU1.3_MODIS": {"prs": (None, None), "tau": (1.3, None)},
    "CLDTOT_TAU1.3_9.4_MODIS": {"prs": (None, None), "tau": (1.3, 9.4)},
    "CLDTOT_TAU9.4_MODIS": {"prs": (None, None), "tau": (9.4, None)},
    "CLDHGH_TAU1.3_MODIS": {"prs": (440, 0), "tau": (1.3, None)},
    "CLDHGH_TAU1.3_9.4_MODIS": {"prs": (440, 0), "tau": (1.3, 9.4)},
    "CLDHGH_TAU9.4_MODIS": {"prs": (440, 0), "tau": (9.4, None)},
    # MISR
    "CLDTOT_TAU1.3_MISR": {"prs": (None, None), "tau": (1.3, None)},
    "CLDTOT_TAU1.3_9.4_MISR": {"prs": (None, None), "tau": (1.3, 9.4)},
    "CLDTOT_TAU9.4_MISR": {"prs": (None, None), "tau": (9.4, None)},
    "CLDLOW_TAU1.3_MISR": {"prs": (0, 3), "tau": (1.3, None)},
    "CLDLOW_TAU1.3_9.4_MISR": {"prs": (0, 3), "tau": (1.3, 9.4)},
    "CLDLOW_TAU9.4_MISR": {"prs": (0, 3), "tau": (9.4, None)},
}

# COSP v2 cosp_pr in units Pa instead of hPa as in v1
# COSP v2 cosp_htmisr in units m instead of km as in v1
PRS_IDS_TO_ADJ_UNITS = {"cosp_prs": 100, "cosp_htmisr": 1000}


def cosp_bin_sum(var: xr.DataArray) -> xr.DataArray:
    # Functions used:  `determine_tau`, `determine_cloud_level`
    # Logic:
    # 1. Get cloud axes
    prs = _get_cloud_axis(var, "prs")
    tau = _get_cloud_axis(var, "tau")

    # 2. Get cloud ranges, lim , and sim
    prs_range = _get_prs_range(prs)
    prs_cld_lvl, prs_sim = _get_prs_cld_lvl_and_sim(var, prs, prs_range)
    tau_range, tau_lim = _get_tau_range_and_lim(tau)

    # 3. Subset the variable on the mask.
    cond = (
        (prs >= prs_range[0])
        & (prs <= prs_range[-1])
        & (tau >= tau_range[0])
        & (tau <= tau_range[-1])
    )
    var_sub = var.where(cond, drop=True)

    # 4. Sum on axis=0 and axis=1 (tau and prs)
    var_sum = var_sub.sum(dims=[prs.name, tau.name])

    # 5. Set the variable's long name.
    if prs_sim is not None and prs_cld_lvl is not None:
        var_sum.long_name = prs_sim + ": " + prs_cld_lvl + " with " + tau_lim

    return var_sum


def _get_prs_range(
    prs: xr.DataArray,
) -> Tuple[float, ...]:
    """Get the pr axis range for subsetting.

    This function loops over the min and max values for the actual and adjusted
    pr range. If the adjusted range is not None, then use that as the range
    value. Otherwise use the actual value from the prs axis.

    Adjust the units if the axis key is in ``PR_IDS_TO_ADJUST_UNITS``.

    Parameters
    ----------
    prs : xr.DataArray
        The prs axis.

    Returns
    -------
    Tuple[float, ...]
    """
    prs_actual_range = (prs[0].item(), prs[-1].item())
    prs_adj_range = CLD_BIN_SUM_RANGE[prs.name]["prs"]  # type: ignore

    prs_final_range: List[float] = []
    for act, adj in zip(prs_actual_range, prs_adj_range):
        if adj is not None:
            if prs.name in PRS_IDS_TO_ADJ_UNITS.keys() and max(prs.item()) > 1000:
                adj = adj * PRS_IDS_TO_ADJ_UNITS[prs.id]

            prs_final_range.append(adj)
        else:
            prs_final_range.append(act)

    return tuple(prs_final_range)


def _get_prs_cld_lvl_and_sim(var: xr.DataArray, prs, prs_range: Tuple[float, ...]):
    prs_low, prs_high = prs_range
    prs_low0, prs_high0 = CLD_BIN_SUM_RANGE[prs.name]["prs"]  # type: ignore

    prs_lim = None
    prs_sim = None

    if prs_low0 is None and prs_high0 is None:
        prs_lim = "total cloud fraction"

    # TODO: This does not cover all cases.
    if var.name == "FISCCP1_COSP":
        prs_sim = "ISCCP"
    elif var.name == "CLMODIS":
        prs_lim = _determine_cloud_level_with_prs(
            prs_low, prs_high, (440, 44000), (680, 68000)
        )
        prs_sim = "MODIS"
    elif var.name == "CLD_MISR":
        prs_lim = _determine_cloud_level_with_prs(
            prs_low, prs_high, (7, 7000), (3, 3000)
        )
        prs_sim = "MISR"

    return prs_lim, prs_sim


def _determine_cloud_level_with_prs(prs_low, prs_high, low_bnds, high_bnds):
    """Determines the cloud type based on prs values and the specified boundaries"""
    # Threshold for cloud top height: high cloud (<440hPa or > 7km), midlevel cloud (440-680hPa, 3-7 km) and low clouds (>680hPa, < 3km)
    # TODO: Refactor htis
    if prs_low in low_bnds and prs_high in high_bnds:
        return "middle cloud fraction"
    elif prs_low in low_bnds:
        return "high cloud fraction"
    elif prs_high in high_bnds:
        return "low cloud fraction"
    else:
        return "total cloud fraction"


def _get_tau_range_and_lim(tau: xr.DataArray) -> Tuple[Tuple[float, float], str]:
    """Get the tau range for subsetting ahd the tau lim.

    Parameters
    ----------
    tau : xr.DataArray
        The tau axis.

    Returns
    -------
    Tuple[Tuple[float, ...], str]
        A tuple consisting of the prs range (min, max) and the tau lim string.
    """
    # The actual values from the tau axis.
    tau_low_act, tau_high_act = tau[0].item(), tau[-1].item()
    # The adjusted values to use if either/both values are not None.
    tau_low_adj, tau_high_adj = CLD_BIN_SUM_RANGE[tau.name]["tau"]  # type: ignore

    # 1. Adjust high.
    if tau_low_adj is None and tau_high_adj is not None:
        return_tuple = (tau_low_act, tau_high_adj)
        tau_lim = "tau <" + str(tau_high_adj)
    # 2. Adjust low.
    elif tau_low_adj is not None and tau_high_adj is None:
        return_tuple = (tau_low_adj, tau_high_act)
        tau_lim = "tau >" + str(tau_low_adj)
    # 3. Adjust low and high.
    elif tau_low_adj is not None and tau_high_adj is not None:
        return_tuple = (tau_low_adj, tau_high_adj)
        tau_lim = str(tau_low_act) + "< tau < " + str(tau_high_act)
    # 4. No adjustments, use actual values.
    else:
        return_tuple = (tau_low_act, tau_high_act)
        tau_lim = str(tau_low_act) + "< tau < " + str(tau_high_act)

    return return_tuple, tau_lim


def cosp_histogram_standardize(
    ds: xr.Dataset, target_var_key: str, var: xr.DataArray
) -> xr.Dataset:
    """Standardize cloud top pressure and cloud thickness bins.

    This standardization makes the cloud variable data suitable for plotting.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset with the cloud variable, cloud axes (prs, tau), and
        cloud axes bounds (optional). If bounds don't exist, they will be added
        using default values (refer to ``CLOUD_HIST_MAP``).
    target_var_key : str
        The key of the target variable used to add to the dataset.
    var : xr.DataArray
        The cloud variable to be standardized.

    Returns
    -------
    xr.Dataset
        The dataset with the standardized cloud variable.
    """
    ds_new = ds.copy()
    prs = _get_cloud_axis(var, "prs")
    tau = _get_cloud_axis(var, "tau")

    # Mask on the min and max of prs and/or tau, then subset by dropping masked
    # values.
    var_std = _subset_cloud_var(var.copy(), prs, tau)

    # Align the dataset dimensions with the standardized variabe. `xr.align`
    # returns both objects and element 0 is the xr.Dataset that is needed.
    ds_final: xr.Dataset = xr.align(ds_new, var_std)[0]  # type: ignore
    ds_final[target_var_key] = var_std
    ds_final = ds_final.drop_vars(var.name)

    ds_final = _add_missing_cloud_bnds(ds_final, prs, "prs")
    ds_final = _add_missing_cloud_bnds(ds_final, tau, "tau")

    return ds_final


def _get_cloud_axis(var: xr.DataArray, axis: CloudAxis) -> xr.DataArray:
    da_axis = None

    keys = CLOUD_HIST_MAP[axis]["keys"]
    for key in keys:
        if key in var.dims:
            da_axis = var[key]

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
    # FIXME: ValueError: conflicting sizes for dimension 'misr_cth': length 7
    # on the data but length 16 on coordinate 'misr_cth'
    # https://github.com/E3SM-Project/e3sm_diags/issues/761
    if bnds_key is None or bnds_key not in ds.data_vars.keys():
        # TODO: Might need to dynamically generate bounds if they don't exist.
        bnds_data = CLOUD_HIST_MAP[axis]["default_bnds"]
        default_bnds = xr.DataArray(
            name=f"{axis_var.name}_bnds",
            data=bnds_data,
            dims=list(axis_var.dims) + ["bnds"],
            coords=axis_var.coords,
        )
        ds[default_bnds.name] = default_bnds

    return ds

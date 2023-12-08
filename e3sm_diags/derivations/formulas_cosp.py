from typing import Dict, List, Literal, Optional, Tuple, TypedDict

import numpy as np
import xarray as xr

from e3sm_diags.derivations.formulas import convert_units

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


# A dictionary mapping the target variable key to the "prs" and "tau" axes
# adjustment ranges. If either value in the (min, max) tuple is None, then
# the actual value from the axis is used instead.
CLOUD_AXES_ADJ_RANGES: Dict[str, Dict[str, Tuple[Optional[float], Optional[float]]]] = {
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

# A dictionary mapping the names of "prs" axes to the unit adjustment value.
# - COSP v2 "cosp_pr" is in units "Pa" instead of "hPa" (v1).
# - COSP v2 "cosp_htmisr" is in units "m" instead of "km" (v1).
PRS_NAMES_TO_ADJ_UNITS = {"cosp_prs": 100, "cosp_htmisr": 1000}


def cosp_bin_sum(target_var_key: str, var: xr.DataArray) -> xr.DataArray:
    # 1. Get cloud axes
    prs = _get_cloud_axis(var, "prs")
    tau = _get_cloud_axis(var, "tau")

    # 2. Get prs range, lim, and sim
    prs_adj_range = CLOUD_AXES_ADJ_RANGES[target_var_key]["prs"]
    prs_act_range = _get_prs_range(prs, prs_adj_range)

    # 3. Get tau range and lim.
    tau_range, tau_lim = _get_tau_range_and_lim(tau)

    # 4. Get the axes mask conditiona and subset the variable on it.
    cond = (
        (prs >= prs_act_range[0])
        & (prs <= prs_act_range[-1])
        & (tau >= tau_range[0])
        & (tau <= tau_range[-1])
    )
    var_sub = var.where(cond, drop=True)

    # 5. Sum on axis=0 and axis=1 (tau and prs)
    var_sum = var_sub.sum(dims=[prs.name, tau.name])

    # 6. Set the variable's long name based on the prs name and ranges.
    prs_name = str(prs.name)
    prs_sim = _get_prs_sim(prs_name)
    prs_cld_lvl = _get_prs_cloud_level(prs_name, prs_act_range, prs_adj_range)

    if prs_sim is not None and prs_cld_lvl is not None:
        var_sum.long_name = prs_sim + ": " + prs_cld_lvl + " with " + tau_lim

    # 7. Conver units to %.
    final_var = convert_units(var_sum, "%")

    return final_var


def _get_prs_range(
    prs: xr.DataArray, adj_range: Tuple[Optional[float], Optional[float]]
) -> Tuple[float, float]:
    """Get the pr axis range for subsetting.

    This function loops over the min and max values for the actual and adjusted
    pr range. If the adjusted range is not None, then use that as the range
    value. Otherwise use the actual value from the prs axis.

    Adjust the units if the axis key is in ``PR_NAMES_TO_ADJUST_UNITS``.

    Parameters
    ----------
    prs : xr.DataArray
        The prs axis.
    target_var_key : str
        The key of the target variable (e.g., "CLDTOT_TAU1.3_ISCCP"). This key
        is used to get the appropriate pr adjustment ranges.

    Returns
    -------
    Tuple[float, float]
        A tuple of the (min, max) for the prs range.
    """
    act_range = (prs[0].item(), prs[-1].item())
    final_range: List[float] = []

    for act_val, adj_val in zip(act_range, adj_range):
        if adj_val is not None:
            if prs.name in PRS_NAMES_TO_ADJ_UNITS.keys() and max(prs.item()) > 1000:
                adj_val = adj_val * PRS_NAMES_TO_ADJ_UNITS[prs.id]

            final_range.append(adj_val)
        else:
            final_range.append(act_val)

    return tuple(final_range)  # type: ignore


def _get_tau_range_and_lim(tau: xr.DataArray) -> Tuple[Tuple[float, float], str]:
    """Get the tau range for subsetting ahd the tau lim.

    Parameters
    ----------
    tau : xr.DataArray
        The tau axis.

    Returns
    -------
    Tuple[Tuple[float, float], str]
        A tuple consisting of the prs range (min, max) and the tau lim string.
    """
    # The actual min and max values from the tau axis.
    act_min, act_max = tau[0].item(), tau[-1].item()

    # The adjusted min and max values to use if one or both values are not None.
    adj_min, adj_max = CLOUD_AXES_ADJ_RANGES[tau.name]["tau"]  # type: ignore

    # 1. Adjust min.
    if adj_min is not None and adj_max is None:
        new_range = (adj_min, act_max)
        lim = "tau >" + str(adj_min)
    # 2. Adjust high.
    elif adj_min is None and adj_max is not None:
        new_range = (act_min, adj_max)
        lim = "tau <" + str(adj_max)
    # 3. Adjust min and max.
    elif adj_min is not None and adj_max is not None:
        new_range = (adj_min, adj_max)
        lim = str(adj_min) + "< tau < " + str(adj_max)
    # 4. No adjustments, use actual axis min and max values.
    else:
        new_range = (act_min, act_max)
        lim = str(act_min) + "< tau < " + str(act_max)

    return new_range, lim


def _get_prs_sim(prs_name: str) -> Optional[str]:
    """Get the prs simulation.

    Parameters
    ----------
    prs_name : str
        The prs axis name.

    Returns
    -------
    Optional[str]
        The optional prs simulation string.
    """
    sim = None

    # NOTE: This does not cover all cases.
    if prs_name == "FISCCP1_COSP":
        sim = "ISCCP"
    elif prs_name == "CLMODIS":
        sim = "MODIS"
    elif prs_name == "CLD_MISR":
        sim = "MISR"

    return sim


def _get_prs_cloud_level(
    prs_name: str,
    prs_act_range: Tuple[float, float],
    prs_adj_range: Tuple[Optional[float], Optional[float]],
) -> Optional[str]:
    """Get the prs cloud level and simulation type.

    Parameters
    ----------
    prs_name : str
        The prs axis name
    prs_act_range : Tuple[float, float]
        The prs actual range.
    prs_adj_range : Tuple[Optional[float], Optional[float]]
        The prs adjusted range.

    Returns
    -------
    Optional[str]
        The optional cloud level string.
    """
    adj_min, adj_max = prs_adj_range
    cloud_level = None

    if adj_min is None and adj_max is None:
        cloud_level = "total cloud fraction"

    # NOTE: This does not cover all cases.
    if prs_name == "CLMODIS":
        cloud_level = _get_cloud_level(prs_act_range, (440, 44000), (680, 68000))
    elif prs_name == "CLD_MISR":
        cloud_level = _get_cloud_level(prs_act_range, (7, 7000), (3, 3000))

    return cloud_level


def _get_cloud_level(prs_act_range, low_bnds, high_bnds) -> str:
    """Determines the cloud type based on prs values and the specified boundaries

    Thresholds for cloud levels:
      - cloud top height: high cloud (<440hPa or > 7km)
       - midlevel cloud (440-680hPa, 3-7 km)
       - low clouds (>680hPa, < 3km)
    """
    min, max = prs_act_range

    if min in low_bnds and max in high_bnds:
        return "middle cloud fraction"
    elif min in low_bnds:
        return "high cloud fraction"
    elif max in high_bnds:
        return "low cloud fraction"

    return "total cloud fraction"


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

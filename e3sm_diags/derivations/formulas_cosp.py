from typing import Literal, TypedDict

import xarray as xr

from e3sm_diags.derivations.formulas import convert_units

# A dictionary storing attributes for "prs" and "tau" cloud axes.
# 1. 'keys': a list of valid variable keys found in the dataset, used for
#     dynamic mapping.
# 2. 'min_mask': the minimum value to standardize the axis on for histogram
#     plotting.
CloudAxis = Literal["prs", "tau"]
CloudHistMapAttrs = TypedDict(
    "CloudHistMapAttrs", {"keys": list[str], "min_mask": float}
)
CLOUD_HIST_MAP: dict[CloudAxis, CloudHistMapAttrs] = {
    "prs": {
        "keys": ["cosp_prs", "cosp_htmisr", "modis_prs", "misr_cth", "isccp_prs"],
        "min_mask": 0.0,
    },
    "tau": {
        "keys": ["cosp_tau", "cosp_tau_modis", "modis_tau", "misr_tau", "isccp_tau"],
        "min_mask": 0.3,
    },
}

# A dictionary mapping the target variable key to the "prs" and "tau" axes
# adjustment ranges. If either value in the (min, max) tuple is None, then
# the actual value from the axis is used instead.
AdjRange = tuple[float | None, float | None]
CLOUD_BIN_SUM_MAP: dict[str, dict[str, AdjRange]] = {
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
PRS_UNIT_ADJ_MAP = {"cosp_prs": 100, "cosp_htmisr": 1000}


def cosp_histogram_standardize(target_var_key: str, var: xr.DataArray) -> xr.DataArray:
    """Standardize cloud top pressure and cloud thickness bins.

    This standardization makes the cloud variable data suitable for plotting.

    Parameters
    ----------
    target_var_key : str
        The target variable key (e.g,. "CLDTOT_TAU1.3_ISCCP").
    var : xr.DataArray
        The variable in the dataset used for deriving the target variable
        (e.g., "CLD_MODIS").

    Returns
    -------
    xr.DataArray
        The target variable, which is the standardized version of the derived
        variable (``var``).
    """
    prs = _get_cloud_axis(var, "prs")
    tau = _get_cloud_axis(var, "tau")

    # Mask on the min and max of prs and/or tau, then subset by dropping masked
    # values.
    var_std = _subset_cloud_var(var.copy(), prs, tau)
    var_std.name = target_var_key

    return var_std


def cosp_bin_sum(target_var_key: str, var: xr.DataArray) -> xr.DataArray:
    """Get the cosp bin sum for the derived variable.

    Parameters
    ----------
    target_var_key : str
        The target variable key (e.g,. "CLDTOT_TAU1.3_ISCCP").
    var : xr.DataArray
        The variable in the dataset used for deriving the target variable
        (e.g., "CLD_MODIS").

    Returns
    -------
    xr.DataArray
        The target variable, which is the cosp bin sum of the derived variable,
        ``var``.
    """
    # 1. Get cloud axes
    prs = _get_cloud_axis(var, "prs")
    tau = _get_cloud_axis(var, "tau")

    # 2. Get the prs and tau axis adjustment ranges if they are set.
    prs_adj_range = CLOUD_BIN_SUM_MAP[target_var_key]["prs"]
    tau_adj_range = CLOUD_BIN_SUM_MAP[target_var_key]["tau"]

    # 3. Get prs range, lim.
    prs_range = _get_prs_subset_range(prs, prs_adj_range)

    # 4. Get tau range and lim.
    tau_range, tau_lim = _get_tau_subset_range_and_str(tau, tau_adj_range)

    # 5. Get the axes mask conditional and subset the variable on it.
    cond = _get_prs_and_tau_cond(prs, tau, prs_range, tau_range)
    var_sub = var.where(cond, drop=True)

    # 7. Sum on axis=0 and axis=1 (tau and prs)
    var_sum = var_sub.sum(
        dim=[prs.name, tau.name], keep_attrs=True, skipna=True, min_count=1
    )

    # 8. Set the variable's long name based on the original variable's name and
    # prs ranges.
    var_key = str(var.name)
    simulation = _get_simulation_str(var_key)
    prs_cloud_level = _get_prs_cloud_level_str(var_key, prs_range, prs_adj_range)

    if simulation is not None and prs_cloud_level is not None:
        var_sum.attrs["long_name"] = f"{simulation}: {prs_cloud_level} with {tau_lim}"

    # 9. Convert units to %.
    final_var = convert_units(var_sum, "%")

    return final_var


def _get_cloud_axis(var: xr.DataArray, axis: CloudAxis) -> xr.DataArray:
    da_axis = None

    keys = CLOUD_HIST_MAP[axis]["keys"]
    for key in keys:
        if key in var.dims:
            da_axis = var[key]

            break

    if da_axis is None:
        raise KeyError(
            f"The {axis!r} axis is not in the {var.name!r} to "
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
        if var.name == "CLD_MISR" and prs.max().item() > 1000:
            var = var.isel({prs.name: slice(1, None)})
            prs_max = 1000.0 * prs_max

        cond = (prs >= prs_min) & (prs <= prs_max) & (tau >= tau_min) & (tau <= tau_max)
    # ISCCP model, # ISCCP obs, and MODIS
    elif var.name in ["FISCCP1_COSP", "CLISCCP", "CLMODIS"]:
        cond = (tau >= tau_min) & (tau <= tau_max)

    var_sub = var.where(cond, drop=True)

    return var_sub


def _get_prs_subset_range(
    prs: xr.DataArray, prs_adj_range: AdjRange
) -> tuple[float, float]:
    """Get the pr axis subset range.

    This function loops over the min and max values for the actual and adjusted
    pr range. If the adjusted range is not None, then use that as the range
    value. Otherwise use the actual value from the prs axis.

    Adjust the units if the axis key is in ``PR_NAMES_TO_ADJUST_UNITS``.

    Parameters
    ----------
    prs : xr.DataArray
        The prs axis.
    prs_adj_range : AdjRange
        The prs axis adjustment range.

    Returns
    -------
    tuple[float, float]
        A tuple of the (min, max) for the prs subset range.
    """
    act_range = (prs[0].item(), prs[-1].item())
    range: list[float] = []
    for act_val, adj_val in zip(act_range, prs_adj_range, strict=False):
        if adj_val is not None:
            if prs.name in PRS_UNIT_ADJ_MAP.keys() and prs.max().item() > 1000:
                adj_val = adj_val * PRS_UNIT_ADJ_MAP[str(prs.name)]

            range.append(adj_val)
        else:
            range.append(act_val)

    return range  # type: ignore


def _get_tau_subset_range_and_str(
    tau: xr.DataArray, tau_adj_range: AdjRange
) -> tuple[tuple[float, float], str]:
    """Get the tau range for subsetting and the tau range string.

    Parameters
    ----------
    tau : xr.DataArray
        The tau axis.
    tau_adj_range : AdjRange
        The tau axis adjustment range.

    Returns
    -------
    tuple[tuple[float, float], str]
        A tuple consisting of the tau range (min, max) and the tau range string.
    """
    # The adjusted min and max values to use if either if them are None.
    adj_min, adj_max = tau_adj_range

    # 1. Adjust min.
    if adj_min is not None and adj_max is None:
        act_max = tau[-1].item()

        range = (adj_min, act_max)
        range_str = f"tau > {adj_min}"
    # 2. Adjust min and max.
    elif adj_min is not None and adj_max is not None:
        range = (adj_min, adj_max)
        range_str = f"{adj_min} < tau < {adj_max}"

    final_range = tuple(range)

    return final_range, range_str


def _get_prs_and_tau_cond(
    prs: xr.DataArray,
    tau: xr.DataArray,
    prs_range: tuple[float, float],
    tau_range: tuple[float, float],
) -> xr.DataArray:
    """Get the prs and tau condition for sub-setting a variable.

    Parameters
    ----------
    prs : xr.DataArray
        The prs axis.
    tau : xr.DataArray
        The tau axis.
    prs_range : tuple[float, float]
        The range of prs values to subset with.
    tau_range : tuple[float, float]
        The range of tau values to subset with.

    Returns
    -------
    xr.DataArray
        The condition dataarray.
    """
    # Values must be sorted in ascending order to correctly subset within a
    # range using Xarray's `.where()` method.
    sorted_prs_range = sorted(list(prs_range))
    sorted_tau_range = sorted(list(tau_range))

    cond = (
        (prs >= sorted_prs_range[0])
        & (prs <= sorted_prs_range[-1])
        & (tau >= sorted_tau_range[0])
        & (tau <= sorted_tau_range[-1])
    )

    return cond


def _get_simulation_str(var_key: str) -> str | None:
    """Get the prs simulation string.

    Parameters
    ----------
    var_key : str
        The key of the variable in the dataset.

    Returns
    -------
    str | None
        The optional prs simulation string.
    """
    sim = None

    # NOTE: This does not cover all variable keys.
    if var_key == "FISCCP1_COSP":
        sim = "ISCCP"
    elif var_key == "CLMODIS":
        sim = "MODIS"
    elif var_key == "CLD_MISR":
        sim = "MISR"

    return sim


def _get_prs_cloud_level_str(
    var_key: str,
    prs_act_range: tuple[float, float],
    prs_adj_range: AdjRange,
) -> str | None:
    """Get the prs cloud level and simulation type.

    Parameters
    ----------
    var_key: str
        The key of the variable in the dataset.
    prs_act_range : tuple[float, float]
        The prs actual range.
    prs_adj_range : AdjRange
        The prs adjusted range.

    Returns
    -------
    str | None
        The optional cloud level string.
    """
    adj_min, adj_max = prs_adj_range
    cloud_level = None

    if adj_min is None and adj_max is None:
        cloud_level = "total cloud fraction"

    # NOTE: This does not cover all variable keys, including "FISCCP1_COSP".
    if var_key == "CLMODIS":
        cloud_level = _get_cloud_level(
            prs_act_range, low_bnds=(440, 44000), high_bnds=(680, 68000)
        )
    elif var_key == "CLD_MISR":
        cloud_level = _get_cloud_level(
            prs_act_range, low_bnds=(7, 7000), high_bnds=(3, 3000)
        )

    return cloud_level


def _get_cloud_level(
    prs_act_range: tuple[float, float],
    low_bnds: tuple[float, float],
    high_bnds: tuple[float, float],
) -> str:
    """Get the cloud type based on prs values and the specified boundaries

    Thresholds for cloud levels:
      1. cloud top height: high cloud (<440hPa or > 7km)
      2. mid-level cloud (440-680hPa, 3-7 km)
      3. low clouds (>680hPa, < 3km)

    Parameters
    ----------
    prs_act_range : tuple[float, float]
        The prs actual range.
    low_bnds : tuple[float, float]
        The low bounds.
    high_bnds : tuple[float, float]
        The high bounds.

    Returns
    -------
    str
        The cloud level string.
    """
    min, max = prs_act_range

    if min in low_bnds and max in high_bnds:
        return "middle cloud fraction"
    elif min in low_bnds:
        return "high cloud fraction"
    elif max in high_bnds:
        return "low cloud fraction"

    return "total cloud fraction"

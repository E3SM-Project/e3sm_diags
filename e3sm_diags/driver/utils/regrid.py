from __future__ import annotations

from typing import TYPE_CHECKING, List, Literal, Tuple

import xarray as xr
import xcdat as xc

from e3sm_diags.derivations.default_regions_xr import REGION_SPECS
from e3sm_diags.driver import FRAC_REGION_KEYS
from e3sm_diags.logger import custom_logger

if TYPE_CHECKING:
    from e3sm_diags.parameter.core_parameter import CoreParameter

logger = custom_logger(__name__)


# Valid hybrid-sigma levels keys that can be found in datasets.
# _get_hybrid_sigma_level()` will search for the variables using the keys
# in the order of elements in the tuple and will return the first match if
# found.
HYBRID_SIGMA_KEYS = {
    "p0": ("p0", "P0"),
    "ps": ("ps", "PS"),
    "hyam": ("hyam", "hyai", "hya", "a"),
    "hyai": ("hyai", "hyam", "hya", "a"),
    "hybm": ("hybm", "hybi", "hyb", "b"),
    "hybi": ("hybi", "hybm", "hyb", "b"),
}

REGRID_TOOLS = Literal["esmf", "xesmf", "regrid2"]


def subset_and_align_datasets(
    parameter: CoreParameter,
    ds_test: xr.Dataset,
    ds_ref: xr.Dataset,
    ds_land_sea_mask: xr.Dataset,
    var_key: str,
    region: str,
) -> Tuple[xr.Dataset, xr.Dataset, xr.Dataset, xr.Dataset, xr.Dataset]:
    """Subset ref and test datasets on a region and regrid to align them.

    Parameters
    ----------
    parameter : CoreParameter
        The parameter for the diagnostic.
    ds_test : xr.Dataset
        The dataset containing the test variable.
    ds_ref : xr.Dataset
        The dataset containing the ref variable.
    ds_land_sea_mask : xr.Dataset
        The land sea mask dataset, which is only used for masking if the region
        is "land" or "ocean".
    var_key : str
        The key of the variable.
    region : str
        The region.

    Returns
    -------
    Tuple[xr.Dataset, xr.Dataset, xr.Dataset, xr.Dataset, xr.Dataset]
        A tuple containing the test dataset, the regridded test
        dataset, the ref dataset, the regridded ref dataset, and the difference
        between regridded datasets.
    """
    ds_test_new = ds_test.copy()
    ds_ref_new = ds_ref.copy()

    logger.info(f"Selected region: {region}")
    parameter.var_region = region

    # Apply a land sea mask.
    if "land" in region or "ocean" in region:
        ds_test_new = _apply_land_sea_mask(
            ds_test_new,
            ds_land_sea_mask,
            var_key,
            region,  # type: ignore
            parameter.regrid_tool,
            parameter.regrid_method,
        )
        ds_ref_new = _apply_land_sea_mask(
            ds_ref_new,
            ds_land_sea_mask,
            var_key,
            region,  # type: ignore
            parameter.regrid_tool,
            parameter.regrid_method,
        )

    # Subset on a specific region.
    if "global" not in region:
        ds_test_new = _subset_on_region(ds_test_new, var_key, region)
        ds_ref_new = _subset_on_region(ds_ref_new, var_key, region)

    ds_test_regrid, ds_ref_regrid = align_grids_to_lower_res(
        ds_test_new,
        ds_ref_new,
        var_key,
        parameter.regrid_tool,
        parameter.regrid_method,
    )

    ds_diff = ds_test_regrid.copy()
    ds_diff[var_key] = ds_test_regrid[var_key] - ds_ref_regrid[var_key]

    return ds_test_new, ds_test_regrid, ds_ref_new, ds_ref_regrid, ds_diff


def has_z_axis(data_var: xr.DataArray) -> bool:
    """Checks whether the data variable has a Z axis.

    Parameters
    ----------
    data_var : xr.DataArray
        The data variable.

    Returns
    -------
    bool
        True if data variable has Z axis, else False.

    Notes
    -----
    Replaces `cdutil.variable.TransientVariable.getLevel()`.
    """
    try:
        get_z_axis(data_var)
        return True
    except KeyError:
        return False


def get_z_axis(obj: xr.Dataset | xr.DataArray) -> xr.DataArray:
    """Gets the Z axis coordinates from an xarray object.

    Returns True if:
    - Data variable has a "Z" axis in the cf-xarray mapping dict
    - A coordinate has a matching "positive" attribute ("up" or "down")
    - A coordinate has a matching "name"
    - # TODO: conditional for valid pressure units with "Pa"

    Parameters
    ----------
    obj : xr.Dataset | xr.DataArray
        The xarray Dataset or DataArray.

    Returns
    -------
    xr.DataArray
        The Z axis coordinates.

    Notes
    -----
    Based on
    - https://cdms.readthedocs.io/en/latest/_modules/cdms2/avariable.html#AbstractVariable.getLevel
    - https://cdms.readthedocs.io/en/latest/_modules/cdms2/axis.html#AbstractAxis.isLevel
    """
    try:
        z_coords = xc.get_dim_coords(obj, axis="Z")

        return z_coords
    except KeyError:
        pass

    for coord in obj.coords.values():
        if coord.name in ["lev", "plev", "depth"]:
            return coord

    raise KeyError(
        f"No Z axis coordinates were found in the {type(obj)}. Make sure the "
        f"{type(obj)} has Z axis coordinates."
    )


def _apply_land_sea_mask(
    ds: xr.Dataset,
    ds_mask: xr.Dataset,
    var_key: str,
    region: Literal["land", "ocean"],
    regrid_tool: str,
    regrid_method: str,
) -> xr.Dataset:
    """Apply a land or sea mask based on the region ("land" or "ocean").

    Parameters
    ----------
    ds: xr.Dataset
        The dataset containing the variable.
    ds_mask : xr.Dataset
        The dataset containing the land sea region mask variable(s).
    var_key : str
        The key the variable
    region : Literal["land", "ocean"]
        The region to mask.
    regrid_tool : {"esmf", "xesmf", "regrid2"}
        The regridding tool to use. Note, "esmf" is accepted for backwards
        compatibility with e3sm_diags and is simply updated to "xesmf".
    regrid_method : str
        The regridding method to use. Refer to [1]_ for more information on
        these options.

        esmf/xesmf options:
          - "bilinear"
          - "conservative"
          - "conservative_normed" -- equivalent to "conservative" in cdms2 ESMF
          - "patch"
          - "nearest_s2d"
          - "nearest_d2s"

        regrid2 options:
          - "conservative"

    Returns
    -------
    xr.Dataset
        The Dataset with the land or sea mask applied to the variable.
    """
    # TODO: Remove this conditional once "esmf" references are updated to
    # "xesmf" throughout the codebase.
    if regrid_tool == "esmf":
        regrid_tool = "xesmf"

    # TODO: Remove this conditional once "conservative" references are updated
    # to "conservative_normed" throughout the codebase.
    # NOTE: this is equivalent to "conservative" in cdms2 ESMF. If
    # "conservative" is chosen, it is updated to "conservative_normed". This
    # logic can be removed once the CoreParameter.regrid_method default
    # value is updated to "conservative_normed" and all sets have been
    # refactored to use this function.
    if regrid_method == "conservative":
        regrid_method = "conservative_normed"

    # A dictionary storing the specifications for this region.
    specs = REGION_SPECS[region]

    # If the region is land or ocean, regrid the land sea mask to the same
    # shape (lat x lon) as the variable then apply the mask to the variable.
    # Land and ocean masks have a region value which is used as the upper limit
    # for masking.
    ds_new = ds.copy()
    ds_new = _drop_unused_ilev_axis(ds)
    output_grid = ds_new.regridder.grid
    mask_var_key = _get_region_mask_var_key(ds_mask, region)

    ds_mask_new = _drop_unused_ilev_axis(ds_mask)
    ds_mask_regrid = ds_mask_new.regridder.horizontal(
        mask_var_key,
        output_grid,
        tool=regrid_tool,
        method=regrid_method,
    )
    # Update the mask variable with a lower limit. All values below the
    # lower limit will be masked.
    land_sea_mask = ds_mask_regrid[mask_var_key]
    lower_limit = specs["value"]  # type: ignore
    cond = land_sea_mask > lower_limit

    # Apply the mask with a condition (`cond`) using `.where()`. Note, the
    # condition matches values to keep, not values to mask out, `drop` is
    # set to False because we want to preserve the masked values (`np.nan`)
    # for plotting purposes.
    masked_var = ds_new[var_key].where(cond=cond, drop=False)

    ds_new[var_key] = masked_var

    return ds_new


def _subset_on_region(ds: xr.Dataset, var_key: str, region: str) -> xr.Dataset:
    """Subset a variable in the dataset based on the region.

    This function makes sure that the axes/axes being subsetted on is in
    ascending order. For example, if the latitude axis is in descending order,
    [90, 0, -90], no matches will be made on the subset if the region has a lat
    slice spec of (-90, -55) (e.g., 'polar_S'). This is because Xarray subsets
    in ascending order. Sorting the axis beforehand will avoid this issue.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset.
    var_key : str
        The variable to subset.
    region : str
        The region

    Returns
    -------
    xr.Dataset
        The dataset with the subsetted variable.

    Notes
    -----
    Replaces `e3sm_diags.utils.general.select_region`.
    """
    specs = REGION_SPECS[region]
    lat, lon = specs.get("lat"), specs.get("lon")  # type: ignore

    ds_new = ds.copy()

    if lon is not None:
        lon_dim = xc.get_dim_keys(ds[var_key], axis="X")
        ds_new = ds_new.sortby(lon_dim)

        # TODO: Add a unit test for this
        is_lon_axis_diff = lon[0] < 0 and ds_new[lon_dim].values[0] >= 0
        if is_lon_axis_diff:
            ds_new = xc.swap_lon_axis(ds_new, to=(-180, 180))

        ds_new = ds_new.sel({f"{lon_dim}": slice(*lon)})

    if lat is not None:
        lat_dim = xc.get_dim_keys(ds[var_key], axis="Y")
        ds_new = ds_new.sortby(lat_dim)
        ds_new = ds_new.sel({f"{lat_dim}": slice(*lat)})

    return ds_new


def _subset_on_arm_coord(ds: xr.Dataset, var_key: str, arm_site: str):
    """Subset a variable in the dataset on the specified ARM site coordinate.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset.
    var_key : str
        The variable to subset.
    arm_site : str
        The ARM site.

    Notes
    -----
    Replaces `e3sm_diags.utils.general.select_point`.
    """
    # TODO: Refactor this method with ARMS diagnostic set.
    pass  # pragma: no cover


def align_grids_to_lower_res(
    ds_a: xr.Dataset,
    ds_b: xr.Dataset,
    var_key: str,
    tool: REGRID_TOOLS,
    method: str,
    axis_to_compare: str = "Y",
) -> Tuple[xr.Dataset, xr.Dataset]:
    """Align the grids of two Dataset using the lower resolution of the two.

    Using the legacy logic, compare the number of latitude coordinates to
    determine if A or B has lower resolution:
      * If A is lower resolution (A <= B), regrid B -> A.
      * If B is lower resolution (A > B), regrid A -> B.

    Parameters
    ----------
    ds_a : xr.Dataset
        The first Dataset containing ``var_key``.
    ds_b : xr.Dataset
        The second Dataset containing ``var_key``.
    var_key : str
        The key of the variable in both datasets to regrid.
    tool : {"esmf", "xesmf", "regrid2"}
        The regridding tool to use. Note, "esmf" is accepted for backwards
        compatibility with e3sm_diags and is simply updated to "xesmf".
    method : str
        The regridding method to use. Refer to [1]_ for more information on
        these options.

        esmf/xesmf options:
          - "bilinear"
          - "conservative"
          - "conservative_normed"
          - "patch"
          - "nearest_s2d"
          - "nearest_d2s"

        regrid2 options:
          - "conservative"
    axis_to_compare : str
        The axis to use to compare resolutions to find the lower resolution
        of both variables, either "X" or "Y". by default "Y".

    Returns
    -------
    Tuple[xr.Dataset, xr.Dataset]
        A tuple of both DataArrays regridded to the lower resolution of the two.

    Notes
    -----
    Replaces `e3sm_diags.driver.utils.general.regrid_to_lower_res`.

    References
    ----------
    .. [1] https://xcdat.readthedocs.io/en/stable/generated/xarray.Dataset.regridder.horizontal.html
    """
    # TODO: Accept "esmf" as `tool` value for now because `CoreParameter`
    # defines `self.regrid_tool="esmf"` by default and
    # `e3sm_diags.driver.utils.general.regrid_to_lower_res()` accepts "esmf".
    # Once this function is deprecated, we can remove "esmf" as an option here
    # and update `CoreParameter.regrid_tool` to "xesmf"`.
    if tool == "esmf":
        tool = "xesmf"

    ds_a_new = _drop_unused_ilev_axis(ds_a)
    ds_b_new = _drop_unused_ilev_axis(ds_b)

    axis_a = xc.get_dim_coords(ds_a_new[var_key], axis=axis_to_compare)
    axis_b = xc.get_dim_coords(ds_b_new[var_key], axis=axis_to_compare)

    is_a_lower_res = len(axis_a) <= len(axis_b)

    if is_a_lower_res:
        output_grid = ds_a_new.regridder.grid
        ds_b_regrid = ds_b_new.regridder.horizontal(
            var_key, output_grid, tool=tool, method=method
        )

        return ds_a_new, ds_b_regrid

    output_grid = ds_b_new.regridder.grid
    ds_a_regrid = ds_a_new.regridder.horizontal(
        var_key, output_grid, tool=tool, method=method
    )

    return ds_a_regrid, ds_b_new


def _drop_unused_ilev_axis(ds: xr.Dataset) -> xr.Dataset:
    """Drop the unused ilev axis in a dataset.

    The ilev axis needs to be dropped prior to regridding with xCDAT. Otherwise,
    this error might be raised: `ValueError: Multiple 'Z' axis dims were found
    in this dataset, ['ilev', 'lev']. Please drop the unused dimension(s) before
    performing grid operations.`

    The ilev axis is usually associated with pressure variables such as "hyam"
    and "hybm".

    Parameters
    ----------
    ds : xr.Dataset
        The dataset with a lev and ilev axes.

    Returns
    -------
    xr.Dataset
        The dataset with a lev axis.
    """
    ds_new = ds.copy()
    if "ilev" in ds_new.dims:
        ds_new = ds_new.drop_dims("ilev")

    return ds_new


def _get_region_mask_var_key(ds_mask: xr.Dataset, region: str):
    """Get the region's mask variable key.

    This variable key can be used to map the the variable data in a dataset.
    Only land and ocean regions are supported.

    Parameters
    ----------
    ds_mask : xr.Dataset
        The dataset containing the land and ocean mask variables.
    region : str
        The region.

    Returns
    -------
    Tuple[str, ...]
        A tuple of valid keys for the land or ocean fraction variable.

    Raises
    ------
    ValueError
        If the region passed is not land or ocean.
    """
    for region_prefix in ["land", "ocean"]:
        if region_prefix in region:
            region_keys = FRAC_REGION_KEYS.get(region_prefix)

    if region_keys is None:
        raise ValueError(f"Only land and ocean regions are supported, not '{region}'.")

    for key in region_keys:
        if key in ds_mask.data_vars:
            return key


def regrid_z_axis_to_plevs(
    dataset: xr.Dataset,
    var_key: str,
    plevs: List[int] | List[float],
) -> xr.Dataset:
    """Regrid a variable's Z axis to the desired pressure levels (mb units).

    The Z axis (e.g., 'lev') must either include hybrid-sigma levels (which
    are converted to pressure coordinates) or pressure coordinates. This is
    determined determined by the "long_name" attribute being set to either
    "hybrid", "isobaric", and "pressure". Afterwards, the pressure coordinates
    are regridded to the specified pressure levels (``plevs``).

    Parameters
    ----------
    dataset : xr.Dataset
        The dataset with the variable on a Z axis.
    var_key : str
        The variable key.
    plevs : List[int] | List[float]
        A 1-D array of floats or integers representing output pressure levels
        in mb units. This parameter is usually set by ``CoreParameter.plevs``
        attribute. For example, ``plevs=[850.0, 200.0]``.

    Returns
    -------
    xr.Dataset
        The dataset with the variables's Z axis regridded to the desired
        pressure levels (mb units).

    Raises
    ------
    KeyError
        If the Z axis has no "long_name" attribute to determine whether it is
        hybrid or pressure.
    ValueError
        If the Z axis "long_name" attribute is not "hybrid", "isobaric",
        or "pressure".

    Notes
    -----
    Replaces `e3sm_diags.driver.utils.general.convert_to_pressure_levels`.
    """
    ds = dataset.copy()

    # Make sure that the input dataset has Z axis bounds, which are required for
    # getting grid positions during vertical regridding.
    try:
        ds.bounds.get_bounds("Z")
    except KeyError:
        ds = ds.bounds.add_bounds("Z")

    z_axis = get_z_axis(ds[var_key])
    z_long_name = z_axis.attrs.get("long_name")
    if z_long_name is None:
        raise KeyError(
            f"The vertical level ({z_axis.name}) for '{var_key}' does "
            "not have a 'long_name' attribute to determine whether it is hybrid "
            "or pressure."
        )
    z_long_name = z_long_name.lower()

    # Hybrid must be the first conditional statement because the long_name attr
    # can be "hybrid sigma pressure coordinate" which includes "hybrid" and
    # "pressure".
    if "hybrid" in z_long_name:
        ds_plevs = _hybrid_to_plevs(ds, var_key, plevs)
    elif "pressure" in z_long_name or "isobaric" in z_long_name:
        ds_plevs = _pressure_to_plevs(ds, var_key, plevs)
    else:
        raise ValueError(
            f"The vertical level ({z_axis.name}) for '{var_key}' is "
            "not hybrid or pressure. Its long name must either include 'hybrid', "
            "'pressure', or 'isobaric'."
        )

    # Add bounds for the new, regridded Z axis if the length is greater than 1.
    # xCDAT does not support adding bounds for singleton coordinates.
    new_z_axis = get_z_axis(ds_plevs[var_key])
    if len(new_z_axis) > 1:
        ds_plevs = ds_plevs.bounds.add_bounds("Z")

    return ds_plevs


def _hybrid_to_plevs(
    ds: xr.Dataset,
    var_key: str,
    plevs: List[int] | List[float],
) -> xr.Dataset:
    """Regrid a variable's hybrid-sigma levels to the desired pressure levels.

    Steps:
        1. Create the output pressure grid using ``plevs``.
        2. Convert hybrid-sigma levels to pressure coordinates.
        3. Regrid the pressure coordinates to the output pressure grid (plevs).

    Parameters
    ----------
    ds : xr.Dataset
        The dataset with the variable using hybrid-sigma levels.
    var_key : var_key.
        The variable key.
    plevs : List[int] | List[float]
        A 1-D array of floats or integers representing output pressure levels
        in mb units. For example, ``plevs=[850.0, 200.0]``. This parameter is
        usually set by the ``CoreParameter.plevs`` attribute.

    Returns
    -------
    xr.Dataset
        The variable with a Z axis regridded to pressure levels (mb units).

    Notes
    -----
    Replaces `e3sm_diags.driver.utils.general.hybrid_to_plevs`.
    """
    # TODO: mb units are always expected, but we should consider checking
    # the units to confirm whether or not unit conversion is needed.
    z_axis, _ = xc.create_axis("lev", plevs, generate_bounds=False)

    pressure_grid = xc.create_grid(z=z_axis)
    pressure_coords = _hybrid_to_pressure(ds, var_key)

    # Keep the "axis" and "coordinate" attributes for CF mapping.
    with xr.set_options(keep_attrs=True):
        result = ds.regridder.vertical(
            var_key,
            output_grid=pressure_grid,
            tool="xgcm",
            method="log",
            target_data=pressure_coords,
        )

    # Vertical regriding sets the units to "mb", but the original units
    # should be preserved.
    result[var_key].attrs["units"] = ds[var_key].attrs["units"]
    result[var_key].attrs["long_name"] = ds[var_key].attrs["long_name"]

    return result


def _hybrid_to_pressure(
    ds: xr.Dataset,
    var_key: str,
    p0: float = 1000.0,
    a_key: Literal["hyam", "hyai"] = "hyam",
    b_key: Literal["hybm", "hybi"] = "hybm",
) -> xr.DataArray:
    """Regrid hybrid-sigma levels to pressure coordinates.

    Formula: p(k) = a(k) * p0 + b(k) * ps
      * p: pressure data (mb).
      * hyam/hyai: 1-D array equal to hybrid A coefficients.
      * p0: Scalar numeric value equal to surface reference pressure with
          the same units as "ps" (mb).
      * hybm/hybi: 1-D array equal to hybrid B coefficients.
      * ps: 2-D array equal to surface pressure data (mb, converted from Pa).

    Parameters
    ----------
    ds : xr.Dataset
        The dataset containing the variable and hybrid levels.
    var_key : str
        The variable key.
    p0 : float
        Scalar numeric value equal to surface reference pressure with
        the same units as "ps", by default 1000.0 (mb).
    a_key : Literal["hyam", "hyai"]
        The key for the hybrid A coefficient variable. By default, "hyam".
    b_key : Literal["hybm", "hybi"]
        The key for the hybrid B coefficient variable. By default, "hybm".

    Returns
    -------
    xr.DataArray
        The variable with a Z axis on pressure coordinates.

    Raises
    ------
    KeyError
        If the dataset does not contain pressure data (ps) or any of the
        hybrid levels (hyam/hyai, hybm/hybi).

    Notes
    -----
    This function is equivalent to `geocat.comp.interp_hybrid_to_pressure()`
    and `cdutil.vertical.reconstructPressureFromHybrid()`.
    """
    # p0 is statically set to mb (1000) instead of retrieved from the dataset
    # because the pressure data should be in mb.
    ps = _get_hybrid_sigma_level(ds, "ps")
    a = _get_hybrid_sigma_level(ds, a_key)
    b = _get_hybrid_sigma_level(ds, b_key)

    if ps is None or a is None or b is None:
        raise KeyError(
            f"The dataset for '{var_key}' does not contain hybrid-sigma level 'ps', "
            f"'{a_key}' and/or '{b_key}' for reconstructing to pressure data."
        )

    if p0 == 1000:
        ps = _convert_dataarray_units_to_mb(ps)

    pressure_coords = a * p0 + b * ps

    if p0 == 1000.0:
        pressure_coords.attrs["units"] = "mb"
    elif p0 == 100000.0:
        pressure_coords.attrs["units"] = "Pa"

    return pressure_coords


def _get_hybrid_sigma_level(
    ds: xr.Dataset, name: Literal["ps", "p0", "hyam", "hyai", "hybm", "hybi"]
) -> xr.DataArray | None:
    """Get the hybrid-sigma level xr.DataArray from the xr.Dataset.

    This function retrieves the valid keys for the specified hybrid-sigma
    level and loops over them. A dictionary look-up is performed and the first
    match is returned. If there are no matches, None is returned.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset.
    name : {"ps", "p0", "hyam", "hyai", "hybm", "hybi"}
        The name of the hybrid-sigma level to get.

    Returns
    -------
    xr.DataArray | None
        The hybrid-sigma level xr.DataArray if found or None.
    """
    keys = HYBRID_SIGMA_KEYS[name]

    for key in keys:
        da = ds.get(key)

        if da is not None:
            return da

    return None


def _pressure_to_plevs(
    ds: xr.Dataset,
    var_key: str,
    plevs: List[int] | List[float],
) -> xr.Dataset:
    """Regrids pressure coordinates to the desired pressure level(s).

    Parameters
    ----------
    ds : xr.Dataset
        The dataset with a variable using pressure data.
    var_key : str
        The variable key.
    plevs : List[int] | List[float]
        A 1-D array of floats or integers representing output pressure levels
        in mb units. This parameter is usually set by ``CoreParameter.plevs``
        attribute. For example, ``plevs=[850.0, 200.0]``.

    Returns
    -------
    xr.Dataset
        The variable with a Z axis on pressure levels (mb).

    Notes
    -----
    Replaces `e3sm_diags.driver.utils.general.pressure_to_plevs`.
    """
    # Convert pressure coordinates and bounds to mb if it is not already in mb.
    ds = _convert_dataset_units_to_mb(ds, var_key)

    # Create the output pressure grid to regrid to using the `plevs` array.
    z_axis, _ = xc.create_axis("lev", plevs, generate_bounds=False)
    pressure_grid = xc.create_grid(z=z_axis)

    # Keep the "axis" and "coordinate" attributes for CF mapping.
    with xr.set_options(keep_attrs=True):
        result = ds.regridder.vertical(
            var_key,
            output_grid=pressure_grid,
            tool="xgcm",
            method="log",
        )

    return result


def _convert_dataset_units_to_mb(ds: xr.Dataset, var_key: str) -> xr.Dataset:
    """Convert a dataset's Z axis and bounds to mb if they are not in mb.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset.
    var_key : str
        The key of the variable.

    Returns
    -------
    xr.Dataset
        The dataset with a Z axis in mb units.

    Raises
    ------
    RuntimeError
        If the Z axis units does not align with the Z bounds units.
    """
    z_axis = xc.get_dim_coords(ds[var_key], axis="Z")
    z_bnds = ds.bounds.get_bounds(axis="Z", var_key=var_key)

    # Make sure that Z and Z bounds units are aligned. If units do not exist
    # assume they are the same because bounds usually don't have a units attr.
    z_axis_units = z_axis.attrs["units"]
    z_bnds_units = z_bnds.attrs.get("units")
    if z_bnds_units is not None and z_bnds_units != z_axis_units:
        raise RuntimeError(
            f"The units for '{z_bnds.name}' ({z_bnds_units}) "
            f"does not align with '{z_axis.name}' ({z_axis_units}). "
        )
    else:
        z_bnds.attrs["units"] = z_axis_units

    # Convert Z and Z bounds and update them in the Dataset.
    z_axis_new = _convert_dataarray_units_to_mb(z_axis)
    ds = ds.assign_coords({z_axis.name: z_axis_new})

    z_bnds_new = _convert_dataarray_units_to_mb(z_bnds)
    z_bnds_new[z_axis.name] = z_axis_new
    ds[z_bnds.name] = z_bnds_new

    return ds


def _convert_dataarray_units_to_mb(da: xr.DataArray) -> xr.DataArray:
    """Convert a dataarray to mb (millibars) if they are not in mb.

    Unit conversion formulas:
      * hPa = mb
      * mb = Pa / 100
      * Pa = (mb * 100)

    The more common unit on weather maps is mb.

    Parameters
    ----------
    da : xr.DataArray
        An xr.DataArray, usually for "lev" or "ps".

    Returns
    -------
    xr.DataArray
        The DataArray in mb units.

    Raises
    ------
    ValueError
        If ``da`` DataArray has no 'units' attribute.
    ValueError
        If ``da`` DataArray has units not in mb or Pa.
    """
    units = da.attrs.get("units")

    if units is None:
        raise ValueError(
            f"'{da.name}' has no 'units' attribute to determine if data is in'mb', "
            "'hPa', or 'Pa' units."
        )

    if units == "Pa":
        with xr.set_options(keep_attrs=True):
            da = da / 100.0

        da.attrs["units"] = "mb"
    elif units == "hPa":
        da.attrs["units"] = "mb"
    elif units == "mb":
        pass
    else:
        raise ValueError(
            f"'{da.name}' should be in 'mb' or 'Pa' (which gets converted to 'mb'), "
            f"not '{units}'."
        )

    return da

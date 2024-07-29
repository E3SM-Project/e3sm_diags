"""
Issue - The slice is excluding points for the ref file.
"""
# %%
import xarray as xr
import xcdat as xc

from e3sm_diags.derivations.default_regions_xr import REGION_SPECS

REGION = "5S5N"
region_slice = (-5.0, 5.0)


test_file = "/global/cfs/cdirs/e3sm/e3sm_diags/postprocessed_e3sm_v2_data_for_e3sm_diags/20210528.v2rc3e.piControl.ne30pg2_EC30to60E2r2.chrysalis/time-series/rgr/U_005101_006012.nc"
ref_file = "/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/time-series/ERA-Interim/ua_197901_201612.nc"


def _subset_on_region(ds: xr.Dataset) -> xr.Dataset:
    """Subset the dataset by the region 5S5N (latitude).

    This function takes into account the CDAT subset flag, "ccb", which can
    add new latitude coordinate points to the beginning and end.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset.

    Returns
    -------
    xr.Dataset
        The dataset subsetted by the region.
    """
    lat_slice = REGION_SPECS[REGION]["lat"]  # type: ignore

    ds_new = ds.copy()
    dim_key = xc.get_dim_keys(ds, "Y")

    # 1. Subset on the region_slice
    # slice = -5.0, 5.0
    ds_new = ds_new.sel({dim_key: slice(*lat_slice)})

    # 2. Add delta to first and last value
    dim_bounds = ds_new.bounds.get_bounds(axis="Y")
    # delta = 1.0 / 2 = 0.5
    delta = (dim_bounds[0][1].item() - dim_bounds[0][0].item()) / 2
    delta_slice = (lat_slice[0] - delta, lat_slice[1] + delta)

    # 3. Check if latitude slice value exists in original latitude.
    # If it exists already, then don't add the coordinate point.
    # If it does not exist, add the coordinate point.
    # delta = 0.5
    # delta slice = -5.5, 5.5
    ds_list = [ds_new]

    try:
        ds.sel({dim_key: delta_slice[0]})
    except KeyError:
        ds_first_pt = ds_new.isel({dim_key: 0})
        ds_first_pt[dim_key] = ds_first_pt[dim_key] - delta

        ds_list.append(ds_first_pt)

    try:
        ds.sel({dim_key: delta_slice[-1]})
    except KeyError:
        ds_last_pt = ds_new.isel({dim_key: -1})
        ds_last_pt[dim_key] = ds_last_pt[dim_key] + delta

        ds_list.append(ds_last_pt)

    ds_new = xr.concat(ds_list, dim=dim_key, data_vars="minimal", coords="minimal")
    ds_new.drop_vars(dim_bounds)

    return ds_new


ds_test = xc.open_dataset(test_file)
ds_ref = xc.open_dataset(ref_file)

# %%
"""
"ccb" flag is adding the bounds delta / 2 to the end and beginning coordinates.

CDAT Expected: 10
    array([-4.5, -3.5, -2.5, -1.5, -0.5,  0.5,  1.5,  2.5,  3.5,  4.5])
xCDAT Result: 10
    array([-4.5, -3.5, -2.5, -1.5, -0.5,  0.5,  1.5,  2.5,  3.5,  4.5])
"""
ds_test_reg = _subset_on_region(ds_test)
ds_test_reg.lat

# %%
"""
"ccb" flag is adding the bounds delta / 2 to the end and beginning coordinates.

CDAT Expected: 15
    array([-4.9375, -4.5   , -3.75  , -3.    , -2.25  , -1.5   , -0.75  ,
        0.    ,  0.75  ,  1.5   ,  2.25  ,  3.    ,  3.75  ,  4.5   ,
        4.9375])
xCDAT Result: 15
    array([-4.875, -4.5  , -3.75 , -3.   , -2.25 , -1.5  , -0.75 ,  0.   ,  0.75 ,
            1.5  ,  2.25 ,  3.   ,  3.75 ,  4.5  ,  4.875])
"""
ds_ref_reg = _subset_on_region(ds_ref)
ds_ref_reg.lat

# %%

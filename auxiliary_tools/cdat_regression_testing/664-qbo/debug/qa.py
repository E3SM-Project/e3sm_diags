"""
Issue - The slice is excluding points for the ref file.
"""
# %%
from typing import Literal

import xarray as xr
import xcdat as xc

region_slice = (-5.0, 5.0)


test_file = "/global/cfs/cdirs/e3sm/e3sm_diags/postprocessed_e3sm_v2_data_for_e3sm_diags/20210528.v2rc3e.piControl.ne30pg2_EC30to60E2r2.chrysalis/time-series/rgr/U_005101_006012.nc"
ref_file = "/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/time-series/ERA-Interim/ua_197901_201612.nc"


def _add_delta_to_dim_slice(
    ds: xr.Dataset, axis: Literal["X", "Y"], region_slice
) -> slice:
    ds_new = ds.copy()
    dim_key = xc.get_dim_keys(ds, axis)

    # 1. Subset on the region_slice
    # slice = -5.0, 5.0
    ds_new = ds_new.sel({dim_key: slice(*region_slice)})

    # 2. Add delta to first and last value
    dim_bounds = ds_new.bounds.get_bounds(axis=axis)
    # delta = 1.0 / 2 = 0.5
    delta = (dim_bounds[0][1].item() - dim_bounds[0][0].item()) / 2
    delta_slice = (region_slice[0] - delta, region_slice[1] + delta)

    # 3. Check if latitude slice value exists in original latitude
    # delta = 0.5
    # delta slice = -5.5, 5.5
    try:
        ds.sel({dim_key: delta_slice[0]})
    # 4. If not, add the latitude value
    except KeyError:
        ds_first_pt = ds_new.isel({dim_key: 0})
        ds_first_pt[dim_key] = ds_first_pt[dim_key] - delta
        ds_new = xr.concat([ds_first_pt, ds_new], dim=dim_key)  # type: ignore

    try:
        ds.sel({dim_key: delta_slice[-1]})
    except KeyError:
        ds_last_pt = ds_new.isel({dim_key: -1})
        ds_last_pt[dim_key] = ds_last_pt[dim_key] + delta
        ds_new = xr.concat([ds_new, ds_last_pt], dim=dim_key)  # type: ignore

    return ds_new


ds

ds_test = xc.open_dataset(test_file)
ds_ref = xc.open_dataset(ref_file)

# %%
"""
"ccb" flag is adding the bounds delta / 2 to the end and beginning coordinates.

CDAT Expected: 15
    array([-4.875, -4.5  , -3.75 , -3.   , -2.25 , -1.5  , -0.75 ,  0.   ,  0.75 ,
            1.5  ,  2.25 ,  3.   ,  3.75 ,  4.5  ,  4.875])
xCDAT Result: 15
    array([-4.5, -3.5, -2.5, -1.5, -0.5,  0.5,  1.5,  2.5,  3.5,  4.5])
"""
ds_test_reg = _add_delta_to_dim_slice(ds_test, "Y", region_slice)
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
ds_ref_reg = _add_delta_to_dim_slice(ds_ref, "Y", region_slice)
ds_ref_reg.lat
# %%


# -5.0, 5.0
# ----------
# %%
# Expected: 10, Result 10
ds_test.sel(lat=slice(*region_slice)).lat
# %%
"""
Expected: 15
    array([-4.9375, -4.5   , -3.75  , -3.    , -2.25  , -1.5   , -0.75  ,
        0.    ,  0.75  ,  1.5   ,  2.25  ,  3.    ,  3.75  ,  4.5   ,
        4.9375])
Result: 13
    array([-4.5 , -3.75, -3.  , -2.25, -1.5 , -0.75,  0.  ,  0.75,  1.5 ,  2.25,
        3.  ,  3.75,  4.5 ])
"""
ds_ref.sel(lat=slice(*region_slice)).lat
# %%
"""
Expected: 10
    array([-4.5, -3.5, -2.5, -1.5, -0.5,  0.5,  1.5,  2.5,  3.5,  4.5])
Result: 12
    array([-5.5, -4.5, -3.5, -2.5, -1.5, -0.5,  0.5,  1.5,  2.5,  3.5,  4.5,  5.5])
"""

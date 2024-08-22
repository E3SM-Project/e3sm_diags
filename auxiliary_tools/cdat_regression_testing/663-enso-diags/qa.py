def _subset_on_region(ds: xr.Dataset, region: str) -> xr.Dataset:
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
    specs = REGION_SPECS[region]
    lat_slice, lon_slice = specs.get("lat"), specs.get("lon")  # type: ignore

    ds_new = ds.copy()
    lat_dim_key = xc.get_dim_keys(ds, "Y")
    lon_dim_key = xc.get_dim_keys(ds, "X")

    # 1. Subset on the region_slice
    # slice = -5.0, 5.0
    slice_dict = {}
    if lat_slice is not None:
        slice_dict[lat_dim_key] = slice(*lat_slice)
    if lon_slice is not None:
        slice_dict[lon_dim_key] = slice(*lon_slice)

    ds_new = ds_new.sel(slice_dict)

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
        ds.sel({lat_dim_key: delta_slice[0]})
    except KeyError:
        ds_first_pt = ds_new.isel({lat_dim_key: 0})
        ds_first_pt[lat_dim_key] = ds_first_pt[lat_dim_key] - delta

        ds_list.append(ds_first_pt)

    try:
        ds.sel({lat_dim_key: delta_slice[-1]})
    except KeyError:
        ds_last_pt = ds_new.isel({lat_dim_key: -1})
        ds_last_pt[lat_dim_key] = ds_last_pt[lat_dim_key] + delta

        ds_list.append(ds_last_pt)

    ds_new = xr.concat(ds_list, dim=lat_dim_key, data_vars="minimal", coords="minimal")
    ds_new = ds_new.sortby(lat_dim_key)

    return ds_new

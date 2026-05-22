import xarray as xr


climo_path ="auxiliary_tools/debug/1048-py314-stall-cont/v3.LR.historical_0051_ANN_198501_201412_climo.nc"

def open_dataset(filepath: str) -> xr.Dataset:
    args = {
        "paths": filepath,
        "decode_times": True,
        "coords": "minimal",
        "compat": "override",
    }

    ds = xr.open_mfdataset(**args)

    return ds


ds = open_dataset(climo_path)
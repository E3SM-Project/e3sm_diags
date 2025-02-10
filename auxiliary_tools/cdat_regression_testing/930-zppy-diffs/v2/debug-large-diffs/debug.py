#%%
import glob
from typing import List

import numpy as np
import xarray as xr

from e3sm_diags.derivations.derivations import DERIVED_VARIABLES
import pandas as pd

DEV_DIR = "25-02-04-branch-930-zppy-diffs"
DEV_PATH = f"/lcrc/group/e3sm/public_html/cdat-migration-fy24/{DEV_DIR}/"

DEV_GLOB_ALL = sorted(glob.glob(DEV_PATH + "**/**/*.nc"))
DEV_NUM_FILES = len(DEV_GLOB_ALL)

MAIN_DIR = "25-02-04-main-zppy-diffs"
MAIN_PATH = f"/lcrc/group/e3sm/public_html/cdat-migration-fy24/{MAIN_DIR}/"
MAIN_GLOB_ALL = sorted(glob.glob(MAIN_PATH + "**/**/*.nc"))
MAIN_NUM_FILES = len(MAIN_GLOB_ALL)

#%%
KEEP_VARS = ["OMI-MLS-TCO-ANN-60S60N"]

DEV_GLOB = [fp for fp in DEV_GLOB_ALL if "diff.nc" not in fp and any(var in fp for var in KEEP_VARS)]
MAIN_GLOB = [fp for fp in MAIN_GLOB_ALL if "diff.nc" not in fp and any(var in fp for var in KEEP_VARS)]

DEV_GLOB_DIFF = [fp for fp in DEV_GLOB_ALL if "diff.nc" in fp and any(var in fp for var in KEEP_VARS)]
MAIN_GLOB_DIFF = [fp for fp in MAIN_GLOB_ALL if "diff.nc" in fp and any(var in fp for var in KEEP_VARS)]

def _get_var_data(ds: xr.Dataset, var_key: str) -> np.ndarray:
    """Get the variable data using a list of matching keys.

    The `main` branch saves the dataset using the original variable name,
    while the dev branch saves the variable with the derived variable name.
    The dev branch is performing the expected behavior here.

    Parameters
    ----------
    ds : xr.Dataset
        _description_
    var_key : str
        _description_

    Returns
    -------
    np.ndarray
        _description_
    """

    data = None

    try:
        data = ds[var_key]
    except KeyError:
        try:
            var_keys = DERIVED_VARIABLES[var_key].keys()
        except KeyError:
            var_keys = DERIVED_VARIABLES[var_key.upper()].keys()

        var_keys = [var_key] + list(sum(var_keys, ()))

        for key in var_keys:
            if key in ds.data_vars.keys():
                data = ds[key]
                break

    return data
#%%
ATOL = 0
RTOL = 1e-5

print(f"Relative tolerance: {RTOL}, Absolute tolerance: {ATOL}")

def compare_files(main_glob: List[str]):
    for fp_main in main_glob:
        var_key = fp_main.split("-")[-3]
        fp_type = fp_main.split("-")[-1].split("_")[-1]

        fp_dev = fp_main.replace(MAIN_DIR, DEV_DIR)

        print(f"{var_key} - {fp_type}")
        print("-" * 50)
        print(f"Main: {fp_main}\nDev: {fp_dev}")

        ds_main = xr.open_dataset(fp_main)
        ds_dev = xr.open_dataset(fp_dev)

        dv_main = _get_var_data(ds_main, var_key)
        dv_dev = _get_var_data(ds_dev, var_key)

        if dv_main is None:
            dv_main = _get_var_data(ds_main, var_key + "_diff")

        try:
            np.testing.assert_allclose(dv_main.values, dv_dev.values, rtol=RTOL, atol=ATOL)
        except AssertionError as e:
            print(f"{e}\n")
        else:
            print(f"Arrays are within relative tolerance.\n")

print("Comparing test.nc and ref.nc files")
print("=" * 50)
compare_files(MAIN_GLOB)

#%%
print("Comparing diff.nc files")
print("=" * 50)
compare_files(MAIN_GLOB_DIFF)

# %%
def compare_stats(main_glob_diff: List[str]):
    for fp_main in main_glob_diff:
        var_key = fp_main.split("-")[-3]
        fp_type = fp_main.split("-")[-1].split("_")[-1]

        fp_dev = fp_main.replace(MAIN_DIR, DEV_DIR)

        print(f"{var_key} - {fp_type}")
        print("-" * 50)
        print(f"Main: {fp_main}\nDev: {fp_dev}")

        ds_main = xr.open_dataset(fp_main)
        ds_dev = xr.open_dataset(fp_dev)

        dv_main = _get_var_data(ds_main, var_key)
        dv_dev = _get_var_data(ds_dev, var_key)

        if dv_main is None:
            dv_main = _get_var_data(ds_main, var_key + "_diff")
            stats_main = {
                "min": dv_main.min().item(),
                "max": dv_main.max().item(),
                "mean": dv_main.mean().item(),
                "sum": dv_main.sum().item(),
                "nan_count": np.isnan(dv_main).sum().item()
            }

            stats_dev = {
                "min": dv_dev.min().item(),
                "max": dv_dev.max().item(),
                "mean": dv_dev.mean().item(),
                "sum": dv_dev.sum().item(),
                "nan_count": np.isnan(dv_dev).sum().item()
            }

            df_stats = pd.DataFrame([stats_main, stats_dev], index=["main", "dev"])
            print(df_stats)

print("Comparing stats of diff.nc files")
print("=" * 50)
compare_stats(MAIN_GLOB_DIFF)


# %%

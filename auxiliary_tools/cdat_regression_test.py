# %%
"""
This script checks for regressions between the refactored and `main` branches
of a diagnostic set.

How it works
------------
  It compares the absolute and relative differences between two sets of
  `.json` files in two separate directories, one for the refactored code
  and the other for the `main` branch. This script will generate an Excel file
  containing:
    1. The raw metrics by each branch for each variable.
    2. The absolute and relative differences of each variable between branches.
    3. The highest relative differences (threshold > 2% difference)

How to use
-----------
    1. mamba env create -f conda/dev-yml -n e3sm_diags_dev_<gh-issue-#>
    2. mamba activate e3sm_diags_dev_<gh-issue-#>
    3. Update `DEV_PATH` and `PROD_PATH` in `/auxiliary_tools/cdat_regression_test.py`
    4. python auxiliary_tools/cdat_regression_test.py
    5. Excel file generated in `/auxiliary_tools`

Tips
-----------
Relative differences should be taken into consideration moreso than absolute
differences.
  - Relative differences show the scale using a percentage unit.
  - Absolute differences is just a raw number that doesn't factor in
    floating point size (e.g., 100.00 vs. 0.0001), which can be misleading.
"""
import glob
import logging
import os
import time
from typing import List

import pandas as pd

log_format = (
    "%(asctime)s [%(levelname)s]: %(filename)s(%(funcName)s:%(lineno)s) >> %(message)s"
)
logging.basicConfig(format=log_format, filemode="w", level=logging.INFO)
logger = logging.getLogger(__name__)

# TODO: Update DEV_RESULTS and PROD_RESULTS.
# ------------------------------------------------------------------------------
DEV_PATH = "/global/cfs/cdirs/e3sm/www/vo13/examples_658/ex1_modTS_vs_modTS_3years/lat_lon/model_vs_model"
PRO_PATH = "/global/cfs/cdirs/e3sm/www/vo13/examples/ex1_modTS_vs_modTS_3years/lat_lon/model_vs_model"
# ------------------------------------------------------------------------------

if not os.path.exists(DEV_PATH):
    raise ValueError(f"DEV_RESULTS path does not exist ({DEV_PATH})")
if not os.path.exists(PRO_PATH):
    raise ValueError(f"PROD_RESULTS path does not exist ({PRO_PATH})")

DEV_GLOB = sorted(glob.glob(DEV_PATH + "/*.json"))
PROD_GLOB = sorted(glob.glob(PRO_PATH + "/*.json"))

TIME_STR = time.strftime("%Y%m%d-%H%M%S")
EXCEL_FILENAME = f"{TIME_STR}-metrics-diffs.xlsx"


def get_metrics(filepaths: List[str]) -> pd.DataFrame:
    """Get the metrics using a glob of `.json` metric files in a directory.

    Parameters
    ----------
    filepaths : List[str]
        The filepaths for metrics `.json` files.

    Returns
    -------
    pd.DataFrame
        The DataFrame containing the metrics for all of the variables in
        the results directory.
    """
    metrics = []

    for filepath in filepaths:
        df = pd.read_json(filepath)

        filename = filepath.split("/")[-1]
        var_key = filename.split("-")[1]

        # Add the variable key to the MultiIndex and update the index
        # before stacking to make the DataFrame easier to parse.
        multiindex = pd.MultiIndex.from_product([[var_key], [*df.index]])
        df = df.set_index(multiindex)
        df.stack()

        metrics.append(df)

    df_final = pd.concat(metrics)

    # Reorder columns and drop "unit" column (string dtype breaks Pandas
    # arithmetic).
    df_final = df_final[["test", "ref", "test_regrid", "ref_regrid", "diff", "misc"]]

    return df_final


def get_diffs(df_a: pd.DataFrame, df_b: pd.DataFrame) -> pd.DataFrame:
    """The metrics differences between two DataFrames.

    Parameters
    ----------
    df_a : pd.DataFrame
        The first DataFrame representing "actual" results (aka development).
    df_b : pd.DataFrame
        The second DataFrame representing "reference" results (aka production).

    Returns
    -------
    pd.DataFrame
        The DataFrame containing absolute and relative differences between
        the metrics DataFrames.
    """
    #  Absolute difference: abs(actual - reference)
    df_abs = abs(df_a - df_b)
    df_abs = df_abs.add_suffix("_abs")

    # Relative difference: abs(actual - reference) / abs(actual)
    df_rel = abs(df_a - df_b) / abs(df_a)
    df_rel = df_rel.add_suffix("_rel")

    # Combine both DataFrames
    df_final = pd.concat([df_abs, df_rel], axis=1, join="outer")

    return df_final


# %% Get the metrics DataFrames.
df_dev = get_metrics(DEV_GLOB)
df_prod = get_metrics(PROD_GLOB)

# %% Combine metrics DataFrames.
df_dev_pref = df_dev.add_prefix("dev_")
df_prod_pref = df_prod.add_prefix("prod_")
df_metrics = pd.concat([df_dev_pref, df_prod_pref], axis=1, join="outer")
#%%
# Sort the columns
df_metrics = df_metrics[
    [
        "dev_test",
        "prod_test",
        "dev_ref",
        "prod_ref",
        "dev_test_regrid",
        "prod_test_regrid",
        "dev_ref_regrid",
        "prod_ref_regrid",
        "dev_diff",
        "prod_diff",
        "dev_misc",
        "prod_misc",
    ]
]

# %% Get differences between metrics.
df_diffs = get_diffs(df_dev, df_prod)


#%%
with pd.ExcelWriter(EXCEL_FILENAME) as writer:
    df_metrics.to_excel(writer, sheet_name="metrics")
    df_diffs.to_excel(writer, sheet_name="metric_diffs")


# %% Only get the metrics where the absolute and relative differences are
# greater than a specific threshold (>1%)

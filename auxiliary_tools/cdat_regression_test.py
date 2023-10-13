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
    1. Make a copy of this script.
    2. Run `mamba create -n cdat_regression_test -y -c conda-forge "python<3.12" pandas matplotlib-base ipykernel`
    3. Run `mamba activate cdat_regression_test`
    4. Update `DEV_PATH` and `PROD_PATH` in the copy of your script.
    5. Run `python auxiliary_tools/<COPY_OF_SCRIPT>.py`
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
import math
import os
import time
from typing import List

import numpy as np
import pandas as pd

log_format = (
    "%(asctime)s [%(levelname)s]: %(filename)s(%(funcName)s:%(lineno)s) >> %(message)s"
)
logging.basicConfig(format=log_format, filemode="w", level=logging.INFO)
logger = logging.getLogger(__name__)

# TODO: Update DEV_RESULTS and PROD_RESULTS.
# ------------------------------------------------------------------------------
DEV_PATH = "/global/cfs/cdirs/e3sm/www/vo13/examples_658/ex1_modTS_vs_modTS_3years/lat_lon/model_vs_model"
PROD_PATH = "/global/cfs/cdirs/e3sm/www/vo13/examples/ex1_modTS_vs_modTS_3years/lat_lon/model_vs_model"
# ------------------------------------------------------------------------------

if not os.path.exists(DEV_PATH):
    raise ValueError(f"DEV_RESULTS path does not exist ({DEV_PATH})")
if not os.path.exists(PROD_PATH):
    raise ValueError(f"PROD_RESULTS path does not exist ({PROD_PATH})")

DEV_GLOB = sorted(glob.glob(DEV_PATH + "/*.json"))
PROD_GLOB = sorted(glob.glob(PROD_PATH + "/*.json"))

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


def _get_rel_diffs(df_a: pd.DataFrame, df_b: pd.DataFrame) -> pd.DataFrame:
    """Get the relative differences between two DataFrames.

    Formula: abs(actual - reference) / abs(actual)

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
    df_diff = abs(df_a - df_b) / abs(df_a)
    df_diff = df_diff.add_suffix(" DIFF (%)")

    return df_diff


# %%
def _sort_columns(df: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "test_dev",
        "test_prod",
        "test DIFF (%)",
        "ref_dev",
        "ref_prod",
        "ref DIFF (%)",
        "test_regrid_dev",
        "test_regrid_prod",
        "test_regrid DIFF (%)",
        "ref_regrid_dev",
        "ref_regrid_prod",
        "ref_regrid DIFF (%)",
        "diff_dev",
        "diff_prod",
        "diff DIFF (%)",
        "misc_dev",
        "misc_prod",
        "misc DIFF (%)",
    ]

    df_new = df[columns]

    return df_new


# %% Get metrics DataFrames and relative differences between both.
df_metrics_dev = get_metrics(DEV_GLOB)
df_metrics_prod = get_metrics(PROD_GLOB)
df_metrics_all = pd.concat(
    [df_metrics_dev.add_suffix("_dev"), df_metrics_prod.add_suffix("_prod")],
    axis=1,
    join="outer",
)

df_metrics_diffs = _get_rel_diffs(df_metrics_dev, df_metrics_prod)

# %%
# Create the final DataFrame by left joining the metrics diffs (>=2%) with
# the related metrics values. NaN values mean it is less than the threshold.
# If all cell in a row are NaN, the entire row is dropped to make the results
# easier to parse.
df_metrics_diffs_thres = df_metrics_diffs[df_metrics_diffs >= 0.02]
df_metrics_diffs_thres = df_metrics_diffs_thres.dropna(
    axis=0, how="all", ignore_index=False
)
df_final = df_metrics_diffs_thres.join(df_metrics_all)
df_final = _sort_columns(df_final)

# Update difference values from float to string percentage.
PERCENTAGE_COLUMNS = [
    "test DIFF (%)",
    "ref DIFF (%)",
    "test_regrid DIFF (%)",
    "ref_regrid DIFF (%)",
    "diff DIFF (%)",
    "misc DIFF (%)",
]
df_final[PERCENTAGE_COLUMNS] = df_final[PERCENTAGE_COLUMNS].map(
    lambda x: "{0:.2f}%".format(x * 100) if not math.isnan(x) else x
)

# %%
# Display the final result.
df_final.reset_index(names=["var_key", "metric"]).style.applymap(
    lambda x: "background-color : red" if isinstance(x, str) else "",
    subset=pd.IndexSlice[:, PERCENTAGE_COLUMNS],
)

# %% Write to an Excel file.
# with pd.ExcelWriter(EXCEL_FILENAME) as writer:
#     df_metrics.to_excel(writer, sheet_name="metrics")
#     df_diffs.to_excel(writer, sheet_name="metric_diffs")


# if __name__ == "__main__":
#     main()

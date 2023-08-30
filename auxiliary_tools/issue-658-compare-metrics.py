"""
This script compares the absolute and relative differences between metrics
generates by the `ex1` script on the `refactor/658-lat-lon-set` and `main`
branch.

Relative differences should be taken into consideration moreso than absolute
differences.
  - Relative differences show the scale using a percentage unit.
  - Absolute differences is just a raw number that doesn't factor in
    floating point size (e.g., 100.00 vs. 0.0001), which can be misleading.
"""
import glob
from typing import List

import pandas as pd

DEV_RESULTS = "/global/cfs/cdirs/e3sm/www/vo13/examples/ex1_modTS_vs_modTS_3years/lat_lon/model_vs_model"
PROD_RESULTS = "/global/cfs/cdirs/e3sm/www/forsyth/examples/ex1_modTS_vs_modTS_3years/lat_lon/model_vs_model"

DEV_GLOB = sorted(glob.glob(DEV_RESULTS + "/*.json"))
PROD_GLOB = sorted(glob.glob(PROD_RESULTS + "/*.json"))


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
    # NOTE: Drop the unit attributes since it breaks Pandas arithmetic.
    df_a = df_a.drop(columns=["unit"])
    df_b = df_b.drop(columns=["unit"])

    #  Absolute difference: abs(actual - reference)
    df_abs = abs(df_a - df_b)
    df_abs = df_abs.add_suffix("_abs")

    # Relative difference: abs(actual - reference) / abs(actual)
    df_rel = abs(df_a - df_b) / abs(df_a)
    df_rel = df_rel.add_suffix("_rel")

    # Combine both DataFrames
    df_final = pd.concat([df_abs, df_rel], axis=1, join="outer")
    df_final = _sort_cols(df_final)

    return df_final


def _sort_cols(df: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "test_abs",
        "test_rel",
        "test_regrid_abs",
        "test_regrid_rel",
        "ref_abs",
        "ref_rel",
        "ref_regrid_abs",
        "ref_regrid_rel",
        "misc_abs",
        "misc_rel",
        "diff_abs",
        "diff_rel",
    ]
    df = df[columns]

    return df


# %%
df_dev = get_metrics(DEV_GLOB)
df_prod = get_metrics(PROD_GLOB)

df_diff = get_diffs(df_dev, df_prod)
df_diff.to_excel("20230830-issue-658-metrics-diff.xlsx")

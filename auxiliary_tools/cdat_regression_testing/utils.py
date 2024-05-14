import math
from typing import List

import pandas as pd
from IPython.display import display
from PIL import Image, ImageChops, ImageDraw

# The names of the columns that store percentage difference values.
PERCENTAGE_COLUMNS = [
    "test DIFF (%)",
    "ref DIFF (%)",
    "test_regrid DIFF (%)",
    "ref_regrid DIFF (%)",
    "misc DIFF (%)",
]


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


def get_rel_diffs(df_actual: pd.DataFrame, df_reference: pd.DataFrame) -> pd.DataFrame:
    """Get the relative differences between two DataFrames.

    Formula: abs(actual - reference) / abs(actual)

    Parameters
    ----------
    df_actual : pd.DataFrame
        The first DataFrame representing "actual" results (dev branch).
    df_reference : pd.DataFrame
        The second DataFrame representing "reference" results (main branch).

    Returns
    -------
    pd.DataFrame
        The DataFrame containing absolute and relative differences between
        the metrics DataFrames.
    """
    df_diff = abs(df_actual - df_reference) / abs(df_actual)
    df_diff = df_diff.add_suffix(" DIFF (%)")

    return df_diff


def sort_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Sorts the order of the columns for the final DataFrame output.

    Parameters
    ----------
    df : pd.DataFrame
        The final DataFrame output.

    Returns
    -------
    pd.DataFrame
        The final DataFrame output with sorted columns.
    """
    columns = [
        "test_dev",
        "test_main",
        "test DIFF (%)",
        "ref_dev",
        "ref_main",
        "ref DIFF (%)",
        "test_regrid_dev",
        "test_regrid_main",
        "test_regrid DIFF (%)",
        "ref_regrid_dev",
        "ref_regrid_main",
        "ref_regrid DIFF (%)",
        "misc_dev",
        "misc_main",
        "misc DIFF (%)",
    ]

    df_new = df.copy()
    df_new = df_new[columns]

    return df_new


def update_diffs_to_pct(df: pd.DataFrame, cols: List[str] = PERCENTAGE_COLUMNS):
    """Update relative diff columns from float to string percentage.

    Parameters
    ----------
    df : pd.DataFrame
        The final DataFrame containing metrics and diffs (floats).

    Returns
    -------
    pd.DataFrame
        The final DataFrame containing metrics and diffs (str percentage).
    """
    df_new = df.copy()
    df_new[cols] = df_new[cols].map(
        lambda x: "{0:.2f}%".format(x * 100) if not math.isnan(x) else x
    )

    return df_new


def highlight_large_diffs(df: pd.DataFrame, cols: List[str] = PERCENTAGE_COLUMNS):
    if "var_key" not in df.columns and "metric" not in df.columns:
        df_new = df.reset_index(names=["var_key", "metric"])
    else:
        df_new = df.copy()

    df_new = df_new.style.map(
        lambda x: "background-color : red" if isinstance(x, str) else "",
        subset=pd.IndexSlice[:, cols],
    )

    display(df_new)


def get_num_metrics_above_diff_thres(
    df_metrics: pd.DataFrame, df_metric_above_thres: pd.DataFrame
):
    var_keys = list(df_metric_above_thres.var_key.unique())
    print(f"* Related variables {var_keys}")

    num_rows = df_metrics.shape[0]
    num_rows_largest_diffs = df_metric_above_thres.shape[0]
    print(
        f"* Number of metrics above 2% max threshold: {num_rows_largest_diffs} / {num_rows}"
    )


def get_image_diffs(actual_path: str, expected_path: str):
    """Get the diffs between two images.

    This function is useful for comparing two datasets that can't be compared
    directly using `np.testing.assert_allclose()` due to `x and y nan location
    mismatch` error. This error might happen after using the land-sea mask
    after regridding, which can differ slightly between xCDAT/xESMF and
    CDAT/ESMF.

    Parameters
    ----------
    actual_path : str
        The path to the actual png (e.g., Xarray/xCDAT).
    expected_path : str
        The path to the expected png (e.g., CDAT).
    """
    actual_png = Image.open(actual_path).convert("RGB")
    expected_png = Image.open(expected_path).convert("RGB")

    diff = ImageChops.difference(actual_png, expected_png)

    draw = ImageDraw.Draw(diff)
    (left, upper, right, lower) = diff.getbbox()
    draw.rectangle(((left, upper), (right, lower)), outline="red")

    diff_path = actual_path.replace("actual", "diff")
    diff.save(diff_path)

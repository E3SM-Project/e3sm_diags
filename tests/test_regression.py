import glob

import numpy as np
import pytest
import xarray as xr

from e3sm_diags.derivations.derivations import DERIVED_VARIABLES
from e3sm_diags.logger import custom_logger
from tests.complete_run_script import params, runner


logger = custom_logger(__name__)


DEV_DIR = "843-migration-phase3-model-vs-obs"
DEV_PATH = f"/global/cfs/cdirs/e3sm/www/cdat-migration-fy24/{DEV_DIR}/"

DEV_GLOB = sorted(glob.glob(DEV_PATH + "**/**/*.nc"))
DEV_NUM_FILES = len(DEV_GLOB)

MAIN_DIR = "main"
MAIN_PATH = f"/global/cfs/cdirs/e3sm/www/cdat-migration-fy24/{MAIN_DIR}/"
MAIN_GLOB = sorted(glob.glob(MAIN_PATH + "**/**/*.nc"))
MAIN_NUM_FILES = len(MAIN_GLOB)

# Absolute and relative tolerance levels for comparison of the data.
# Absolute is in floating point terms, relative is in percentage terms.
ATOL = 0
RTOL = 1e-5


@pytest.fixture(scope="module")
def run_diags_and_get_results_dir() -> str:
    """Run the diagnostics and get the results directory containing the images.

    The scope of this fixture is at the module level so that it only runs
    once, then each individual test can reference the result directory.

    Returns
    -------
    str
        The path to the results directory.
    """
    results = runner.run_diags(params)

    if results is not None:
        results_dir = results[0].results_dir
    else:
        results_dir = params[0].results_dir

    return results_dir


class TestRegression:
    @pytest.fixture(autouse=True)
    def setup(self, run_diags_and_get_results_dir):
        # TODO: We need to store `main` results on a data container
        self.results_dir = run_diags_and_get_results_dir

    def test_check_if_files_found(self):
        if DEV_NUM_FILES == 0 or MAIN_NUM_FILES == 0:
            raise IOError(
                "No files found at DEV_PATH and/or MAIN_PATH. "
                f"Please check {DEV_PATH} and {MAIN_PATH}."
            )

    def test_check_if_matching_filecount(self):
        if DEV_NUM_FILES != MAIN_NUM_FILES:
            raise IOError(
                "Number of files do not match at DEV_PATH and MAIN_PATH "
                f"({DEV_NUM_FILES} vs. {MAIN_NUM_FILES})."
            )

        logger.info(f"Matching file count ({DEV_NUM_FILES} and {MAIN_NUM_FILES}).")

    def test_check_if_missing_files(self):
        missing_dev_files, missing_main_files = _check_if_missing_files()

        assert len(missing_dev_files) == 0
        assert len(missing_main_files) == 0

    def test_get_relative_diffs(self):
        results = _get_relative_diffs()

        assert len(results["missing_files"]) == 0
        assert len(results["missing_vars"]) == 0
        assert len(results["matching_files"]) > 0
        assert len(results["mismatch_errors"]) == 0
        assert len(results["not_equal_errors"]) == 0
        assert len(results["key_errors"]) == 0


def _get_relative_diffs():
    results = {
        "missing_files": [],
        "missing_vars": [],
        "matching_files": [],
        "mismatch_errors": [],
        "not_equal_errors": [],
        "key_errors": [],
    }

    for fp_main in MAIN_GLOB:
        fp_dev = fp_main.replace(MAIN_DIR, DEV_DIR)

        logger.info("Comparing:")
        logger.info(f"    * {fp_dev}")
        logger.info(f"    * {fp_main}")

        try:
            ds1 = xr.open_dataset(fp_dev)
            ds2 = xr.open_dataset(fp_main)
        except FileNotFoundError as e:
            logger.info(f"    {e}")

            if isinstance(e, FileNotFoundError) or isinstance(e, OSError):
                results["missing_files"].append(fp_dev)

            continue

        var_key = fp_main.split("-")[-3]

        # for 3d vars such as T-200
        var_key.isdigit()
        if var_key.isdigit():
            var_key = fp_main.split("-")[-4]

        dev_data = _get_var_data(ds1, var_key)
        main_data = _get_var_data(ds2, var_key)

        logger.info(f"    * var_key: {var_key}")

        if dev_data is None or main_data is None:
            if dev_data is None:
                results["missing_vars"].append(fp_dev)
            elif main_data is None:
                results["missing_vars"].append(fp_main)

            logger.error("    * Could not find variable key in the dataset(s)")

            continue

        try:
            np.testing.assert_allclose(
                dev_data,
                main_data,
                atol=ATOL,
                rtol=RTOL,
            )
            results["matching_files"].append(fp_main)
        except (KeyError, AssertionError) as e:
            msg = str(e)

            logger.info(f"    {msg}")

            if "mismatch" in msg:
                results["mismatch_errors"].append(fp_dev)
            elif "Not equal to tolerance" in msg:
                results["not_equal_errors"].append(fp_dev)
        else:
            logger.info(f"    * All close and within relative tolerance ({RTOL})")

    return results


def _get_var_data(ds: xr.Dataset, var_key: str) -> np.ndarray | None:
    """Retrieve variable data from an xarray Dataset.

    Parameters
    ----------
    ds : xr.Dataset
        The xarray Dataset from which to retrieve the variable data.
    var_key : str
        The key of the variable to retrieve.

    Returns
    -------
    np.ndarray
        The data of the specified variable as a NumPy array. If the variable is
        not found, returns None.

    Raises
    ------
    KeyError
        If the variable key is not found in the Dataset and is not a derived
        variable.
    """
    data = None

    try:
        var_keys = DERIVED_VARIABLES[var_key].keys()
    except KeyError:
        var_keys = DERIVED_VARIABLES[var_key.upper()].keys()

    var_keys = [var_key] + list(sum(var_keys, ()))

    for key in var_keys:
        if key in ds.data_vars.keys():
            data = ds[key].values

            break

    return data


def _check_if_missing_files():
    missing_dev_files = []
    missing_main_files = []

    for fp_main in MAIN_GLOB:
        fp_dev = fp_main.replace(MAIN_DIR, DEV_DIR)

        if fp_dev not in DEV_GLOB:
            missing_dev_files.append(fp_dev)

    for fp_dev in DEV_GLOB:
        fp_main = fp_dev.replace(DEV_DIR, MAIN_DIR)

        if fp_main not in MAIN_GLOB:
            missing_main_files.append(fp_main)

    return missing_dev_files, missing_main_files

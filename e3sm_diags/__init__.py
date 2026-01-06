import os
import sys

__version__ = "v3.1.1"
INSTALL_PATH = os.path.join(sys.prefix, "share/e3sm_diags/")

# Needed for when using hdf5 >= 1.10.0,
# without this, errors are thrown on Edison compute nodes.
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
# Used by numpy, causes too many threads to spawn otherwise.
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

# Settings to preserve legacy Xarray behavior when merging datasets.
# See https://xarray.pydata.org/en/stable/user-guide/io.html#combining-multiple-files
# and https://xarray.pydata.org/en/stable/whats-new.html#id14
LEGACY_XARRAY_MERGE_KWARGS = {
    # "override", "exact" are the new defaults as of Xarray v2025.08.0
    "compat": "no_conflicts",
    "join": "outer",
}

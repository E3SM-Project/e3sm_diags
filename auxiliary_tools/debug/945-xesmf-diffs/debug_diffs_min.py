# %%
"""
This script demonstrates the difference in regridding results when using
datasets with and without pre-existing lat/lon bounds.

Two datasets are used:
1. `ds_a`: Contains pre-existing lat/lon bounds.
2. `ds_b`: Does not have sorted lat/lon bounds initially.

The script highlights that regridded arrays differ even when the input arrays
are the same, due to differences in how bounds are handled during regridding.

Key Steps:
1. Load datasets.
2. Regrid with pre-existing bounds.
3. Regrid without pre-existing bounds (by adding missing bounds).
4. Compare statistics and floating-point differences.
"""

# %%
import xcdat as xc
import numpy as np
import xesmf as xe


def print_stats(*arrays, labels=None):
    """Prints statistical comparison of multiple arrays."""
    if labels is None:
        labels = [f"Array {i + 1}" for i in range(len(arrays))]
    elif len(labels) != len(arrays):
        raise ValueError("Number of labels must match the number of arrays.")

    stats = {
        "Min": [np.min(arr) for arr in arrays],
        "Max": [np.max(arr) for arr in arrays],
        "Mean": [np.mean(arr) for arr in arrays],
        "Sum": [np.sum(arr) for arr in arrays],
        "Std": [np.std(arr) for arr in arrays],
    }

    header = f"\n{'Stat':<10} " + " ".join([f"{label:<20}" for label in labels])
    print(header)
    print("-" * (10 + 20 * len(labels)))
    for stat, values in stats.items():
        row = f"{stat:<10} " + " ".join([f"{value:<20.6f}" for value in values])
        print(row)


# %%
# Load datasets
ds_a = xc.open_dataset(
    "/lcrc/group/e3sm/public_html/cdat-migration-fy24/25-02-14-branch-930-polar-after/polar/GPCP_v3.2/test.nc"
)
ds_b = xc.open_dataset(
    "/lcrc/group/e3sm/public_html/cdat-migration-fy24/25-02-14-branch-930-polar-after/polar/GPCP_v3.2/ref.nc"
)

# %%
# 1. Regrid with unsorted, descending lat_bnds in reference dataset
# ----------------------------------------------------
# Regrid with pre-existing lat/lon bounds
# NOTE: Need to copy dataset becase it seems to affect xCDAT's regridder to xESMF.
ds_b_cp = ds_b.copy(deep=True)
output_grid1 = ds_a.regridder.grid
regridder1 = xe.Regridder(ds_b_cp.copy(), output_grid1, method="conservative_normed")

# Perform regridding and get the weights
prect1 = regridder1(ds_b_cp["PRECT"])
weights1 = regridder1.weights

# Regrid using xcdat's horizontal regrid method
prect_xcdat1 = ds_b.copy().regridder.horizontal(
    "PRECT", output_grid1, tool="xesmf", method="conservative_normed"
)["PRECT"]


# %%
# 2. Regrid with sorted, ascending lat_bnds in reference dataset
# ----------------------------------------------------
# Fix unsorted lat_bnds in ds_b
ds_b2 = ds_b.copy(deep=True)
ds_b2["lat_bnds"].values = np.sort(ds_b2["lat_bnds"], axis=-1)

# Regrid with pre-existing lat/lon bounds
output_grid2 = ds_a.regridder.grid
regridder2 = xe.Regridder(ds_b2, output_grid2, method="conservative_normed")

# Perform regridding and get the weights
prect2 = regridder2(ds_b2["PRECT"])
weights2 = regridder2.weights

# Regrid using xcdat's horizontal regrid method
prect_xcdat2 = ds_b2.regridder.horizontal(
    "PRECT", output_grid2, method="conservative_normed"
)["PRECT"]

# %%
# 3. Regrid with generated generated lat_bnds in reference dataset
# ----------------------------------------------------
ds_a3 = ds_a.drop_vars(["lat_bnds"]).bounds.add_missing_bounds(axes=["Y"])
ds_b3 = ds_b.drop_vars(["lat_bnds"]).bounds.add_missing_bounds(axes=["Y"])

output_grid3 = ds_a3.regridder.grid

# Create the regridder using xesmf directly
regridder3 = xe.Regridder(ds_b3, output_grid3, method="conservative_normed")

# Perform regridding and get the weights
prect3 = regridder3(ds_b3["PRECT"])
weights3 = regridder3.weights

# Regrid using xcdat's horizontal regrid method
prect_xcdat3 = ds_b3.regridder.horizontal(
    "PRECT", output_grid3, method="conservative_normed"
)["PRECT"]


# %%
# Convert sparse weights to dense arrays for comparison
weights1_dense = weights1.data.todense()
weights2_dense = weights2.data.todense()
weights3_dense = weights3.data.todense()

# Compare dense weights arrays for floating-point differences
print_stats(
    weights1_dense,
    weights2_dense,
    weights3_dense,
    labels=[
        "Weights (unsorted bounds)",
        "Weights (sorted bounds)",
        "Weights (generated bounds)",
    ],
)

# %%

# Compare dense weights arrays for floating-point differences
try:
    np.testing.assert_allclose(weights1_dense, weights2_dense, rtol=1e-5, atol=0)
    print(
        "\nDense Weights (unsorted bounds) and Dense Weights (sorted bounds) are within rtol=1e-5"
    )
except AssertionError as e:
    print(
        "\nDense Weights (unsorted bounds) and Dense Weights (sorted bounds) are not within rtol=1e-5"
    )
    print(e)

try:
    np.testing.assert_allclose(weights2_dense, weights3_dense, rtol=1e-5, atol=0)
    print(
        "\nDense Weights (sorted bounds) and Dense Weights (generated bounds) are within rtol=1e-5"
    )
except AssertionError as e:
    print(
        "\nDense Weights (sorted bounds) and Dense Weights (generated bounds) are not within rtol=1e-5"
    )
    print(e)

try:
    np.testing.assert_allclose(weights1_dense, weights3_dense, rtol=1e-5, atol=0)
    print(
        "\nDense Weights (unsorted bounds) and Dense Weights (generated bounds) are within rtol=1e-5"
    )
except AssertionError as e:
    print(
        "\nDense Weights (unsorted bounds) and Dense Weights (generated bounds) are not within rtol=1e-5"
    )
    print(e)

# %%
print_stats(
    prect1.values,
    prect2.values,
    prect3.values,
    labels=[
        "PRECT (unsorted bounds)",
        "PRECT (sorted bounds)",
        "PRECT (generated bounds)",
    ],
)

# %%
print_stats(
    prect_xcdat1.values,
    prect_xcdat2.values,
    prect_xcdat3.values,
    labels=[
        "PRECT (unsorted bounds)",
        "PRECT (sorted bounds)",
        "PRECT (generated bounds)",
    ],
)


# %%
# Compare arrays for floating-point differences

# %%
# Compare weights arrays for floating-point differences
try:
    np.testing.assert_allclose(weights1, weights2, rtol=1e-5, atol=0)
    print(
        "\nWeights (unsorted bounds) and Weights (sorted bounds) are within rtol=1e-5"
    )
except AssertionError as e:
    print(
        "\nWeights (unsorted bounds) and Weights (sorted bounds) are not within rtol=1e-5"
    )
    print(e)

try:
    np.testing.assert_allclose(weights2, weights3, rtol=1e-5, atol=0)
    print(
        "\nWeights (sorted bounds) and Weights (generated bounds) are within rtol=1e-5"
    )
except AssertionError as e:
    print(
        "\nWeights (sorted bounds) and Weights (generated bounds) are not within rtol=1e-5"
    )
    print(e)

try:
    np.testing.assert_allclose(weights1, weights3, rtol=1e-5, atol=0)
    print(
        "\nWeights (unsorted bounds) and Weights (generated bounds) are within rtol=1e-5"
    )
except AssertionError as e:
    print(
        "\nWeights (unsorted bounds) and Weights (generated bounds) are not within rtol=1e-5"
    )
    print(e)

# %%

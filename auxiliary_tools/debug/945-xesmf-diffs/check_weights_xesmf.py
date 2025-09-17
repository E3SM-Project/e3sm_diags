# %%
"""
This script demonstrates the difference in regridding weights when using
datasets with and without pre-existing lat/lon bounds.

Key Steps:
1. Load datasets.
2. Regrid with pre-existing bounds.
3. Regrid with sorted bounds.
4. Regrid with generated bounds.
5. Compare weights for floating-point differences.

Findings:
- The order of latitude bounds (ascending vs. descending) affects the generated
    regridding weights. xESMF considers the order of bounds when calculating weights
    to ensure consistency with the conservative regridding method, which relies on
    the spatial orientation of grid cells.
- Regridding weights generated with unsorted bounds differ significantly from those
    generated with sorted or newly created bounds, as shown by the statistical and
    floating-point comparison results.
- Weights generated with sorted bounds and generated bounds are consistent within
    the specified tolerance, indicating that proper bounds handling is critical for
    accurate regridding.

This highlights the importance of verifying and standardizing bounds in datasets
before performing regridding operations to avoid unintended discrepancies.
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
ds_b_cp = ds_b.copy(deep=True)
output_grid1 = ds_a.regridder.grid
regridder1 = xe.Regridder(ds_b_cp.copy(), output_grid1, method="conservative_normed")
weights1 = regridder1.weights

# %%
# 2. Regrid with sorted, ascending lat_bnds in reference dataset
ds_b2 = ds_b.copy(deep=True)
ds_b2["lat_bnds"].values = np.sort(ds_b2["lat_bnds"], axis=-1)
regridder2 = xe.Regridder(ds_b2, output_grid1, method="conservative_normed")
weights2 = regridder2.weights

# %%
# 3. Regrid with generated lat_bnds in reference dataset
ds_a3 = ds_a.drop_vars(["lat_bnds"]).bounds.add_missing_bounds(axes=["Y"])
ds_b3 = ds_b.drop_vars(["lat_bnds"]).bounds.add_missing_bounds(axes=["Y"])
output_grid3 = ds_a3.regridder.grid
regridder3 = xe.Regridder(ds_b3, output_grid3, method="conservative_normed")
weights3 = regridder3.weights

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

"""
Stat       Weights (unsorted bounds) Weights (sorted bounds) Weights (generated bounds)
----------------------------------------------------------------------
Min        0.000000             0.000000             0.000000
Max        0.499470             0.374993             0.374993
Mean       0.000017             0.000017             0.000017
Sum        14400.000000         14400.000000         14400.000000
Std        0.002107             0.002091             0.002091
"""


# %%
# Compare weights arrays for floating-point differences
def compare_weights(weights_list, labels, rtol=1e-5, atol=0):
    """Compares dense weights arrays for floating-point differences."""
    results = []
    for i in range(len(weights_list)):
        for j in range(i + 1, len(weights_list)):
            try:
                np.testing.assert_allclose(
                    weights_list[i], weights_list[j], rtol=rtol, atol=atol
                )
                results.append(
                    f"✅ Dense Weights ({labels[i]}) and Dense Weights ({labels[j]}) are within rtol={rtol}"
                )
            except AssertionError as e:
                results.append(
                    f"❌ Dense Weights ({labels[i]}) and Dense Weights ({labels[j]}) are NOT within rtol={rtol}\n{e}"
                )
    print("\nComparison Results:")
    print("\n".join(results))


compare_weights(
    [weights1_dense, weights2_dense, weights3_dense],
    [
        "unsorted bounds",
        "sorted bounds",
        "generated bounds",
    ],
)

# %%
"""
Comparison Results:
❌ Dense Weights (unsorted bounds) and Dense Weights (sorted bounds) are NOT within rtol=1e-05

Not equal to tolerance rtol=1e-05, atol=0

Mismatched elements: 113760 / 829440000 (0.0137%)
Max absolute difference among violations: 0.29999086
Max relative difference among violations: 1.
 ACTUAL: array([[0.49947 , 0.49947 , 0.      , ..., 0.      , 0.      , 0.      ],
       [0.      , 0.      , 0.49947 , ..., 0.      , 0.      , 0.      ],
       [0.      , 0.      , 0.      , ..., 0.      , 0.      , 0.      ],...
 DESIRED: array([[0.251056, 0.251056, 0.      , ..., 0.      , 0.      , 0.      ],
       [0.      , 0.      , 0.251056, ..., 0.      , 0.      , 0.      ],
       [0.      , 0.      , 0.      , ..., 0.      , 0.      , 0.      ],...
❌ Dense Weights (unsorted bounds) and Dense Weights (generated bounds) are NOT within rtol=1e-05

Not equal to tolerance rtol=1e-05, atol=0

Mismatched elements: 113760 / 829440000 (0.0137%)
Max absolute difference among violations: 0.29999086
Max relative difference among violations: 1.
 ACTUAL: array([[0.49947 , 0.49947 , 0.      , ..., 0.      , 0.      , 0.      ],
       [0.      , 0.      , 0.49947 , ..., 0.      , 0.      , 0.      ],
       [0.      , 0.      , 0.      , ..., 0.      , 0.      , 0.      ],...
 DESIRED: array([[0.251056, 0.251056, 0.      , ..., 0.      , 0.      , 0.      ],
       [0.      , 0.      , 0.251056, ..., 0.      , 0.      , 0.      ],
       [0.      , 0.      , 0.      , ..., 0.      , 0.      , 0.      ],...

✅ Dense Weights (sorted bounds) and Dense Weights (generated bounds) are within rtol=1e-05
"""

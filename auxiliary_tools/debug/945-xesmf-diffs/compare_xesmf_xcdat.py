import cdms2
import numpy as np
import pandas as pd
import xarray as xr
import xcdat as xc
import xesmf as xe

from regrid2 import Regridder

def print_stats(*arrays, labels=None):
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

    df = pd.DataFrame(stats, index=labels)
    print("\nStatistical Comparison:")
    print(df)

# Load datasets
ds_a = xr.open_dataset(
    "/lcrc/group/e3sm/public_html/cdat-migration-fy24/25-02-14-branch-930-polar-after/polar/GPCP_v3.2/test.nc"
)
ds_b = xr.open_dataset(
    "/lcrc/group/e3sm/public_html/cdat-migration-fy24/25-02-14-branch-930-polar-after/polar/GPCP_v3.2/ref.nc"
)

# Regridding tests
output_grid_xesmf = ds_a.regridder.grid

# 1. Asc lat, Desc lat_bnds
#unaligned1 = regrid_with_xesmf(ds_b, output_grid_xesmf)
ds_b1 = ds_b.copy(deep=True)
unaligned1 = ds_b1.regridder.horizontal(
            "PRECT", output_grid_xesmf, tool='xesmf', method='conservative_normed'
        )

# 2. Desc lat, Asc lat_bnds
ds_b2 = ds_b.copy(deep=True).sortby("lat", ascending=False)
ds_b2["lat_bnds"].values = np.sort(ds_b2["lat_bnds"], axis=-1)
unaligned2 = ds_b2.regridder.horizontal(
            "PRECT", output_grid_xesmf, tool='xesmf', method='conservative_normed'
        )

# 3. Asc lat, Asc lat_bnds
ds_b3 = ds_b.copy(deep=True)
ds_b3["lat_bnds"].values = np.sort(ds_b3["lat_bnds"], axis=-1)
#aligned1 = regrid_with_xesmf(ds_b3, output_grid_xesmf)
aligned1 = ds_b3.regridder.horizontal(
            "PRECT", output_grid_xesmf, tool='xesmf', method='conservative_normed'
        )

# 4. Desc lat, Desc lat_bnds
ds_b4 = ds_b.copy(deep=True).sortby("lat", ascending=False)
ds_b4["lat_bnds"].values = np.sort(ds_b4["lat_bnds"], axis=-1)[:, ::-1]
#aligned2 = regrid_with_xesmf(ds_b4, output_grid_xesmf)
aligned2 = ds_b4.regridder.horizontal(
            "PRECT", output_grid_xesmf, tool='xesmf', method='conservative_normed'
        )

# 5. Desc lat, no lat_bnds
ds_b5 = ds_b.copy(deep=True).sortby("lat", ascending=False)
ds_b5 = ds_b5.drop("lat_bnds")
nobounds1 = ds_b5.regridder.horizontal(
            "PRECT", output_grid_xesmf, tool='xesmf', method='conservative_normed'
        )

# 6. Asc lat, no lat_bnds
ds_b6 = ds_b.copy(deep=True).sortby("lat", ascending=True)
ds_b6 = ds_b6.drop("lat_bnds")
nobounds2 = ds_b6.regridder.horizontal(
            "PRECT", output_grid_xesmf, tool='xesmf', method='conservative_normed'
        )

print_stats(
    unaligned1.PRECT.values,
    unaligned2.PRECT.values,
    aligned1.PRECT.values,
    aligned2.PRECT.values,
    nobounds1.PRECT.values,
    nobounds2.PRECT.values,
    labels=[
        "asc lat, desc lat_bnds",
        "desc lat, asc lat_bnds",
        "asc lat, asc lat_bnds",
        "desc lat, desc lat_bnds",
        "desc lat, no lat_bnds",
        "asc lat, no lat_bnds",
    ],
)

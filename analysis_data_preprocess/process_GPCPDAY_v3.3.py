#!/usr/bin/env python3
"""
Process GPCPDAY v3.3 (0.5 Degree 360x720) into a single time_series file for e3sm_diags.

Replaces:
- ncks --mk_rec_dmn time (not needed; xarray writes an unlimited time dimension)
- ncrcat (xarray concat along time)

Author: Jill Zhang (zhang40@llnl.gov) - Python rewrite
Dataset: https://disc.gsfc.nasa.gov/datasets/GPCPDAY_3.3/summary
"""
from pathlib import Path
import xarray as xr

# paths + years
path = Path("/global/cfs/cdirs/e3sm/zhang40/analysis_data_e3sm_diags/GPCPDAY_0.5D_3.3/")
orig = path / "original"
outdir = path / "time_series_py"
outdir.mkdir(parents=True, exist_ok=True)

start_yr, end_yr = 1998, 2023
out = outdir / f"precip_{start_yr}01_{end_yr}12.nc"

# list all daily files in year range
files = []
for y in range(start_yr, end_yr + 1):
    files += sorted(orig.glob(f"GPCPDAY_L3_{y:04d}????_V3.3.nc4"))

# open + concat along time (each file has time=1)
ds = xr.open_mfdataset(
    [str(f) for f in files],
    combine="by_coords",
    engine="netcdf4",
    decode_cf=True,
    mask_and_scale=True,
)[["precip"]]  # keep only precip

# write one time series file (time unlimited)
ds.to_netcdf(out, format="NETCDF4", unlimited_dims=("time",))
ds.close()

print("Wrote:", out)

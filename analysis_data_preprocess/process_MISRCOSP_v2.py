import glob
import os
from datetime import datetime, timezone

import dask
import numpy as np
import pandas as pd
import xarray as xr
from pcmdi_metrics.io import xcdat_open

SCRIPT_NAME = os.path.basename(__file__)


def prepend_history(ds_out, message):
    """Prepend a processing-history line (CF convention)."""
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    new_line = f"{timestamp}: {SCRIPT_NAME}: {message}"
    existing = ds_out.attrs.get("history", "")
    ds_out.attrs["history"] = f"{new_line}\n{existing}" if existing else new_line
    return ds_out

"""
Process MISR clMISR obs4MIPs monthly cloud joint-histograms into a monthly
time series and ANN/DJF/MAM/JJA/SON climatology files for the COSP
cloud-fraction diags in e3sm_diags.

Ported from misr_diags_{ANN,DJF,JJA}.pro (originally written in IDL by
Yuying Zhang).

MISR related URL: https://misr.jpl.nasa.gov/
Note: the IDL source multiplies clMISR by 100 even though the obs4MIPs
input is already in percent; that scaling is dropped here so output
values are correct percent (0..100).
"""

# settings
original_path = "/global/cfs/cdirs/e3sm/zhang40/analysis_data_e3sm_diags/COSP/obs_sat/MISR/clMISR/"
output_path = "/global/cfs/cdirs/e3sm/zhang40/analysis_data_e3sm_diags/COSP/obs_sat/MISRCOSP_v2/climatology/"
output_path_time_series = "/global/cfs/cdirs/e3sm/zhang40/analysis_data_e3sm_diags/COSP/obs_sat/MISRCOSP_v2/time_series/"
start_yr = 2000
end_yr = 2020

os.makedirs(output_path, exist_ok=True)
os.makedirs(output_path_time_series, exist_ok=True)

case_id = "MISRCOSP_v2"


def preprocess_one_month(file_path):
    """Read one clMISR obs4MIPs monthly file and reshape to e3sm_diags layout.

    Source: clMISR(time=1, tau=8, cth=16, lat=180, lon=360), units=%,
    lat 89.5..-89.5, lon -179.5..179.5, cth in meters with -1 sentinel as
    bin 0, _FillValue=1e+20 (xarray auto-masks to NaN).

    Output: CLMISR(misr_cth, misr_tau, lat, lon), lat -89.5..89.5,
    lon 0.5..359.5, misr_cth in km. Values are percent (0..100).
    """
    # chunks={} keeps the whole file lazy as one dask block (no splits) — lets
    # xr.concat build a lazy 8 GB time series instead of loading it all in RAM.
    ds = xr.open_dataset(file_path, decode_times=False, chunks={})

    ds = ds.rename(
        {
            "clMISR": "CLMISR",
            "tau": "misr_tau",
            "tau_bnds": "misr_tau_bnds",
            "cth": "misr_cth",
            "cth_bnds": "misr_cth_bnds",
        }
    )
    ds = ds[["CLMISR", "misr_cth", "misr_cth_bnds", "misr_tau", "misr_tau_bnds"]]

    ds["CLMISR"] = ds["CLMISR"].squeeze("time", drop=True)

    ds["CLMISR"] = ds["CLMISR"].transpose("misr_cth", "misr_tau", "lat", "lon")
    ds["misr_cth_bnds"] = ds["misr_cth_bnds"].transpose("misr_cth", "bnds")
    ds["misr_tau_bnds"] = ds["misr_tau_bnds"].transpose("misr_tau", "bnds")

    # convert cloud-top height from m to km
    ds["misr_cth"] = ds["misr_cth"] / 1000.0
    ds["misr_cth_bnds"] = ds["misr_cth_bnds"] / 1000.0
    ds["misr_cth"].attrs.update(
        {"long_name": "cloud top height", "units": "km", "bounds": "misr_cth_bnds"}
    )
    ds["misr_cth_bnds"].attrs["units"] = "km"
    ds["misr_tau"].attrs.update(
        {"long_name": "cloud optical depth", "units": "unitless", "bounds": "misr_tau_bnds"}
    )
    ds["misr_tau_bnds"].attrs["units"] = "unitless"

    # reverse latitude to ascending
    ds = ds.isel(lat=slice(None, None, -1))

    # rotate longitude from [-180, 180) to [0, 360)
    ds = ds.assign_coords(lon=((ds["lon"] + 360) % 360))
    ds = ds.sortby("lon")

    ds["CLMISR"].attrs.update(
        {"long_name": "Cloud Occurrence Histogram", "units": "percent"}
    )
    return ds


# build monthly time series
monthly_list = []
time_stamps = []
missing_months = []
for iy in range(start_yr, end_yr + 1):
    for im in range(1, 13):
        patt = os.path.join(
            original_path, f"clMISR_obs4MIPs_MISR_V7_{iy:04d}{im:02d}*.nc"
        )
        matches = sorted(glob.glob(patt))
        if not matches:
            missing_months.append(f"{iy:04d}-{im:02d}")
            continue
        monthly_list.append(preprocess_one_month(matches[0]))
        time_stamps.append(f"{iy:04d}-{im:02d}-15")

# report only gaps inside the actual coverage window
if time_stamps:
    first, last = time_stamps[0][:7], time_stamps[-1][:7]
    interior_gaps = [m for m in missing_months if first <= m <= last]
    if interior_gaps:
        print(f"WARNING: {len(interior_gaps)} month(s) missing inside coverage "
              f"{first}..{last}: {interior_gaps}")

time_index = pd.to_datetime(time_stamps)
# data_vars=["CLMISR"] keeps the small bounds/coord variables un-broadcast over
# the 243-month time axis — concat="all" (the default) replicates them per file
# and inflates both the dask graph and peak RAM during to_netcdf.
time_series = xr.concat(monthly_list, dim="time", data_vars=["CLMISR"])
time_series = time_series.assign_coords(time=("time", time_index))
time_series["time"].attrs["axis"] = "T"

# actual coverage YYYYMM (data may not start in Jan or end in Dec)
actual_start = f"{time_index[0].year:04d}{time_index[0].month:02d}"
actual_end = f"{time_index[-1].year:04d}{time_index[-1].month:02d}"

time_series.attrs["Description"] = "MISR cloud observations for climate model evaluation"
prepend_history(
    time_series,
    f"concatenated monthly source files into time series ({actual_start}-{actual_end})",
)

output_filename = f"{case_id}_{actual_start}_{actual_end}.nc"
# single-threaded scheduler keeps peak RAM near one chunk; the threaded default
# fans out 243 concurrent chunk reads and exceeds the 30 GB cgroup cap on the
# Perlmutter login node, killing the process silently mid-write.
with dask.config.set(scheduler="single-threaded"):
    time_series.to_netcdf(output_path_time_series + output_filename)
print(f"Saved monthly time series: {output_path_time_series + output_filename}")

# compute seasonal/annual climatologies via xcdat. Chunk by time so the
# 1.6 GB time series is processed in a streaming fashion (avoids OOM).
ds = xcdat_open(output_path_time_series + output_filename, chunks={"time": 12})
ds = ds.bounds.add_missing_bounds("T")

# midpoint of the averaging period — used as the singleton time coord
period_midpoint = time_index[len(time_index) // 2]

seasons = ["ANN", "DJF", "MAM", "JJA", "SON"]
season_idx = {"DJF": 0, "MAM": 1, "JJA": 2, "SON": 3}


def _finalize(ds_out, time_value, season):
    """Cast CLMISR to float32, ensure a singleton time dim, and stamp
    climatology metadata (yrs_averaged + a CF-style history line)."""
    if "CLMISR" in ds_out:
        ds_out["CLMISR"] = ds_out["CLMISR"].astype("float32")
    if "time" not in ds_out.dims:
        ds_out = ds_out.expand_dims({"time": [time_value]})
    ds_out.attrs["yrs_averaged"] = f"{actual_start}-{actual_end}"
    label = "annual mean" if season == "ANN" else f"{season} seasonal mean"
    prepend_history(
        ds_out, f"computed {label} over {actual_start}-{actual_end}"
    )
    return ds_out


for season in seasons:
    if season == "ANN":
        season_ds = ds.temporal.average("CLMISR", weighted=True)
    else:
        clim = ds.temporal.climatology(
            "CLMISR",
            freq="season",
            weighted=True,
            season_config={"dec_mode": "DJF"},
        )
        season_ds = clim.isel(time=season_idx[season])
        # xcdat folds time to a placeholder year-1 cftime stamp; replace with
        # the period midpoint so the time coord is a real, serializable date
        season_ds = season_ds.drop_vars(
            [v for v in ("time", "time_bnds") if v in season_ds.variables]
        )

    season_ds = _finalize(season_ds, period_midpoint, season)

    season_filename = f"{case_id}_{actual_start}_{actual_end}_{season}.nc"
    with dask.config.set(scheduler="single-threaded"):
        season_ds.to_netcdf(os.path.join(output_path, season_filename))
    print(f"Saved {season_filename}")

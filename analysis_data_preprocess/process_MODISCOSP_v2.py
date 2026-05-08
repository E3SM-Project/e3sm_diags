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
Process MODIS L3 (MCD08_M3) monthly cloud joint-histograms into a monthly
time series and ANN/DJF/MAM/JJA/SON climatology files for the COSP
cloud-fraction diags in e3sm_diags.

Ported from modis_diags_{ANN,DJF,JJA}.pro (originally written in IDL by
Yuying Zhang).

MODIS related URL: https://modis-atmos.gsfc.nasa.gov/products/cloud
"""

# settings
original_path = "/global/cfs/cdirs/e3sm/zhang40/analysis_data_e3sm_diags/COSP/obs_sat/MCD_2002-2016/"
output_path = "/global/cfs/cdirs/e3sm/zhang40/analysis_data_e3sm_diags/COSP/obs_sat/MODISCOSP_v2/climatology/"
output_path_time_series = "/global/cfs/cdirs/e3sm/zhang40/analysis_data_e3sm_diags/COSP/obs_sat/MODISCOSP_v2/time_series/"
start_yr = 2002
end_yr = 2016

os.makedirs(output_path, exist_ok=True)
os.makedirs(output_path_time_series, exist_ok=True)

case_id = "MODISCOSP_v2"


def preprocess_one_month(file_path):
    """Read one MCD08_M3 monthly file and reshape to e3sm_diags layout.

    Source: Optical_Thickness_vs_Cloud_Top_Pressure(lat, lon, ctp, cot),
    lat 89.5..-89.5, lon -179.5..179.5, ctp ascending Pa, fraction (0..1)
    with -999 fill.

    Output: CLMODIS(modis_prs, modis_tau, lat, lon), lat -89.5..89.5,
    lon 0.5..359.5, modis_prs descending hPa, percent.
    """
    ds = xr.open_dataset(file_path, decode_times=False)

    ds = ds.rename(
        {
            "Optical_Thickness_vs_Cloud_Top_Pressure": "CLMODIS",
            "Cloud_Top_Pressure": "modis_prs",
            "Cloud_Top_Pressure_bnds": "modis_prs_bnds",
            "Cloud_Optical_Thickness": "modis_tau",
            "Cloud_Optical_Thickness_bnds": "modis_tau_bnds",
            "nbnds": "bnds",
        }
    )
    ds = ds[["CLMODIS", "modis_prs", "modis_prs_bnds", "modis_tau", "modis_tau_bnds"]]

    ds["CLMODIS"] = xr.where(ds["CLMODIS"] == -999.0, np.nan, ds["CLMODIS"]) * 100.0

    ds["CLMODIS"] = ds["CLMODIS"].transpose("modis_prs", "modis_tau", "lat", "lon")
    ds["modis_prs_bnds"] = ds["modis_prs_bnds"].transpose("modis_prs", "bnds")
    ds["modis_tau_bnds"] = ds["modis_tau_bnds"].transpose("modis_tau", "bnds")

    # reverse pressure axis (highest pressure first), Pa -> hPa
    ds = ds.isel(modis_prs=slice(None, None, -1))
    ds["modis_prs"] = ds["modis_prs"] / 100.0
    ds["modis_prs_bnds"] = ds["modis_prs_bnds"] / 100.0
    ds["modis_prs"].attrs.update({"units": "hPa", "long_name": "Cloud Top Pressure", "bounds": "modis_prs_bnds"})
    ds["modis_prs_bnds"].attrs["units"] = "hPa"
    ds["modis_tau"].attrs.update({"units": "unitless", "long_name": "cloud optical depth", "bounds": "modis_tau_bnds"})
    ds["modis_tau_bnds"].attrs["units"] = "unitless"

    # reverse latitude to ascending
    ds = ds.isel(lat=slice(None, None, -1))

    # rotate longitude from [-180, 180) to [0, 360)
    ds = ds.assign_coords(lon=((ds["lon"] + 360) % 360))
    ds = ds.sortby("lon")

    ds["CLMODIS"].attrs.update(
        {"long_name": "MODIS histogram of cloud occurrence", "units": "percent"}
    )
    return ds


# build monthly time series
monthly_list = []
time_stamps = []
missing_months = []
for iy in range(start_yr, end_yr + 1):
    for im in range(1, 13):
        patt = os.path.join(original_path, f"MCD08_M3_NC.{iy:04d}.{im:02d}*.C051.V02.nc")
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
time_series = xr.concat(monthly_list, dim="time")
time_series = time_series.assign_coords(time=("time", time_index))
time_series["time"].attrs["axis"] = "T"

# actual coverage YYYYMM (data may not start in Jan or end in Dec)
actual_start = f"{time_index[0].year:04d}{time_index[0].month:02d}"
actual_end = f"{time_index[-1].year:04d}{time_index[-1].month:02d}"

time_series.attrs["Retrieval_Version"] = "Collection 51 v2"
time_series.attrs["Description"] = "MODIS cloud observations for climate model evaluation"
prepend_history(
    time_series,
    f"concatenated monthly source files into time series ({actual_start}-{actual_end})",
)

output_filename = f"{case_id}_{actual_start}_{actual_end}.nc"
time_series.to_netcdf(output_path_time_series + output_filename)
print(f"Saved monthly time series: {output_path_time_series + output_filename}")

# compute seasonal/annual climatologies via xcdat. Chunk by time so the
# 2.2 GB time series is processed in a streaming fashion (avoids the 30 GB
# Perlmutter login-node cgroup cap when xcdat builds intermediate seasonal
# arrays).
ds = xcdat_open(output_path_time_series + output_filename, chunks={"time": 12})
ds = ds.bounds.add_missing_bounds("T")

# midpoint of the averaging period — used as the singleton time coord
period_midpoint = time_index[len(time_index) // 2]

seasons = ["ANN", "DJF", "MAM", "JJA", "SON"]
season_idx = {"DJF": 0, "MAM": 1, "JJA": 2, "SON": 3}


def _finalize(ds_out, time_value, season):
    """Cast CLMODIS to float32, ensure a singleton time dim, and stamp
    climatology metadata (yrs_averaged + a CF-style history line)."""
    if "CLMODIS" in ds_out:
        ds_out["CLMODIS"] = ds_out["CLMODIS"].astype("float32")
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
        season_ds = ds.temporal.average("CLMODIS", weighted=True)
    else:
        clim = ds.temporal.climatology(
            "CLMODIS",
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

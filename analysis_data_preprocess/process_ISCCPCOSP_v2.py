import os
from datetime import datetime, timezone

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
Process the single ISCCP-H Series file (1983-01..2017-12, 420 months) into a
monthly time series and ANN/DJF/MAM/JJA/SON climatology files for the COSP
cloud-fraction diags in e3sm_diags.

Ported from isccp_diags_{ANN,DJF,JJA}.pro (originally written in IDL by
Yuying Zhang).

ISCCP-H Series related URL: https://www.ncei.noaa.gov/products/international-satellite-cloud-climatology
"""

# settings
original_file = "/global/cfs/cdirs/e3sm/zhang40/analysis_data_e3sm_diags/COSP/obs_sat/ISCCP-H/clisccp_ISCCP_HGG_198301-201712.nc"
output_path = "/global/cfs/cdirs/e3sm/zhang40/analysis_data_e3sm_diags/COSP/obs_sat/ISCCPCOSP_v2/climatology/"
output_path_time_series = "/global/cfs/cdirs/e3sm/zhang40/analysis_data_e3sm_diags/COSP/obs_sat/ISCCPCOSP_v2/time_series/"
start_yr = 1983
end_yr = 2017

os.makedirs(output_path, exist_ok=True)
os.makedirs(output_path_time_series, exist_ok=True)

case_id = "ISCCPCOSP_v2"

# ISCCP pressure-bin bounds (hPa) ordered to match p_midpt = [900, 740, 620, 500, 375, 245, 105]
# Set in IDL via altb[1,*] = [1000,800,680,560,440,310,180]; altb[0,*] = [800,680,560,440,310,180,1]
isccp_prs_bnds = np.array(
    [
        [800.0, 1000.0],
        [680.0, 800.0],
        [560.0, 680.0],
        [440.0, 560.0],
        [310.0, 440.0],
        [180.0, 310.0],
        [1.0, 180.0],
    ],
    dtype=np.float32,
)


# --- read source ---
# source clisccp is float64 9GB — chunk by time to avoid OOM. Default
# mask_and_scale=True converts the 1e+20 fill to NaN automatically.
src = xr.open_dataset(original_file, decode_times=False, chunks={"time": 12})
# clisccp dims (time, levtau, p_midpt, lat, lon); values in percent already (0..100)
clisccp = src["clisccp"].astype(np.float32)

# rename to e3sm_diags conventions
clisccp = clisccp.rename({"levtau": "isccp_tau", "p_midpt": "isccp_prs"})
clisccp = clisccp.rename("CLISCCP")
clisccp = clisccp.transpose("time", "isccp_prs", "isccp_tau", "lat", "lon")
clisccp.attrs.update(
    {"long_name": "Cloud Occurrence Histogram", "units": "percent"}
)

# build CF time coord (the file's nominal range is 1983-01 .. 2017-12)
n_time = clisccp.sizes["time"]
time_stamps = []
for ti in range(n_time):
    yr = start_yr + ti // 12
    mo = ti % 12 + 1
    time_stamps.append(f"{yr:04d}-{mo:02d}-15")
time_index = pd.to_datetime(time_stamps)
clisccp = clisccp.assign_coords(time=("time", time_index))

# trim all-fill (entirely-NaN) leading/trailing months — ISCCP-H's actual
# coverage is narrower than the file's nominal range
valid_per_time = (
    clisccp.notnull()
    .any(dim=("isccp_prs", "isccp_tau", "lat", "lon"))
    .compute()
)
valid_idx = np.flatnonzero(valid_per_time.values)
if valid_idx.size == 0:
    raise RuntimeError("No valid time slices found in source")
first_idx, last_idx = int(valid_idx[0]), int(valid_idx[-1])
clisccp = clisccp.isel(time=slice(first_idx, last_idx + 1))
trimmed_index = time_index[first_idx : last_idx + 1]

# any all-fill months interior to the trimmed range are interior gaps
interior_gaps = []
for ti, ts in enumerate(trimmed_index):
    if not bool(valid_per_time.values[first_idx + ti]):
        interior_gaps.append(f"{ts.year:04d}-{ts.month:02d}")
if interior_gaps:
    print(f"WARNING: {len(interior_gaps)} all-fill month(s) interior to coverage: "
          f"{interior_gaps}")

# coordinate variables (already correct order/orientation in source)
ds_out = xr.Dataset(
    {
        "CLISCCP": clisccp,
        "isccp_prs_bnds": (("isccp_prs", "bnds"), isccp_prs_bnds),
        "isccp_tau_bnds": (("isccp_tau", "bnds"), src["tau_bounds"].values.astype(np.float32)),
    }
)
ds_out["isccp_prs"].attrs.update(
    {"long_name": "Cloud Top Pressure", "units": "hPa", "bounds": "isccp_prs_bnds"}
)
ds_out["isccp_tau"].attrs.update(
    {"long_name": "Cloud Optical thickness", "units": "unitless", "bounds": "isccp_tau_bnds"}
)
ds_out["lat"].attrs.update({"long_name": "latitude", "units": "degrees_north"})
ds_out["lon"].attrs.update({"long_name": "longitude", "units": "degrees_east"})
ds_out["time"].attrs["axis"] = "T"

# actual coverage YYYYMM (ISCCP-H starts/ends inside the file's nominal range)
actual_start = f"{trimmed_index[0].year:04d}{trimmed_index[0].month:02d}"
actual_end = f"{trimmed_index[-1].year:04d}{trimmed_index[-1].month:02d}"

ds_out.attrs["Comments"] = "This data was generated from the ISCCP-H Series dataset."
prepend_history(
    ds_out,
    f"trimmed all-fill leading/trailing months and built time series ({actual_start}-{actual_end})",
)

output_filename = f"{case_id}_{actual_start}_{actual_end}.nc"
ds_out.to_netcdf(output_path_time_series + output_filename)
print(f"Saved monthly time series: {output_path_time_series + output_filename}")

# compute seasonal/annual climatologies via xcdat. Chunk by time so the
# 4 GB time series is processed in a streaming fashion (avoids OOM).
ds = xcdat_open(output_path_time_series + output_filename, chunks={"time": 12})
ds = ds.bounds.add_missing_bounds("T")

# midpoint of the averaging period — used as the singleton time coord
period_midpoint = trimmed_index[len(trimmed_index) // 2]

seasons = ["ANN", "DJF", "MAM", "JJA", "SON"]
season_idx = {"DJF": 0, "MAM": 1, "JJA": 2, "SON": 3}


def _finalize(ds_out, time_value, season):
    """Cast CLISCCP to float32, ensure a singleton time dim, and stamp
    climatology metadata (yrs_averaged + a CF-style history line)."""
    if "CLISCCP" in ds_out:
        ds_out["CLISCCP"] = ds_out["CLISCCP"].astype("float32")
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
        season_ds = ds.temporal.average("CLISCCP", weighted=True)
    else:
        clim = ds.temporal.climatology(
            "CLISCCP",
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
    season_ds.to_netcdf(os.path.join(output_path, season_filename))
    print(f"Saved {season_filename}")

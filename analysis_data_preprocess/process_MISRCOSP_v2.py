import glob
import os

import numpy as np
import pandas as pd
import xarray as xr
from pcmdi_metrics.io import xcdat_open

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
    ds = xr.open_dataset(file_path, decode_times=False)

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
for iy in range(start_yr, end_yr + 1):
    for im in range(1, 13):
        patt = os.path.join(
            original_path, f"clMISR_obs4MIPs_MISR_V7_{iy:04d}{im:02d}*.nc"
        )
        matches = sorted(glob.glob(patt))
        if not matches:
            continue
        monthly_list.append(preprocess_one_month(matches[0]))
        time_stamps.append(f"{iy:04d}-{im:02d}-15")

time_index = pd.to_datetime(time_stamps)
time_series = xr.concat(monthly_list, dim="time")
time_series = time_series.assign_coords(time=("time", time_index))
time_series["time"].attrs["axis"] = "T"

time_series.attrs["Description"] = "MISR cloud observations for climate model evaluation"
time_series.attrs["Coverage"] = (
    f"{time_stamps[0][:4]}{time_stamps[0][5:7]}-{time_stamps[-1][:4]}{time_stamps[-1][5:7]}"
)

output_filename = f"{case_id}_{start_yr}01_{end_yr}12.nc"
time_series.to_netcdf(output_path_time_series + output_filename)
print(f"Saved monthly time series: {output_path_time_series + output_filename}")

# compute seasonal/annual climatologies via xcdat
ds = xcdat_open(output_path_time_series + output_filename)
ds = ds.bounds.add_missing_bounds("T")

seasons = ["ANN", "DJF", "MAM", "JJA", "SON"]
season_idx = {"DJF": 0, "MAM": 1, "JJA": 2, "SON": 3}
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

    season_filename = f"{case_id}_{start_yr}01_{end_yr}12_{season}.nc"
    season_ds.to_netcdf(os.path.join(output_path, season_filename))
    print(f"Saved {season_filename}")

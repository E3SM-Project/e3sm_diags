import glob
import os

import numpy as np
import pandas as pd
import xarray as xr
from pcmdi_metrics.io import xcdat_open

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
        {"long_name": "MODIS histogram of cloud occurrence", "units": "precent"}
    )
    return ds


# build monthly time series
monthly_list = []
time_stamps = []
for iy in range(start_yr, end_yr + 1):
    for im in range(1, 13):
        patt = os.path.join(original_path, f"MCD08_M3_NC.{iy:04d}.{im:02d}*.C051.V02.nc")
        matches = sorted(glob.glob(patt))
        if not matches:
            continue
        monthly_list.append(preprocess_one_month(matches[0]))
        time_stamps.append(f"{iy:04d}-{im:02d}-15")

time_index = pd.to_datetime(time_stamps)
time_series = xr.concat(monthly_list, dim="time")
time_series = time_series.assign_coords(time=("time", time_index))
time_series["time"].attrs["axis"] = "T"

time_series.attrs["Retrieval_Version"] = "Collection 51 v2"
time_series.attrs["Description"] = "MODIS cloud observations for climate model evaluation"
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
        season_ds = ds.temporal.average("CLMODIS", weighted=True)
    else:
        clim = ds.temporal.climatology(
            "CLMODIS",
            freq="season",
            weighted=True,
            season_config={"dec_mode": "DJF"},
        )
        season_ds = clim.isel(time=season_idx[season])

    season_filename = f"{case_id}_{start_yr}01_{end_yr}12_{season}.nc"
    season_ds.to_netcdf(os.path.join(output_path, season_filename))
    print(f"Saved {season_filename}")

import glob
import os

import numpy as np
import pandas as pd
import xarray as xr
from pcmdi_metrics.io import xcdat_open

"""
Process CALIPSO GOCCP monthly cloud-fraction maps (high/mid/low/total) into a
monthly time series and ANN/DJF/MAM/JJA/SON climatology files for the COSP
cloud-fraction diags in e3sm_diags.

Ported from cal_hml_diags_{ANN,DJF,JJA}.pro (originally written in IDL by
Yuying Zhang).

CALIPSO GOCCP related URL: https://climserv.ipsl.polytechnique.fr/cfmip-obs/
"""

# settings
original_path = "/global/cfs/cdirs/e3sm/zhang40/analysis_data_e3sm_diags/COSP/obs_sat/GOCCP/2D_map/"
output_path = "/global/cfs/cdirs/e3sm/zhang40/analysis_data_e3sm_diags/COSP/obs_sat/CALIPSOCOSP_v2/climatology/"
output_path_time_series = "/global/cfs/cdirs/e3sm/zhang40/analysis_data_e3sm_diags/COSP/obs_sat/CALIPSOCOSP_v2/time_series/"
start_yr = 2006
end_yr = 2020

os.makedirs(output_path, exist_ok=True)
os.makedirs(output_path_time_series, exist_ok=True)

case_id = "CALIPSOCOSP_v2"

# source variable -> output variable name and long_name
VAR_MAP = {
    "clhcalipso": ("CLDHGH_CAL", "CALIPSO High-level Cloud Fraction"),
    "clmcalipso": ("CLDMED_CAL", "CALIPSO Mid-level Cloud Fraction"),
    "cllcalipso": ("CLDLOW_CAL", "CALIPSO Low-level Cloud Fraction"),
    "cltcalipso": ("CLDTOT_CAL", "CALIPSO Total Cloud Fraction"),
}


def preprocess_one_month(file_path):
    """Read one GOCCP monthly file and reshape to e3sm_diags layout.

    Source: cl{h,m,l,t}calipso(time=1, latitude=90, longitude=180), values 0..1,
    lat -89..89 (ascending), lon -179..179.

    Output: CLD{HGH,MED,LOW,TOT}_CAL(lat, lon), lat -89..89, lon 1..359, percent.
    """
    ds = xr.open_dataset(file_path, decode_times=False)
    ds = ds.rename({"latitude": "lat", "longitude": "lon"})

    out = xr.Dataset()
    for src_name, (out_name, long_name) in VAR_MAP.items():
        v = ds[src_name].squeeze("time", drop=True)
        # mask out anything outside [0, 1] (matches IDL filter 0..100 on the 0..1 source)
        v = xr.where((v < 0.0) | (v > 1.0), np.nan, v) * 100.0
        v.attrs.update({"long_name": long_name, "units": "percent"})
        out[out_name] = v

    # rotate longitude from [-180, 180) to [0, 360)
    out = out.assign_coords(lon=((out["lon"] + 360) % 360))
    out = out.sortby("lon")

    out["lat"].attrs.update({"long_name": "latitude", "units": "degrees_north", "axis": "Y"})
    out["lon"].attrs.update({"long_name": "longitude", "units": "degrees_east", "axis": "X"})
    return out


# build monthly time series
monthly_list = []
time_stamps = []
for iy in range(start_yr, end_yr + 1):
    for im in range(1, 13):
        patt = os.path.join(
            original_path, f"MapLowMidHigh330m_{iy:04d}{im:02d}_avg_CFMIP2_sat_3.1.2.nc"
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

time_series.attrs["Description"] = "CALIPSO GOCCP v3.1.2"
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
out_vars = list(v[0] for v in VAR_MAP.values())

for season in seasons:
    if season == "ANN":
        merged = xr.merge(
            [ds.temporal.average(v, weighted=True) for v in out_vars]
        )
    else:
        season_dsets = []
        for v in out_vars:
            clim = ds.temporal.climatology(
                v,
                freq="season",
                weighted=True,
                season_config={"dec_mode": "DJF"},
            )
            season_dsets.append(clim.isel(time=season_idx[season]))
        merged = xr.merge(season_dsets)

    season_filename = f"{case_id}_{start_yr}01_{end_yr}12_{season}.nc"
    merged.to_netcdf(os.path.join(output_path, season_filename))
    print(f"Saved {season_filename}")

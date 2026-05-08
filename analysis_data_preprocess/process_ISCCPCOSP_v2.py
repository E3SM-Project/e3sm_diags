import os

import numpy as np
import pandas as pd
import xarray as xr
from pcmdi_metrics.io import xcdat_open

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
src = xr.open_dataset(original_file, decode_times=False)
# clisccp dims (time, levtau, p_midpt, lat, lon); values in percent already; -1e34 fill
clisccp = src["clisccp"]
fill = clisccp.attrs.get("_FillValue", -1e34)
clisccp = xr.where(clisccp == fill, np.nan, clisccp)

# rename to e3sm_diags conventions
clisccp = clisccp.rename({"levtau": "isccp_tau", "p_midpt": "isccp_prs"})
clisccp = clisccp.rename("CLISCCP")
clisccp = clisccp.transpose("time", "isccp_prs", "isccp_tau", "lat", "lon")
clisccp.attrs.update(
    {"long_name": "Cloud Occurrence Histogram", "units": "percent"}
)

# build CF time coord (months since 1983-01)
n_time = clisccp.sizes["time"]
time_stamps = []
for ti in range(n_time):
    yr = start_yr + ti // 12
    mo = ti % 12 + 1
    time_stamps.append(f"{yr:04d}-{mo:02d}-15")
time_index = pd.to_datetime(time_stamps)

# coordinate variables (already correct order/orientation in source)
ds_out = xr.Dataset(
    {
        "CLISCCP": clisccp.assign_coords(time=("time", time_index)),
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

ds_out.attrs["Coverage"] = f"{start_yr}-{end_yr}"
ds_out.attrs["Comments"] = "This data was generated from the ISCCP-H Series dataset."

output_filename = f"{case_id}_{start_yr}01_{end_yr}12.nc"
ds_out.to_netcdf(output_path_time_series + output_filename)
print(f"Saved monthly time series: {output_path_time_series + output_filename}")

# compute seasonal/annual climatologies via xcdat
ds = xcdat_open(output_path_time_series + output_filename)
ds = ds.bounds.add_missing_bounds("T")

seasons = ["ANN", "DJF", "MAM", "JJA", "SON"]
season_idx = {"DJF": 0, "MAM": 1, "JJA": 2, "SON": 3}
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

    season_filename = f"{case_id}_{start_yr}01_{end_yr}12_{season}.nc"
    season_ds.to_netcdf(os.path.join(output_path, season_filename))
    print(f"Saved {season_filename}")

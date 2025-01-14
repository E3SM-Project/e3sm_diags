"""
A script to convert high frequency single point E3SM output to per-variable
per-site netcdf files as input for ARM diagostics.

subprocess.run('source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_cori-haswell.sh', shell=True)

In this example 3 hourly output at ARM sites are saved on h4 tape using namelist as follows:
fincl2 = 'PS', 'Q', 'T', 'Z3', 'CLOUD', 'CONCLD', 'CLDICE', 'CLDLIQ', 'LS_FLXPRC', 'LS_FLXSNW', 'ZMFLXPRC', 'ZMFLXSNW', 'FREQR', 'REI', 'REL', 'CV_REFFICE', 'CV_REFFLIQ', 'LS_REFFRAIN', 'LS_REFFSNOW', 'PRECT', 'TMQ', 'PRECC', 'TREFHT', 'QREFHT', 'OMEGA','CLDTOT', 'LHFLX', 'SHFLX', 'FLDS', 'FSDS', 'FLNS', 'FSNS', 'FLNSC', 'FSDSC', 'FSNSC', 'AODVIS', 'AODABS'

fincl2lonlat = '262.5e_36.6n','203.4e_71.3n','147.4e_2.0s','166.9e_0.5s','130.9e_12.4s','331.97e_39.09n'
"""
import glob
import os
from pathlib import Path

import xarray as xr

# Parameters to define.
DATA_PATH = (
    "/global/cfs/cdirs/e3sm/www/Tutorials/2024/simulations/"
    "extendedOutput.v3.LR.historical_0101/archive/atm/hist/"
)
OUTPUT_PATH = (
    "/global/cfs/cdirs/e3sm/vo13/tutorial2024/v3.LR.historical_0101/post/atm/site/"
)
Path(OUTPUT_PATH).mkdir(parents=True, exist_ok=True)

H_NUM = "h6"
START_YEAR = "2000"
END_YEAR = "2014"
TIME_RANGE = f"{START_YEAR}01_{END_YEAR}12"

print(f"Opening datasets from input path: {os.path.join(DATA_PATH, H_NUM)}")

filepaths = glob.glob(DATA_PATH + f"*{H_NUM}.20*")
ds_sub = xr.open_mfdataset(filepaths, parallel=True, chunks="auto")
ds = ds_sub.sel(time=slice(f"{START_YEAR}-01-01", f"{END_YEAR}-12-31"))

print(f"Xarray Dataset object: {ds_sub}")
print("Create time-series (it can take long depends on number of years being analyzed)")

VARIABLES = [
    "CLOUD",
    "CONCLD",
    "CLDICE",
    "CLDLIQ",
    "PRECT",
    "TMQ",
    "TREFHT",
    "QREFHT",
    "CLDTOT",
    "LHFLX",
    "SHFLX",
    "FLDS",
    "FSDS",
    "FLNS",
    "FSNS",
    "FLNSC",
    "FSDSC",
    "FSNSC",
    "AODVIS",
    "AODABS",
    "PS",
    "num_a1",  # Accumu mode aerosol concentration (1/kg) at lowest level
    "num_a2",  # Aitken mode aerosol concentration (1/kg) at lowest level
    "num_a3",  # Coarse mode aerosol concentration (1/kg) at lowest level
    "so4_a1",  # Accumu mode SO4 mass conc. (kg/kg) at lowest level
    "so4_a2",  # Aitken mode SO4 mass conc. (kg/kg) at lowest level
    "CCN3",  # CCN 0.1%SS concentration (1/CC) at lowest level
    "CCN4",  # CCN 0.2%SS concentration (1/CC) at lowest level
    "CCN5",  # CCN 0.5%SS concentration (1/CC) at lowest level
]

SITE_INFO = {
    "sgpc1": [262.5, 36.6],
    "nsac1": [203.4, 71.3],
    "twpc1": [147.4, -2.0],
    "twpc2": [166.9, -0.5],
    "twpc3": [130.9, -12.4],
    "enac1": [331.97, 39.09],
}
SITES = ["sgpc1", "nsac1", "twpc1", "twpc2", "twpc3", "enac1"]

datasets = {}

for site in SITES:
    if SITE_INFO[site][1] > 0:
        lon_lat = str(SITE_INFO[site][0]) + "e_" + str(SITE_INFO[site][1]) + "n"
    else:
        lon_lat = str(SITE_INFO[site][0]) + "e_" + str(abs(SITE_INFO[site][1])) + "s"

    for variable in VARIABLES:
        fname_out = f"{OUTPUT_PATH}{variable}_{site}_{TIME_RANGE}.nc"
        print(fname_out)

        var_name = variable + "_" + lon_lat
        print(var_name)

        da_var = ds_sub[var_name].copy()
        ds_var = da_var.rename(variable).squeeze().to_dataset()

        ds_var["lat"] = xr.DataArray(
            data=ds_sub["lat_" + lon_lat].values[0],
            attrs=dict(
                units="degrees_north",
                long_name="latitude",
            ),
        )
        ds_var["lon"] = xr.DataArray(
            data=ds_sub["lon_" + lon_lat].values[0],
            attrs=dict(
                units="degrees_east",
                long_name="longitude",
            ),
        )

        if "lev" in ds_sub[var_name].coords:
            ds_var["PS"] = ds_sub["PS" + "_" + lon_lat].squeeze()
            ds_var["P0"] = ds_sub["P0"]
            ds_var["hyam"] = ds_sub["hyam"]
            ds_var["hybm"] = ds_sub["hybm"]

        datasets[fname_out] = ds_var

for k, v in datasets.items():
    v.to_netcdf(k)
    v.close()

    print(f"Saved {k}")

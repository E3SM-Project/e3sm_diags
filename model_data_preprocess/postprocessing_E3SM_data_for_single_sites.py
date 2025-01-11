from pathlib import Path
import glob
import xarray as xr

# subprocess.run('source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_cori-haswell.sh', shell=True)
# A script to convert high frequency single point E3SM output to per-variable per-site netcdf files as input for ARM diagostics.
# In this example 3 hourly output at ARM sites are saved on h4 tape using namelist as follows:
# fincl2 = 'PS', 'Q', 'T', 'Z3', 'CLOUD', 'CONCLD', 'CLDICE', 'CLDLIQ', 'LS_FLXPRC', 'LS_FLXSNW', 'ZMFLXPRC', 'ZMFLXSNW', 'FREQR', 'REI', 'REL', 'CV_REFFICE', 'CV_REFFLIQ', 'LS_REFFRAIN', 'LS_REFFSNOW', 'PRECT', 'TMQ', 'PRECC', 'TREFHT', 'QREFHT', 'OMEGA','CLDTOT', 'LHFLX', 'SHFLX', 'FLDS', 'FSDS', 'FLNS', 'FSNS', 'FLNSC', 'FSDSC', 'FSNSC', 'AODVIS', 'AODABS'
# fincl2lonlat = '262.5e_36.6n','203.4e_71.3n','147.4e_2.0s','166.9e_0.5s','130.9e_12.4s','331.97e_39.09n'
data_path = "/global/cfs/cdirs/e3sm/www/Tutorials/2024/simulations/extendedOutput.v3.LR.historical_0101/archive/atm/hist/"
out_path = "/global/cfs/cdirs/e3sm/chengzhu/tutorial2024/v3.LR.historical_0101/post/atm/site/"
Path(out_path).mkdir(parents=True, exist_ok=True)
h_num = "h6"
start="2000"
end="2014"
time_range = f"{start}01_{end}12"
print(data_path,h_num)

ds = xr.open_mfdataset(sorted(glob.glob(data_path + f'*{h_num}.20*'))).sel(
    time=slice(f'{start}-01-01', f'{end}-12-31'))
print(ds)
print("Create time-series (it can take long depends on number of years being analyzed)")

variables = [
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

sites_info = {
    "sgpc1": [262.5, 36.6],
    "nsac1": [203.4, 71.3],
    "twpc1": [147.4, -2.0],
    "twpc2": [166.9, -0.5],
    "twpc3": [130.9, -12.4],
    "enac1": [331.97, 39.09],
}
sites = ["sgpc1", "nsac1", "twpc1", "twpc2", "twpc3", "enac1"]

for site in sites:
    if sites_info[site][1] > 0:
        lon_lat = str(sites_info[site][0]) + "e_" + str(sites_info[site][1]) + "n"
    else:
        lon_lat = str(sites_info[site][0]) + "e_" + str(abs(sites_info[site][1])) + "s"

    for variable in variables:
   
        fname_out = f"{out_path}{variable}_{site}_{time_range}.nc"
        print(fname_out)
        var_name = variable + "_" + lon_lat
        print(var_name)
        ds_new = ds[var_name].rename(variable).squeeze().to_dataset()
        lat = xr.DataArray(data = ds["lat_" + lon_lat].values[0],attrs=dict(units = "degrees_north", long_name = "latitude",),)
        lon = xr.DataArray(data = ds["lon_" + lon_lat].values[0],attrs=dict(units = "degrees_east", long_name = "longitude",),)
        ds_new["lat"] = lat
        ds_new["lon"] = lon
        if 'lev' in ds[var_name].coords:
            ds_new[
                'PS'
            ] = ds["PS" + "_" + lon_lat].squeeze()
            ds_new['P0'] = ds['P0']
            ds_new['hyam'] = ds['hyam']
            ds_new['hybm'] = ds['hybm']
        ds_new.to_netcdf(fname_out)
        ds_new.close()

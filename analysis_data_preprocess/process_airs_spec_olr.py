import numpy as np
import xarray as xr
import os
from pcmdi_metrics.io import xcdat_open

"""
This script converts AIRS spectral OLR to band OLR (E3SM band02 and band06) 
and saves monthly (Broadband, band02, and band06), seasonal and annual mean climatology.

AIRS related URL: https://cmr.earthdata.nasa.gov/search/concepts/C1697372449-GES_DISC.html."

Data and scripts provided by Professor Xianglei Huang's group at UMICH
"""

# settings
airs_ceres_path = '/global/cfs/cdirs/e3sm/zhang40/analysis_data_e3sm_diags/AIRS_specOLR/original/AIRS2CERES/'
airs_modis_path = '/global/cfs/cdirs/e3sm/zhang40/analysis_data_e3sm_diags/AIRS_specOLR/original/AIRS2MODIS/'
output_path = "/global/cfs/cdirs/e3sm/zhang40/analysis_data_e3sm_diags/AIRS_specOLR/climatology/"
output_path_time_series = "/global/cfs/cdirs/e3sm/zhang40/analysis_data_e3sm_diags/AIRS_specOLR/time_series/"
start_yr = 2003
end_yr = 2021
bands = ["broadband", "band02", "band06"]
output_filename = "AIRS_specOLR_"+str(start_yr)+"01_"+str(end_yr)+"12.nc"

# E3SM band information
wnum_edge = {"band02":[350,500], "band06":[820,980]}

# AIRS specOLR fill value
fill_value = 9.96921e+36

# make input file lists
airs_ceres_files = os.listdir(airs_ceres_path)
airs_modis_files = os.listdir(airs_modis_path)


for ib, band in enumerate(bands):
    for im in range(1,13):
        for iy in range(start_yr, end_yr+1):
            date_sel = "."+str(iy)+str(im).zfill(2)+"01."
            airs_ceres_file = [item for item in airs_ceres_files if date_sel in item]
            airs_modis_file = [item for item in airs_modis_files if date_sel in item]

            # read data
            airs_ceres_data = xr.open_dataset(airs_ceres_path+airs_ceres_file[0])
            airs_modis_data = xr.open_dataset(airs_modis_path+airs_modis_file[0])

            if band == "broadband":
                olr = airs_ceres_data.olr
                olr_clr = airs_modis_data.olr_clr
            else:
                olr = airs_ceres_data.olr_spectral
                olr_clr = airs_modis_data.olr_clr_spectral
                # select wnum
                olr = olr.sel(wnum=slice(wnum_edge[band][0],wnum_edge[band][1]))
                olr_clr = olr_clr.sel(wnum=slice(wnum_edge[band][0],wnum_edge[band][1]))
                
            # fill number to nan
            olr = xr.where(olr == fill_value, np.nan, olr)
            olr_clr = xr.where(olr_clr == fill_value, np.nan, olr_clr)
            
            # average along orbit_pass
            olr = olr.mean(dim="orbit_pass")
            olr_clr = olr_clr.mean(dim="orbit_pass")
          
            if iy == start_yr:
                monthly_olr = olr
                monthly_olr_clr = olr_clr
            else:
                monthly_olr = xr.concat([monthly_olr, olr], "time")
                monthly_olr_clr = xr.concat([monthly_olr_clr, olr_clr], "time")
 
        if im == 1:
            output_olr = monthly_olr.mean(dim="time")
            output_olr_clr = monthly_olr_clr.mean(dim="time")
        else:
            output_olr = xr.concat([output_olr, monthly_olr.mean(dim="time")], "time")
            output_olr_clr = xr.concat([output_olr_clr, monthly_olr_clr.mean(dim="time")], "time")
            
    if band != "broadband":
        # sum along wnum 
        output_olr = output_olr.sum(dim="wnum", skipna=False)
        output_olr_clr = output_olr_clr.sum(dim="wnum", skipna=False)

        # rename variables
        output_olr = output_olr.rename("olr_"+band)
        output_olr_clr = output_olr_clr.rename("olr_clr_"+band)

    # make time series for output
    dates = xr.date_range(start="0001", periods=12, freq="MS", calendar="noleap", use_cftime=True)

    # add time corrdinate value
    output_olr = output_olr.assign_coords({"time":dates})
    output_olr_clr = output_olr_clr.assign_coords({"time":dates})

    # merge to single dataset
    if ib == 0:
        output_olr_vars = xr.merge([output_olr, output_olr_clr])
    else:
        output_olr_vars = xr.merge([output_olr_vars, output_olr, output_olr_clr])

# save monthly climatology netcdf
output_olr_vars.to_netcdf(output_path_time_series+output_filename)

# Read the saved file using xcdat_open for proper time handling
airs_vars = xcdat_open(output_path_time_series+output_filename)
airs_vars = airs_vars.bounds.add_missing_bounds("T")
print("Loaded dataset for seasonal/annual calculations")

# Calculate seasonal and annual means
seasons = ["ANN", "DJF", "MAM", "JJA", "SON"]

for season in seasons:
    season_output = None

    for band in bands:
        if band == "broadband":
            airs_varname_as = "olr"
            airs_varname_cs = "olr_clr"
        else:
            airs_varname_as = "olr_"+band
            airs_varname_cs = "olr_clr_"+band

        # Calculate seasonal or annual mean values
        if season == "ANN":
            airs_var_as = airs_vars.temporal.average(airs_varname_as, weighted=True)
            airs_var_cs = airs_vars.temporal.average(airs_varname_cs, weighted=True)
        else:
            airs_var_as = airs_vars.temporal.climatology(
                airs_varname_as,
                freq="season",
                weighted=True,
                season_config={"dec_mode": "DJF"}
            )
            airs_var_cs = airs_vars.temporal.climatology(
                airs_varname_cs,
                freq="season",
                weighted=True,
                season_config={"dec_mode": "DJF"}
            )

            # Select the appropriate season
            if season == "DJF":
                airs_var_as = airs_var_as.isel(time=0)
                airs_var_cs = airs_var_cs.isel(time=0)
            elif season == "MAM":
                airs_var_as = airs_var_as.isel(time=1)
                airs_var_cs = airs_var_cs.isel(time=1)
            elif season == "JJA":
                airs_var_as = airs_var_as.isel(time=2)
                airs_var_cs = airs_var_cs.isel(time=2)
            elif season == "SON":
                airs_var_as = airs_var_as.isel(time=3)
                airs_var_cs = airs_var_cs.isel(time=3)

        # Merge variables for this band
        if season_output is None:
            season_output = xr.merge([airs_var_as, airs_var_cs])
        else:
            season_output = xr.merge([season_output, airs_var_as, airs_var_cs])

    # Save seasonal/annual mean netcdf
    season_filename = f"AIRS_specOLR_{start_yr}01_{end_yr}12_{season}.nc"
    season_output.to_netcdf(output_path + season_filename)
    print(f"Saved {season_filename}")

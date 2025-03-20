#%%

import xcdat as xc

test = '/lcrc/group/e3sm/ac.forsyth2/zppy_weekly_comprehensive_v2_output/test_pr651_both_commits_20250117/v2.LR.historical_0201/post/atm/180x360_aave/clim/2yr/v2.LR.historical_0201_SON_198009_198111_climo.nc'
ref = '/lcrc/group/e3sm/diagnostics/observations/Atm/climatology/MERRA2/MERRA2_SON_198009_201611_climo.nc'

# Open without Dask and load (e3sm_diags method)
ds_a = xc.open_mfdataset(test).sel("time", slice("1980-09", "1981-11"))
ds_b = xc.open_mfdataset(ref)

ds_a.load()
ds_b.load()

# "T" is derived from "ta"
ds_b = ds_b.rename({"ta": "T"})

# Check if the data is contiguous
print(ds_a["T"].data.flags["C_CONTIGUOUS"]) # True
print(ds_b["T"].data.flags["C_CONTIGUOUS"]) # True

# Drop unused ilev dimension.
ds_a = ds_a.drop_dims("ilev")

#%%
output_grid = ds_a.regridder.grid
ds_b_regrid = ds_b.regridder.horizontal(
    "T", output_grid, tool="xesmf", method="conservative_normed"
)

#%%
# Open without Dask (just for testing)
ds_c = xc.open_dataset(test)
ds_d = xc.open_dataset(ref)

# "T" is derived from "ta"
ds_d = ds_d.rename({"ta": "T"})

# Check if the data is contiguous
print(ds_c["T"].data.flags["C_CONTIGUOUS"]) # True
print(ds_d["T"].data.flags["C_CONTIGUOUS"]) # True


#%%
output_grid = ds_a.regridder.grid
ds_b_regrid = ds_b.regridder.horizontal(
    "T", output_grid, tool="xesmf", method="conservative_normed"
)

# %%
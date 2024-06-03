# %%
import os

import cdms2
import numpy as np
import xarray as xr
import xcdat as xc  # noqa: F401

var_key = "TREFHT"
regrid_tool = "esmf"
regrid_method = "bilinear"

dir_path = "/global/u2/v/vo13/E3SM-Project/e3sm_diags/auxiliary_tools/cdat_regression_testing/759-slice-flag/debug/"
fp_a = os.path.join(dir_path, "ds_a.nc")
fp_b = os.path.join(dir_path, "ds_b.nc")
fp_mv1 = os.path.join(dir_path, "mv1.nc")
fp_mv2 = os.path.join(dir_path, "mv2.nc")


# %%
# Regridding with xCDAT
ds_a = xr.open_dataset(fp_a)
ds_b = xr.open_dataset(fp_b)

output_grid = ds_a.regridder.grid
ds_b_regrid = ds_b.regridder.horizontal(
    var_key, output_grid, tool=f"x{regrid_tool}", method=regrid_method
)

# Write out to netCDF for visualization
# ds_b_regrid.to_netcdf("ds_b_regrid.nc")

# %%
# Regridding with CDAT
ds_mv1 = cdms2.open(fp_mv1)
ds_mv2 = cdms2.open(fp_mv2)

mv1 = ds_mv1("variable_21")
mv2 = ds_mv2("variable_36")
mv_grid = mv1.getGrid()
mv2_reg = mv2.regrid(mv_grid, regridTool=regrid_tool, regridMethod=regrid_method)

# Write out to netCDF for visualization
# f1 = cdms2.open("mv2_regrid.nc", "w")
# f1.write(mv2)
# f1.close()

# %%
"""
1. Compare CDAT original data with regridded data

Result: There is no difference between the original and regridded reference
data.
"""
# Test closeness (same)
np.testing.assert_allclose(mv2, mv2_reg, atol=0, rtol=0)

mv2_filled = mv2.filled(np.nan)
mv2_reg_filled = mv2_reg.filled(np.nan)

#  Count of np.nan: 22170 vs. 22170 (same)
np.count_nonzero(np.isnan(mv2_filled.flatten())) == np.count_nonzero(
    np.isnan(mv2_reg_filled.flatten())
)
# Sum: -68837.01023016105 vs. -68837.01023016105 (same)
np.nansum(mv2_filled) == np.nansum(mv2_reg_filled)
# Mean: -6.342086809486 == -6.342086809486 (same)
np.nanmean(mv2_filled) == np.nanmean(mv2_reg_filled)

# %%
"""
2. Compare regridded data between xCDAT and CDAT

Result: The regridded reference data produced by xCDAT results in more `np.nan`,
which subsequently results in different np.nan locations, sum, and mean.

In https://github.com/E3SM-Project/e3sm_diags/pull/794, I noted that there
seems tobe a minor difference in how the bilinear regridding works with xCDAT + xESMF vs.
CDAT + ESMF.
"""
# Test closeness (x and y nan location mismatch)
np.testing.assert_allclose(ds_b_regrid[var_key].values, mv2_reg, atol=0, rtol=0)

# Count of np.nan: 23150 vs. 22170 (different)
np.count_nonzero(np.isnan(ds_b_regrid[var_key].values.flatten())) == np.count_nonzero(
    np.isnan(mv2_reg_filled.flatten())
)
# Sum: -77674.06388742107 vs. -68837.01023016105 (different)
np.nansum(ds_b_regrid[var_key].values) == np.nansum(mv2_reg_filled.data)
# Mean: -7.866524598685545 vs. -6.342086809486 (different)
np.nanmean(ds_b_regrid[var_key].values) == np.nanmean(mv2_reg_filled.data)

# %%

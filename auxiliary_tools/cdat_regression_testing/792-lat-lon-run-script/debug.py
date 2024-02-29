import xarray as xr

from e3sm_diags.derivations.derivations import DERIVED_VARIABLES

ds1 = xr.open_dataset(
    "/global/cfs/cdirs/e3sm/www/cdat-migration-fy24/792-lat-lon/lat_lon/CERES-EBAF-TOA-v4.1/ceres_ebaf_toa_v4.1-FLUTC-JJA-global_ref.nc"
)
ds2 = xr.open_dataset(
    "/global/cfs/cdirs/e3sm/www/cdat-migration-fy24/main/lat_lon/CERES-EBAF-TOA-v4.1/ceres_ebaf_toa_v4.1-FLUTC-JJA-global_ref.nc"
)

DERIVED_VARIABLES["FLUTC"]

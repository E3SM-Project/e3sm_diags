# %%
import xcdat as xc

fp1 = "/global/cfs/cdirs/e3sm/www/cdat-migration-fy24/666-diurnal-cycle/diurnal_cycle/TRMM-3B43v-7_3hr/TRMM-3B43v-7_3hr-PRECT-ANN-20S20N_test.nc"


fp2 = "/global/cfs/cdirs/e3sm/www/cdat-migration-fy24/666-diurnal-cycle/diurnal_cycle/TRMM-3B43v-7_3hr/TRMM-3B43v-7_3hr-PRECT-ANN-20S20N_test.nc"


ds1 = xc.open_dataset(fp1)
ds2 = xc.open_dataset(fp2)

import numpy as np
import xcdat as xc

dev_path = "/global/cfs/cdirs/e3sm/www/cdat-migration-fy24/669-annual_cycle_zonal_mean-debug/annual_cycle_zonal_mean/AOD_550/AOD_550-AODVIS-ANNUALCYCLE-global_ref.nc"
main_path = "/global/cfs/cdirs/e3sm/www/cdat-migration-fy24/main/annual_cycle_zonal_mean/AOD_550/AOD_550-AODVIS-Annual-Cycle_test.nc"


var_a = xc.open_dataset(dev_path)["AODVIS"]
var_b = xc.open_dataset(main_path)["AODVIS"]

"""
Floating point comparison

AssertionError:
Not equal to tolerance rtol=1e-07, atol=0

Mismatched elements: 1808 / 2160 (83.7%)
Max absolute difference: 0.12250582
Max relative difference: 91.14554689
 x: array([[0., 0., 0., ..., 0., 0., 0.],
       [0., 0., 0., ..., 0., 0., 0.],
       [0., 0., 0., ..., 0., 0., 0.],...
 y: array([[0., 0., 0., ..., 0., 0., 0.],
       [0., 0., 0., ..., 0., 0., 0.],
       [0., 0., 0., ..., 0., 0., 0.],...
"""
np.testing.assert_allclose(var_a, var_b)

# Get the max of all values
# -------------------------
# 0.28664299845695496
print(var_a.max().item())
# 0.2866430557436412
print(var_b.max().item())

# Get the min of all values
# -------------------------
# 0.0
print(var_a.min().item())
# 0.0
print(var_b.min().item())

# Get the sum of all values
# -------------------------
# 224.2569122314453
print(var_a.sum().item())
# 224.25691348856003
print(var_b.sum().item())

# Get the mean of all values
# -------------------------
# 0.10382264107465744
print(var_a.mean().item())
# 0.1038226451335926
print(var_b.mean().item())


# %%
# Get the max absolute diff
# -------------------------
# 0.12250582128763199
print((var_a - var_b).max().item())

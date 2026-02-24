#!/usr/bin/env python
"""Inspect the structure of EAMxx COSP histogram data."""

import xarray as xr

data_file = "/pscratch/sd/t/terai/EAMxx/ne256pg2_ne256pg2.F20TR-SCREAMv1.July-1.spanc800.2xauto.acc150.n0032.test2.1/rgr/climo/1ma_ne30pg2.AVERAGE.nmonths_x1_ANN_199501_200412_climo.nc"

print("=" * 70)
print("Inspecting EAMxx COSP Data Structure")
print("=" * 70)

ds = xr.open_dataset(data_file)

# Check for COSP histogram variables
cosp_vars = ["isccp_ctptau", "modis_ctptau", "misr_cthtau"]

for var_name in cosp_vars:
    print(f"\n{var_name}:")
    print("-" * 70)

    if var_name in ds:
        var = ds[var_name]
        print(f"  Dimensions: {var.dims}")
        print(f"  Shape: {var.shape}")
        print(f"  Coords:")
        for dim in var.dims:
            if dim in ds.coords:
                coord = ds[dim]
                print(f"    {dim}: shape={coord.shape}, values={coord.values if coord.size <= 10 else f'[{coord.size} values]'}")
            else:
                print(f"    {dim}: NOT A COORDINATE (dimension only)")
    else:
        print(f"  NOT FOUND in dataset")

# List all dimension coordinates
print("\n" + "=" * 70)
print("All COSP-related dimensions and coordinates:")
print("-" * 70)
for name in ds.dims:
    if 'cosp' in name.lower() or 'misr' in name.lower() or 'isccp' in name.lower():
        size = ds.dims[name]
        is_coord = name in ds.coords
        print(f"  {name}: size={size}, is_coordinate={is_coord}")
        if is_coord:
            print(f"    values: {ds[name].values}")

ds.close()

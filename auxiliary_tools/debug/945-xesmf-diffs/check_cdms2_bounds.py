"""
This script demonstrates the behavior of the `cdms2` library in handling axis bounds,
specifically for latitude axes. It includes tests to verify the automatic generation
and preservation of bounds under different scenarios. The findings are summarized below:

1. **Automatic Bounds Generation**:
    - When a latitude axis is created without explicitly setting bounds,
      `cdms2` automatically generates bounds for the axis.

2. **Preservation of Existing Bounds**:
    - If bounds are explicitly set for a latitude axis, `cdms2` preserves
      these bounds without overwriting them.

3. **Bounds Preservation in NetCDF Files**:
    - When a latitude axis with explicitly set bounds is written to a NetCDF file
      and subsequently read back, the bounds are preserved as expected.

The script includes three tests:
- **Test 1**: Verifies that bounds are auto-generated for a latitude axis.
- **Test 2**: Confirms that explicitly set bounds for a latitude axis are preserved.
- **Test 3**: Ensures that bounds for a latitude axis are preserved when written
  to and read from a NetCDF file.

cdms2.setAutoBounds(mode) controls automatic generation of bounds for axes:
0 or 'off'   = No bounds are generated automatically.
1 or 'on'    = Bounds are auto-generated for all 1D axes (if missing).
2 or 'grid'  = Bounds are auto-generated only for lat/lon axes (default).
Use cdms2.getAutoBounds() to check the current setting.
"""

# %%x
import cdms2
import numpy as np

# Check the current setting for auto bounds generation
print(cdms2.getAutoBounds())
# Expected output: 2 (default setting for lat/lon axes)

# %%
# Test 1: Check if lat bounds are auto-generated with no bounds existing (True)
# -----------------------------------------------------------------------------
axis = cdms2.createAxis([0, 1, 2])
axis.id = "lat"  # Identify the axis as latitude

print(axis.getBounds())  # These are the auto-generated ones, not manually set.
# [[0.  0.5]
#  [0.5 1.5]
#  [1.5 2.5]]

# %%
# Test 2: Check if existing lat bounds are preserved (True)
# -----------------------------------------------------------------------------
axis = cdms2.createAxis([0, 1, 2])
axis.id = "lat"  # Identify the axis as latitude
axis.setBounds(np.array([[0, 0.2], [0.5, 1.5], [1.5, 2.5]]))

print(axis.getBounds())  # These are the manually set ones.
# [[0.  0.2]
#  [0.5 1.5]
#  [1.5 2.5]]

# %%
# Test 3: Create a dummy netCDF file with lat axis, read it, and check bounds are preserved (True)
# -----------------------------------------------------------------------------
# Create a dummy netCDF file
with cdms2.open("dummy.nc", "w") as f:
    lat = cdms2.createAxis([0, 1, 2])
    lat.id = "lat"
    lat.setBounds(np.array([[0, 0.2], [0.5, 1.5], [1.5, 2.5]]))

    var = cdms2.createVariable([0, 1, 2], id="var", axes=[lat])
    f.write(var)

obj = cdms2.open("dummy.nc")["var"]

# %%
print(obj.getLatitude().getBounds())  # These are the manually set ones.
# [[0.  0.2]
#  [0.5 1.5]
#  [1.5 2.5]]

# %%

# Source: https://github.com/pydata/xarray/issues/11205

import xarray as xr
xr.show_versions()
# your reproducer code ...

import numpy

size = 512  # 511 works, 512 segfaults
data = numpy.random.rand(size, size).astype(numpy.float32)
data[0:5, 0:5] = 65534.0

da = xr.DataArray(data, dims=["rows", "columns"])
da.encoding["_FillValue"] = numpy.float32(65534.0)
da.to_netcdf("/tmp/test_fill.nc")

ds = xr.open_dataset("/tmp/test_fill.nc")
var = ds["__xarray_dataarray_variable__"]
print(var.values)  # segfaults on Python 3.14 Linux
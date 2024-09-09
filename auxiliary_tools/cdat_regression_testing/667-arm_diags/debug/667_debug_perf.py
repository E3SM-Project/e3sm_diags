# %%
import timeit

setup_code = """
import xarray as xr

AIR_DENS = 1.225  # standard air density 1.225kg/m3

a1 = xr.open_dataarray("qa/667-arms-diags/a1.nc")
a2 = xr.open_dataarray("qa/667-arms-diags/a2.nc")
a3 = xr.open_dataarray("qa/667-arms-diags/a3.nc")
"""

setup_code2 = """
import xarray as xr

AIR_DENS = 1.225  # standard air density 1.225kg/m3

a1 = xr.open_dataarray("qa/667-arms-diags/a1.nc")
a2 = xr.open_dataarray("qa/667-arms-diags/a2.nc")
a3 = xr.open_dataarray("qa/667-arms-diags/a3.nc")

a1.load(scheduler="sync")
a2.load(scheduler="sync")
a3.load(scheduler="sync")
"""

setup_code3 = """
import xarray as xr

AIR_DENS = 1.225  # standard air density 1.225kg/m3

a1_chunked = xr.open_dataarray("qa/667-arms-diags/a1.nc", chunks={"time": "auto"})
a2_chunked = xr.open_dataarray("qa/667-arms-diags/a2.nc", chunks={"time": "auto"})
a3_chunked = xr.open_dataarray("qa/667-arms-diags/a3.nc", chunks={"time": "auto"})
"""

code_statement1 = """
with xr.set_options(keep_attrs=True):
    var = (a1 + a2 + a3) * AIR_DENS / 1e6
"""


code_statement2 = """
with xr.set_options(keep_attrs=True):
    var = (a1_chunked + a2_chunked + a3_chunked) * AIR_DENS / 1e6
"""

code_statement3 = """
var_data = (a1.values + a2.values + a3.values) * AIR_DENS / 1e6
var_new = xr.DataArray(
    var_data,
    dims=a1.dims,
    coords=a1.coords,
    name="a_num",
    attrs={"units": "/cm3", "long_name": "aerosol number concentration"},
)
"""

code_statement4 = """
var_data2 = (a1.data + a2.data + a3.data) * AIR_DENS / 1e6
var_new2 = xr.DataArray(
    name="a_num", data=var_data2, dims=a1.dims, coords=a1.coords, attrs=a1.attrs
)
var_new2.attrs.update(
    {"units": "/cm3", "long_name": "aerosol number concentration"}
)
"""


def run_timeit(code_statement: str, setup_code: str) -> float:
    elapsed_time = timeit.repeat(
        code_statement, setup=setup_code, globals=globals(), repeat=3, number=1
    )

    return min(elapsed_time)


elapsed_time_xarray = run_timeit(code_statement1, setup_code)
print(f"1. Elapsed time (Xarray non-chunked): {elapsed_time_xarray} seconds")

elapsed_time_xarray_load = run_timeit(code_statement1, setup_code2)
print(
    f"2. Elapsed time (Xarray non-chunked with .load()): {elapsed_time_xarray_load} seconds"
)
elapsed_time_xarray_chunked = run_timeit(code_statement2, setup_code3)
print(f"3. Elapsed time (Xarray chunked): {elapsed_time_xarray_chunked} seconds")

elapsed_time_numpy_1 = run_timeit(code_statement3, setup_code)
print(f"4. Elapsed time (numpy .values): {elapsed_time_numpy_1} seconds")

elapsed_time_numpy_2 = run_timeit(code_statement4, setup_code)
print(f"5. Elapsed time (numpy .data): {elapsed_time_numpy_2} seconds")


"""
Results
----------
1. Elapsed time (Xarray non-chunked): 6.540755605790764 seconds
2. Elapsed time (Xarray non-chunked with .load()): 0.17097265785560012 seconds
3. Elapsed time (Xarray chunked): 0.1452920027077198 seconds
4. Elapsed time (numpy .values): 6.418793010059744 seconds
5. Elapsed time (numpy .data): 7.334999438840896 seconds
"""

# %%
import numpy as np
import xarray as xr
from PIL import Image, ImageChops, ImageDraw

ds1 = xr.open_dataset(
    "/global/cfs/cdirs/e3sm/www/cdat-migration-fy24/792-lat-lon/lat_lon/CRU_IPCC/CRU-TREFHT-ANN-land_60S90N_ref.nc"
)
ds2 = xr.open_dataset(
    "/global/cfs/cdirs/e3sm/www/cdat-migration-fy24/main/lat_lon/CRU_IPCC/CRU-TREFHT-ANN-land_60S90N_ref.nc"
)

var_key = "TREFHT"
# %%
np.testing.assert_allclose(ds1[var_key], ds2[var_key])

# %%
# Check the sum values -- close
# 124810.94
ds1[var_key].sum()
# 119496.11
ds2[var_key].sum()

# %%
# Check the mean values -- close
# 8.036763
ds1[var_key].mean()

# 8.352283
ds2[var_key].mean()

# Check the plots
ax1 = ds1[var_key].plot()
ax1.figure.savefig("plot1.png")

# %%
ax2 = ds2[var_key].plot()
ax2.figure.savefig("plot2.png")

# %%
actual_png = Image.open("plot1.png").convert("RGB")
expected_png = Image.open("plot2.png").convert("RGB")
diff = ImageChops.difference(actual_png, expected_png)

# %%
bbox = diff.getbbox()
nonzero_pixels = (
    diff.crop(bbox).point(lambda x: 255 if x else 0).convert("L").point(bool).getdata()
)

draw = ImageDraw.Draw(diff)
(left, upper, right, lower) = diff.getbbox()
draw.rectangle(((left, upper), (right, lower)), outline="red")
diff.save("diff.png")

# %%

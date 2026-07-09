"""A minimum, self-contained example demonstrating Cartopy contourf() plotting
differences between cartopy 0.24.0 and 0.25.0.

Issue: The contourf() plot appears inverted (flipped vertically) on the Y-axis
when using cartopy 0.25.0 compared to cartopy 0.24.0

Cause:
- Commit https://github.com/SciTools/cartopy/commit/eb4d0423ee9940067fa18c71cb38e178471ed888
introduces changes to "FIX: create a single inverted polygon when no exteriors found".
- This commit is where the issue first appears (11/4/2024).

Variables:
- x: longitude (0.5 o 360.5)
- y: latitude (-89.5 to 89.5)
- var: variable on (y, x) grid with masked values present at the bottom edge.

Workaround:
- Converting x and y to 2D meshgrid using numpy.meshgrid() and using transform_first=True.
- However, this workaround causes issues with many other plots.

Setup Environment 1 (Without Plot Issues):
----------------------------------------------------------
- Commit: https://github.com/SciTools/cartopy/commit/5c4ef99f7093f43e387b975391dd77179a8df9ac
- 11/3/2024: DOC: add deprecation note for clip_path [skip actions]
- Hash: 5c4ef99f7093f43e387b975391dd77179a8df9ac

Commands:
    git clone https://github.com/SciTools/cartopy.git

    conda create -n cartopy_5c4e4f -c conda-forge xarray=2025.12.0 netcdf4 matplotlib-base=3.10.8 ipykernel
    conda activate cartopy_5c4e4f

    cd cartopy
    git checkout 5c4ef99f7093f43e387b975391dd77179a8df9ac
    python -m pip install .


Setup Environment 2 (With Plot Issues):
----------------------------------------------------------
- Commit: https://github.com/SciTools/cartopy/commit/eb4d0423ee9940067fa18c71cb38e178471ed888
- 11/4/2024: FIX: create a single inverted polygon when no exteriors found
- Hash: eb4d0423ee9940067fa18c71cb38e178471ed888

Commands:
    git clone https://github.com/SciTools/cartopy.git

    conda create -n cartopy_eb4d04 -c conda-forge xarray=2025.12.0 netcdf4 matplotlib-base=3.10.8 ipykernel
    conda activate cartopy_eb4d04

    cd cartopy
    git checkout eb4d0423ee9940067fa18c71cb38e178471ed888
    python -m pip install .
"""
import urllib.request
import os

import cartopy
import matplotlib
import numpy as np
import xarray as xr
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # isort:skip  # noqa: E402

# Download the data
# --------------------------------------------------------------------------
# Source: /lcrc/group/e3sm/public_html/diagnostic_output/ac.tvo/tests/1026_cartopy_0250_mvce/
public_path = "https://web.lcrc.anl.gov/public/e3sm/diagnostic_output/ac.tvo/tests/1026_cartopy_0250_mvce/"
local_data_dir = "./data_cartopy_0250_mvce"
os.makedirs(local_data_dir, exist_ok=True)

def download_data():
    files = ["x.nc", "y.nc", "var.nc"]

    for fname in files:
        url = public_path + fname
        local_path = os.path.join(local_data_dir, fname)

        if not os.path.exists(local_path):
            print(f"Downloading {url} to {local_path} ...")
            # Set a User-Agent header to avoid HTTP 403 errors
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req) as response, open(local_path, 'wb') as out_file:
                out_file.write(response.read())
        else:
            print(f"File {local_path} already exists.")

download_data()

# Open the data
# --------------------------------------------------------------------------
x = xr.open_dataarray(os.path.join(local_data_dir, "x.nc"))
y = xr.open_dataarray(os.path.join(local_data_dir, "y.nc"))
var = xr.open_dataarray(os.path.join(local_data_dir, "var.nc"))

# Create the figure and axis with Cartopy projection
# --------------------------------------------------------------------------
fig = plt.figure(figsize=(8.5, 11.0), dpi=150)
projection = cartopy.crs.PlateCarree(central_longitude=180)
ax = fig.add_axes((0.1691, 0.1112, 0.6465, 0.2258), projection=projection)

# Set the longitude and latitude limits.
lon_west = 0
lon_east=360
lat_south = -90
lat_north = 90

ax.set_extent([lon_west, lon_east, lat_south, lat_north], crs=projection)

# Create the contour plot
# --------------------------------------------------------------------------
transform = cartopy.crs.PlateCarree()
contour_plot = ax.contourf(
    x,
    y,
    var,
    cmap="BrBG",
    transform=transform,
    norm=None,
    levels=None,
    extend="both",
)

# Configure aspect ratio and coast lines.
# --------------------------------------------------------------------------
ax.set_aspect((lon_east - lon_west) / (2 * (lat_north - lat_south)))
ax.coastlines(lw=0.3)

# Configure x and y axes:
# --------------------------------------------------------------------------
x_ticks = np.array([  0. ,  60. , 120. , 180. , 240. , 300. , 359.5])
y_ticks = np.array([-90, -60, -30,   0,  30,  60,  90])
ax.set_xticks(x_ticks, crs=transform)
ax.set_yticks(y_ticks, crs=transform)
ax.tick_params(labelsize=8.0, direction="out", width=1)
ax.tick_params(labelsize=8.0, direction="out", width=1)
ax.xaxis.set_ticks_position("bottom")
ax.yaxis.set_ticks_position("left")

# Add longitude and latitude formatters
lon_formatter = LongitudeFormatter(
    zero_direction_label=True, number_format=".0f"
)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)

# Add color bar
# --------------------------------------------------------------------------
cbax_rect = (0.8326, 0.13269999999999998, 0.0326, 0.1792)
cbax = fig.add_axes(cbax_rect)
cbar = fig.colorbar(contour_plot, cax=cbax)
cbar.ax.tick_params(labelsize=9.0, length=0)

# Save the figure
# --------------------------------------------------------------------------
# Get the name of the current conda environment
conda_env = os.environ.get("CONDA_DEFAULT_ENV", "unknown_env")

output_path = os.path.abspath(f"/lcrc/group/e3sm/public_html/diagnostic_output/ac.tvo/tests/1026_cartopy_0250_mvce/{conda_env}_mvce.png")
fig.savefig(output_path)

print(f"Saved figure at {output_path}")

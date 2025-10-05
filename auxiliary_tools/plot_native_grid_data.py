import uxarray as ux
import matplotlib.pyplot as pl
import cartopy.crs as ccrs
import cartopy.feature as cfeature

base_path = "/Users/zhang40/Documents/ACME_simulations/E3SM_v2/native_grid_data/"
grid_info = "ne30pg2"
grid_path = base_path + f"{grid_info}.nc"
data_path = base_path + f"PRECC.{grid_info}.nc"

uxds = ux.open_dataset(grid_path, data_path)

pc = uxds["PRECT"].to_polycollection()
#pc = uxds["PRECT"].to_polycollection(periodic_elements="split")

# disables grid lines
pc.set_antialiased(False)

pc.set_cmap("plasma")

fig, ax = plt.subplots(
    1,
    1,
    figsize=(10, 5),
    facecolor="w",
    constrained_layout=True,
    subplot_kw=dict(projection=ccrs.PlateCarree()),
)

ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS)

ax.add_collection(pc)
ax.set_global()
plt.title("PolyCollection Plot with Projection & Features")

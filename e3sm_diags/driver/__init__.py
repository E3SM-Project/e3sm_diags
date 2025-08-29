import os

from e3sm_diags import INSTALL_PATH

# The path to the land ocean mask file, which is bundled with the installation
# of e3sm_diags in the conda environment.
LAND_OCEAN_MASK_PATH = os.path.join(INSTALL_PATH, "acme_ne30_ocean_land_mask.nc")

# The keys for the land and ocean fraction variables in the
# `LAND_OCEAN_MASK_PATH` file.
FRAC_REGION_KEYS = {"land": ("LANDFRAC", "landfrac"), "ocean": ("OCNFRAC", "ocnfrac")}

# The default value for metrics if it is not calculated. This value was
# preserved from the legacy CDAT codebase because the viewer expects this
# value for metrics that aren't calculated.
# TODO: Update `lat_lon_viewer.py` to handle missing metrics with None value.
METRICS_DEFAULT_VALUE = 999.999

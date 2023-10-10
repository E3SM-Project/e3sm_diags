import os

from e3sm_diags import INSTALL_PATH

# The path to the land ocean mask file, which is bundled with the installation
# of e3sm_diags in the conda environment.
LAND_OCEAN_MASK_PATH = os.path.join(INSTALL_PATH, "acme_ne30_ocean_land_mask.nc")

# The keys for the land and ocean fraction variables in the
# `LAND_OCEAN_MASK_PATH` file.
LAND_FRAC_KEY = "LANDFRAC"
OCEAN_FRAC_KEY = "OCNFRAC"

MASK_REGION_TO_VAR_KEY = {"land": LAND_FRAC_KEY, "ocean": OCEAN_FRAC_KEY}

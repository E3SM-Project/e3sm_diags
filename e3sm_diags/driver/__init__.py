import os

from e3sm_diags import INSTALL_PATH

# The path to the land ocean mask file, which is bundled with the installation
# of e3sm_diags in the conda environment.
LAND_OCEAN_MASK_PATH = os.path.join(INSTALL_PATH, "acme_ne30_ocean_land_mask.nc")

# The keys for the land and ocean fraction variables in the
# `LAND_OCEAN_MASK_PATH` file.
LAND_FRAC_KEY = "LANDFRAC"
OCEAN_FRAC_KEY = "OCNFRAC"


def _get_region_mask_var_key(region: str):
    """Get the region's mask variable key.

    This variable key can be used to map the the variable data in a sdataset.
    Only land and ocean regions are supported.

    Parameters
    ----------
    region : str
        The region.

    Returns
    -------
    str
       The variable key, either "LANDFRAC" or "OCNFRAC".

    Raises
    ------
    ValueError
        If the region passed is not land or ocean.
    """
    if "land" in region:
        return LAND_FRAC_KEY
    elif "ocean" in region:
        return OCEAN_FRAC_KEY

    raise ValueError(f"Only land and ocean regions are supported, not '{region}'.")

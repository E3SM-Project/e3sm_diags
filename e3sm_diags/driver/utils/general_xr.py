import cf_xarray as cfxr  # noqa: F401
import xarray as xr


def has_z_axis(data_var: xr.DataArray) -> bool:
    """Checks whether the data variable has a Z axis.

    Conditionals are based on:
    - https://cdms.readthedocs.io/en/latest/_modules/cdms2/avariable.html#AbstractVariable.getLevel
    - https://cdms.readthedocs.io/en/latest/_modules/cdms2/axis.html#AbstractAxis.isLevel

    Returns True if:
    - Data variable has a "Z" axis in the cf-xarray mapping dict
    - A coordinate in the data variable has a matching "positive" attribute or
      "name"
        # TODO: conditional for valid pressure units with "Pa"

    Parameters
    ----------
    data_var : xr.DataArray
        The data variable.

    Returns
    -------
    bool
        True if data variable has Z axis, else False.
    """

    if "Z" in data_var.cf.axes:
        return True

    for coord in data_var.coords.values():
        positive = coord.attrs.get("positive", "")
        positive = positive.strip()

        if positive in ["up", "down"]:
            return True

        if coord.name in ["lev", "plev", "depth"]:
            return True

    return False


def select_region(
    region: str,
    var: xr.DataArray,
    land_frac: xr.DataArray,
    ocean_frac: xr.DataArray,
    regrid_tool: str,
    regrid_method: str,
) -> xr.DataArray:
    """Select desired regions for the variable.

    Parameters
    ----------
    region : str.
        TODO: The region, options include "global"....
    var : xr.DataArray
        The variable.
    land_frac : xr.DataArray
        The land mask.
    ocean_frac : xr.DataArray
        The ocean mask.
    regrid_tool : str
        The regridding tool to use. Options include "xesmf" and "regrid2".
    regrid_method: str
        TODO: The regridding method to use. Options include: "conservation",
        "bilinear", ....

    Returns
    -------
    xr.DataArray
        The variable subsetted by region(s).
    """

    pass

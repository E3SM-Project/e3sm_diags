import numpy as np
import pytest
import xarray as xr

from e3sm_diags.derivations.utils import convert_units


class TestConvertUnitsCarbonFluxPerYear:
    @pytest.mark.parametrize(
        "units, factor",
        [
            ("gC/m^2/s", 365 * 86400),
            ("gC/m2/s", 365 * 86400),
            ("kgC/m^2/s", 1000 * 365 * 86400),
            ("kgC/m2/s", 1000 * 365 * 86400),
            ("gC/m^2/day", 365),
        ],
    )
    def test_converts_carbon_flux_to_gC_per_m2_per_year(self, units, factor):
        var = xr.DataArray(
            data=np.array([1.0e-8, 2.0e-8]),
            dims=["x"],
            attrs={"units": units},
            name="GPP",
        )

        result = convert_units(var, target_units="g*/m^2/yr")

        np.testing.assert_allclose(result.values, var.values * factor)
        assert result.attrs["units"] == "gC/m^2/yr"


class TestConvertUnitsAerosolFluxPerYear:
    @pytest.mark.parametrize("units", ["kg/m^2/s", "kg/m2/s"])
    def test_converts_aerosol_flux_to_g_per_m2_per_year(self, units):
        var = xr.DataArray(
            data=np.array([1.0e-12, 1.0e-8]),
            dims=["x"],
            attrs={"units": units},
            name="DSTDEP",
        )

        result = convert_units(var, target_units="g/m^2/yr")

        np.testing.assert_allclose(result.values, var.values * 1000 * 365 * 86400)
        assert result.attrs["units"] == "g/m^2/yr"

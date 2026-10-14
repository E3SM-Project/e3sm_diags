import numpy as np
import pytest
import xarray as xr

from e3sm_diags.derivations.formulas import qflxconvert_units


class TestQflxconvertUnits:
    @pytest.mark.parametrize("units", ["kg/m2/s", "kg m-2 s-1", "mm/s", "mm H2O/s"])
    def test_converts_per_second_water_flux_to_mm_per_day(self, units):
        var = xr.DataArray(
            data=np.array([1.0e-5, 2.0e-5]), dims=["x"], attrs={"units": units}
        )

        result = qflxconvert_units(var)

        np.testing.assert_allclose(result.values, [0.864, 1.728])
        assert result.attrs["units"] == "mm/day"

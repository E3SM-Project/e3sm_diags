import pytest
import xarray as xr


class TestClimo:
    @pytest.fixture(autouse=True)
    def setup(self):
        self.var = xr.DataArray()

    def returns_annual_cycle_climatology(self):
        assert 0

    def returns_seasonal_cycle_climatology(self):
        assert 0

    def returns_DJF_season_climatology(self):
        assert 0

    def returns_MAM_season_climatology(self):
        assert 0

    def returns_SON_season_climatology(self):
        assert 0

    def returns_climatology_for_derived_variable(self):
        assert 0

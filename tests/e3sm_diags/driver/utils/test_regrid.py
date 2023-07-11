import numpy as np
import pytest
import xarray as xr

from e3sm_diags.driver.utils.regrid import (
    convert_z_axis_to_pressure_levels,
    get_z_axis,
    has_z_axis,
)


class TestHasZAxis:
    def test_returns_true_if_data_array_has_have_z_axis(self):
        # Has Z axis
        z_axis1 = xr.DataArray(
            dims="height",
            data=np.array([0]),
            coords={"height": np.array([0])},
            attrs={"axis": "Z"},
        )
        dv1 = xr.DataArray(data=[0], coords=[z_axis1])

        dv_has_z_axis = has_z_axis(dv1)
        assert dv_has_z_axis

    def test_returns_true_if_data_array_has_z_coords_with_matching_positive_attr(self):
        # Has "positive" attribute equal to "up"
        z_axis1 = xr.DataArray(data=np.array([0]), attrs={"positive": "up"})
        dv1 = xr.DataArray(data=[0], coords=[z_axis1])

        dv_has_z_axis = has_z_axis(dv1)
        assert dv_has_z_axis

        # Has "positive" attribute equal to "down"
        z_axis2 = xr.DataArray(data=np.array([0]), attrs={"positive": "down"})
        dv2 = xr.DataArray(data=[0], coords=[z_axis2])

        dv_has_z_axis = has_z_axis(dv2)
        assert dv_has_z_axis

    def test_returns_true_if_data_array_has_z_coords_with_matching_name(self):
        # Has name equal to "lev"
        z_axis1 = xr.DataArray(name="lev", dims=["lev"], data=np.array([0]))
        dv1 = xr.DataArray(data=[0], coords={"lev": z_axis1})

        dv_has_z_axis = has_z_axis(dv1)
        assert dv_has_z_axis

        # Has name equal to "plev"
        z_axis2 = xr.DataArray(name="plev", dims=["plev"], data=np.array([0]))
        dv2 = xr.DataArray(data=[0], coords=[z_axis2])

        dv_has_z_axis = has_z_axis(dv2)
        assert dv_has_z_axis

        # Has name equal to "depth"
        z_axis3 = xr.DataArray(name="depth", dims=["depth"], data=np.array([0]))
        dv3 = xr.DataArray(data=[0], coords=[z_axis3])

        dv_has_z_axis = has_z_axis(dv3)
        assert dv_has_z_axis

    def test_raises_error_if_data_array_does_not_have_z_axis(self):
        dv1 = xr.DataArray(data=[0])

        dv_has_z_axis = has_z_axis(dv1)
        assert not dv_has_z_axis


class TestGetZAxis:
    def test_returns_true_if_data_array_has_have_z_axis(self):
        # Has Z axis
        z_axis1 = xr.DataArray(
            dims="height",
            data=np.array([0]),
            coords={"height": np.array([0])},
            attrs={"axis": "Z"},
        )
        dv1 = xr.DataArray(data=[0], coords=[z_axis1])

        result = get_z_axis(dv1)
        assert result.identical(dv1["height"])

    def test_returns_true_if_data_array_has_z_coords_with_matching_positive_attr(self):
        # Has "positive" attribute equal to "up"
        z_axis1 = xr.DataArray(data=np.array([0]), attrs={"positive": "up"})
        dv1 = xr.DataArray(data=[0], coords=[z_axis1])

        result1 = get_z_axis(dv1)
        assert result1.identical(dv1["dim_0"])

        # Has "positive" attribute equal to "down"
        z_axis2 = xr.DataArray(data=np.array([0]), attrs={"positive": "down"})
        dv2 = xr.DataArray(data=[0], coords=[z_axis2])

        result2 = get_z_axis(dv2)
        assert result2.identical(dv2["dim_0"])

    def test_returns_true_if_data_array_has_z_coords_with_matching_name(self):
        # Has name equal to "lev"
        z_axis1 = xr.DataArray(name="lev", dims=["lev"], data=np.array([0]))
        dv1 = xr.DataArray(data=[0], coords={"lev": z_axis1})

        result = get_z_axis(dv1)
        assert result.identical(dv1["lev"])

        # Has name equal to "plev"
        z_axis2 = xr.DataArray(name="plev", dims=["plev"], data=np.array([0]))
        dv2 = xr.DataArray(data=[0], coords=[z_axis2])

        result2 = get_z_axis(dv2)
        assert result2.identical(dv2["plev"])

        # Has name equal to "depth"
        z_axis3 = xr.DataArray(name="depth", dims=["depth"], data=np.array([0]))
        dv3 = xr.DataArray(data=[0], coords=[z_axis3])

        result = get_z_axis(dv3)
        assert result.identical(dv3["depth"])

    def test_raises_error_if_data_array_does_not_have_z_axis(self):
        dv1 = xr.DataArray(data=[0])

        with pytest.raises(KeyError):
            get_z_axis(dv1)


class TestConvertZAxisToPressureLevels:
    @pytest.fixture(autouse=True)
    def setup(self):
        z_axis = xr.DataArray(name="lev", dims=["lev"], data=np.array([0]))

        self.ds = xr.Dataset(
            data_vars={
                "lev": xr.DataArray(
                    data=[0], coords=[z_axis], attrs={"standard_name": "height"}
                )
            }
        )

    def test_raises_error_if_long_name_attr_is_None(self):
        ds = self.ds.copy()

        with pytest.raises(KeyError):
            convert_z_axis_to_pressure_levels(ds, ds["lev"], [1, 2, 3])

    def test_raises_error_if_long_name_attr_is_not_hybrid_or_pressure(self):
        ds = self.ds.copy()
        ds["lev"].attrs["long_name"] = "invalid"

        with pytest.raises(ValueError):
            convert_z_axis_to_pressure_levels(ds, ds["lev"], [1, 2, 3])

    def test_raises_error_if_dataset_does_not_contain_ps_hya_or_hyb_vars(self):
        ds = self.ds.copy()
        ds["lev"].attrs["long_name"] = "hybrid"

        with pytest.raises(KeyError):
            convert_z_axis_to_pressure_levels(ds, ds["lev"], [1, 2, 3])

    def test_converts_pressure_coordinates_to_pressure_levels(self):
        ds = self.ds.copy()
        ds["lev"].attrs["long_name"] = "pressure"

        # expected = xr.DataArray()
        # result = convert_z_axis_to_pressure_levels(ds, ds["lev"], [1, 2, 3])
        assert 0

        # ds["lev"].attrs["long_name"] = "isobaric"

    def test_converts_hybrid_levels_to_pressure_levels(self):
        ds = self.ds.copy()
        ds["lev"].attrs["long_name"] = "hybrid"
        assert 0

import numpy as np
import pytest
import xarray as xr
from xarray.testing import assert_identical

from e3sm_diags.driver.utils.regrid import (
    get_z_axis,
    has_z_axis,
    regrid_z_axis_to_plevs,
)
from tests.e3sm_diags.fixtures import generate_lev_dataset


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


class TestRegridZAxisToPlevs:
    @pytest.fixture(autouse=True)
    def setup(self):
        self.ds = generate_lev_dataset()
        self.plevs = [800, 200]

    def test_raises_error_if_long_name_attr_is_not_set(self):
        ds = generate_lev_dataset("hybrid")
        del ds["lev"].attrs["long_name"]

        with pytest.raises(KeyError):
            regrid_z_axis_to_plevs(ds, "so", self.plevs)

    def test_raises_error_if_long_name_attr_is_not_hybrid_or_pressure(self):
        ds = generate_lev_dataset("hybrid")
        ds["lev"].attrs["long_name"] = "invalid"

        with pytest.raises(ValueError):
            regrid_z_axis_to_plevs(ds, "so", self.plevs)

    def test_raises_error_if_dataset_does_not_contain_ps_hya_or_hyb_vars(self):
        ds = generate_lev_dataset("hybrid")
        ds = ds.drop_vars(["ps", "hyam", "hybm"])

        with pytest.raises(KeyError):
            regrid_z_axis_to_plevs(ds, "so", self.plevs)

    def test_raises_error_if_ps_variable_units_attr_is_None(self):
        ds = generate_lev_dataset("hybrid")
        ds.ps.attrs["units"] = None

        with pytest.raises(ValueError):
            regrid_z_axis_to_plevs(ds, "so", self.plevs)

    def test_raises_error_if_ps_variable_units_attr_is_not_mb_or_pa(self):
        ds = generate_lev_dataset("hybrid")
        ds.ps.attrs["units"] = "invalid"

        with pytest.raises(ValueError):
            regrid_z_axis_to_plevs(ds, "so", self.plevs)

    def test_regrids_hybrid_levels_to_pressure_levels(self):
        ds = generate_lev_dataset("hybrid")

        # Create the expected dataset using the original dataset. This involves
        # updating the arrays and attributes of data variables and coordinates.
        expected = ds.sel(lev=[800, 200]).drop_vars(["ps", "hyam", "hybm"])
        expected["so"].data[:] = np.nan
        expected["so"].attrs["units"] = "mb"
        expected["lev"].attrs = {
            "axis": "Z",
            "coordinate": "vertical",
            "bounds": "lev_bnds",
        }
        expected["lev_bnds"] = xr.DataArray(
            name="lev_bnds",
            data=np.array([[1100.0, 500.0], [500.0, -100.0]]),
            dims=["lev", "bnds"],
            attrs={"xcdat_bounds": "True"},
        )

        result = regrid_z_axis_to_plevs(ds, "so", self.plevs)

        assert_identical(expected, result)

    def test_regrids_hybrid_levels_to_pressure_levels_with_Pa_units(self):
        ds = generate_lev_dataset("hybrid")

        # Create the expected dataset using the original dataset. This involves
        # updating the arrays and attributes of data variables and coordinates.
        expected = ds.sel(lev=[800, 200]).drop_vars(["ps", "hyam", "hybm"])
        expected["so"].data[:] = np.nan
        expected["so"].attrs["units"] = "mb"
        expected["lev"].attrs = {
            "axis": "Z",
            "coordinate": "vertical",
            "bounds": "lev_bnds",
        }
        expected["lev_bnds"] = xr.DataArray(
            name="lev_bnds",
            data=np.array([[1100.0, 500.0], [500.0, -100.0]]),
            dims=["lev", "bnds"],
            attrs={"xcdat_bounds": "True"},
        )

        # Update from Pa to mb.
        ds_pa = ds.copy()
        with xr.set_options(keep_attrs=True):
            ds_pa["ps"] = ds_pa.ps * 100
        ds_pa.ps.attrs["units"] = "Pa"

        result = regrid_z_axis_to_plevs(ds_pa, "so", self.plevs)

        assert_identical(expected, result)

    @pytest.mark.parametrize("long_name", ("pressure", "isobaric"))
    def test_regrids_pressure_coordinates_to_pressure_levels(self, long_name):
        ds = generate_lev_dataset(long_name)

        # Create the expected dataset using the original dataset. This involves
        # updating the arrays and attributes of data variables and coordinates.
        expected = ds.sel(lev=[800, 200]).drop_vars("ps")
        expected["lev"].attrs = {
            "axis": "Z",
            "coordinate": "vertical",
            "bounds": "lev_bnds",
        }
        expected["lev_bnds"] = xr.DataArray(
            name="lev_bnds",
            data=np.array([[1100.0, 500.0], [500.0, -100.0]]),
            dims=["lev", "bnds"],
            attrs={"xcdat_bounds": "True"},
        )
        result = regrid_z_axis_to_plevs(ds, "so", self.plevs)

        assert_identical(expected, result)

    @pytest.mark.parametrize("long_name", ("pressure", "isobaric"))
    def test_regrids_pressure_coordinates_to_pressure_levels_with_Pa_units(
        self, long_name
    ):
        ds = generate_lev_dataset(long_name)

        expected = ds.sel(lev=[800, 200]).drop_vars("ps")
        expected["lev"].attrs = {
            "axis": "Z",
            "coordinate": "vertical",
            "bounds": "lev_bnds",
        }
        expected["lev_bnds"] = xr.DataArray(
            name="lev_bnds",
            data=np.array([[1100.0, 500.0], [500.0, -100.0]]),
            dims=["lev", "bnds"],
            attrs={"xcdat_bounds": "True"},
        )

        # Update from Pa to mb.
        ds_pa = ds.copy()
        with xr.set_options(keep_attrs=True):
            ds_pa["lev"] = ds_pa.lev * 100
        ds_pa.lev.attrs["units"] = "Pa"

        result = regrid_z_axis_to_plevs(ds_pa, "so", self.plevs)

        assert_identical(expected, result)

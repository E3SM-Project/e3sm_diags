import numpy as np
import pytest
import xarray as xr
from xarray.testing import assert_identical

from e3sm_diags.driver.utils.regrid import (
    _apply_land_sea_mask,
    _subset_on_region,
    align_grids_to_lower_res,
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


class Test_ApplyLandSeaMask:
    @pytest.fixture(autouse=True)
    def setup(self):
        self.lat = xr.DataArray(
            data=np.array([-90, -88.75]),
            dims=["lat"],
            attrs={"units": "degrees_north", "axis": "Y", "standard_name": "latitude"},
        )

        self.lon = xr.DataArray(
            data=np.array([0, 1.875]),
            dims=["lon"],
            attrs={"units": "degrees_east", "axis": "X", "standard_name": "longitude"},
        )

    @pytest.mark.filterwarnings("ignore:.*Latitude is outside of.*:UserWarning")
    @pytest.mark.parametrize("regrid_tool", ("xesmf",))
    def test_applies_land_mask_on_variable(self, regrid_tool):
        ds = generate_lev_dataset("pressure").isel(time=1)

        # Create the land mask with different grid.
        land_frac = xr.DataArray(
            name="LANDFRAC",
            data=[[np.nan, 1.0], [1.0, 1.0]],
            dims=["lat", "lon"],
            coords={"lat": self.lat, "lon": self.lon},
        )
        ds_mask = land_frac.to_dataset()

        # Create the expected array for the "so" variable after masking.
        # Updating specific indexes is somewhat hacky but it gets the job done
        # here.
        # TODO: Consider making this part of the test more robust.
        expected_arr = np.empty((4, 4, 4))
        expected_arr[:] = np.nan
        for idx in range(len(ds.lev)):
            expected_arr[idx, 0, 1] = 1

        expected = ds.copy()
        expected.so[:] = expected_arr
        expected["mask"] = xr.where(~np.isnan(expected["so"]), 1, 0)

        result = _apply_land_sea_mask(
            ds, ds_mask, "so", "land", regrid_tool, "conservative_normed"
        )

        assert_identical(expected, result)

    @pytest.mark.filterwarnings("ignore:.*Latitude is outside of.*:UserWarning")
    @pytest.mark.parametrize("regrid_tool", ("xesmf",))
    def test_applies_sea_mask_on_variable(self, regrid_tool):
        ds = generate_lev_dataset("pressure").isel(time=1)

        # Create the land mask with different grid.
        ocean_frac = xr.DataArray(
            name="OCNFRAC",
            data=[[np.nan, 1.0], [1.0, 1.0]],
            dims=["lat", "lon"],
            coords={"lat": self.lat, "lon": self.lon},
        )
        ds_mask = ocean_frac.to_dataset()

        # Create the expected array for the "so" variable after masking.
        # Updating specific indexes is somewhat hacky but it gets the job done
        # here. TODO: Consider making this part of the test more robust.
        expected_arr = np.empty((4, 4, 4))
        expected_arr[:] = np.nan
        for idx in range(len(ds.lev)):
            expected_arr[idx, 0, 1] = 1

        expected = ds.copy()
        expected.so[:] = expected_arr
        expected["mask"] = xr.where(~np.isnan(expected["so"]), 1, 0)

        result = _apply_land_sea_mask(
            ds, ds_mask, "so", "ocean", regrid_tool, "conservative_normed"
        )

        assert_identical(expected, result)

    def test_raises_error_if_region_is_not_valid_land_or_ocean_region(self):
        ds = generate_lev_dataset("pressure").isel(time=1)
        ds_mask = generate_lev_dataset("pressure").isel(time=1)

        with pytest.raises(
            ValueError,
            match="Only land and ocean regions are supported, not 'invalid_region'.",
        ):
            _apply_land_sea_mask(
                ds,
                ds_mask,
                "so",
                "invalid_region",  # type: ignore
                "xesmf",
                "conservative_normed",
            )

    def test_raises_error_if_region_is_does_not_have_region_specs_defined(self):
        ds = generate_lev_dataset("pressure").isel(time=1)
        ds_mask = generate_lev_dataset("pressure").isel(time=1)

        with pytest.raises(
            ValueError, match="No region specifications found for 'land_invalid'."
        ):
            _apply_land_sea_mask(
                ds,
                ds_mask,
                "so",
                "land_invalid",  # type: ignore
                "xesmf",
                "conservative_normed",
            )


class Test_SubsetOnDomain:
    def test_subsets_on_domain_if_region_specs_has_domain_defined(self):
        ds = generate_lev_dataset("pressure").isel(time=1)
        expected = ds.sel(lat=slice(0.0, 45.0), lon=slice(210.0, 310.0))

        result = _subset_on_region(ds, "so", "NAMM")

        assert_identical(expected, result)


class TestAlignGridstoLowerRes:
    @pytest.mark.parametrize("tool", ("xesmf", "regrid2"))
    def test_regrids_to_first_dataset_with_equal_latitude_points(self, tool):
        ds_a = generate_lev_dataset("pressure", pressure_vars=False)
        ds_b = generate_lev_dataset("pressure", pressure_vars=False)

        result_a, result_b = align_grids_to_lower_res(
            ds_a, ds_b, "so", tool, "conservative_normed"
        )

        expected_a = ds_a.copy()
        expected_b = ds_a.copy()

        if tool == "xesmf":
            expected_b.so.attrs["regrid_method"] = "conservative_normed"

        # A has lower resolution (A = B), regrid B -> A.
        assert_identical(result_a, expected_a)

        # NOTE: xesmf regridding changes the order of the dimensions, resulting
        # in lon being before lat. We use assert_equal instead of
        # assert_identical for this case.
        if tool == "xesmf":
            np.testing.assert_equal(result_b, expected_b)
        else:
            assert_identical(result_b, expected_b)

    @pytest.mark.parametrize("tool", ("xesmf", "regrid2"))
    def test_regrids_to_first_dataset_with_conservative_normed_method(self, tool):
        ds_a = generate_lev_dataset("pressure", pressure_vars=False)
        ds_b = generate_lev_dataset("pressure", pressure_vars=False)

        # Subset the first dataset's latitude to make it "lower resolution".
        ds_a = ds_a.isel(lat=slice(0, 3, 1))

        result_a, result_b = align_grids_to_lower_res(
            ds_a, ds_b, "so", tool, "conservative_normed"
        )

        expected_a = ds_a.copy()
        expected_b = ds_a.copy()
        # regrid2 only supports conservative and does not set "regrid_method".
        if tool == "xesmf":
            expected_b.so.attrs["regrid_method"] = "conservative_normed"

        # A has lower resolution (A < B), regrid B -> A.
        assert_identical(result_a, expected_a)

        # NOTE: xesmf regridding changes the order of the dimensions, resulting
        # in lon being before lat. We use assert_equal instead of
        # assert_identical for this case.
        if tool == "xesmf":
            np.testing.assert_equal(result_b, expected_b)
        else:
            assert_identical(result_b, expected_b)

    @pytest.mark.parametrize("tool", ("xesmf", "regrid2"))
    def test_regrids_to_first_dataset_with_conservative_normed_method_and_drops_ilev(
        self, tool
    ):
        ds_a = generate_lev_dataset("pressure", pressure_vars=False)
        ds_b = generate_lev_dataset("pressure", pressure_vars=False)

        # Add an unused ilev axis to ds_a and ds_b
        ds_a["ilev"] = xr.DataArray(data=np.arange(5), dims=["ilev"])
        ds_b["ilev"] = xr.DataArray(data=np.arange(5), dims=["ilev"])

        # Subset the first dataset's latitude to make it "lower resolution".
        ds_a = ds_a.isel(lat=slice(0, 3, 1))

        result_a, result_b = align_grids_to_lower_res(
            ds_a, ds_b, "so", tool, "conservative_normed"
        )

        expected_a = ds_a.copy().drop_vars("ilev")
        expected_b = ds_a.copy().drop_vars("ilev")
        # regrid2 only supports conservative and does not set "regrid_method".
        if tool == "xesmf":
            expected_b.so.attrs["regrid_method"] = "conservative_normed"

        # A has lower resolution (A < B), regrid B -> A.
        assert_identical(result_a, expected_a)

        # NOTE: xesmf regridding changes the order of the dimensions, resulting
        # in lon being before lat. We use assert_equal instead of
        # assert_identical for this case.
        if tool == "xesmf":
            np.testing.assert_equal(result_b, expected_b)
        else:
            assert_identical(result_b, expected_b)

    @pytest.mark.parametrize("tool", ("xesmf", "regrid2"))
    def test_regrids_to_second_dataset_with_conservative_normed_method(self, tool):
        ds_a = generate_lev_dataset("pressure", pressure_vars=False)
        ds_b = generate_lev_dataset("pressure", pressure_vars=False)

        # Subset the second dataset's latitude to make it "lower resolution".
        ds_b = ds_b.isel(lat=slice(0, 3, 1))
        result_a, result_b = align_grids_to_lower_res(
            ds_a, ds_b, "so", tool, "conservative_normed"
        )

        expected_a = ds_b.copy()
        expected_b = ds_b.copy()
        # regrid2 only supports conservative and does not set "regrid_method".
        if tool == "xesmf":
            expected_a.so.attrs["regrid_method"] = "conservative_normed"

        # B has lower resolution (A > B), regrid A -> B.
        # NOTE: xesmf regridding changes the order of the dimensions, resulting
        # in lon being before lat. We use assert_equal instead of
        # assert_identical for this case.
        if tool == "xesmf":
            np.testing.assert_equal(result_a, expected_a)
        else:
            assert_identical(result_a, expected_a)

        assert_identical(result_b, expected_b)


class TestRegridZAxisToPlevs:
    @pytest.fixture(autouse=True)
    def setup(self):
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

    @pytest.mark.filterwarnings(
        "ignore:.*From version 0.8.0 the Axis computation methods will be removed.*:FutureWarning",
        "ignore:.*The `xgcm.Axis` class will be deprecated.*:DeprecationWarning",
    )
    def test_regrids_hybrid_levels_to_pressure_levels_with_existing_z_bounds(self):
        ds = generate_lev_dataset("hybrid")
        del ds.lev_bnds.attrs["xcdat_bounds"]

        # Create the expected dataset using the original dataset. This involves
        # updating the arrays and attributes of data variables and coordinates.
        expected = ds.sel(lev=[800, 200]).drop_vars(["ps", "hyam", "hybm"])
        expected["so"].data[:] = np.nan
        expected["so"].attrs["units"] = "ppt"
        expected["lev"].attrs = {
            "axis": "Z",
            "coordinate": "vertical",
            "bounds": "lev_bnds",
        }
        # New Z bounds are generated for the updated Z axis.
        expected["lev_bnds"] = xr.DataArray(
            name="lev_bnds",
            data=np.array([[1100.0, 500.0], [500.0, -100.0]]),
            dims=["lev", "bnds"],
            attrs={"xcdat_bounds": "True"},
        )

        result = regrid_z_axis_to_plevs(ds, "so", self.plevs)

        assert_identical(expected, result)

    @pytest.mark.filterwarnings(
        "ignore:.*From version 0.8.0 the Axis computation methods will be removed.*:FutureWarning",
        "ignore:.*The `xgcm.Axis` class will be deprecated.*:DeprecationWarning",
    )
    def test_regrids_hybrid_levels_to_pressure_levels_with_generated_z_bounds(self):
        ds = generate_lev_dataset("hybrid")
        ds = ds.drop_vars("lev_bnds")

        # Create the expected dataset using the original dataset. This involves
        # updating the arrays and attributes of data variables and coordinates.
        expected = ds.sel(lev=[800, 200]).drop_vars(["ps", "hyam", "hybm"])
        expected["so"].data[:] = np.nan
        expected["so"].attrs["units"] = "ppt"
        expected["lev"].attrs = {
            "axis": "Z",
            "coordinate": "vertical",
            "bounds": "lev_bnds",
        }
        # New Z bounds are generated for the updated Z axis.
        expected["lev_bnds"] = xr.DataArray(
            name="lev_bnds",
            data=np.array([[1100.0, 500.0], [500.0, -100.0]]),
            dims=["lev", "bnds"],
            attrs={"xcdat_bounds": "True"},
        )

        result = regrid_z_axis_to_plevs(ds, "so", self.plevs)

        assert_identical(expected, result)

    @pytest.mark.filterwarnings(
        "ignore:.*From version 0.8.0 the Axis computation methods will be removed.*:FutureWarning",
        "ignore:.*The `xgcm.Axis` class will be deprecated.*:DeprecationWarning",
    )
    def test_regrids_hybrid_levels_to_pressure_levels_with_Pa_units(self):
        ds = generate_lev_dataset("hybrid")

        # Create the expected dataset using the original dataset. This involves
        # updating the arrays and attributes of data variables and coordinates.
        expected = ds.sel(lev=[800, 200]).drop_vars(["ps", "hyam", "hybm"])
        expected["so"].data[:] = np.nan
        expected["so"].attrs["units"] = "ppt"
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

    @pytest.mark.filterwarnings(
        "ignore:.*From version 0.8.0 the Axis computation methods will be removed.*:FutureWarning",
        "ignore:.*The `xgcm.Axis` class will be deprecated.*:DeprecationWarning",
    )
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

    @pytest.mark.filterwarnings(
        "ignore:.*From version 0.8.0 the Axis computation methods will be removed.*:FutureWarning",
        "ignore:.*The `xgcm.Axis` class will be deprecated.*:DeprecationWarning",
    )
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

        # Update mb to Pa so this test can make sure conversions to mb are done.
        ds_pa = ds.copy()
        with xr.set_options(keep_attrs=True):
            ds_pa["lev"] = ds_pa.lev * 100
            ds_pa["lev_bnds"] = ds_pa.lev_bnds * 100
        ds_pa.lev.attrs["units"] = "Pa"

        result = regrid_z_axis_to_plevs(ds_pa, "so", self.plevs)

        assert_identical(expected, result)

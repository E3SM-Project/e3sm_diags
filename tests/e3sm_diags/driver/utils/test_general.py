import numpy as np
import pytest
import xarray as xr

from e3sm_diags.driver.utils.general import pad_year, subtract_dataarrays


class TestPadYear:
    @pytest.mark.parametrize(
        "input_year, expected_output",
        [
            (42, "0042"),
            ("42", "0042"),
            (2023, "2023"),
            ("2023", "2023"),
            (0, "0000"),
            ("0", "0000"),
            (9999, "9999"),
            ("9999", "9999"),
        ],
    )
    def test_valid_years(self, input_year, expected_output):
        assert pad_year(input_year) == expected_output

    @pytest.mark.parametrize(
        "invalid_year",
        [
            -1,
            "abcd",
            None,
            10000,
            "10000",
            1.5,
            "1.5",
            float("nan"),
            float("inf"),
            "-42",
        ],
    )
    def test_invalid_years(self, invalid_year):
        with pytest.raises(ValueError):
            pad_year(invalid_year)


class TestSubtractDataArrays:
    def test_subtract_basic(self):
        a = xr.DataArray([1, 2, 3], attrs={"units": "meters"})
        b = xr.DataArray([0, 1, 2], attrs={"units": "seconds"})

        result = subtract_dataarrays(a, b)
        expected = xr.DataArray([1, 1, 1], attrs={"units": "meters"})

        assert np.all(result.values == expected.values)
        assert result.attrs == expected.attrs

    def test_subtract_preserves_left_attrs(self):
        a = xr.DataArray([10, 20], attrs={"foo": "bar"})
        b = xr.DataArray([1, 2], attrs={"foo": "baz"})

        result = subtract_dataarrays(a, b)
        assert result.attrs == {"foo": "bar"}

    def test_subtract_shape_mismatch(self):
        a = xr.DataArray([1, 2, 3])
        b = xr.DataArray([1, 2])

        with pytest.raises(
            ValueError, match="Input DataArrays must have the same shape."
        ):
            subtract_dataarrays(a, b)

    def test_subtract_multidimensional(self):
        a = xr.DataArray([[1, 2], [3, 4]], attrs={"desc": "test"})
        b = xr.DataArray([[0, 1], [1, 2]], attrs={"desc": "other"})

        result = subtract_dataarrays(a, b)
        expected = xr.DataArray([[1, 1], [2, 2]], attrs={"desc": "test"})

        assert np.all(result.values == expected.values)
        assert result.attrs == expected.attrs

    def test_subtract_empty_arrays(self):
        a = xr.DataArray([], attrs={"empty": "yes"})
        b = xr.DataArray([], attrs={"empty": "no"})

        result = subtract_dataarrays(a, b)

        assert result.size == 0
        assert result.attrs == {"empty": "yes"}

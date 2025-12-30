import pytest

from e3sm_diags.driver.utils.general import pad_year


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

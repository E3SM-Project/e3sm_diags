import numpy as np
import pytest
import xarray as xr

from e3sm_diags.derivations.formulas_cosp import (
    cosp_bin_sum,
    cosp_histogram_standardize,
)


class TestCospHistogramStandardize:
    @pytest.fixture(autouse=True)
    def setup(self):
        self.target_var_key = "CLDTOT_TAU1.3_MISR"
        self.ds = xr.Dataset(
            data_vars={
                "cld_var_dummy": xr.DataArray(
                    data=np.array(
                        [
                            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
                            [2.0, 2.0, 2.0, 2.0, 2.0, 2.0],
                            [3.0, 3.0, 3.0, 3.0, 3.0, 3.0],
                            [4.0, 4.0, 4.0, 4.0, 4.0, 4.0],
                            [5.0, 5.0, 5.0, 5.0, 5.0, 5.0],
                            [6.0, 6.0, 6.0, 6.0, 6.0, 6.0],
                            [7.0, 7.0, 7.0, 7.0, 7.0, 7.0],
                        ],
                        dtype="float64",
                    ),
                    dims=["cosp_htmisr", "cosp_tau"],
                    attrs={"test_attr": "test"},
                ),
                "cosp_htmisr_bnds": xr.DataArray(
                    data=np.array(
                        [
                            [-1.0, 0.0],
                            [0.0, 1.0],
                            [0.5, 1.5],
                            [1.5, 2.5],
                            [2.5, 3.5],
                            [3.5, 4.5],
                            [4.5, 5.5],
                        ]
                    ),
                    dims=["cosp_htmisr", "bnds"],
                ),
                "cosp_tau_bnds": xr.DataArray(
                    data=np.array(
                        [
                            [-1.0, 0.0],
                            [0.0, 1.0],
                            [0.5, 1.5],
                            [1.5, 2.5],
                            [2.5, 3.5],
                            [3.5, 4.5],
                        ]
                    ),
                    dims=["cosp_tau", "bnds"],
                ),
            },
            coords={
                "cosp_htmisr": xr.DataArray(
                    data=np.array([-0.5, 0.5, 1, 2, 3, 4, 5]),
                    dims=["cosp_htmisr"],
                    attrs={
                        "bounds": "cosp_htmisr_bnds",
                        "units": "km",
                        "long_name": "COSP MISR height",
                        "realtopology": "linear",
                    },
                ),
                "cosp_tau": xr.DataArray(
                    data=np.array([0.0, 0.5, 1, 2, 3, 4]),
                    dims=["cosp_tau"],
                    attrs={
                        "bounds": "cosp_tau_bnds",
                        "units": "1",
                        "long_name": "COSP MEAN ISCCP optional depth",
                        "realtopology": "linear",
                    },
                ),
            },
        )

    def test_raises_error_if_dataset_is_missing_prs_or_tau_axis(self):
        ds1 = self.ds.copy()
        ds1 = ds1.rename({"cld_var_dummy": "CLD_MISR", "cosp_htmisr": "invalid_name"})

        with pytest.raises(KeyError):
            cosp_histogram_standardize(self.target_var_key, ds1["CLD_MISR"])

        ds2 = self.ds.copy()
        ds2 = ds2.rename({"cld_var_dummy": "CLD_MISR", "cosp_tau": "invalid_name"})

        with pytest.raises(KeyError):
            cosp_histogram_standardize(self.target_var_key, ds2["CLD_MISR"])

    @pytest.mark.parametrize("var_key", ("CLD_MISR", "CLMISR"))
    def test_standardizes_cld_misr_and_clmisr(self, var_key):
        ds = self.ds.copy()

        ds = ds.rename({"cld_var_dummy": var_key})

        result = cosp_histogram_standardize(self.target_var_key, ds[var_key])
        expected = xr.DataArray(
            name=self.target_var_key,
            data=np.array(
                [
                    [2.0, 2.0, 2.0, 2.0, 2.0],
                    [3.0, 3.0, 3.0, 3.0, 3.0],
                    [4.0, 4.0, 4.0, 4.0, 4.0],
                    [5.0, 5.0, 5.0, 5.0, 5.0],
                    [6.0, 6.0, 6.0, 6.0, 6.0],
                    [7.0, 7.0, 7.0, 7.0, 7.0],
                ],
                dtype="float64",
            ),
            dims=["cosp_htmisr", "cosp_tau"],
            coords={
                "cosp_htmisr": xr.DataArray(
                    data=np.array([0.5, 1, 2, 3, 4, 5]),
                    dims=["cosp_htmisr"],
                    attrs={
                        "bounds": "cosp_htmisr_bnds",
                        "units": "km",
                        "long_name": "COSP MISR height",
                        "realtopology": "linear",
                    },
                ),
                "cosp_tau": xr.DataArray(
                    data=np.array([0.5, 1, 2, 3, 4]),
                    dims=["cosp_tau"],
                    attrs={
                        "bounds": "cosp_tau_bnds",
                        "units": "1",
                        "long_name": "COSP MEAN ISCCP optional depth",
                        "realtopology": "linear",
                    },
                ),
            },
            attrs={"test_attr": "test"},
        )

        xr.testing.assert_identical(result, expected)

    def test_returns_sum_with_cld_misr_with_unit_adjustment(self):
        target_var_key = "CLDLOW_TAU1.3_MISR"

        ds1 = xr.Dataset(
            data_vars={
                "CLD_MISR": xr.DataArray(
                    data=np.array(
                        [
                            [1.0, 1.0, 1.0],
                            [2.0, 2.0, 2.0],
                            [3.0, 3.0, 3.0],
                        ],
                        dtype="float64",
                    ),
                    dims=["cosp_htmisr", "cosp_tau"],
                    attrs={"test_attr": "test"},
                ),
            },
            coords={
                "cosp_htmisr": xr.DataArray(
                    data=np.array([0, 1, 2000]),
                    dims=["cosp_htmisr"],
                    attrs={
                        "bounds": "cosp_htmisr_bnds",
                        "units": "km",
                        "long_name": "COSP MISR height",
                        "realtopology": "linear",
                    },
                ),
                "cosp_tau": xr.DataArray(
                    data=np.array([-0.5, 1.5, 2.0]),
                    dims=["cosp_tau"],
                    attrs={
                        "bounds": "cosp_tau_bnds",
                        "units": "1",
                        "long_name": "COSP MEAN ISCCP optional depth",
                        "realtopology": "linear",
                    },
                ),
            },
        )

        result = cosp_histogram_standardize(target_var_key, ds1["CLD_MISR"])

        expected = xr.DataArray(
            name=target_var_key,
            data=np.array(
                [[2.0, 2.0], [3.0, 3.0]],
                dtype="float64",
            ),
            dims=["cosp_htmisr", "cosp_tau"],
            coords={
                "cosp_htmisr": xr.DataArray(
                    data=np.array([1, 2000]),
                    dims=["cosp_htmisr"],
                    attrs={
                        "bounds": "cosp_htmisr_bnds",
                        "units": "km",
                        "long_name": "COSP MISR height",
                        "realtopology": "linear",
                    },
                ),
                "cosp_tau": xr.DataArray(
                    data=np.array([1.5, 2.0]),
                    dims=["cosp_tau"],
                    attrs={
                        "bounds": "cosp_tau_bnds",
                        "units": "1",
                        "long_name": "COSP MEAN ISCCP optional depth",
                        "realtopology": "linear",
                    },
                ),
            },
            attrs={"test_attr": "test"},
        )

        xr.testing.assert_identical(result, expected)

    @pytest.mark.xfail
    @pytest.mark.parametrize("var_key", ("CLD_MISR", "CLMISR"))
    def test_standardizes_cld_misr_and_cldmisr_and_adds_default_bnds_if_bnds_are_missing(
        self, var_key
    ):
        # TODO: Update this test if missing bounds are dynamically generated.
        ds = self.ds.copy()
        ds = ds.drop_vars(["cosp_htmisr_bnds", "cosp_tau_bnds"])

        # Rename the cloud variable placeholder to the variable to be tested.
        ds = ds.rename({"cld_var_dummy": var_key})
        result = cosp_histogram_standardize(self.target_var_key, ds[var_key])

        expected = xr.DataArray(
            name=self.target_var_key,
            data=np.array(
                [
                    [2.0, 2.0, 2.0, 2.0, 2.0],
                    [3.0, 3.0, 3.0, 3.0, 3.0],
                    [4.0, 4.0, 4.0, 4.0, 4.0],
                    [5.0, 5.0, 5.0, 5.0, 5.0],
                    [6.0, 6.0, 6.0, 6.0, 6.0],
                    [7.0, 7.0, 7.0, 7.0, 7.0],
                ],
                dtype="float64",
            ),
            dims=["cosp_htmisr", "cosp_tau"],
            coords={
                "cosp_htmisr": xr.DataArray(
                    data=np.array([0.5, 1, 2, 3, 4, 5]),
                    dims=["cosp_htmisr"],
                    attrs={
                        "bounds": "cosp_htmisr_bnds",
                        "units": "km",
                        "long_name": "COSP MISR height",
                        "realtopology": "linear",
                    },
                ),
                "cosp_tau": xr.DataArray(
                    data=np.array([0.5, 1, 2, 3, 4]),
                    dims=["cosp_tau"],
                    attrs={
                        "bounds": "cosp_tau_bnds",
                        "units": "1",
                        "long_name": "COSP MEAN ISCCP optional depth",
                        "realtopology": "linear",
                    },
                ),
            },
            attrs={"test_attr": "test"},
        )

        xr.testing.assert_identical(result, expected)

    @pytest.mark.parametrize(
        "var_key",
        ("FISCCP1_COSP", "CLISCCP", "CLMODIS", "isccp_ctptau", "modis_ctptau"),
    )
    def test_standardizes_fisccp1_cosp_clisccp_and_clmodis(self, var_key):
        ds = self.ds.copy()

        # Rename the cloud variable placeholder to the variable to be tested
        # and also rename the "cosp_tau" dimension to test other dimension keys.
        ds = ds.rename({"cld_var_dummy": var_key})
        result = cosp_histogram_standardize(self.target_var_key, ds[var_key])
        expected = xr.DataArray(
            name=self.target_var_key,
            data=np.array(
                [
                    [1.0, 1.0, 1.0, 1.0, 1.0],
                    [2.0, 2.0, 2.0, 2.0, 2.0],
                    [3.0, 3.0, 3.0, 3.0, 3.0],
                    [4.0, 4.0, 4.0, 4.0, 4.0],
                    [5.0, 5.0, 5.0, 5.0, 5.0],
                    [6.0, 6.0, 6.0, 6.0, 6.0],
                    [7.0, 7.0, 7.0, 7.0, 7.0],
                ],
                dtype="float64",
            ),
            dims=["cosp_htmisr", "cosp_tau"],
            coords={
                "cosp_htmisr": xr.DataArray(
                    data=np.array([-0.5, 0.5, 1, 2, 3, 4, 5]),
                    dims=["cosp_htmisr"],
                    attrs={
                        "bounds": "cosp_htmisr_bnds",
                        "units": "km",
                        "long_name": "COSP MISR height",
                        "realtopology": "linear",
                    },
                ),
                "cosp_tau": xr.DataArray(
                    data=np.array([0.5, 1, 2, 3, 4]),
                    dims=["cosp_tau"],
                    attrs={
                        "bounds": "cosp_tau_bnds",
                        "units": "1",
                        "long_name": "COSP MEAN ISCCP optional depth",
                        "realtopology": "linear",
                    },
                ),
            },
            attrs={"test_attr": "test"},
        )

        xr.testing.assert_identical(result, expected)


class TestCospBinSum:
    @pytest.fixture(autouse=True)
    def setup(self):
        self.ds = xr.Dataset(
            data_vars={
                "cld_var_dummy": xr.DataArray(
                    data=np.array(
                        [
                            [1.0, 1.0, 1.0],
                            [2.0, 2.0, 2.0],
                            [3.0, 3.0, 3.0],
                        ],
                        dtype="float64",
                    ),
                    dims=["cosp_htmisr", "cosp_tau"],
                    attrs={"test_attr": "test"},
                ),
            },
            coords={
                "cosp_htmisr": xr.DataArray(
                    data=np.array([-0.5, 0.5, 1]),
                    dims=["cosp_htmisr"],
                    attrs={
                        "bounds": "cosp_htmisr_bnds",
                        "units": "km",
                        "long_name": "COSP MISR height",
                        "realtopology": "linear",
                    },
                ),
                "cosp_tau": xr.DataArray(
                    data=np.array([0.0, 0.5, 1]),
                    dims=["cosp_tau"],
                    attrs={
                        "bounds": "cosp_tau_bnds",
                        "units": "1",
                        "long_name": "COSP MEAN ISCCP optional depth",
                        "realtopology": "linear",
                    },
                ),
            },
        )

    def test_raises_error_if_dataset_is_missing_prs_or_tau_axis(self):
        target_var_key = "CLDTOT_TAU1.3_MISR"

        ds1 = self.ds.copy()
        ds1 = ds1.rename({"cosp_htmisr": "invalid_name"})

        with pytest.raises(KeyError):
            cosp_bin_sum(target_var_key, ds1["CLD_MISR"])

        ds2 = self.ds.copy()
        ds2 = ds2.rename({"cosp_tau": "invalid_name"})

        with pytest.raises(KeyError):
            cosp_bin_sum(target_var_key, ds1["CLD_MISR"])

    def test_returns_sum(self):
        target_var_key = "CLDTOT_TAU1.3_MISR"

        ds1 = self.ds.copy()
        ds1 = ds1.rename({"cld_var_dummy": "CLD_MISR"})

        result = cosp_bin_sum(target_var_key, ds1["CLD_MISR"])

        expected = xr.DataArray(
            name="CLD_MISR",
            data=np.array(6.0),
            attrs={
                "test_attr": "test",
                "long_name": "MISR: total cloud fraction with tau > 1.3",
            },
        )

        xr.testing.assert_identical(result, expected)

    def test_returns_sum_using_prs_subset(self):
        target_var_key = "CLDLOW_TAU1.3_MISR"

        ds1 = xr.Dataset(
            data_vars={
                "CLD_MISR": xr.DataArray(
                    data=np.array(
                        [
                            [1.0, 1.0, 1.0],
                            [2.0, 2.0, 2.0],
                            [3.0, 3.0, 3.0],
                        ],
                        dtype="float64",
                    ),
                    dims=["cosp_htmisr", "cosp_tau"],
                    attrs={"test_attr": "test"},
                ),
            },
            coords={
                "cosp_htmisr": xr.DataArray(
                    data=np.array([0, 1, 2]),
                    dims=["cosp_htmisr"],
                    attrs={
                        "bounds": "cosp_htmisr_bnds",
                        "units": "km",
                        "long_name": "COSP MISR height",
                        "realtopology": "linear",
                    },
                ),
                "cosp_tau": xr.DataArray(
                    data=np.array([-0.5, 1.5, 2.0]),
                    dims=["cosp_tau"],
                    attrs={
                        "bounds": "cosp_tau_bnds",
                        "units": "1",
                        "long_name": "COSP MEAN ISCCP optional depth",
                        "realtopology": "linear",
                    },
                ),
            },
        )

        result = cosp_bin_sum(target_var_key, ds1["CLD_MISR"])

        expected = xr.DataArray(
            name="CLD_MISR",
            data=np.array(12),
            attrs={
                "test_attr": "test",
                "long_name": "MISR: low cloud fraction with tau > 1.3",
            },
        )

        xr.testing.assert_identical(result, expected)

    def test_returns_sum_using_prs_subset_with_unit_adjustment(self):
        target_var_key = "CLDLOW_TAU1.3_MISR"

        ds1 = xr.Dataset(
            data_vars={
                "CLD_MISR": xr.DataArray(
                    data=np.array(
                        [
                            [1.0, 1.0, 1.0],
                            [2.0, 2.0, 2.0],
                            [3.0, 3.0, 3.0],
                        ],
                        dtype="float64",
                    ),
                    dims=["cosp_htmisr", "cosp_tau"],
                    attrs={"test_attr": "test"},
                ),
            },
            coords={
                "cosp_htmisr": xr.DataArray(
                    data=np.array([0, 1, 2000]),
                    dims=["cosp_htmisr"],
                    attrs={
                        "bounds": "cosp_htmisr_bnds",
                        "units": "km",
                        "long_name": "COSP MISR height",
                        "realtopology": "linear",
                    },
                ),
                "cosp_tau": xr.DataArray(
                    data=np.array([-0.5, 1.5, 2.0]),
                    dims=["cosp_tau"],
                    attrs={
                        "bounds": "cosp_tau_bnds",
                        "units": "1",
                        "long_name": "COSP MEAN ISCCP optional depth",
                        "realtopology": "linear",
                    },
                ),
            },
        )

        result = cosp_bin_sum(target_var_key, ds1["CLD_MISR"])

        expected = xr.DataArray(
            name="CLD_MISR",
            data=np.array(12),
            attrs={
                "test_attr": "test",
                "long_name": "MISR: low cloud fraction with tau > 1.3",
            },
        )

        xr.testing.assert_identical(result, expected)

    def test_returns_sum_using_tau_subset_with_adjusted_min_and_max(self):
        target_var_key = "CLDTOT_TAU1.3_9.4_ISCCP"

        ds1 = xr.Dataset(
            data_vars={
                "CLD_MISR": xr.DataArray(
                    data=np.array(
                        [
                            [1.0, 1.0, 1.0],
                            [2.0, 2.0, 2.0],
                            [3.0, 3.0, 3.0],
                        ],
                        dtype="float64",
                    ),
                    dims=["cosp_htmisr", "cosp_tau"],
                    attrs={"test_attr": "test"},
                ),
            },
            coords={
                "cosp_htmisr": xr.DataArray(
                    data=np.array([0, 1, 2]),
                    dims=["cosp_htmisr"],
                    attrs={
                        "bounds": "cosp_htmisr_bnds",
                        "units": "km",
                        "long_name": "COSP MISR height",
                        "realtopology": "linear",
                    },
                ),
                "cosp_tau": xr.DataArray(
                    data=np.array([-0.5, 1.5, 10.0]),
                    dims=["cosp_tau"],
                    attrs={
                        "bounds": "cosp_tau_bnds",
                        "units": "1",
                        "long_name": "COSP MEAN ISCCP optional depth",
                        "realtopology": "linear",
                    },
                ),
            },
        )

        result = cosp_bin_sum(target_var_key, ds1["CLD_MISR"])

        expected = xr.DataArray(
            name="CLD_MISR",
            data=np.array(6),
            attrs={
                "test_attr": "test",
                "long_name": "MISR: total cloud fraction with 1.3 < tau < 9.4",
            },
        )

        xr.testing.assert_identical(result, expected)

    def test_returns_sum_using_tau_subset_with_adjusted_min_only(self):
        target_var_key = "CLDTOT_TAU1.3_ISCCP"

        ds1 = xr.Dataset(
            data_vars={
                "CLD_MISR": xr.DataArray(
                    data=np.array(
                        [
                            [1.0, 1.0, 1.0],
                            [2.0, 2.0, 2.0],
                            [3.0, 3.0, 3.0],
                        ],
                        dtype="float64",
                    ),
                    dims=["cosp_htmisr", "cosp_tau"],
                    attrs={"test_attr": "test"},
                ),
            },
            coords={
                "cosp_htmisr": xr.DataArray(
                    data=np.array([0, 1, 2]),
                    dims=["cosp_htmisr"],
                    attrs={
                        "bounds": "cosp_htmisr_bnds",
                        "units": "km",
                        "long_name": "COSP MISR height",
                        "realtopology": "linear",
                    },
                ),
                "cosp_tau": xr.DataArray(
                    data=np.array([-0.5, 1.5, 2.0]),
                    dims=["cosp_tau"],
                    attrs={
                        "bounds": "cosp_tau_bnds",
                        "units": "1",
                        "long_name": "COSP MEAN ISCCP optional depth",
                        "realtopology": "linear",
                    },
                ),
            },
        )

        result = cosp_bin_sum(target_var_key, ds1["CLD_MISR"])

        expected = xr.DataArray(
            name="CLD_MISR",
            data=np.array(12),
            attrs={
                "test_attr": "test",
                "long_name": "MISR: total cloud fraction with tau > 1.3",
            },
        )

        xr.testing.assert_identical(result, expected)

    @pytest.mark.parametrize(
        "var_key,expected",
        [
            ("FISCCP1_COSP", "ISCCP: total cloud fraction with tau > 1.3"),
            (
                "CLMODIS",
                "MODIS: total cloud fraction with tau > 1.3",
            ),
            ("CLD_MISR", "MISR: total cloud fraction with tau > 1.3"),
        ],
    )
    def test_sets_variable_long_name_attr_if_matching_simulation_and_cloud_levels_are_set(
        self, var_key, expected
    ):
        target_var_key = "CLDTOT_TAU1.3_MODIS"

        ds1 = self.ds.copy()
        ds1 = ds1.rename({"cld_var_dummy": var_key})

        result_var = cosp_bin_sum(target_var_key, ds1[var_key])
        result = result_var.attrs["long_name"]
        assert result == expected

    @pytest.mark.parametrize(
        "var_key,expected,prs_axis_data",
        [
            (
                "CLMODIS",
                "MODIS: middle cloud fraction with tau > 1.3",
                np.array([440, 500, 680]),
            ),
            (
                "CLD_MISR",
                "MISR: middle cloud fraction with tau > 1.3",
                np.array([7, 5, 3]),
            ),
        ],
    )
    def test_sets_variable_long_name_attr_with_middle_cloud_fraction(
        self, var_key, expected, prs_axis_data
    ):
        # Min in low_bnds and max in high_bnds
        target_var_key = "CLDTOT_TAU1.3_MODIS"

        ds1 = xr.Dataset(
            data_vars={
                var_key: xr.DataArray(
                    data=np.array(
                        [
                            [1.0, 1.0, 1.0],
                            [2.0, 2.0, 2.0],
                            [3.0, 3.0, 3.0],
                        ],
                        dtype="float64",
                    ),
                    dims=["cosp_htmisr", "cosp_tau"],
                    attrs={"test_attr": "test"},
                ),
            },
            coords={
                "cosp_htmisr": xr.DataArray(
                    data=prs_axis_data,
                    dims=["cosp_htmisr"],
                    attrs={
                        "bounds": "cosp_htmisr_bnds",
                        "units": "km",
                        "long_name": "COSP MISR height",
                        "realtopology": "linear",
                    },
                ),
                "cosp_tau": xr.DataArray(
                    data=np.array([-0.5, 1.5, 2.0]),
                    dims=["cosp_tau"],
                    attrs={
                        "bounds": "cosp_tau_bnds",
                        "units": "1",
                        "long_name": "COSP MEAN ISCCP optional depth",
                        "realtopology": "linear",
                    },
                ),
            },
        )

        result_var = cosp_bin_sum(target_var_key, ds1[var_key])
        result = result_var.attrs["long_name"]
        assert result == expected

    @pytest.mark.parametrize(
        "var_key,expected,prs_axis_data",
        [
            (
                "CLMODIS",
                "MODIS: high cloud fraction with tau > 1.3",
                np.array([440, 500, 600]),
            ),
            (
                "CLD_MISR",
                "MISR: high cloud fraction with tau > 1.3",
                np.array([7, 5, 6]),
            ),
        ],
    )
    def test_sets_variable_long_name_attr_with_high_cloud_fraction(
        self, var_key, expected, prs_axis_data
    ):
        target_var_key = "CLDTOT_TAU1.3_MODIS"

        ds1 = xr.Dataset(
            data_vars={
                var_key: xr.DataArray(
                    data=np.array(
                        [
                            [1.0, 1.0, 1.0],
                            [2.0, 2.0, 2.0],
                            [3.0, 3.0, 3.0],
                        ],
                        dtype="float64",
                    ),
                    dims=["cosp_htmisr", "cosp_tau"],
                    attrs={"test_attr": "test"},
                ),
            },
            coords={
                "cosp_htmisr": xr.DataArray(
                    data=prs_axis_data,
                    dims=["cosp_htmisr"],
                    attrs={
                        "bounds": "cosp_htmisr_bnds",
                        "units": "km",
                        "long_name": "COSP MISR height",
                        "realtopology": "linear",
                    },
                ),
                "cosp_tau": xr.DataArray(
                    data=np.array([-0.5, 1.5, 2.0]),
                    dims=["cosp_tau"],
                    attrs={
                        "bounds": "cosp_tau_bnds",
                        "units": "1",
                        "long_name": "COSP MEAN ISCCP optional depth",
                        "realtopology": "linear",
                    },
                ),
            },
        )

        result_var = cosp_bin_sum(target_var_key, ds1[var_key])
        result = result_var.attrs["long_name"]
        assert result == expected

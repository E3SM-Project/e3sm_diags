import numpy as np
import pytest
import xarray as xr

from e3sm_diags.derivations.formulas_cosp import (
    CLOUD_HIST_MAP,
    cosp_histogram_standardize,
)


class TestCospHistogramStandardizeXr:
    @pytest.fixture(autouse=True)
    def setup(self):
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
        ds1 = ds1.rename({"cosp_htmisr": "invalid_name"})

        with pytest.raises(KeyError):
            cosp_histogram_standardize(ds1, ds1["CLD_MISR"])

        ds2 = self.ds.copy()
        ds2 = ds2.rename({"cosp_tau": "invalid_name"})

        with pytest.raises(KeyError):
            cosp_histogram_standardize(ds2, ds2["CLD_MISR"])

    @pytest.mark.parametrize("var_key", ("CLD_MISR", "CLMISR"))
    def test_standardizes_cld_misr_and_clmisr(self, var_key):
        ds = self.ds.copy()

        ds = ds.rename({"cld_var_dummy": var_key})
        result = cosp_histogram_standardize(ds, ds[var_key])
        expected = xr.Dataset(
            data_vars={
                var_key: xr.DataArray(
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
                    attrs={"test_attr": "test"},
                ),
                "cosp_htmisr_bnds": xr.DataArray(
                    data=np.array(
                        [
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
        )

        xr.testing.assert_allclose(result, expected)

    @pytest.mark.parametrize("var_key", ("CLD_MISR", "CLMISR"))
    def test_standardizes_cld_misr_and_cldmisr_and_adds_default_bnds_if_bnds_are_missing(
        self, var_key
    ):
        ds = self.ds.copy()
        ds = ds.drop_vars(["cosp_htmisr_bnds", "cosp_tau_bnds"])
        default_prs_bnds = CLOUD_HIST_MAP["prs"]["default_bnds"][1:]
        default_tau_bnds = CLOUD_HIST_MAP["tau"]["default_bnds"][1:]

        # Rename the cloud variable placeholder to the variable to be tested.
        ds = ds.rename({"cld_var_dummy": var_key})
        result = cosp_histogram_standardize(ds, ds[var_key])

        expected = xr.Dataset(
            data_vars={
                var_key: xr.DataArray(
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
                    attrs={"test_attr": "test"},
                ),
                "cosp_htmisr_bnds": xr.DataArray(
                    data=default_prs_bnds,
                    dims=["cosp_htmisr", "bnds"],
                ),
                "cosp_tau_bnds": xr.DataArray(
                    data=default_tau_bnds,
                    dims=["cosp_tau", "bnds"],
                ),
            },
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
        )

        xr.testing.assert_allclose(result, expected)

    @pytest.mark.parametrize("var_key", ("FISCCP1_COSP", "CLISCCP", "CLMODIS"))
    def test_standardizes_fisccp1_cosp_clisccp_and_clmodis(self, var_key):
        ds = self.ds.copy()

        # Rename the cloud variable placeholder to the variable to be tested
        # and also rename the "cosp_tau" dimension to test other dimension keys.
        ds = ds.rename({"cld_var_dummy": var_key})
        result = cosp_histogram_standardize(ds, ds[var_key])
        expected = xr.Dataset(
            data_vars={
                var_key: xr.DataArray(
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
        )

        xr.testing.assert_allclose(result, expected)

from __future__ import annotations

from unittest.mock import MagicMock

import numpy as np
import pytest
import xarray as xr

from e3sm_diags.driver import enso_diags_driver
from e3sm_diags.parameter.enso_diags_parameter import EnsoDiagsParameter


def _create_parameter(run_type: str = "model_vs_model") -> EnsoDiagsParameter:
    parameter = EnsoDiagsParameter()
    parameter.run_type = run_type
    parameter.nino_regions = ["NINO3", "NINO34"]
    parameter.variables = ["TS"]
    parameter.viewer_descr = {}
    parameter.save_netcdf = False

    return parameter


class TestCalculateNinoIndices:
    def test_model_vs_model_uses_model_indices_for_test_and_ref(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        parameter = _create_parameter()
        test_ds = MagicMock()
        ref_ds = MagicMock()
        calls: list[tuple[object, str]] = []

        def _calculate_model(
            ds: object, _parameter: EnsoDiagsParameter, region: str
        ) -> xr.DataArray:
            calls.append((ds, region))
            return xr.DataArray([len(calls)])

        calculate_obs = MagicMock()
        monkeypatch.setattr(
            enso_diags_driver, "calculate_nino_index_model", _calculate_model
        )
        monkeypatch.setattr(
            enso_diags_driver, "calculate_nino_index_obs", calculate_obs
        )

        test_indices, ref_indices = enso_diags_driver._calculate_nino_indices(
            parameter, test_ds, ref_ds, parameter.nino_regions
        )

        assert list(test_indices) == parameter.nino_regions
        assert list(ref_indices) == parameter.nino_regions
        assert calls == [
            (test_ds, "NINO3"),
            (ref_ds, "NINO3"),
            (test_ds, "NINO34"),
            (ref_ds, "NINO34"),
        ]
        calculate_obs.assert_not_called()

    def test_model_vs_obs_uses_observational_indices_for_ref(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        parameter = _create_parameter(run_type="model_vs_obs")
        test_ds = MagicMock()
        ref_ds = MagicMock()
        calculate_model = MagicMock(side_effect=[xr.DataArray([1]), xr.DataArray([2])])
        calculate_obs = MagicMock(side_effect=[xr.DataArray([3]), xr.DataArray([4])])
        monkeypatch.setattr(
            enso_diags_driver, "calculate_nino_index_model", calculate_model
        )
        monkeypatch.setattr(
            enso_diags_driver, "calculate_nino_index_obs", calculate_obs
        )

        test_indices, ref_indices = enso_diags_driver._calculate_nino_indices(
            parameter, test_ds, ref_ds, parameter.nino_regions
        )

        assert list(test_indices) == parameter.nino_regions
        assert list(ref_indices) == parameter.nino_regions
        assert [call.args for call in calculate_model.call_args_list] == [
            (test_ds, parameter, "NINO3"),
            (test_ds, parameter, "NINO34"),
        ]
        assert [call.args for call in calculate_obs.call_args_list] == [
            (parameter, "NINO3", "ref"),
            (parameter, "NINO34", "ref"),
        ]

    def test_invalid_run_type_raises_exception(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        parameter = _create_parameter(run_type="invalid")
        calculate_model = MagicMock(return_value=xr.DataArray([1]))
        monkeypatch.setattr(
            enso_diags_driver, "calculate_nino_index_model", calculate_model
        )

        with pytest.raises(Exception, match="Invalid run_type=invalid"):  # noqa: B017
            enso_diags_driver._calculate_nino_indices(
                parameter, MagicMock(), MagicMock(), ["NINO3"]
            )


class TestNinoIndexRunners:
    def test_timeseries_passes_indices_to_plot(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        parameter = _create_parameter()
        parameter.variables = []
        set_name_yrs_attrs = MagicMock()
        monkeypatch.setattr(parameter, "_set_name_yrs_attrs", set_name_yrs_attrs)
        test_ds = MagicMock()
        ref_ds = MagicMock()
        monkeypatch.setattr(
            enso_diags_driver,
            "Dataset",
            MagicMock(side_effect=[test_ds, ref_ds]),
        )
        test_indices = {"NINO3": xr.DataArray([1])}
        ref_indices = {"NINO3": xr.DataArray([2])}
        calculate_indices = MagicMock(return_value=(test_indices, ref_indices))
        plot = MagicMock()
        monkeypatch.setattr(
            enso_diags_driver, "_calculate_nino_indices", calculate_indices
        )
        monkeypatch.setattr(enso_diags_driver, "plot_nino_index_timeseries", plot)

        result = enso_diags_driver.run_diag_nino_index_timeseries(parameter)

        set_name_yrs_attrs.assert_called_once_with(test_ds, ref_ds, None)
        calculate_indices.assert_called_once_with(
            parameter, test_ds, ref_ds, parameter.nino_regions
        )
        plot.assert_called_once_with(parameter, test_indices, ref_indices)
        assert parameter.var_id == "NINO-index"
        assert parameter.main_title == "Nino index time series"
        assert parameter.output_file == "nino-index-timeseries"
        assert parameter.viewer_descr["TS"] == parameter.main_title
        assert result is parameter

    def test_seasonality_transforms_indices_before_plot(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        parameter = _create_parameter()
        set_name_yrs_attrs = MagicMock()
        monkeypatch.setattr(parameter, "_set_name_yrs_attrs", set_name_yrs_attrs)
        test_ds = MagicMock()
        ref_ds = MagicMock()
        monkeypatch.setattr(
            enso_diags_driver,
            "Dataset",
            MagicMock(side_effect=[test_ds, ref_ds]),
        )
        test_indices = {
            "NINO3": xr.DataArray([1]),
            "NINO34": xr.DataArray([2]),
        }
        ref_indices = {
            "NINO3": xr.DataArray([3]),
            "NINO34": xr.DataArray([4]),
        }
        calculate_indices = MagicMock(return_value=(test_indices, ref_indices))
        monkeypatch.setattr(
            enso_diags_driver, "_calculate_nino_indices", calculate_indices
        )
        calculate_seasonality = MagicMock(
            side_effect=[
                np.array([10]),
                np.array([20]),
                np.array([30]),
                np.array([40]),
            ]
        )
        monkeypatch.setattr(
            enso_diags_driver, "_calculate_seasonality", calculate_seasonality
        )
        plot = MagicMock()
        monkeypatch.setattr(enso_diags_driver, "plot_seasonality", plot)

        result = enso_diags_driver.run_diag_seasonality(parameter)

        set_name_yrs_attrs.assert_called_once_with(test_ds, ref_ds, None)
        calculate_indices.assert_called_once_with(
            parameter, test_ds, ref_ds, parameter.nino_regions
        )
        plot.assert_called_once()
        test_seasonality = plot.call_args.args[1]
        ref_seasonality = plot.call_args.args[2]
        np.testing.assert_array_equal(test_seasonality["NINO3"], [10])
        np.testing.assert_array_equal(test_seasonality["NINO34"], [20])
        np.testing.assert_array_equal(ref_seasonality["NINO3"], [30])
        np.testing.assert_array_equal(ref_seasonality["NINO34"], [40])
        assert parameter.var_id == "NINO-index"
        assert parameter.main_title == "ENSO Seasonality: NINO index"
        assert parameter.output_file == "nino-index-seasonality"
        assert parameter.viewer_descr["TS"] == parameter.main_title
        assert result is parameter

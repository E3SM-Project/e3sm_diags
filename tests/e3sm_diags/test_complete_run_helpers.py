"""Test pure helpers used by the manual complete-run comparison flow.

This module contains focused unit coverage for helper logic that does not
require HPC data or a full diagnostics run, such as variable-key inference,
derived-variable fallback, file-tree matching, and summary classification.
It exists to keep a small amount of automated regression protection around
the refactored manual workflow while leaving the actual complete run as a
developer-driven validation path.

Usage
-----
Run this module with pytest as part of the normal unit-test suite, for
example with ``pytest tests/e3sm_diags/test_complete_run_helpers.py``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import xarray as xr

from tests.complete_run.helpers import (
    classify_array_difference,
    expand_candidate_var_keys,
    get_var_data,
    infer_variable_key_from_path,
    match_netcdf_files,
)


class TestInferVariableKeyFromPath:
    def test_returns_standard_variable_key(self):
        result = infer_variable_key_from_path("lat_lon/plot-ALBEDO-ANN-global.nc")

        assert result == "ALBEDO"

    def test_returns_3d_variable_key_when_pressure_suffix_present(self):
        result = infer_variable_key_from_path("zonal_mean_2d/plot-T-200-ANN-global.nc")

        assert result == "T"


class TestExpandCandidateVarKeys:
    def test_includes_derived_source_variables(self):
        result = expand_candidate_var_keys("ALBEDO")

        assert result[0] == "ALBEDO"
        assert "SOLIN" in result
        assert "FSNTOA" in result


class TestGetVarData:
    def test_resolves_first_available_candidate_variable(self):
        ds = xr.Dataset(
            {
                "SOLIN": xr.DataArray(np.array([1.0, 2.0])),
                "FSNTOA": xr.DataArray(np.array([0.5, 1.5])),
            }
        )

        data, matched_key = get_var_data(ds, "ALBEDO")

        assert matched_key == "SOLIN"
        np.testing.assert_array_equal(data, np.array([1.0, 2.0]))


class TestMatchNetcdfFiles:
    def test_reports_missing_paths_on_both_sides(self, tmp_path: Path):
        dev_root = tmp_path / "dev"
        baseline_root = tmp_path / "baseline"
        (dev_root / "lat_lon").mkdir(parents=True)
        (baseline_root / "lat_lon").mkdir(parents=True)

        (dev_root / "lat_lon" / "shared-ALBEDO-ANN-global.nc").touch()
        (dev_root / "lat_lon" / "dev-only-ALBEDO-ANN-global.nc").touch()
        (baseline_root / "lat_lon" / "shared-ALBEDO-ANN-global.nc").touch()
        (baseline_root / "lat_lon" / "baseline-only-ALBEDO-ANN-global.nc").touch()

        result = match_netcdf_files(dev_root, baseline_root)

        assert result.shared_paths == [Path("lat_lon/shared-ALBEDO-ANN-global.nc")]
        assert result.missing_dev_paths == [
            Path("lat_lon/baseline-only-ALBEDO-ANN-global.nc")
        ]
        assert result.missing_baseline_paths == [
            Path("lat_lon/dev-only-ALBEDO-ANN-global.nc")
        ]


class TestClassifyArrayDifference:
    def test_reports_nan_location_mismatch(self):
        result = classify_array_difference(
            np.array([1.0, np.nan]),
            np.array([1.0, 2.0]),
            atol=0.0,
            rtol=1e-5,
        )

        assert result[0] == "nan_location_mismatch"
        assert result[1] is not None
        assert "NaN locations differ" in result[1]

    def test_reports_tolerance_failure(self):
        result = classify_array_difference(
            np.array([1.0, 2.0]),
            np.array([1.0, 2.1]),
            atol=0.0,
            rtol=1e-6,
        )

        assert result[0] == "tolerance_failure"
        assert result[1] is not None
        assert "Not equal to tolerance" in result[1]

    def test_reports_match(self):
        result = classify_array_difference(
            np.array([1.0, np.nan]),
            np.array([1.0, np.nan]),
            atol=0.0,
            rtol=1e-5,
        )

        assert result == ("matching", None)

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

from argparse import Namespace
from pathlib import Path

import numpy as np
import pytest
import xarray as xr
from PIL import Image

from tests.complete_run import baseline, run
from tests.complete_run.helpers import (
    classify_array_difference,
    compare_dataset_pair,
    compare_png_trees,
    expand_candidate_var_keys,
    get_var_data,
    infer_variable_key_from_path,
    match_netcdf_files,
    match_png_files,
)
from tests.complete_run.params import (
    CompleteRunConfig,
    CompleteRunPaths,
    build_complete_run_config,
)


class TestInferVariableKeyFromPath:
    def test_returns_standard_variable_key(self):
        result = infer_variable_key_from_path("lat_lon/plot-ALBEDO-ANN-global.nc")

        assert result == "ALBEDO"

    def test_returns_3d_variable_key_when_pressure_suffix_present(self):
        result = infer_variable_key_from_path("zonal_mean_2d/plot-T-200-ANN-global.nc")

        assert result == "T"

    @pytest.mark.parametrize("data_type", ["test", "ref"])
    def test_returns_qbo_zonal_wind_key(self, data_type: str):
        result = infer_variable_key_from_path(f"qbo/qbo_diags_qbo_{data_type}.nc")

        assert result == "U"


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


class TestCompleteRunImageComparison:
    @staticmethod
    def _write_png(path: Path, color: str = "black") -> None:
        Image.new("RGB", (10, 10), color).save(path)

    def test_matches_png_files_by_relative_path(self, tmp_path: Path):
        dev_root = tmp_path / "dev"
        baseline_root = tmp_path / "baseline"
        (dev_root / "lat_lon").mkdir(parents=True)
        (baseline_root / "lat_lon").mkdir(parents=True)
        self._write_png(dev_root / "lat_lon" / "shared.png")
        self._write_png(dev_root / "lat_lon" / "dev-only.png")
        self._write_png(baseline_root / "lat_lon" / "shared.png")
        self._write_png(baseline_root / "lat_lon" / "baseline-only.png")

        result = match_png_files(dev_root, baseline_root)

        assert result.shared_paths == [Path("lat_lon/shared.png")]
        assert result.missing_dev_paths == [Path("lat_lon/baseline-only.png")]
        assert result.missing_baseline_paths == [Path("lat_lon/dev-only.png")]

    def test_reports_pixel_mismatch_and_writes_artifacts(self, tmp_path: Path):
        dev_root = tmp_path / "dev"
        baseline_root = tmp_path / "baseline"
        (dev_root / "lat_lon").mkdir(parents=True)
        (baseline_root / "lat_lon").mkdir(parents=True)
        self._write_png(dev_root / "lat_lon" / "plot.png", "white")
        self._write_png(baseline_root / "lat_lon" / "plot.png", "black")

        summary = compare_png_trees(
            dev_root,
            baseline_root,
            mismatch_threshold=0.0002,
            diff_artifact_dir=tmp_path / "artifacts",
        )

        assert summary.matching_images == []
        assert len(summary.image_mismatches) == 1
        assert "Mismatched pixel fraction: 1" in summary.image_mismatches[0].detail
        assert summary.image_mismatches[0].artifact_path == (
            tmp_path / "artifacts" / "image-diffs" / "lat_lon" / "plot_diff.png"
        )
        assert summary.image_mismatches[0].artifact_path.exists()


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


class TestCompareDatasetPair:
    def test_compares_all_shared_variables(self):
        dev_ds = xr.Dataset(
            {
                "NINO3": xr.DataArray(np.array([1.0, 2.0])),
                "NINO34": xr.DataArray(np.array([3.0, 4.0])),
            }
        )
        baseline_ds = dev_ds.copy(deep=True)

        outcomes = compare_dataset_pair(
            dev_ds,
            baseline_ds,
            relative_path="enso_diags/nino-index-timeseries_test.nc",
            atol=0.0,
            rtol=1e-5,
        )

        assert [(outcome.var_key, outcome.status) for outcome in outcomes] == [
            ("NINO3", "matching"),
            ("NINO34", "matching"),
        ]

    def test_reports_actual_variable_missing_from_one_dataset(self):
        outcomes = compare_dataset_pair(
            xr.Dataset({"FLNS": xr.DataArray(np.array([1.0]))}),
            xr.Dataset({"FSNS": xr.DataArray(np.array([1.0]))}),
            relative_path="enso_diags/feedback-FLNS-NINO3-TS-NINO3_test.nc",
            atol=0.0,
            rtol=1e-5,
        )

        assert [(outcome.var_key, outcome.detail) for outcome in outcomes] == [
            ("FLNS", "Variable 'FLNS' is missing from the baseline dataset."),
            ("FSNS", "Variable 'FSNS' is missing from the dev dataset."),
        ]


def _complete_run_config(results_dir: Path) -> CompleteRunConfig:
    paths = CompleteRunPaths(
        results_dir=str(results_dir),
        test_climo="test-climo",
        test_ts="test-ts",
        test_ts_daily_dir="test-ts-daily",
        test_diurnal_climo="test-diurnal",
        test_streamflow_ts="test-streamflow",
        test_tc_analysis="test-tc",
        test_arm_site="test-arm",
        ref_climo="ref-climo",
        ref_ts="ref-ts",
        ref_tc_analysis="ref-tc",
        ref_arm="ref-arm",
    )
    return build_complete_run_config(case="case", short_name="short", paths=paths)


class TestCompleteRunManifest:
    def test_collect_package_versions_uses_not_installed_fallback(
        self, monkeypatch: pytest.MonkeyPatch
    ):
        def missing_distribution(package: str) -> str:
            raise baseline.metadata.PackageNotFoundError(package)

        monkeypatch.setattr(baseline.metadata, "version", missing_distribution)

        assert baseline._collect_package_versions() == dict.fromkeys(
            baseline._CURATED_PACKAGES, baseline._NOT_INSTALLED
        )

    def test_manifest_includes_and_validates_environment(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        monkeypatch.setattr(
            baseline,
            "_collect_package_versions",
            lambda: dict.fromkeys(baseline._CURATED_PACKAGES, "1.0"),
        )
        manifest = baseline._build_manifest(_complete_run_config(tmp_path), [])

        assert manifest["environment"]["packages"]["xarray"] == "1.0"
        assert manifest["environment"]["python_version"]
        del manifest["environment"]

        with pytest.raises(baseline._ManifestError, match="incomplete schema"):
            baseline._validate_manifest(manifest, tmp_path / "manifest.json")

    def test_run_writes_manifest_after_runner_completes(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        results_dir = tmp_path / "run"
        results_dir.mkdir()
        config = _complete_run_config(results_dir)
        args = Namespace(
            case=config.case,
            short_name=config.short_name,
            start_yr=config.start_yr,
            end_yr=config.end_yr,
            num_workers=config.num_workers,
            save_netcdf=config.save_netcdf,
            sets_to_run=["lat_lon"],
            workflow_revision=None,
            **config.paths.__dict__,
        )
        monkeypatch.setattr(run, "build_complete_run_config", lambda **_: config)
        monkeypatch.setattr(run, "_validate_input_paths", lambda _: None)
        monkeypatch.setattr(run, "build_complete_run_params", lambda _: [])
        monkeypatch.setattr(run.runner, "run_diags", lambda _: [])
        monkeypatch.setattr(
            baseline, "_get_git_metadata", lambda: {"branch": "main", "sha": "abc123"}
        )

        run._run_complete_run(args)

        manifest = baseline._load_manifest(results_dir)
        assert manifest["git"] == {"branch": "main", "sha": "abc123"}
        assert manifest["workflow_revision"] == "abc123"
        assert manifest["result_dir"] == str(results_dir)
        assert manifest["test_paths"]["climo"] == "test-climo"
        assert manifest["reference_paths"]["arm"] == "ref-arm"
        assert manifest["config"]["selected_sets"] == ["lat_lon"]
        assert manifest["config"]["save_netcdf"] is True
        assert (
            results_dir / baseline._MANIFEST_FILENAME
        ).stat().st_mode & 0o777 == 0o644

    def test_run_records_explicit_workflow_revision(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        results_dir = tmp_path / "run"
        results_dir.mkdir()
        config = _complete_run_config(results_dir)
        args = Namespace(
            case=config.case,
            short_name=config.short_name,
            start_yr=config.start_yr,
            end_yr=config.end_yr,
            num_workers=config.num_workers,
            save_netcdf=config.save_netcdf,
            sets_to_run=["lat_lon"],
            workflow_revision="workflow-sha",
            **config.paths.__dict__,
        )
        monkeypatch.setattr(run, "build_complete_run_config", lambda **_: config)
        monkeypatch.setattr(run, "_validate_input_paths", lambda _: None)
        monkeypatch.setattr(run, "build_complete_run_params", lambda _: [])
        monkeypatch.setattr(run.runner, "run_diags", lambda _: [])
        monkeypatch.setattr(
            baseline, "_get_git_metadata", lambda: {"branch": "main", "sha": "code-sha"}
        )

        run._run_complete_run(args)

        assert (
            baseline._load_manifest(results_dir)["workflow_revision"] == "workflow-sha"
        )

    def test_manifest_write_cleans_temporary_file_on_serialization_failure(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        run_dir = tmp_path / "run"
        run_dir.mkdir()
        manifest = baseline._build_manifest(_complete_run_config(run_dir), ["lat_lon"])

        def fail_dump(*args: object, **kwargs: object) -> None:
            raise TypeError("not serializable")

        monkeypatch.setattr(baseline.json, "dump", fail_dump)

        with pytest.raises(TypeError, match="not serializable"):
            baseline._write_manifest(run_dir, manifest)

        assert not (run_dir / baseline._MANIFEST_FILENAME).exists()
        assert list(run_dir.glob(f".{baseline._MANIFEST_FILENAME}.*.tmp")) == []

    def test_manifest_write_rejects_and_preserves_existing_manifest(
        self, tmp_path: Path
    ):
        run_dir = tmp_path / "run"
        run_dir.mkdir()
        manifest_path = run_dir / baseline._MANIFEST_FILENAME
        original_contents = '{"existing": true}\n'
        manifest_path.write_text(original_contents, encoding="utf-8")

        with pytest.raises(FileExistsError, match="immutable manifest"):
            baseline._write_manifest(
                run_dir, baseline._build_manifest(_complete_run_config(run_dir), [])
            )

        assert manifest_path.read_text(encoding="utf-8") == original_contents

    def test_run_rejects_existing_manifest_before_runner(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        results_dir = tmp_path / "run"
        results_dir.mkdir()
        config = _complete_run_config(results_dir)
        baseline._write_manifest(results_dir, baseline._build_manifest(config, []))
        args = Namespace(
            case=config.case,
            short_name=config.short_name,
            start_yr=config.start_yr,
            end_yr=config.end_yr,
            num_workers=config.num_workers,
            save_netcdf=config.save_netcdf,
            sets_to_run=["lat_lon"],
            workflow_revision=None,
            **config.paths.__dict__,
        )
        monkeypatch.setattr(run, "build_complete_run_config", lambda **_: config)
        monkeypatch.setattr(run, "_validate_input_paths", lambda _: None)
        monkeypatch.setattr(run, "build_complete_run_params", lambda _: [])
        monkeypatch.setattr(
            run.runner,
            "run_diags",
            lambda _: pytest.fail("runner must not run for immutable results"),
        )

        with pytest.raises(FileExistsError, match="immutable results directory"):
            run._run_complete_run(args)

    def test_promote_rejects_non_main_manifest(self, tmp_path: Path):
        root = tmp_path / "baselines"
        run_dir = tmp_path / "run"
        root.mkdir()
        run_dir.mkdir()
        manifest = baseline._build_manifest(_complete_run_config(run_dir), ["lat_lon"])
        manifest["git"] = {"branch": "feature", "sha": "abc123"}
        baseline._write_manifest(run_dir, manifest)

        with pytest.raises(baseline._ManifestError, match="non-main"):
            baseline.promote_baseline(run_dir, "main", results_root=root)

        link = baseline.promote_baseline(
            run_dir, "main", allow_non_main=True, results_root=root
        )
        assert link.resolve() == run_dir.resolve()

    def test_show_resolves_latest_target_and_manifest(self, tmp_path: Path):
        root = tmp_path / "baselines"
        run_dir = tmp_path / "run"
        root.mkdir()
        run_dir.mkdir()
        manifest = baseline._build_manifest(_complete_run_config(run_dir), ["lat_lon"])
        manifest["git"] = {"branch": "main", "sha": "abc123"}
        baseline._write_manifest(run_dir, manifest)
        baseline.promote_baseline(run_dir, "main", results_root=root)

        target, shown_manifest = baseline.show_baseline("main", results_root=root)

        assert target == run_dir.resolve()
        assert shown_manifest == manifest

"""Tests for the manual complete-run comparison CLI."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
import xarray as xr
from PIL import Image

from tests.complete_run import baseline, compare
from tests.complete_run.helpers import ComparisonSummary
from tests.complete_run.params import DEFAULT_RESULTS_DIR


def _write_manifest_with_environment(
    directory: Path, environment: dict[str, object]
) -> None:
    manifest = {
        "schema_version": baseline._MANIFEST_SCHEMA_VERSION,
        "created_at_utc": "2026-01-01T00:00:00+00:00",
        "result_dir": str(directory),
        "git": {"branch": "main", "sha": "abc123"},
        "workflow_revision": "abc123",
        "environment": environment,
        "test_paths": {},
        "reference_paths": {},
        "config": {"selected_sets": []},
    }
    baseline._write_manifest(directory, manifest)


def _environment() -> dict[str, object]:
    return {
        "python_version": "3.11.0",
        "python_implementation": "CPython",
        "platform": "Linux-test",
        "conda_environment": "test",
        "packages": dict.fromkeys(baseline._CURATED_PACKAGES, "1.0"),
    }


def test_parser_defaults_to_latest_main_baseline():
    args = compare._build_parser().parse_args(["--dev-dir", "dev-results"])

    assert compare.DEFAULT_BASELINE_DIR == Path(DEFAULT_RESULTS_DIR) / "latest-main"
    assert args.baseline_dir == compare.DEFAULT_BASELINE_DIR


def test_parser_requires_dev_dir():
    with pytest.raises(SystemExit):
        compare._build_parser().parse_args([])


def test_missing_latest_main_pointer_explains_promotion(tmp_path: Path, monkeypatch):
    dev_dir = tmp_path / "dev"
    dev_dir.mkdir()
    latest_main = tmp_path / "latest-main"
    monkeypatch.setattr(compare, "DEFAULT_BASELINE_DIR", latest_main)

    with pytest.raises(
        FileNotFoundError, match="No accepted main baseline is promoted"
    ) as error:
        compare._validate_compare_dirs(dev_dir, latest_main)

    assert "baseline promote --run-dir <run-dir> --channel main" in str(error.value)


def test_environment_differences_emit_warning(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    dev_dir = tmp_path / "dev"
    baseline_dir = tmp_path / "baseline"
    dev_dir.mkdir()
    baseline_dir.mkdir()
    dev_environment = _environment()
    dev_environment["conda_environment"] = "newer"
    dev_environment["packages"]["xarray"] = "2.0"  # type: ignore[index]
    _write_manifest_with_environment(dev_dir, dev_environment)
    _write_manifest_with_environment(baseline_dir, _environment())
    for directory, environment_name, xarray_version in (
        (baseline_dir, "baseline", "1.0"),
        (dev_dir, "development", "2.0"),
    ):
        provenance_dir = directory / "prov"
        provenance_dir.mkdir()
        (provenance_dir / "environment.yml").write_text(
            f"name: {environment_name}\ndependencies:\n  - xarray={xarray_version}\n",
            encoding="utf-8",
        )
    warnings: list[str] = []
    monkeypatch.setattr(
        compare.logger,
        "warning",
        lambda message, *args: warnings.append(message % args),
    )

    environment_comparison = compare._warn_environment_differences(
        dev_dir, baseline_dir
    )

    assert len(warnings) == 1
    assert "conda_environment" in warnings[0]
    assert "xarray" in warnings[0]
    assert str(baseline_dir / "prov" / "environment.yml") in warnings[0]
    assert str(dev_dir / "prov" / "environment.yml") in warnings[0]
    assert "-  - xarray=1.0" in warnings[0]
    assert "+  - xarray=2.0" in warnings[0]
    assert "name: baseline" not in warnings[0]
    assert "name: development" not in warnings[0]
    assert environment_comparison["environment_file_diff"] == {
        "available": True,
        "changes": [
            {
                "operation": "replace",
                "baseline_lines": ["  - xarray=1.0"],
                "dev_lines": ["  - xarray=2.0"],
            }
        ],
    }


def test_missing_manifests_only_log_info_and_do_not_fail_comparison(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    dev_dir = tmp_path / "dev"
    baseline_dir = tmp_path / "baseline"
    dev_dir.mkdir()
    baseline_dir.mkdir()
    warnings: list[str] = []
    monkeypatch.setattr(
        compare.logger,
        "warning",
        lambda message, *args: warnings.append(message % args),
    )
    monkeypatch.setattr(
        compare, "compare_netcdf_trees", lambda **_: ComparisonSummary()
    )

    assert (
        compare.main(["--dev-dir", str(dev_dir), "--baseline-dir", str(baseline_dir)])
        == 0
    )
    assert warnings == []


@pytest.mark.parametrize(
    ("summary", "expected_exit_code"),
    [
        (ComparisonSummary(matching_files=[Path("matching.nc")]), 0),
        (ComparisonSummary(missing_dev_files=[Path("missing.nc")]), 1),
    ],
)
def test_main_returns_comparison_status(
    tmp_path: Path,
    monkeypatch,
    summary: ComparisonSummary,
    expected_exit_code: int,
):
    dev_dir = tmp_path / "dev"
    baseline_dir = tmp_path / "baseline"
    dev_dir.mkdir()
    baseline_dir.mkdir()
    monkeypatch.setattr(compare, "compare_netcdf_trees", lambda **_: summary)

    result = compare.main(
        ["--dev-dir", str(dev_dir), "--baseline-dir", str(baseline_dir)]
    )

    assert result == expected_exit_code
    report_path = tmp_path / "comparison" / "dev-vs-baseline" / "comparison-report.json"
    report = json.loads(report_path.read_text(encoding="utf-8"))
    assert report["exit_code"] == expected_exit_code
    assert report["summary"]["failure_count"] == summary.failure_count


def test_images_mode_skips_netcdf_checks(tmp_path: Path):
    """``--mode images`` must narrow the comparison, not add to the default."""
    dev_dir = tmp_path / "dev"
    baseline_dir = tmp_path / "baseline"
    (dev_dir / "lat_lon").mkdir(parents=True)
    (baseline_dir / "lat_lon").mkdir(parents=True)
    Image.new("RGB", (10, 10), "white").save(dev_dir / "lat_lon" / "plot.png")
    Image.new("RGB", (10, 10), "white").save(baseline_dir / "lat_lon" / "plot.png")
    xr.Dataset({"ts": ("x", [1.0])}).to_netcdf(baseline_dir / "lat_lon" / "ts.nc")

    result = compare.main(
        [
            "--dev-dir",
            str(dev_dir),
            "--baseline-dir",
            str(baseline_dir),
            "--mode",
            "images",
        ]
    )

    assert result == 0
    report_path = tmp_path / "comparison" / "dev-vs-baseline" / "comparison-report.json"
    report = json.loads(report_path.read_text(encoding="utf-8"))
    assert report["comparison_settings"]["modes"] == ["images"]
    assert report["summary"]["missing_dev_files"] == []


def test_images_mode_fails_and_reports_png_mismatches(tmp_path: Path):
    dev_dir = tmp_path / "dev"
    baseline_dir = tmp_path / "baseline"
    (dev_dir / "lat_lon").mkdir(parents=True)
    (baseline_dir / "lat_lon").mkdir(parents=True)
    Image.new("RGB", (10, 10), "white").save(dev_dir / "lat_lon" / "plot.png")
    Image.new("RGB", (10, 10), "black").save(baseline_dir / "lat_lon" / "plot.png")

    result = compare.main(
        [
            "--dev-dir",
            str(dev_dir),
            "--baseline-dir",
            str(baseline_dir),
            "--mode",
            "images",
            "--write-diff-pngs",
        ]
    )

    assert result == 1
    report_path = tmp_path / "comparison" / "dev-vs-baseline" / "comparison-report.json"
    report = json.loads(report_path.read_text(encoding="utf-8"))
    assert report["summary"]["image_mismatches"][0]["relative_path"] == (
        "lat_lon/plot.png"
    )
    assert (
        tmp_path
        / "comparison"
        / "dev-vs-baseline"
        / "diff-pngs"
        / "image-diffs"
        / "lat_lon"
        / "plot_diff.png"
    ).exists()

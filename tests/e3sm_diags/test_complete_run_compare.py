"""Tests for the manual complete-run comparison CLI."""

from __future__ import annotations

import json
import re
from datetime import datetime, timezone
from pathlib import Path

import pytest
import xarray as xr
from PIL import Image

from tests.complete_run import baseline, compare, diff_html
from tests.complete_run.helpers import ComparisonIssue, ComparisonSummary
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


def _find_comparison_report(tmp_path: Path) -> Path:
    reports = list(tmp_path.glob("comparison/*/comparison-report.json"))
    assert len(reports) == 1
    return reports[0]


def test_parser_defaults_to_latest_main_baseline():
    args = compare._build_parser().parse_args(["--dev-dir", "dev-results"])

    assert compare.DEFAULT_BASELINE_DIR == Path(DEFAULT_RESULTS_DIR) / "latest-main"
    assert args.baseline_dir == compare.DEFAULT_BASELINE_DIR


def test_parser_requires_dev_dir():
    with pytest.raises(SystemExit):
        compare._build_parser().parse_args([])


def test_parser_rejects_obsolete_raw_image_threshold():
    with pytest.raises(SystemExit):
        compare._build_parser().parse_args(
            ["--dev-dir", "dev-results", "--image-mismatch-threshold", "0.1"]
        )


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
    publicized_paths: list[Path] = []
    monkeypatch.setattr(compare, "make_tree_public", publicized_paths.append)

    result = compare.main(
        ["--dev-dir", str(dev_dir), "--baseline-dir", str(baseline_dir)]
    )

    assert result == expected_exit_code
    report_path = _find_comparison_report(tmp_path)
    report = json.loads(report_path.read_text(encoding="utf-8"))
    assert report["exit_code"] == expected_exit_code
    assert report["summary"]["failure_count"] == summary.failure_count
    assert report["summary"]["coverage"]["netcdf"] == {
        "compared": summary.compared_file_count,
        "identical": len(summary.matching_files),
        "cosmetic": 0,
        "different": summary.compared_file_count - len(summary.matching_files),
        "missing_dev": len(summary.missing_dev_files),
        "missing_baseline": len(summary.missing_baseline_files),
    }
    assert publicized_paths == [report_path.parent]
    assert re.fullmatch(r"dev-vs-baseline-\d{8}-\d{6}", report_path.parent.name)


def test_comparison_coverage_counts_netcdf_and_png_artifacts():
    summary = ComparisonSummary(
        matching_files=[Path("matching.nc")],
        tolerance_failures=[ComparisonIssue(Path("different.nc"))],
        missing_dev_files=[Path("missing-dev.nc")],
        missing_baseline_files=[Path("missing-baseline.nc")],
        matching_images=[Path("same.png"), Path("cosmetic.png")],
        identical_images=[Path("same.png")],
        cosmetic_images=[Path("cosmetic.png")],
        image_mismatches=[ComparisonIssue(Path("different.png"))],
        missing_dev_images=[Path("missing-dev.png")],
        missing_baseline_images=[Path("missing-baseline.png")],
    )

    assert compare._comparison_coverage(summary) == {
        "netcdf": {
            "compared": 2,
            "identical": 1,
            "cosmetic": 0,
            "different": 1,
            "missing_dev": 1,
            "missing_baseline": 1,
        },
        "png": {
            "compared": 3,
            "identical": 1,
            "cosmetic": 1,
            "different": 1,
            "missing_dev": 1,
            "missing_baseline": 1,
        },
    }


@pytest.mark.parametrize("artifact_flag", ["--write-diff-pngs", "--write-diff-html"])
def test_clean_diff_artifact_directory_is_not_publicized(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, artifact_flag: str
):
    dev_dir = tmp_path / "dev"
    baseline_dir = tmp_path / "baseline"
    dev_dir.mkdir()
    baseline_dir.mkdir()
    monkeypatch.setattr(
        compare, "compare_netcdf_trees", lambda **_: ComparisonSummary()
    )
    publicized_paths: list[Path] = []
    monkeypatch.setattr(compare, "make_tree_public", publicized_paths.append)

    assert (
        compare.main(
            [
                "--dev-dir",
                str(dev_dir),
                "--baseline-dir",
                str(baseline_dir),
                artifact_flag,
            ]
        )
        == 0
    )

    report_path = _find_comparison_report(tmp_path)
    assert publicized_paths == [report_path.parent]


def test_comparison_report_path_reserves_unique_directory(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    class FixedDatetime:
        @classmethod
        def now(cls, tz: timezone) -> datetime:
            return datetime(2026, 9, 16, 12, 0, tzinfo=tz)

    monkeypatch.setattr(compare, "datetime", FixedDatetime)
    dev_dir = tmp_path / "dev"
    baseline_dir = tmp_path / "baseline"

    first = compare._comparison_report_path(dev_dir, baseline_dir, tmp_path / "reports")
    second = compare._comparison_report_path(
        dev_dir, baseline_dir, tmp_path / "reports"
    )

    assert first.parent.name == "dev-vs-baseline-20260916-120000"
    assert second.parent.name == "dev-vs-baseline-20260916-120000-2"
    assert first.parent.is_dir()
    assert second.parent.is_dir()


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
    report_path = _find_comparison_report(tmp_path)
    report = json.loads(report_path.read_text(encoding="utf-8"))
    assert report["comparison_settings"]["modes"] == ["images"]
    assert report["comparison_settings"]["image_comparison"] == "severity-v1"
    assert report["summary"]["missing_dev_files"] == []
    assert report["summary"]["identical_images"] == ["lat_lon/plot.png"]
    assert report["summary"]["cosmetic_images"] == []


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
    report_path = _find_comparison_report(tmp_path)
    report = json.loads(report_path.read_text(encoding="utf-8"))
    assert report["summary"]["image_mismatches"][0]["relative_path"] == (
        "lat_lon/plot.png"
    )
    assert report["summary"]["image_mismatches"][0]["severity"] == "MAJOR"
    assert report["summary"]["image_mismatches"][0]["content_fraction"] == 1.0
    assert "content fraction: 1" in report["summary"]["image_mismatches"][0]["detail"]
    assert (
        report_path.parent / "diff-pngs" / "image-diffs" / "lat_lon" / "plot_diff.png"
    ).exists()


class TestDiffHtml:
    def _report(self, tmp_path: Path, mismatches: list[dict]) -> dict:
        return {
            "paths": {
                "dev_dir": str(tmp_path / "dev"),
                "baseline_dir": str(tmp_path / "baseline"),
            },
            "environment": {"differences": ["matplotlib (3.10.9 -> 3.11.1)"]},
            "summary": {
                "matching_files": ["a.nc"],
                "matching_images": [],
                "identical_images": [],
                "cosmetic_images": [],
                "compared_file_count": 1,
                "failure_count": len(mismatches),
                "image_mismatches": mismatches,
                "missing_baseline_files": [],
                "missing_baseline_images": [],
            },
        }

    def test_returns_none_without_image_mismatches(self, tmp_path: Path):
        report_path = tmp_path / "comparison-report.json"

        assert (
            diff_html.write_diff_html(self._report(tmp_path, []), report_path) is None
        )
        assert not (tmp_path / "index.html").exists()

    def test_sorts_by_severity_then_content_fraction_and_links_triptych(
        self, tmp_path: Path
    ):
        diffs = tmp_path / "diff-pngs"
        diffs.mkdir()
        mismatches = [
            {
                "relative_path": "polar/minor.png",
                "severity": "MINOR",
                "content_fraction": 0.9,
                "cause": "small isolated difference",
                "artifact_path": str(diffs / "small_diff.png"),
            },
            {
                "relative_path": "lat_lon/major.png",
                "severity": "MAJOR",
                "content_fraction": 0.05,
                "cause": "same size, content differs",
                "artifact_path": str(diffs / "big_diff.png"),
            },
            {
                "relative_path": "polar/structural.png",
                "severity": "STRUCTURAL",
                "content_fraction": 0.01,
                "cause": "figure size changed a lot",
                "artifact_path": str(diffs / "structural_diff.png"),
            },
        ]
        report_path = tmp_path / "comparison-report.json"

        out = diff_html.write_diff_html(self._report(tmp_path, mismatches), report_path)

        assert out == tmp_path / "index.html"
        page = out.read_text(encoding="utf-8")
        match = re.search(r"const ROWS = (\[.*?\]);", page, re.S)
        assert match is not None
        rows = json.loads(match.group(1))
        assert [row["path"] for row in rows] == [
            "polar/structural.png",
            "lat_lon/major.png",
            "polar/minor.png",
        ]
        assert rows[0]["severity"] == "STRUCTURAL"
        assert rows[0]["level"] == 6
        assert rows[0]["cause"] == "figure size changed a lot"
        # Paths are relative to the report, and the baseline/current panels are
        # derived from the diff artifact's name.
        assert rows[1]["diff"] == "diff-pngs/big_diff.png"
        assert rows[1]["expected"] == "diff-pngs/big_expected.png"
        assert rows[1]["actual"] == "diff-pngs/big_actual.png"

    def test_keeps_phase_one_content_fraction_visible_in_viewer(self, tmp_path: Path):
        diffs = tmp_path / "diff-pngs"
        diffs.mkdir()
        mismatches = [
            {
                "relative_path": "lat_lon/plot.png",
                "severity": "MAJOR",
                "content_fraction": 0.25,
                "cause": "same size, content differs",
                "artifact_path": str(diffs / "plot_diff.png"),
            }
        ]

        page = diff_html.write_diff_html(
            self._report(tmp_path, mismatches), tmp_path / "comparison-report.json"
        )

        assert page is not None
        assert '"frac": 0.25' in page.read_text(encoding="utf-8")

    def test_renders_severity_controls_and_cosmetic_counts(self, tmp_path: Path):
        report = self._report(
            tmp_path,
            [
                {
                    "relative_path": "lat_lon/plot.png",
                    "severity": "MAJOR",
                    "content_fraction": 0.25,
                    "cause": "same size, content differs",
                    "artifact_path": str(tmp_path / "plot_diff.png"),
                }
            ],
        )
        report["summary"]["identical_images"] = ["lat_lon/exact.png"]
        report["summary"]["cosmetic_images"] = ["lat_lon/shifted.png"]
        report["summary"]["cosmetic_samples"] = [
            {
                "relative_path": "lat_lon/shifted.png",
                "severity": "NEGLIGIBLE",
                "content_fraction": 0.0,
                "raw_fraction": 0.1,
                "cause": "same size, content differs",
                "artifact_path": str(tmp_path / "shifted_diff.png"),
            }
        ]

        page = diff_html.write_diff_html(report, tmp_path / "comparison-report.json")

        assert page is not None
        content = page.read_text(encoding="utf-8")
        assert 'data-severity="MAJOR"' in content
        assert 'data-severity="NEGLIGIBLE"' in content
        assert "Severity guide" in content
        assert '<table class="severity-table">' in content
        assert "Severity, highest first" in content
        assert "review levels (1)" in content
        assert "6. structural (0)" in content
        assert "5. major (1)" in content
        assert "4. moderate (0)" in content
        assert "3. minor (0)" in content
        assert "1. Identical" in content
        assert "2. Negligible" in content
        assert "Passed; sample available" in content
        assert "1 total, 1 sampled" in content
        assert "unmatched content" in content
        assert "pixels differ" in content
        assert "images passing" in content
        assert "cosmetic" in content

    def test_html_flag_implies_diff_artifacts(self, tmp_path: Path):
        """The index links to diff PNGs, so requesting it must produce them."""
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
                "--write-diff-html",
            ]
        )

        assert result == 1
        index = _find_comparison_report(tmp_path).parent / "index.html"
        assert index.exists()
        assert "lat_lon/plot.png" in index.read_text(encoding="utf-8")

    def test_html_includes_a_cosmetic_sample_without_reviewable_images(
        self, tmp_path: Path
    ):
        dev_dir = tmp_path / "dev"
        baseline_dir = tmp_path / "baseline"
        (dev_dir / "lat_lon").mkdir(parents=True)
        (baseline_dir / "lat_lon").mkdir(parents=True)
        baseline = Image.new("RGB", (100, 100), "white")
        baseline.paste("black", (30, 30, 70, 70))
        baseline.save(baseline_dir / "lat_lon" / "shifted.png")
        actual = Image.new("RGB", (100, 100), "white")
        actual.paste("black", (30, 31, 70, 71))
        actual.save(dev_dir / "lat_lon" / "shifted.png")

        assert (
            compare.main(
                [
                    "--dev-dir",
                    str(dev_dir),
                    "--baseline-dir",
                    str(baseline_dir),
                    "--write-diff-html",
                ]
            )
            == 0
        )

        index = _find_comparison_report(tmp_path).parent / "index.html"
        content = index.read_text(encoding="utf-8")
        assert "lat_lon/shifted.png" in content
        assert '"severity": "NEGLIGIBLE"' in content
        assert "1 total, 1 sampled" in content
        assert (
            '<div class="stat warn"><b>0</b><span>images needing review</span>'
            in content
        )
        assert (
            '<div class="stat"><b>0</b><span>other comparison findings</span>'
            in content
        )

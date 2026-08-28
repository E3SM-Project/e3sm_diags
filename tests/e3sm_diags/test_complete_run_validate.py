"""Tests for the combined complete-run validation entry point."""

from __future__ import annotations

from argparse import Namespace
from pathlib import Path

from tests.complete_run import compare, run, validate


def test_validate_runs_candidate_then_compares_with_forwarded_options(
    tmp_path: Path, monkeypatch
):
    results_dir = tmp_path / "candidate"
    baseline_dir = tmp_path / "baseline"
    captured_run_args: list[Namespace] = []
    captured_compare_argv: list[list[str]] = []

    def fake_run(args: Namespace):
        captured_run_args.append(args)
        results_dir.mkdir()
        return []

    def fake_compare(argv: list[str]) -> int:
        captured_compare_argv.append(argv)
        return 0

    monkeypatch.setattr(run, "_run_complete_run", fake_run)
    monkeypatch.setattr(compare, "main", fake_compare)

    result = validate.main(
        [
            "--results-dir",
            str(results_dir),
            "--baseline-dir",
            str(baseline_dir),
            "--set",
            "lat_lon",
            "--atol",
            "0.25",
            "--rtol",
            "0.5",
            "--mode",
            "images",
            "--image-mismatch-threshold",
            "0.75",
            "--diff-artifact-dir",
            str(tmp_path / "artifacts"),
            "--report-dir",
            str(tmp_path / "reports"),
        ]
    )

    assert result == 0
    assert captured_run_args[0].sets_to_run == ["lat_lon"]
    assert captured_compare_argv == [
        [
            "--dev-dir",
            str(results_dir),
            "--baseline-dir",
            str(baseline_dir),
            "--atol",
            "0.25",
            "--rtol",
            "0.5",
            "--image-mismatch-threshold",
            "0.75",
            "--write-diff-pngs",
            "--mode",
            "images",
            "--show",
            "all",
            "--diff-artifact-dir",
            str(tmp_path / "artifacts"),
            "--report-dir",
            str(tmp_path / "reports"),
        ]
    ]


def test_validate_returns_comparison_failure_and_preserves_candidate(
    tmp_path: Path, monkeypatch
):
    results_dir = tmp_path / "candidate"

    def fake_run(_: Namespace):
        results_dir.mkdir()
        (results_dir / "completed-output.nc").touch()
        return []

    monkeypatch.setattr(run, "_run_complete_run", fake_run)
    monkeypatch.setattr(compare, "main", lambda _: 1)

    result = validate.main(["--results-dir", str(results_dir)])

    assert result == 1
    assert (results_dir / "completed-output.nc").exists()

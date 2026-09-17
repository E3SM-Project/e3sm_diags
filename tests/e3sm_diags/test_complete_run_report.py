"""Tests for automated complete-run report rendering."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from tests.complete_run import report


def _status(tmp_path: Path, stage: str) -> Path:
    result_dir = tmp_path / "www" / "complete" / "run"
    result_dir.mkdir(parents=True)
    path = tmp_path / "status.json"
    path.write_text(
        json.dumps(
            {
                "stage": stage,
                "git_sha": "abc",
                "selected_sets": ["lat_lon"],
                "result_dir": str(result_dir),
                "environment_name": "ci_abc",
                "environment_prefix": "/scratch/ci_abc",
            }
        ),
        encoding="utf-8",
    )
    return path


def _comparison(tmp_path: Path, exit_code: int) -> Path:
    path = tmp_path / "www" / "complete" / "comparison-report.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            {
                "exit_code": exit_code,
                "status": "passed" if exit_code == 0 else "failed",
                "summary": {"missing_dev_files": ["a"]},
            }
        ),
        encoding="utf-8",
    )
    return path


@pytest.mark.parametrize(
    ("stage", "comparison_exit", "expected"),
    [
        ("passed", 0, "passed"),
        ("comparison_failed", 1, "comparison_failed"),
        ("diagnostics_failed", None, "diagnostics_failed"),
        ("timed_out", None, "incomplete"),
    ],
)
def test_render_report_terminal_statuses(
    tmp_path: Path, stage: str, comparison_exit: int | None, expected: str
):
    comparison = (
        _comparison(tmp_path, comparison_exit) if comparison_exit is not None else None
    )

    rendered = report.render_report(
        _status(tmp_path, stage),
        comparison,
        cfs_root=tmp_path / "www",
        portal_root="https://portal.example",
    )

    assert rendered["status"] == expected
    assert rendered["paths"]["result_url"] == "https://portal.example/complete/run"
    assert list(rendered["comparison"]["failure_counts"]) == list(
        report.FAILURE_CATEGORIES
    )


def test_missing_comparison_is_incomplete_and_output_is_stable(tmp_path: Path):
    rendered = report.render_report(
        _status(tmp_path, "passed"),
        None,
        cfs_root=tmp_path / "www",
        portal_root="https://portal.example",
    )
    json_path, markdown_path = report.write_report(rendered, tmp_path / "output")

    assert rendered["status"] == "incomplete"
    assert json.loads(json_path.read_text(encoding="utf-8"))["status"] == "incomplete"
    assert "Failures require human review" in markdown_path.read_text(encoding="utf-8")


def test_public_url_returns_none_outside_configured_cfs_root(tmp_path: Path):
    assert (
        report.public_url(
            tmp_path / "private", tmp_path / "www", "https://portal.example"
        )
        is None
    )

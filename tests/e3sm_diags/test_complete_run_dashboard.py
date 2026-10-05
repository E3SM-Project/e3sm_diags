"""Regression tests for the readable, read-only operations dashboard."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

from tests.complete_run import ops, scrontab
from tests.e3sm_diags.test_complete_run_scrontab import _config


def _managed(entries: str) -> str:
    """Wrap sample entries in the actual managed-block markers."""
    return f"{scrontab._MANAGED_BEGIN}\n{entries}\n{scrontab._MANAGED_END}\n"


def test_dashboard_summarizes_schedule_jobs_and_metadata(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    config = _config(tmp_path)
    schedule = scrontab.validate_config(config)
    monkeypatch.setattr(
        scrontab,
        "_read_scrontab",
        lambda: "0 0 * * * /unrelated/job\n" + _managed(schedule),
    )
    calls: list[list[str]] = []

    def inspect(arguments: list[str]) -> str:
        calls.append(arguments)
        if arguments[0] == "squeue":
            return (
                "59173844|/a/long/path/complete-run-controller.sh|PENDING|2026-10-11T06:00:00\n"
                "59173847|/a/long/path/complete-run-reporter.sh|PENDING|N/A\n"
                "123|other-cron-job|RUNNING|2026-10-05T09:00:00"
            )
        return "42 COMPLETED 0:0"

    monkeypatch.setattr(ops, "_inspect_command", inspect)
    run = tmp_path / "results" / "automation" / "abc123-20261001-010000"
    run.mkdir(parents=True)
    status = {
        "stage": "comparison_failed",
        "git_sha": "abc123",
        "job_id": "42",
        "result_dir": "/results/abc123",
        "environment_name": "ci-env",
        "selected_sets": ["lat_lon", "qbo"],
        "error": "first line\nsecond line",
    }
    (run / "status.json").write_text(json.dumps(status), encoding="utf-8")
    (run / "automation-report.json").write_text(
        json.dumps(
            {
                "status": "comparison_failed",
                "publication": {
                    "status": "published",
                    "discussion_url": "https://example.org/discussion/1",
                },
            }
        ),
        encoding="utf-8",
    )
    before = sorted(tmp_path.rglob("*"))
    ops.dashboard(config)
    output = capsys.readouterr().out
    assert "Deployment\n----------" in output
    assert "Configuration health: valid" in output
    assert "Repository: [missing]" in output
    assert "controller" in output and "Sunday 13:00" in output
    assert "reporter" in output and "Monday 17:00" in output
    assert "#SCRON" not in output and "# BEGIN" not in output
    assert "/unrelated/job" not in output
    assert "/a/long/path/" not in output
    assert "controller  " in output and "PENDING  " in output
    assert "other-cron-job" in output
    assert "Eligible time (Slurm)" in output
    assert "Recorded outcome: comparison_failed" in output
    assert "Diagnostics: 2 selected sets" in output
    assert "    second line" in output
    assert "42 COMPLETED 0:0" in output
    assert "Publication: published" in output
    assert "Discussion: https://example.org/discussion/1" in output
    assert "Not present: publication-receipt.json" not in output
    assert '{"' not in output
    assert sorted(tmp_path.rglob("*")) == before
    assert calls[0][-1] == "JobID:0|,Name:0|,State:0|,EligibleTime:0"


@pytest.mark.parametrize(
    "schedule,expected",
    [
        ("", "No managed complete-run schedule installed"),
        ("0 0 * * * /unrelated/job", "No managed complete-run schedule installed"),
        (scrontab._MANAGED_BEGIN, "Invalid managed schedule markers"),
        (
            f"{scrontab._MANAGED_END}\n{scrontab._MANAGED_BEGIN}",
            "Invalid managed schedule markers",
        ),
        (_managed("# comments only"), "No recognized cron entries"),
        (_managed("@daily /custom/job"), "Unrecognized schedule entry"),
        (_managed("0 0 * * * 'unterminated"), "Unrecognized schedule entry"),
        (_managed("30 3 * * 2 /custom/job --flag"), "/custom/job --flag"),
    ],
)
def test_schedule_edge_cases(
    schedule: str,
    expected: str,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    monkeypatch.setattr(scrontab, "_read_scrontab", lambda: schedule)
    ops._installed_schedule()
    assert expected in capsys.readouterr().out


def test_schedule_tool_failure(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    def fail() -> str:
        raise subprocess.CalledProcessError(1, ["scrontab", "-l"])

    monkeypatch.setattr(scrontab, "_read_scrontab", fail)
    ops._installed_schedule()
    assert "Schedule unavailable" in capsys.readouterr().out


@pytest.mark.parametrize(
    "result,expected",
    [
        ("(none)", "No cron jobs in the queue"),
        (
            "Unavailable: squeue is not installed",
            "Unavailable: squeue is not installed",
        ),
        ("unexpected output", "Unrecognized Slurm output: unexpected output"),
    ],
)
def test_queue_empty_or_unavailable(
    result: str,
    expected: str,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    monkeypatch.setattr(ops, "_inspect_command", lambda arguments: result)
    ops._cron_jobs()
    assert expected in capsys.readouterr().out


def test_publication_markers_remain_visible(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    (tmp_path / "publication-receipt.json").write_text("{bad", encoding="utf-8")
    (tmp_path / "publication-failure.json").write_text(
        '{"status": "failed", "discussion_url": null}', encoding="utf-8"
    )
    ops._publication_summary(tmp_path)
    output = capsys.readouterr().out
    assert "Report: Not present: automation-report.json" in output
    assert "Publication receipt: Malformed or unreadable metadata" in output
    assert "Publication failure: failed" in output
    assert "a successful run need not have a Discussion" in output


def test_invalid_publication_field(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    (tmp_path / "automation-report.json").write_text(
        '{"status": "passed", "publication": []}', encoding="utf-8"
    )
    ops._publication_summary(tmp_path)
    assert "missing or invalid report publication field" in capsys.readouterr().out


@pytest.mark.parametrize("configured", [False, True])
def test_no_latest_run_uses_labeled_fields(
    configured: bool, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    ops._latest_run({"RESULTS_ROOT": str(tmp_path)} if configured else {})
    output = capsys.readouterr().out
    expected = (
        f"none in {tmp_path / 'automation'}"
        if configured
        else "none (RESULTS_ROOT is missing)"
    )
    assert output == f"  Run: {expected}\n"

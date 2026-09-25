"""Finalize and publish reports for completed complete-run Slurm jobs."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path
from typing import Any, Sequence

from e3sm_diags.logger import _setup_child_logger
from tests.complete_run.report import (
    DEFAULT_CFS_ROOT,
    DEFAULT_PORTAL_ROOT,
    publish_discussion,
    render_report,
    write_report,
)

logger = _setup_child_logger(__name__)


def report_runs(args: argparse.Namespace) -> int:
    """Process completed automation runs, leaving active and unready runs alone."""
    automation_root = args.results_root / "automation"
    if not automation_root.is_dir():
        return 0
    for run_root in automation_root.iterdir():
        if not run_root.is_dir():
            continue
        try:
            if (run_root / "automation-report.json").is_file():
                _retry_publication(run_root, args)
            else:
                _report_run(run_root, args)
        except (OSError, ValueError, json.JSONDecodeError) as error:
            logger.warning(
                "Unable to report complete-run artifacts in %s: %s", run_root, error
            )
    return 0


def _report_run(run_root: Path, args: argparse.Namespace) -> None:
    status_path = run_root / "status.json"
    status = _load_status(status_path)
    if status is None or status.get("stage") == "submitted":
        if status is None or not _finalize_submitted(status_path, status):
            return
    comparison_report = _comparison_report(run_root / "comparison")
    report = render_report(
        status_path,
        comparison_report,
        cfs_root=args.cfs_root,
        portal_root=args.portal_root,
    )
    write_report(report, run_root)
    if comparison_report is None or not _has_comparison_failures(comparison_report):
        return
    _publish_comparison_failure(run_root, comparison_report, args)
    write_report(
        render_report(
            status_path,
            comparison_report,
            cfs_root=args.cfs_root,
            portal_root=args.portal_root,
        ),
        run_root,
    )


def _retry_publication(run_root: Path, args: argparse.Namespace) -> None:
    """Retry only an unpublished comparison failure from an existing report."""
    comparison_report = _comparison_report(run_root / "comparison")
    if comparison_report is not None and _has_comparison_failures(comparison_report):
        _publish_comparison_failure(run_root, comparison_report, args)
        write_report(
            render_report(
                run_root / "status.json",
                comparison_report,
                cfs_root=args.cfs_root,
                portal_root=args.portal_root,
            ),
            run_root,
        )


def _publish_comparison_failure(
    run_root: Path, comparison_report: Path, args: argparse.Namespace
) -> None:
    """Publish once; the receipt makes later reporting invocations idempotent."""
    receipt = comparison_report.parent / "publication-receipt.json"
    if receipt.is_file():
        return
    try:
        publish_discussion(
            run_root / "automation-report.md",
            receipt,
            repository_id=args.repository_id,
            category_id=args.category_id,
            token_path=args.token_file,
            title=_report_title(run_root),
        )
    except (OSError, RuntimeError, ValueError):
        logger.warning("Unable to publish complete-run Discussion for %s", run_root)
        return


def _finalize_submitted(status_path: Path, status: dict[str, Any]) -> bool:
    """Classify a departed submitted job, deferring while accounting catches up."""
    job_id = status.get("job_id")
    if not isinstance(job_id, str) or not job_id:
        return False
    try:
        queued = _command(["squeue", "-h", "-j", job_id])
    except (OSError, subprocess.CalledProcessError) as error:
        logger.warning(
            "Unable to query Slurm queue for complete-run job %s: %s", job_id, error
        )
        return False
    if queued:
        return False
    try:
        state = (
            _command(
                ["sacct", "-j", job_id, "--format=State", "--noheader", "--parsable2"]
            )
            .splitlines()[0]
            .split("|")[0]
        )
    except (OSError, subprocess.CalledProcessError, IndexError):
        logger.warning("Slurm accounting is not ready for complete-run job %s", job_id)
        return False
    status["stage"] = _terminal_stage(state)
    _write_json(status_path, status)
    return True


def _terminal_stage(state: str) -> str:
    """Map Slurm accounting state to a final automation stage."""
    if state.startswith("CANCELLED"):
        return "cancelled"
    if state.startswith("TIMEOUT"):
        return "timed_out"
    if state.startswith("COMPLETED"):
        return "job_completed_without_status"
    return "slurm_failed"


def _comparison_report(directory: Path) -> Path | None:
    return (
        next(directory.glob("*/comparison-report.json"), None)
        if directory.is_dir()
        else None
    )


def _has_comparison_failures(path: Path) -> bool:
    payload = json.loads(path.read_text(encoding="utf-8"))
    return payload.get("summary", {}).get("failure_count", 0) > 0


def _load_status(path: Path) -> dict[str, Any] | None:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    return payload if isinstance(payload, dict) else None


def _report_title(run_root: Path) -> str:
    """Read the deterministic Discussion title from an automation report."""
    report = _load_status(run_root / "automation-report.json")
    title = report.get("title") if report is not None else None
    return title if isinstance(title, str) else "E3SM Diags complete-run report"


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _command(command: list[str]) -> str:
    return subprocess.run(
        command, check=True, capture_output=True, text=True
    ).stdout.strip()


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-root", required=True, type=Path)
    parser.add_argument("--repository-id", required=True)
    parser.add_argument("--category-id", required=True)
    parser.add_argument("--token-file", required=True, type=Path)
    parser.add_argument("--cfs-root", type=Path, default=DEFAULT_CFS_ROOT)
    parser.add_argument("--portal-root", default=DEFAULT_PORTAL_ROOT)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    return report_runs(_build_parser().parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())

"""Render durable, deterministic reports for automated complete runs."""

from __future__ import annotations

import argparse
import json
from collections import OrderedDict
from pathlib import Path
from typing import Any, Sequence

from tests.complete_run.baseline import _MANIFEST_FILENAME

DEFAULT_CFS_ROOT = Path("/global/cfs/cdirs/e3sm/www")
DEFAULT_PORTAL_ROOT = "https://portal.nersc.gov/cfs/e3sm"
FAILURE_CATEGORIES = (
    "missing_dev_files",
    "missing_baseline_files",
    "missing_variables",
    "nan_location_mismatches",
    "shape_mismatches",
    "tolerance_failures",
    "missing_dev_images",
    "missing_baseline_images",
    "image_mismatches",
)


def public_url(path: str | Path, cfs_root: Path, portal_root: str) -> str | None:
    """Map a CFS artifact path to its NERSC Portal URL, when public."""
    try:
        relative_path = Path(path).resolve().relative_to(cfs_root.resolve())
    except ValueError:
        return None
    return f"{portal_root.rstrip('/')}/{relative_path.as_posix()}"


def render_report(
    status_path: Path,
    comparison_report_path: Path | None,
    *,
    cfs_root: Path = DEFAULT_CFS_ROOT,
    portal_root: str = DEFAULT_PORTAL_ROOT,
) -> dict[str, Any]:
    """Build a report even if the diagnostic or comparison artifact is absent."""
    status = _load_json(status_path, "orchestration status")
    result_dir = Path(status["result_dir"])
    comparison = (
        _load_json(comparison_report_path, "comparison report")
        if comparison_report_path is not None and comparison_report_path.is_file()
        else None
    )
    manifest_path = result_dir / _MANIFEST_FILENAME
    manifest = (
        _load_json(manifest_path, "run manifest") if manifest_path.is_file() else None
    )
    stage = status.get("stage", "unknown")
    status_value = _report_status(stage, comparison)
    comparison_summary = comparison.get("summary", {}) if comparison else {}
    failure_counts = OrderedDict(
        (category, len(comparison_summary.get(category, [])))
        for category in FAILURE_CATEGORIES
    )
    paths = {
        "result_dir": str(result_dir),
        "result_url": public_url(result_dir, cfs_root, portal_root),
        "comparison_report": str(comparison_report_path)
        if comparison_report_path is not None
        else None,
        "comparison_url": public_url(comparison_report_path, cfs_root, portal_root)
        if comparison_report_path is not None
        else None,
        "status": str(status_path),
        "status_url": public_url(status_path, cfs_root, portal_root),
    }
    return {
        "schema_version": 1,
        "status": status_value,
        "git_sha": status.get("git_sha"),
        "selected_sets": status.get("selected_sets", []),
        "environment": {
            "name": status.get("environment_name"),
            "prefix": status.get("environment_prefix"),
            "provenance": str(result_dir / "prov" / "environment.yml"),
            "manifest_environment": manifest.get("environment") if manifest else None,
        },
        "orchestration": status,
        "comparison": {
            "status": comparison.get("status") if comparison else "not-produced",
            "exit_code": comparison.get("exit_code") if comparison else None,
            "failure_counts": failure_counts,
        },
        "paths": paths,
    }


def _report_status(stage: str, comparison: dict[str, Any] | None) -> str:
    if stage in {"submission_failed", "cancelled", "timed_out"}:
        return "incomplete"
    if stage == "diagnostics_failed":
        return "diagnostics_failed"
    if comparison is None:
        return "incomplete"
    return "passed" if comparison.get("exit_code") == 0 else "comparison_failed"


def write_report(report: dict[str, Any], output_dir: Path) -> tuple[Path, Path]:
    """Write sorted JSON and stable Markdown report artifacts."""
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "automation-report.json"
    markdown_path = output_dir / "automation-report.md"
    json_path.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    markdown_path.write_text(_render_markdown(report), encoding="utf-8")
    return json_path, markdown_path


def _render_markdown(report: dict[str, Any]) -> str:
    paths = report["paths"]
    comparison = report["comparison"]
    lines = [
        "# E3SM Diags complete-run report",
        "",
        f"- Status: **{report['status']}**",
        f"- Git SHA: `{report['git_sha']}`",
        f"- Selected sets: {', '.join(report['selected_sets']) or 'all'}",
        f"- Result: {_link(paths['result_dir'], paths['result_url'])}",
        f"- Comparison: {_link(paths['comparison_report'], paths['comparison_url'])}",
        "",
        "## Comparison failure counts",
        "",
    ]
    lines.extend(
        f"- {name}: {count}" for name, count in comparison["failure_counts"].items()
    )
    lines.extend(
        ["", "Failures require human review and never promote a baseline.", ""]
    )
    return "\n".join(lines)


def _link(label: str | None, url: str | None) -> str:
    if label is None:
        return "not produced"
    return f"[{label}]({url})" if url else label


def _load_json(path: Path, label: str) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"Invalid {label}: {path}") from error
    if not isinstance(payload, dict):
        raise ValueError(f"Invalid {label}: {path}")
    return payload


def main(argv: Sequence[str] | None = None) -> int:
    """Render an automation report from an orchestration status file."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--status", required=True, type=Path)
    parser.add_argument("--comparison-report", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--cfs-root", type=Path, default=DEFAULT_CFS_ROOT)
    parser.add_argument("--portal-root", default=DEFAULT_PORTAL_ROOT)
    args = parser.parse_args(argv)
    report = render_report(
        args.status,
        args.comparison_report,
        cfs_root=args.cfs_root,
        portal_root=args.portal_root,
    )
    write_report(report, args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

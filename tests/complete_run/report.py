"""Render durable, deterministic reports for automated complete runs."""

from __future__ import annotations

import argparse
import json
import os
from collections import OrderedDict
from pathlib import Path
from typing import Any, Sequence
from urllib import error, request

from tests.complete_run.baseline import _MANIFEST_FILENAME

DEFAULT_CFS_ROOT = Path("/global/cfs/cdirs/e3sm/www")
DEFAULT_PORTAL_ROOT = "https://portal.nersc.gov/cfs/e3sm"
GITHUB_GRAPHQL_URL = "https://api.github.com/graphql"
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
    """Map a public CFS artifact path to its NERSC Portal URL."""
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
    """Build a report even when diagnostic or comparison artifacts are absent."""
    status = _load_json(status_path, "orchestration status")
    comparison = _load_optional_json(comparison_report_path, "comparison report")
    result_dir = Path(status["result_dir"])
    manifest = _load_optional_json(result_dir / _MANIFEST_FILENAME, "run manifest")
    receipt_path = _publication_receipt_path(comparison_report_path)
    receipt = _load_optional_json(receipt_path, "publication receipt")
    publication_failure = _load_optional_json(
        _publication_failure_path(receipt_path), "publication failure"
    )

    return {
        "schema_version": 1,
        "status": _report_status(status.get("stage", "unknown"), comparison),
        "git_sha": status.get("git_sha"),
        "selected_sets": status.get("selected_sets", []),
        "environment": {
            "name": status.get("environment_name"),
            "prefix": status.get("environment_prefix"),
            "provenance": str(result_dir / "prov" / "environment.yml"),
            "manifest_environment": manifest.get("environment") if manifest else None,
        },
        "orchestration": status,
        "comparison": _comparison_summary(comparison),
        "publication": receipt
        or publication_failure
        or {"status": "not-published", "discussion_url": None},
        "paths": _report_paths(
            status_path, comparison_report_path, result_dir, cfs_root, portal_root
        ),
    }


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


def publish_discussion(
    markdown_path: Path,
    receipt_path: Path,
    *,
    repository_id: str,
    category_id: str,
    token_path: Path,
) -> dict[str, str]:
    """Create one Discussion and atomically retain its receipt for retries."""
    if receipt_path.is_file():
        return _load_receipt(receipt_path)

    token = token_path.read_text(encoding="utf-8").strip()
    if not token:
        raise ValueError("Discussion token file is empty.")

    payload = _discussion_payload(
        markdown_path.read_text(encoding="utf-8"), repository_id, category_id
    )
    http_request = request.Request(
        GITHUB_GRAPHQL_URL,
        data=json.dumps(payload).encode("utf-8"),
        headers={
            "Authorization": f"Bearer {token}",
            "Content-Type": "application/json",
        },
        method="POST",
    )
    response_payload = _github_response(http_request)
    receipt = _discussion_receipt(response_payload)
    _write_receipt(receipt_path, receipt)
    return receipt


def main(argv: Sequence[str] | None = None) -> int:
    """Render a report or publish its Markdown to an E3SM Diags Discussion."""
    args = _build_parser().parse_args(argv)
    if args.command == "publish":
        return _publish_command(args)

    report = render_report(
        args.status,
        args.comparison_report,
        cfs_root=args.cfs_root,
        portal_root=args.portal_root,
    )
    write_report(report, args.output_dir)
    return 0


def _comparison_summary(comparison: dict[str, Any] | None) -> dict[str, Any]:
    summary = comparison.get("summary", {}) if comparison else {}
    failure_counts = OrderedDict(
        (category, len(summary.get(category, []))) for category in FAILURE_CATEGORIES
    )
    return {
        "status": comparison.get("status") if comparison else "not-produced",
        "exit_code": comparison.get("exit_code") if comparison else None,
        "failure_counts": failure_counts,
    }


def _report_paths(
    status_path: Path,
    comparison_report_path: Path | None,
    result_dir: Path,
    cfs_root: Path,
    portal_root: str,
) -> dict[str, str | None]:
    """Build local artifact paths with optional public Portal URLs."""
    return {
        "result_dir": str(result_dir),
        "result_url": public_url(result_dir, cfs_root, portal_root),
        "comparison_report": (
            str(comparison_report_path) if comparison_report_path is not None else None
        ),
        "comparison_url": (
            public_url(comparison_report_path, cfs_root, portal_root)
            if comparison_report_path is not None
            else None
        ),
        "status": str(status_path),
        "status_url": public_url(status_path, cfs_root, portal_root),
    }


def _report_status(stage: str, comparison: dict[str, Any] | None) -> str:
    """Map orchestration and comparison outcomes to the report status."""
    if stage in {"submission_failed", "cancelled", "timed_out"}:
        return "incomplete"
    if stage == "diagnostics_failed":
        return "diagnostics_failed"
    if comparison is None:
        return "incomplete"
    if comparison.get("exit_code") == 0:
        return "passed"

    return "comparison_failed"


def _publication_receipt_path(comparison_report_path: Path | None) -> Path | None:
    if comparison_report_path is None:
        return None

    return comparison_report_path.parent / "publication-receipt.json"


def _publication_failure_path(receipt_path: Path | None) -> Path | None:
    if receipt_path is None:
        return None

    return receipt_path.parent / "publication-failure.json"


def _load_optional_json(path: Path | None, label: str) -> dict[str, Any] | None:
    if path is None or not path.is_file():
        return None

    return _load_json(path, label)


def _load_json(path: Path, label: str) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"Invalid {label}: {path}") from error
    if not isinstance(payload, dict):
        raise ValueError(f"Invalid {label}: {path}")

    return payload


def _render_markdown(report: dict[str, Any]) -> str:
    paths = report["paths"]
    comparison = report["comparison"]
    publication_url = report["publication"]["discussion_url"]
    lines = [
        "# E3SM Diags complete-run report",
        "",
        f"- Status: **{report['status']}**",
        f"- Git SHA: `{report['git_sha']}`",
        f"- Selected sets: {', '.join(report['selected_sets']) or 'all'}",
        f"- Result: {_link(paths['result_dir'], paths['result_url'])}",
        f"- Comparison: {_link(paths['comparison_report'], paths['comparison_url'])}",
        f"- Discussion: {_link(publication_url, publication_url)}",
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
    if url is None:
        return label

    return f"[{label}]({url})"


def _discussion_payload(
    body: str, repository_id: str, category_id: str
) -> dict[str, Any]:
    """Build the GraphQL mutation without exposing the authentication token."""
    return {
        "query": (
            "mutation CreateDiscussion($repositoryId: ID!, $categoryId: ID!, "
            "$title: String!, $body: String!) { createDiscussion(input: {repositoryId: "
            "$repositoryId, categoryId: $categoryId, title: $title, body: $body}) "
            "{ discussion { id url } } }"
        ),
        "variables": {
            "repositoryId": repository_id,
            "categoryId": category_id,
            "title": "E3SM Diags complete-run report",
            "body": body,
        },
    }


def _github_response(http_request: request.Request) -> dict[str, Any]:
    """Send a GraphQL request without retaining authentication details."""
    try:
        with request.urlopen(http_request, timeout=30) as response:  # noqa: S310
            return json.loads(response.read().decode("utf-8"))
    except (OSError, error.URLError, json.JSONDecodeError) as exception:
        raise RuntimeError("Unable to publish complete-run Discussion.") from exception


def _discussion_receipt(response_payload: dict[str, Any]) -> dict[str, str]:
    """Extract the immutable Discussion receipt from a GraphQL response."""
    try:
        discussion = response_payload["data"]["createDiscussion"]["discussion"]
        return {
            "status": "published",
            "discussion_id": discussion["id"],
            "discussion_url": discussion["url"],
        }
    except (KeyError, TypeError):
        raise RuntimeError("GitHub returned an invalid Discussion response.") from None


def _load_receipt(receipt_path: Path) -> dict[str, str]:
    receipt = _load_json(receipt_path, "publication receipt")
    keys = ("status", "discussion_id", "discussion_url")
    if not all(isinstance(receipt.get(key), str) and receipt[key] for key in keys):
        raise ValueError(f"Invalid publication receipt: {receipt_path}")

    return {key: receipt[key] for key in keys}


def _write_receipt(receipt_path: Path, receipt: dict[str, str]) -> None:
    """Publish a receipt without replacing a concurrent publisher's record."""
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(receipt_path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
        json.dump(receipt, stream, indent=2, sort_keys=True)
        stream.write("\n")


def _write_publication_failure(receipt_path: Path) -> None:
    """Record a retryable publication failure without persisting API details."""
    failure_path = receipt_path.parent / "publication-failure.json"
    failure_path.write_text(
        json.dumps({"status": "failed", "discussion_url": None}, indent=2) + "\n",
        encoding="utf-8",
    )


def _publish_command(args: argparse.Namespace) -> int:
    """Publish a Discussion and retain a generic retryable failure marker."""
    try:
        publish_discussion(
            args.markdown,
            args.receipt,
            repository_id=args.repository_id,
            category_id=args.category_id,
            token_path=args.token_file,
        )
    except (OSError, RuntimeError, ValueError):
        _write_publication_failure(args.receipt)
        return 1

    return 0


def _build_parser() -> argparse.ArgumentParser:
    """Build rendering and publication subcommands."""
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    render = subparsers.add_parser("render")
    render.add_argument("--status", required=True, type=Path)
    render.add_argument("--comparison-report", type=Path)
    render.add_argument("--output-dir", required=True, type=Path)
    render.add_argument("--cfs-root", type=Path, default=DEFAULT_CFS_ROOT)
    render.add_argument("--portal-root", default=DEFAULT_PORTAL_ROOT)

    publish = subparsers.add_parser("publish")
    publish.add_argument("--markdown", required=True, type=Path)
    publish.add_argument("--receipt", required=True, type=Path)
    publish.add_argument("--repository-id", required=True)
    publish.add_argument("--category-id", required=True)
    publish.add_argument("--token-file", required=True, type=Path)
    return parser


if __name__ == "__main__":
    raise SystemExit(main())

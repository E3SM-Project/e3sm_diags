"""Manage immutable manifests and promoted complete-run baselines.

Use ``python -m tests.complete_run.baseline promote --run-dir PATH --channel
main`` to select a completed run as a baseline.  A non-main run requires the
explicit ``--allow-non-main`` escape hatch.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import subprocess
import tempfile
import uuid
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path
from typing import Any, Sequence

from e3sm_diags.logger import _setup_child_logger, _setup_root_logger
from tests.complete_run.params import DEFAULT_RESULTS_DIR, CompleteRunConfig

logger = _setup_child_logger(__name__)

_MANIFEST_FILENAME = "baseline_manifest.json"
_MANIFEST_SCHEMA_VERSION = 2
_REPO_ROOT = Path(__file__).resolve().parents[2]
_CURATED_PACKAGES = (
    "cartopy",
    "dask",
    "esmf",
    "esmpy",
    "matplotlib",
    "numpy",
    "pandas",
    "proj",
    "pyproj",
    "xarray",
    "xcdat",
    "xesmf",
    "xgcm",
)
_NOT_INSTALLED = "not-installed"


class _ManifestError(ValueError):
    """Raised when a baseline manifest cannot be used."""


def _get_git_metadata() -> dict[str, str]:
    """Return the current Git branch and full commit SHA, with safe fallbacks."""
    return {
        "branch": _get_git_output(["rev-parse", "--abbrev-ref", "HEAD"], "unknown"),
        "sha": _get_git_output(["rev-parse", "HEAD"], "unknown"),
    }


def _collect_package_versions() -> dict[str, str]:
    """Collect curated distribution versions without importing their packages."""
    versions: dict[str, str] = {}
    for package in _CURATED_PACKAGES:
        try:
            versions[package] = metadata.version(package)
        except metadata.PackageNotFoundError:
            versions[package] = _NOT_INSTALLED
    return versions


def _get_environment_metadata() -> dict[str, object]:
    """Return runtime and curated dependency provenance for a complete run."""
    return {
        "python_version": platform.python_version(),
        "python_implementation": platform.python_implementation(),
        "platform": platform.platform(),
        "conda_environment": os.environ.get("CONDA_DEFAULT_ENV", "not-set"),
        "packages": _collect_package_versions(),
    }


def _build_manifest(
    config: CompleteRunConfig,
    selected_sets: Sequence[str],
    workflow_revision: str | None = None,
) -> dict[str, Any]:
    """Build the JSON-serializable manifest for a completed run."""
    paths = config.paths
    git_metadata = _get_git_metadata()
    return {
        "schema_version": _MANIFEST_SCHEMA_VERSION,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "result_dir": paths.results_dir,
        "git": git_metadata,
        "workflow_revision": workflow_revision or git_metadata["sha"],
        "environment": _get_environment_metadata(),
        "test_paths": {
            "climo": paths.test_climo,
            "ts": paths.test_ts,
            "ts_daily_dir": paths.test_ts_daily_dir,
            "diurnal_climo": paths.test_diurnal_climo,
            "streamflow_ts": paths.test_streamflow_ts,
            "tc_analysis": paths.test_tc_analysis,
            "arm_site": paths.test_arm_site,
        },
        "reference_paths": {
            "climo": paths.ref_climo,
            "ts": paths.ref_ts,
            "tc_analysis": paths.ref_tc_analysis,
            "arm": paths.ref_arm,
        },
        "config": {
            "case": config.case,
            "short_name": config.short_name,
            "start_year": config.start_yr,
            "end_year": config.end_yr,
            "num_workers": config.num_workers,
            "save_netcdf": config.save_netcdf,
            "selected_sets": list(selected_sets),
        },
    }


def _write_manifest(run_dir: str | Path, manifest: dict[str, Any]) -> Path:
    """Write a manifest exactly once, preserving completed-run provenance.

    Raises
    ------
    FileNotFoundError
        If the completed run directory does not exist.
    FileExistsError
        If a manifest is already present.
    """
    directory = Path(run_dir)
    if not directory.is_dir():
        raise FileNotFoundError(
            f"Complete-run results directory does not exist: {directory}"
        )

    manifest_path = directory / _MANIFEST_FILENAME
    if os.path.lexists(manifest_path):
        raise FileExistsError(
            f"Refusing to replace immutable manifest: {manifest_path}"
        )

    file_descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{_MANIFEST_FILENAME}.", suffix=".tmp", dir=directory
    )
    temporary_path = Path(temporary_name)
    try:
        with os.fdopen(file_descriptor, "w", encoding="utf-8") as stream:
            json.dump(manifest, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        # Results are shared on CFS, so the immutable provenance record must
        # remain readable after the temporary file is published.
        os.chmod(temporary_path, 0o644)
        try:
            # ``link`` publishes the fully written file without replacing a
            # concurrent writer's manifest. Unlike ``os.replace``, it fails
            # with FileExistsError when the immutable destination exists.
            os.link(temporary_path, manifest_path)
        except FileExistsError:
            raise FileExistsError(
                f"Refusing to replace immutable manifest: {manifest_path}"
            ) from None
    except BaseException:
        temporary_path.unlink(missing_ok=True)
        raise
    temporary_path.unlink()
    return manifest_path


def _load_manifest(run_dir: str | Path) -> dict[str, Any]:
    """Load and validate a complete-run manifest."""
    directory = Path(run_dir)
    if not directory.is_dir():
        raise FileNotFoundError(
            f"Run directory does not exist or is not a directory: {directory}"
        )
    manifest_path = directory / _MANIFEST_FILENAME
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Run manifest does not exist: {manifest_path}")
    try:
        with manifest_path.open(encoding="utf-8") as stream:
            manifest = json.load(stream)
    except json.JSONDecodeError as error:
        raise _ManifestError(
            f"Run manifest is not valid JSON: {manifest_path}"
        ) from error
    _validate_manifest(manifest, manifest_path)
    return manifest


def promote_baseline(
    run_dir: str | Path,
    channel: str,
    *,
    allow_non_main: bool = False,
    results_root: str | Path = DEFAULT_RESULTS_DIR,
) -> Path:
    """Atomically promote a valid run directory to a channel's latest link."""
    manifest = _load_manifest(run_dir)
    branch = manifest["git"]["branch"]
    if channel == "main" and branch != "main" and not allow_non_main:
        raise _ManifestError(
            "Refusing to promote a non-main manifest to channel 'main' "
            f"(manifest branch: {branch!r}). Use --allow-non-main to override."
        )

    root = Path(results_root)
    if not root.is_dir():
        raise FileNotFoundError(f"Baseline results root does not exist: {root}")
    target = Path(run_dir).resolve()
    link = root / f"latest-{channel}"
    temporary_link = root / f".{link.name}.{os.getpid()}.{uuid.uuid4().hex}.tmp"
    try:
        temporary_link.symlink_to(target, target_is_directory=True)
        os.replace(temporary_link, link)
    finally:
        if temporary_link.is_symlink():
            temporary_link.unlink()
    logger.info("Promoted %s to %s", target, link)
    return link


def show_baseline(
    channel: str, *, results_root: str | Path = DEFAULT_RESULTS_DIR
) -> tuple[Path, dict[str, Any]]:
    """Resolve a promoted baseline and return its target and manifest."""
    link = Path(results_root) / f"latest-{channel}"
    if not link.is_symlink():
        raise FileNotFoundError(f"Baseline link does not exist: {link}")
    target = link.resolve(strict=False)
    if not target.is_dir():
        raise FileNotFoundError(f"Baseline link is broken: {link} -> {target}")
    return target, _load_manifest(target)


def _validate_manifest(manifest: object, manifest_path: Path) -> None:
    if not isinstance(manifest, dict):
        raise _ManifestError(
            f"Run manifest must contain a JSON object: {manifest_path}"
        )
    required = {
        "schema_version",
        "created_at_utc",
        "result_dir",
        "git",
        "workflow_revision",
        "environment",
        "test_paths",
        "reference_paths",
        "config",
    }
    missing = required - manifest.keys()
    if missing or manifest.get("schema_version") != _MANIFEST_SCHEMA_VERSION:
        raise _ManifestError(
            f"Run manifest has an unsupported or incomplete schema: {manifest_path}"
        )
    git = manifest["git"]
    config = manifest["config"]
    if (
        not isinstance(git, dict)
        or not isinstance(git.get("branch"), str)
        or not isinstance(git.get("sha"), str)
    ):
        raise _ManifestError(f"Run manifest has invalid Git metadata: {manifest_path}")
    if (
        not isinstance(manifest["workflow_revision"], str)
        or not manifest["workflow_revision"].strip()
    ):
        raise _ManifestError(
            f"Run manifest has invalid workflow revision: {manifest_path}"
        )
    environment = manifest["environment"]
    if not isinstance(environment, dict):
        raise _ManifestError(
            f"Run manifest has invalid environment metadata: {manifest_path}"
        )
    required_environment = {
        "python_version",
        "python_implementation",
        "platform",
        "conda_environment",
        "packages",
    }
    packages = environment.get("packages")
    if (
        required_environment - environment.keys()
        or any(
            not isinstance(environment.get(key), str)
            for key in required_environment - {"packages"}
        )
        or not isinstance(packages, dict)
        or set(packages) != set(_CURATED_PACKAGES)
        or not all(
            isinstance(version, str) and version for version in packages.values()
        )
    ):
        raise _ManifestError(
            f"Run manifest has invalid environment metadata: {manifest_path}"
        )
    if not isinstance(config, dict) or not isinstance(
        config.get("selected_sets"), list
    ):
        raise _ManifestError(
            f"Run manifest has invalid configuration metadata: {manifest_path}"
        )


def _get_git_output(args: list[str], fallback: str) -> str:
    try:
        completed = subprocess.run(
            ["git", *args], check=True, capture_output=True, text=True, cwd=_REPO_ROOT
        )
    except (OSError, subprocess.CalledProcessError):
        return fallback
    return completed.stdout.strip() or fallback


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Manage complete-run baselines.")
    subparsers = parser.add_subparsers(dest="command", required=True)
    promote = subparsers.add_parser("promote", help="Promote a completed run.")
    promote.add_argument("--run-dir", required=True, help="Completed run directory.")
    promote.add_argument(
        "--channel", required=True, help="Baseline channel (for example main)."
    )
    promote.add_argument(
        "--allow-non-main",
        action="store_true",
        help="Allow a non-main manifest on the main channel.",
    )
    show = subparsers.add_parser("show", help="Show the promoted baseline.")
    show.add_argument(
        "--channel", required=True, help="Baseline channel (for example main)."
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entrypoint for baseline promotion and inspection."""
    _setup_root_logger()
    args = _build_parser().parse_args(argv)
    try:
        if args.command == "promote":
            promote_baseline(
                args.run_dir, args.channel, allow_non_main=args.allow_non_main
            )
        else:
            target, manifest = show_baseline(args.channel)
            logger.info("Baseline target: %s", target)
            logger.info("Manifest:\n%s", json.dumps(manifest, indent=2, sort_keys=True))
    except (FileNotFoundError, _ManifestError, OSError) as error:
        logger.error("%s", error)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

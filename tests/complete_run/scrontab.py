"""Create, validate, and install complete-run NERSC scrontab configuration."""

from __future__ import annotations

import argparse
import os
import subprocess
from pathlib import Path
from typing import Sequence

_ROOT = Path(__file__).parent
_CONFIG_TEMPLATE = _ROOT / "complete-run-controller.env.template"
_SCRONTAB_TEMPLATE = _ROOT / "complete-run.scrontab.template"
_REQUIRED_CONFIG_KEYS = (
    "REPOSITORY",
    "LOG_DIR",
    "CONDA_BASE",
    "CONTROLLER_ENV",
    "RESULTS_ROOT",
    "SLURM_ACCOUNT",
    "SIMBOARD_REPOSITORY_ID",
    "SIMBOARD_CATEGORY_ID",
    "SIMBOARD_TOKEN_FILE",
)


def create_config(config_path: Path) -> Path:
    """Create a mode-0600 controller configuration from the tracked template.

    Raises
    ------
    FileExistsError
        If the target configuration file already exists.
    """
    if config_path.exists() or config_path.is_symlink():
        raise FileExistsError(
            f"Refusing to replace controller configuration: {config_path}"
        )

    config_path.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text(
        _CONFIG_TEMPLATE.read_text(encoding="utf-8"), encoding="utf-8"
    )
    os.chmod(config_path, 0o600)
    return config_path


def validate_config(config_path: Path) -> str:
    """Validate controller configuration and return a rendered scrontab.

    Raises
    ------
    ValueError
        If required values are absent, unresolved, or not absolute paths.
    """
    config = _read_config(config_path)
    _validate_config_values(config)
    rendered = _render_scrontab(config, config_path)
    if "{{" in rendered or "}}" in rendered:
        raise ValueError("Scrontab template contains unresolved placeholders.")
    return rendered


def install_scrontab(config_path: Path) -> None:
    """Validate and install the complete-run schedule through ``scrontab``."""
    rendered = validate_config(config_path)
    subprocess.run(["scrontab"], input=rendered, text=True, check=True)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the configuration creation, validation, or installation CLI."""
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    for command in ("create-config", "validate", "install"):
        subparser = subparsers.add_parser(command)
        subparser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args(argv)

    try:
        if args.command == "create-config":
            create_config(args.config)
        elif args.command == "validate":
            validate_config(args.config)
        else:
            install_scrontab(args.config)
    except (OSError, subprocess.CalledProcessError, ValueError) as error:
        parser.error(str(error))

    return 0


def _read_config(config_path: Path) -> dict[str, str]:
    """Read a deliberately simple ``KEY=VALUE`` controller configuration."""
    config: dict[str, str] = {}
    for line_number, line in enumerate(
        config_path.read_text(encoding="utf-8").splitlines(), start=1
    ):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if "=" not in stripped:
            raise ValueError(f"Invalid configuration line {line_number}: {config_path}")
        key, value = stripped.split("=", maxsplit=1)
        if not key.isidentifier():
            raise ValueError(f"Invalid configuration key on line {line_number}: {key}")
        config[key] = value
    return config


def _validate_config_values(config: dict[str, str]) -> None:
    """Reject missing template values and relative operational paths."""
    missing = [key for key in _REQUIRED_CONFIG_KEYS if not config.get(key)]
    unresolved = [
        key
        for key in _REQUIRED_CONFIG_KEYS
        if config.get(key, "").startswith("replace-with-")
    ]
    if missing or unresolved:
        values = ", ".join([*missing, *unresolved])
        raise ValueError(f"Configuration has missing or unresolved values: {values}")

    for key in (
        "REPOSITORY",
        "LOG_DIR",
        "CONDA_BASE",
        "RESULTS_ROOT",
        "SIMBOARD_TOKEN_FILE",
    ):
        if not Path(config[key]).is_absolute():
            raise ValueError(f"Configuration value must be an absolute path: {key}")


def _render_scrontab(config: dict[str, str], config_path: Path) -> str:
    """Replace the explicit deployment placeholders in the versioned template."""
    replacements = {
        "{{ACCOUNT}}": config["SLURM_ACCOUNT"],
        "{{LOG_DIR}}": config["LOG_DIR"],
        "{{REPOSITORY}}": config["REPOSITORY"],
        "{{CONFIG_FILE}}": str(config_path.resolve()),
    }
    rendered = _SCRONTAB_TEMPLATE.read_text(encoding="utf-8")
    for placeholder, value in replacements.items():
        rendered = rendered.replace(placeholder, value)
    return rendered


if __name__ == "__main__":
    raise SystemExit(main())

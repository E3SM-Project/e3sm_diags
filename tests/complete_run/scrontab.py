"""Create, validate, and install complete-run NERSC scrontab configuration."""

from __future__ import annotations

import argparse
import fcntl
import json
import os
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence, TextIO

_ROOT = Path(__file__).parent
_CONFIG_TEMPLATE = _ROOT / "complete-run-controller.env.template"
_SCRONTAB_TEMPLATE = _ROOT / "complete-run.scrontab.template"
_DEFAULT_CONTROLLER_ENV_PREFIX = (
    "/global/cfs/projectdirs/e3sm/e3sm_diags/operations/controller-env"
)
_REQUIRED_CONFIG_KEYS = (
    "REPOSITORY",
    "LOG_DIR",
    "CONDA_BASE",
    "CONTROLLER_ENV_PREFIX",
    "RESULTS_ROOT",
    "SLURM_ACCOUNT",
    "SCRON_CPUS",
    "SCRON_MEMORY_PER_CPU",
    "E3SM_DIAGS_REPOSITORY_ID",
    "E3SM_DIAGS_CATEGORY_ID",
    "E3SM_DIAGS_TOKEN_FILE",
)


def create_config(config_path: Path, controller_env_prefix: Path | None = None) -> Path:
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
    content = _CONFIG_TEMPLATE.read_text(encoding="utf-8")
    if controller_env_prefix is not None:
        _validate_absolute_path("CONTROLLER_ENV_PREFIX", controller_env_prefix)
        content = content.replace(
            f"CONTROLLER_ENV_PREFIX={_DEFAULT_CONTROLLER_ENV_PREFIX}",
            f"CONTROLLER_ENV_PREFIX={controller_env_prefix}",
        )
    config_path.write_text(content, encoding="utf-8")
    os.chmod(config_path, 0o600)
    return config_path


def initialize_operations(
    operations_dir: Path, repository_url: str, branch: str
) -> tuple[Path, Path]:
    """Create an operations layout without changing an existing checkout.

    Parameters
    ----------
    operations_dir : Path
        Persistent, non-public directory for controller administration.
    repository_url : str
        Git remote used only when the controller checkout is absent.
    branch : str
        Branch checked out only in a newly cloned controller checkout.

    Returns
    -------
    tuple[Path, Path]
        The controller checkout and external configuration paths.

    Raises
    ------
    FileExistsError
        If the operations checkout or configuration exists but is not valid for
        a non-destructive bootstrap.
    """
    operations_dir.mkdir(parents=True, exist_ok=True)
    checkout_path = operations_dir / "e3sm_diags"
    config_path = operations_dir / "controller.env"
    (operations_dir / "logs").mkdir(exist_ok=True)
    _create_checkout(checkout_path, repository_url, branch)
    if config_path.exists() or config_path.is_symlink():
        raise FileExistsError(
            f"Refusing to replace controller configuration: {config_path}"
        )
    create_config(config_path, operations_dir / "controller-env")
    return checkout_path, config_path


def create_controller_environment(config_path: Path) -> None:
    """Create the persistent controller prefix from the checkout's CI spec."""
    environment = _controller_environment(config_path)
    if environment.prefix.exists():
        raise FileExistsError(
            f"Refusing to replace controller environment: {environment.prefix}"
        )
    _run_conda(
        environment,
        "env",
        "create",
        "--prefix",
        str(environment.prefix),
        "--file",
        str(environment.specification),
    )
    _install_controller_checkout(environment)


def update_controller_environment(config_path: Path, *, confirmed: bool) -> Path:
    """Export and deliberately update the persistent controller prefix.

    Raises
    ------
    ValueError
        If the required explicit confirmation is absent.
    RuntimeError
        If a scheduled controller currently holds the shared lock.
    """
    if not confirmed:
        raise ValueError("Refusing controller environment update without confirmation.")
    environment = _controller_environment(config_path)
    if not environment.prefix.is_dir():
        raise FileNotFoundError(
            f"Controller environment does not exist: {environment.prefix}"
        )
    with _controller_lock(environment.results_root):
        export_path = _export_controller_environment(environment)
        _run_conda(
            environment,
            "env",
            "update",
            "--prune",
            "--prefix",
            str(environment.prefix),
            "--file",
            str(environment.specification),
        )
        _install_controller_checkout(environment)
        _verify_controller_environment(environment)
    return export_path


def show_controller_environment(config_path: Path) -> dict[str, str]:
    """Return immutable metadata for the configured controller prefix."""
    environment = _controller_environment(config_path)
    completed = _run_conda(
        environment,
        "run",
        "--prefix",
        str(environment.prefix),
        "python",
        "--version",
        capture_output=True,
    )
    return {
        "prefix": str(environment.prefix),
        "specification": str(environment.specification),
        "python_version": completed.stdout.strip(),
    }


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
    initialize = subparsers.add_parser("initialize-operations")
    initialize.add_argument("--operations-dir", required=True, type=Path)
    initialize.add_argument("--repository-url", required=True)
    initialize.add_argument("--branch", required=True)
    create_environment = subparsers.add_parser("create-controller-env")
    create_environment.add_argument("--config", required=True, type=Path)
    update_environment = subparsers.add_parser("update-controller-env")
    update_environment.add_argument("--config", required=True, type=Path)
    update_environment.add_argument("--confirm", action="store_true")
    show_environment = subparsers.add_parser("show-controller-env")
    show_environment.add_argument("--config", required=True, type=Path)
    for command in ("create-config", "validate", "install"):
        subparser = subparsers.add_parser(command)
        subparser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args(argv)

    try:
        if args.command == "initialize-operations":
            initialize_operations(args.operations_dir, args.repository_url, args.branch)
        elif args.command == "create-config":
            create_config(args.config)
        elif args.command == "create-controller-env":
            create_controller_environment(args.config)
        elif args.command == "update-controller-env":
            update_controller_environment(args.config, confirmed=args.confirm)
        elif args.command == "show-controller-env":
            print(json.dumps(show_controller_environment(args.config), indent=2))
        elif args.command == "validate":
            validate_config(args.config)
        else:
            install_scrontab(args.config)
    except (OSError, subprocess.CalledProcessError, ValueError) as error:
        parser.error(str(error))

    return 0


@dataclass(frozen=True)
class _ControllerEnvironment:
    """Resolved configuration for the persistent login-node controller prefix."""

    config_path: Path
    conda_base: Path
    prefix: Path
    repository: Path
    results_root: Path

    @property
    def conda_executable(self) -> Path:
        """Return the configured Conda executable."""
        return self.conda_base / "bin" / "conda"

    @property
    def specification(self) -> Path:
        """Return the controller dependency specification in its checkout."""
        return self.repository / "conda-env" / "ci.yml"


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


def _controller_environment(config_path: Path) -> _ControllerEnvironment:
    """Resolve the configuration required for controller-environment actions."""
    config = _read_config(config_path)
    keys = ("CONDA_BASE", "CONTROLLER_ENV_PREFIX", "REPOSITORY", "RESULTS_ROOT")
    missing = [key for key in keys if not config.get(key)]
    if missing:
        raise ValueError(
            f"Controller configuration is missing values: {', '.join(missing)}"
        )
    unresolved = [key for key in keys if config[key].startswith("replace-with-")]
    if unresolved:
        raise ValueError(
            f"Controller configuration has unresolved values: {', '.join(unresolved)}"
        )
    paths = {key: Path(config[key]) for key in keys}
    for key, path in paths.items():
        _validate_absolute_path(key, path)
    return _ControllerEnvironment(
        config_path=config_path,
        conda_base=paths["CONDA_BASE"],
        prefix=paths["CONTROLLER_ENV_PREFIX"],
        repository=paths["REPOSITORY"],
        results_root=paths["RESULTS_ROOT"],
    )


def _validate_absolute_path(key: str, path: Path) -> None:
    """Reject relative configured paths before they reach shell or Conda."""
    if not path.is_absolute():
        raise ValueError(f"Configuration value must be an absolute path: {key}")


def _run_conda(
    environment: _ControllerEnvironment,
    *arguments: str,
    capture_output: bool = False,
) -> subprocess.CompletedProcess[str]:
    """Run Conda through the configured base installation."""
    return subprocess.run(
        [str(environment.conda_executable), *arguments],
        check=True,
        text=True,
        capture_output=capture_output,
        cwd=environment.repository,
    )


def _install_controller_checkout(environment: _ControllerEnvironment) -> None:
    """Install the configured checkout into its persistent controller prefix."""
    _run_conda(
        environment,
        "run",
        "--prefix",
        str(environment.prefix),
        "python",
        "-m",
        "pip",
        "install",
        "--editable",
        str(environment.repository),
    )


def _export_controller_environment(environment: _ControllerEnvironment) -> Path:
    """Write a timestamped pre-update controller-environment export to CFS."""
    completed = _run_conda(
        environment,
        "env",
        "export",
        "--prefix",
        str(environment.prefix),
        capture_output=True,
    )
    provenance_dir = environment.config_path.parent / "provenance"
    provenance_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
    export_path = provenance_dir / f"controller-env-{timestamp}.yml"
    export_path.write_text(completed.stdout, encoding="utf-8")
    return export_path


def _verify_controller_environment(environment: _ControllerEnvironment) -> None:
    """Ensure the updated prefix can import and invoke the controller CLI."""
    _run_conda(
        environment,
        "run",
        "--prefix",
        str(environment.prefix),
        "python",
        "-m",
        "tests.complete_run.automation",
        "--help",
    )


class _ControllerLock:
    """Hold the scheduler lock while a manual controller update is in progress."""

    def __init__(self, results_root: Path):
        self._path = results_root / "automation" / "controller.lock"
        self._stream: TextIO | None = None

    def __enter__(self) -> _ControllerLock:
        self._path.parent.mkdir(parents=True, exist_ok=True)
        self._stream = self._path.open("a", encoding="utf-8")
        try:
            fcntl.flock(self._stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            self._stream.close()
            raise RuntimeError(
                "A complete-run controller is active; refusing environment update."
            ) from error
        return self

    def __exit__(self, *_: object) -> None:
        if self._stream is not None:
            fcntl.flock(self._stream.fileno(), fcntl.LOCK_UN)
            self._stream.close()


def _controller_lock(results_root: Path) -> _ControllerLock:
    """Return the lock shared by scheduled controllers and manual updates."""
    return _ControllerLock(results_root)


def _create_checkout(checkout_path: Path, repository_url: str, branch: str) -> None:
    """Clone a requested controller branch only when no checkout exists."""
    if checkout_path.exists():
        if (checkout_path / ".git").exists():
            return
        raise FileExistsError(
            "Operations checkout path exists but is not a Git checkout: "
            f"{checkout_path}"
        )
    subprocess.run(
        [
            "git",
            "clone",
            "--branch",
            branch,
            "--single-branch",
            repository_url,
            str(checkout_path),
        ],
        check=True,
    )


def _validate_config_values(config: dict[str, str]) -> None:
    """Reject missing template values and relative operational paths."""
    missing = [key for key in _REQUIRED_CONFIG_KEYS if not config.get(key)]
    unresolved = [
        key
        for key in _REQUIRED_CONFIG_KEYS
        if config.get(key, "").startswith(("replace-with-", "/absolute/path"))
    ]
    if missing or unresolved:
        values = ", ".join([*missing, *unresolved])
        raise ValueError(f"Configuration has missing or unresolved values: {values}")

    for key in (
        "REPOSITORY",
        "LOG_DIR",
        "CONDA_BASE",
        "CONTROLLER_ENV_PREFIX",
        "RESULTS_ROOT",
        "E3SM_DIAGS_TOKEN_FILE",
    ):
        if not Path(config[key]).is_absolute():
            raise ValueError(f"Configuration value must be an absolute path: {key}")


def _render_scrontab(config: dict[str, str], config_path: Path) -> str:
    """Replace the explicit deployment placeholders in the versioned template."""
    replacements = {
        "{{ACCOUNT}}": config["SLURM_ACCOUNT"],
        "{{SCRON_CPUS}}": config["SCRON_CPUS"],
        "{{SCRON_MEMORY_PER_CPU}}": config["SCRON_MEMORY_PER_CPU"],
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

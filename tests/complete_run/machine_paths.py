"""Known deployment paths for complete-run tooling (no filesystem search)."""

from __future__ import annotations

from configparser import Error as ConfigParserError
from dataclasses import dataclass
from pathlib import Path

try:
    from mache import MachineInfo
except ImportError:
    MachineInfo = None


@dataclass(frozen=True)
class MachinePaths:
    """Standard private operations and public results directories."""

    operations_dir: Path
    results_root: Path


NERSC_PATHS = MachinePaths(
    operations_dir=Path("/global/cfs/projectdirs/e3sm/e3sm_diags/operations"),
    results_root=Path("/global/cfs/cdirs/e3sm/www/e3sm_diags/complete-run-test"),
)
MACHINE_PATHS = {"pm-cpu": NERSC_PATHS, "pm-gpu": NERSC_PATHS}


def detect_machine_paths() -> MachinePaths | None:
    """Return known defaults, or None if mache cannot identify a mapped machine.

    The optional import permits explicit operations configuration when mache
    is unavailable; discovery failure likewise leaves overrides usable.
    """
    try:
        if MachineInfo is None:
            return None
        machine = MachineInfo(quiet=True).machine
    except (
        ImportError,
        OSError,
        ValueError,
        RuntimeError,
        KeyError,
        ConfigParserError,
    ):
        return None

    return MACHINE_PATHS.get(machine)


def default_results_root() -> Path:
    """Select machine defaults, retaining the historical results fallback.

    Keeping the fallback permits importing comparison and parameter modules on
    other machines; callers can still supply explicit results paths.
    """
    paths = detect_machine_paths()

    return paths.results_root if paths is not None else NERSC_PATHS.results_root

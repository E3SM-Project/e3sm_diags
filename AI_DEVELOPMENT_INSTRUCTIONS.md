# AI Development Instructions for E3SM Diagnostics

This is the canonical, tool-agnostic source of AI development rules for the
`e3sm_diags` repository. Derived files (`AGENTS.md`,
`.github/copilot-instructions.md`, `.claude/CLAUDE.md`) must reflect this
document and must not introduce rules absent from it.

## Project Overview

E3SM Diagnostics (`e3sm_diags`) is a Python-based scientific diagnostics
package for the DOE Energy Exascale Earth System Model (E3SM). It evaluates
climate model performance by comparing model output against observations and
reanalysis products.

- **License:** BSD 3-Clause
- **Python:** ≥3.11 (see `pyproject.toml` for the full support matrix)
- **Backend:** Xarray and xCDAT
- **Entry points:** `e3sm_diags`, `e3sm_diags_vars`

## Architecture

The package follows a **driver–parameter** architecture:

| Directory               | Purpose                              |
| ----------------------- | ------------------------------------ |
| `e3sm_diags/driver/`    | One driver per diagnostic set        |
| `e3sm_diags/parameter/` | Paired parameter class per driver    |
| `e3sm_diags/derivations/` | Derived variable calculations      |
| `e3sm_diags/metrics/`   | Statistical metric functions         |
| `e3sm_diags/plot/`      | Matplotlib/Cartopy plotting modules  |
| `e3sm_diags/viewer/`    | HTML results viewer                  |
| `e3sm_diags/parser/`    | CLI and configuration parsing        |
| `tests/e3sm_diags/`     | Unit tests (pytest)                  |
| `tests/integration/`    | Integration tests                    |
| `docs/`                 | Sphinx documentation                 |

When adding a new diagnostic set, create both a driver module in `driver/` and
a corresponding parameter class in `parameter/`.

## Coding Standards

### Style and Formatting

- **Formatter/Linter:** Ruff (see `pyproject.toml` for rule configuration)
- **Import sorting:** Ruff isort
- **Line length:** `E501` is ignored, but keep lines reasonable
- **Complexity:** McCabe complexity is enforced (see `pyproject.toml`)
- **Docstrings:** NumPy convention
- Use `from __future__ import annotations` at the top of every module

### Type Annotations

- **Type checker:** mypy (see `pyproject.toml` for configuration)
- Add type hints to all function signatures
- Use `TYPE_CHECKING` guard for imports used only in annotations:
  ```python
  from __future__ import annotations

  from typing import TYPE_CHECKING

  if TYPE_CHECKING:
      import xarray as xr
  ```

### Logging

- Use the project's custom logger, not `print()` or bare `logging`:
  ```python
  from e3sm_diags.logger import _setup_child_logger

  logger = _setup_child_logger(__name__)
  ```

### Error Handling

- Raise specific exceptions (`ValueError`, `RuntimeError`, `FileNotFoundError`)
  with descriptive messages
- Document raised exceptions in the `Raises` section of NumPy docstrings

## Dependencies

- See `pyproject.toml` for the full dependency list and version constraints
- Do not add new dependencies without justification
- Prefer the existing scientific stack

## Concurrency and Thread Safety

- **ESMF/xESMF is not thread-safe.** It uses global state and must not be
  called concurrently across threads.
- When using Dask distributed, workers must use `processes=True`,
  `threads_per_worker=1`, and `resources={"ESMF": 1}`.
- See `CoreParameter.dask_scheduler_type` for supported Dask scheduler
  configurations.

## Testing

- **Framework:** pytest with pytest-cov
- **Unit tests:** `tests/e3sm_diags/` (run by default)
- **Integration tests:** `tests/integration/` (run manually)
- Run unit tests: `pytest` (from repo root)
- Do not modify or remove existing tests unless directly related to your change

## Pre-commit Hooks

See `.pre-commit-config.yaml` for the full list of hooks.
Run `pre-commit run --all-files` to validate before committing.

## Documentation

- Sphinx-based; source files in `docs/`
- Build with `make -C docs html`
- Update documentation when changing user-facing behavior

## Pull Request Conventions

Follow the PR template in `.github/pull_request_template.md`:

- Reference the issue being closed
- Ensure code follows style guidelines
- Self-review changes before submitting
- Verify no new warnings are introduced
- Add tests for new functionality
- Update documentation if applicable

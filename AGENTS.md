# AI Development Instructions for E3SM Diagnostics

This is the canonical, tool-agnostic source of AI development rules for the
`e3sm_diags` repository. Tool-specific files (`.github/copilot-instructions.md`,
`.claude/CLAUDE.md`) derive their guidance from this document.

## Project Overview

E3SM Diagnostics (`e3sm_diags`) is a Python-based scientific diagnostics
package for the DOE Energy Exascale Earth System Model (E3SM). It evaluates
climate model performance by comparing model output against observations and
reanalysis products.

- **License:** BSD 3-Clause
- **Python:** ≥3.11 (supports 3.11–3.14)
- **Backend:** Xarray and xCDAT (replaced CDAT in v3.0.0)
- **Entry points:** `e3sm_diags` (main driver), `e3sm_diags_vars`

## Architecture

The package follows a **driver–parameter** architecture:

| Directory             | Purpose                                        |
| --------------------- | ---------------------------------------------- |
| `e3sm_diags/driver/`  | One driver per diagnostic set (~25 drivers)     |
| `e3sm_diags/parameter/` | Paired parameter class per driver             |
| `e3sm_diags/derivations/` | Derived variable calculations               |
| `e3sm_diags/metrics/`    | Statistical metric functions                 |
| `e3sm_diags/plot/`       | Matplotlib/Cartopy plotting modules          |
| `e3sm_diags/viewer/`    | HTML results viewer                           |
| `e3sm_diags/parser/`    | CLI and configuration parsing                 |
| `tests/e3sm_diags/`     | Unit tests (pytest)                           |
| `tests/integration/`    | Integration tests                             |
| `docs/`                 | Sphinx documentation                          |

When adding a new diagnostic set, create both a driver module in `driver/` and
a corresponding parameter class in `parameter/`.

## Coding Standards

### Style and Formatting

- **Formatter/Linter:** Ruff (configured in `pyproject.toml`)
- **Import sorting:** Ruff isort (`I` rules)
- **Line length:** No hard limit enforced (`E501` ignored), but keep lines
  reasonable
- **Complexity:** McCabe max complexity is 18 (`C901`)
- **Docstrings:** NumPy convention (`ruff.lint.pydocstyle.convention = "numpy"`)
- Use `from __future__ import annotations` at the top of every module

### Type Annotations

- **Type checker:** mypy (configured in `pyproject.toml`, target Python 3.13)
- Add type hints to all function signatures
- Use `TYPE_CHECKING` guard for imports used only in annotations:
  ```python
  from __future__ import annotations

  from typing import TYPE_CHECKING

  if TYPE_CHECKING:
      import xarray as xr
  ```
- `check_untyped_defs = true` — untyped functions are still checked
- `ignore_missing_imports = true` — third-party stubs are not required

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

### Core Scientific Stack

xarray, xcdat, uxarray, cartopy, matplotlib, dask, numpy, scipy, xesmf,
xskillscore, netcdf4, cf-units

### Key Constraints

- `cartopy >=0.17.0,<0.25.0` (plot regression constraint)
- `numpy >=2.0.0,<3.0.0`
- `xcdat >=0.11.1,<1.0.0`
- `dask !=2024.12.0,!=2024.12.1`

Do not add new dependencies without justification. Prefer the existing stack.

## Concurrency and Thread Safety

- **ESMF/xESMF is not thread-safe.** It uses global state and must not be
  called concurrently across threads.
- When using Dask distributed, workers must use `processes=True`,
  `threads_per_worker=1`, and `resources={"ESMF": 1}`.
- The project supports two Dask scheduler types via
  `CoreParameter.dask_scheduler_type`: `"processes"` (default, uses
  `dask.bag`) and `"distributed"` (experimental, uses
  `dask.distributed.Client`).

## Testing

- **Framework:** pytest with pytest-cov
- **Unit tests:** `tests/e3sm_diags/` (run by default)
- **Integration tests:** `tests/integration/` (run manually via `test.sh`)
- Run unit tests: `pytest` (from repo root)
- Coverage reports are generated as HTML and XML under
  `tests_coverage_reports/`
- Do not modify or remove existing tests unless directly related to your change

## Pre-commit Hooks

The repository uses pre-commit with the following hooks:

1. Trailing whitespace removal
2. End-of-file fixer
3. YAML validation
4. Ruff — import sorting, linting, formatting
5. mypy — type checking

Run `pre-commit run --all-files` to validate before committing.

## Documentation

- Sphinx with `sphinx_rtd_theme`
- Source files in `docs/source/`
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

# Claude Instructions for E3SM Diagnostics

> Full AI development rules: see [`AGENTS.md`](../AGENTS.md) at the repository
> root. This file provides concise, enforceable guidance for Claude-based tools.

## Project Context

E3SM Diagnostics is a Python scientific diagnostics package for E3SM
evaluation, built on Xarray/xCDAT. It follows a **driver–parameter**
architecture: each diagnostic set has a paired driver in `e3sm_diags/driver/`
and parameter class in `e3sm_diags/parameter/`.

## Code Generation Rules

- Use `from __future__ import annotations` at the top of every module.
- Add type hints to all function signatures.
- Guard type-only imports with `TYPE_CHECKING`:
  ```python
  from typing import TYPE_CHECKING

  if TYPE_CHECKING:
      import xarray as xr
  ```
- Write NumPy-style docstrings for all public functions and classes. Include
  `Parameters`, `Returns`, and `Raises` sections.
- Use the project logger, not `print()`:
  ```python
  from e3sm_diags.logger import _setup_child_logger
  logger = _setup_child_logger(__name__)
  ```
- Raise specific exceptions (`ValueError`, `RuntimeError`) with descriptive
  messages.

## Style

- **Linter/Formatter:** Ruff (see `pyproject.toml` for configuration).
- **Type checker:** mypy (see `pyproject.toml` for configuration).
- `E501` (line-too-long) is ignored but keep lines reasonable.
- Import order is enforced by Ruff isort.

## Dependencies

Do not introduce new dependencies without justification. Prefer the existing
scientific stack. See `pyproject.toml` for the dependency list and version
constraints.

## Concurrency

ESMF/xESMF is **not thread-safe**. Dask workers must use `processes=True`,
`threads_per_worker=1`, and `resources={"ESMF": 1}` when using the distributed
scheduler.

## Testing

- Framework: pytest with pytest-cov.
- Unit tests live in `tests/e3sm_diags/`.
- Do not remove or modify unrelated tests.
- Run `pytest` from the repo root to execute unit tests.

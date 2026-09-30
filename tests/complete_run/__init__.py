"""Manual complete-run utilities for HPC validation workflows.

This package groups the manual-only modules used to run and compare
large complete-run diagnostics outside the normal unit-test and CI paths.
It exists so the workflow can live under `tests/` for developer access
without behaving like an import-time script or an automatically collected
pytest suite.

Usage
-----
Import the package modules through their explicit entrypoints, such as
``tests.complete_run.run`` and ``tests.complete_run.compare``, or invoke
them from the command line with ``python -m tests.complete_run.run`` and
``python -m tests.complete_run.compare``.
"""

from __future__ import annotations

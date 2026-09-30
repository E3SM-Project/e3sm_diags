# Plan: Rank complete-run image failures (#1083)

## Scope

Implement severity-based PNG comparison only for the manual complete-run
workflow. Do not change the targeted CI image-regression suite or the legacy
broad-integration image checker in this issue.

This follows [zppy #865](https://github.com/E3SM-Project/zppy/pull/865), which
addresses false positives from matplotlib text-metric and anti-aliasing shifts.
The complete-run workflow is the analogous high-volume, dependency-regression
comparison and already provides JSON reports, PNG artifacts, and an HTML
reviewer.

## Phase 1: Severity scoring and reports

1. Add `tests/complete_run/image_severity.py`, using existing NumPy, SciPy, and
   Pillow dependencies. Port the focused zppy scorer:
   - identify byte- and pixel-identical images;
   - trim uniform borders before comparison;
   - measure geometry change for layout/panel loss;
   - compare colors within a small neighborhood in both directions, forgiving
     small shifts while detecting added or deleted features;
   - detect compact high-contrast changes separately so changed statistics or
     labels remain reviewable;
   - classify each pair as `IDENTICAL`, `NEGLIGIBLE`, `MINOR`, `MODERATE`,
     `MAJOR`, or `STRUCTURAL` and return numeric scores plus a measured cause.

2. Integrate the scorer into `tests/complete_run/helpers.py`:
   - retain missing images as unconditional failures;
   - treat `IDENTICAL` and `NEGLIGIBLE` as passing and `MINOR` or worse as
     image mismatches;
   - add a structured image result for every shared pair, along with distinct
     identical and cosmetic collections in `ComparisonSummary`; retain
     `matching_images` as the combined passing collection for compatibility;
   - attach optional severity, score, geometry, and cause fields to the
     `ComparisonIssue` records used for reviewable failures;
   - create triptych artifacts only for reviewable images when artifacts are
     requested.

3. Update `tests/complete_run/compare.py` and `tests/complete_run/validate.py`:
   - report identical, cosmetic, and reviewable image counts in the JSON and
     command summary;
   - include severity metadata for each image mismatch;
   - remove the raw `--image-mismatch-threshold` option and its forwarding
     through the combined validate command, since raw pixel fraction no longer
     determines pass/fail status.

4. Update `docs/source/dev_guide/testing.rst` with the complete-run review
   behavior, including that `NEGLIGIBLE` means a cosmetic rendering difference
   rather than byte-for-byte identity.

## Phase 2: Severity-aware comparison viewer

After Phase 1 establishes the scored JSON report, update
`tests/complete_run/diff_html.py` to consume the severity metadata. The
existing comparison viewer should:

- sort triptychs by severity and then score, worst first;
- display each image's severity and measured cause;
- add severity filter controls and counts; and
- retain the existing diagnostic-set filter, search, and lazy loading.

The viewer lists only reviewable failures because cosmetic images have no diff
artifacts. Its summary must still show the cosmetic and identical counts from
the Phase 1 report.

This is a follow-up phase because it depends only on the Phase 1 report schema
and does not alter the scorer's pass/fail decisions.

## Tests and validation

1. Add Phase 1 synthetic-image unit tests covering identity, a tolerated small shift,
   added/deleted thin features, recoloring, structural size changes, compact
   changed-number-like differences, severity ordering, summary counts, JSON
   output, and removal of the obsolete CLI option.
2. Add Phase 2 viewer tests covering severity ordering, labels, filters, and
   the retained existing viewer controls.
3. Run the focused complete-run tests and applicable pre-commit hooks.
4. On NERSC, recompare an existing complete-run result with its accepted
   baseline. Review all images newly classified as cosmetic before relying on
   the auto-pass behavior; do not promote a baseline as part of this feature.

## Why targeted image regression is excluded

`tests/integration/test_plot_image_regressions.py` compares a small, targeted
set of committed plotting baselines in a fixed CI environment and already has
narrow, case-specific tolerances and dependency metadata. Ranking brings no
material triage benefit there and auto-passing cosmetic classifications would
weaken its strict regression signal.

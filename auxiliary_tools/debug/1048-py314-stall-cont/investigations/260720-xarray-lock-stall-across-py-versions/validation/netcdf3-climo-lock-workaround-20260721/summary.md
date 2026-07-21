# NetCDF3 Climatology Lock Workaround Validation

**Overall result: PASS**

This validation compares workaround-disabled and workaround-enabled full
diagnostic outputs within Python 3.13.14 and Python 3.14.6. Incomplete
Python 3.14.1–3.14.4 outputs are explicitly excluded.

## Criteria

- Numeric tolerance: `atol=1e-12`, `rtol=1e-12`.
- Missing/additional paths, structural or metadata differences, dtype
  changes, NaN-mask changes, nonnumeric value changes, and numeric values
  outside tolerance fail validation.
- JSON is canonicalized by sorting object keys, replacing each configured
  absolute results root with `<RESULTS_DIR>`, and replacing the generated
  `viewer/index.json` timestamp with `<VIEWER_GENERATED_AT>`.
- PNG dimensions and modes must match, and decoded RGBA pixels must be
  exactly equal. PNG container metadata and compression are not compared.
- Raw byte/JSON equality is reported before canonical and tolerance results.

## Results

| Python | Inventory | Files | Byte exact | NetCDF exact | JSON canonical exact | PNG pixel exact | Max abs | Max rel | Result |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 3.13.14 | 3300 common, 0 missing, 0 additional | 3300 | 3299/3300 | 5/5 | 1428/1428 | 1867/1867 | 0 | 0 | PASS |
| 3.14.6 | 3300 common, 0 missing, 0 additional | 3300 | 3299/3300 | 5/5 | 1428/1428 | 1867/1867 | 0 | 0 | PASS |

For each Python version, the sole byte-level mismatch is
`viewer/index.json`: its raw JSON has four expected provenance differences
(three embedded results-directory strings and one generation timestamp).
All four disappear under the documented canonicalization rules.

All five NetCDF files in each comparison were checked for data model,
groups, dimensions, variable names/dimensions/shapes, metadata, dtypes,
NaN masks, exact values, and tolerance. All JSON files were parsed and
compared recursively after canonicalization. All PNG files were decoded
to RGBA and checked for matching dimensions, modes, and exact pixels.

## Artifacts and reproduction

- [Inventory differences](inventory.tsv)
- [Concise per-file results](file-results.tsv)
- [Complete machine-readable details](comparison-details.json)
- [Standalone validator](validate_outputs.py)

From this directory, rerun with:

```bash
conda run -n ed_1048_xr_2026070_py31314 python validate_outputs.py
```

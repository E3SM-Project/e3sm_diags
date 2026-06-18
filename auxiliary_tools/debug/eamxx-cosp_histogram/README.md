# EAMxx COSP Histogram Testing

This directory contains tools for testing EAMxx COSP histogram support.

## Background

EAMxx output uses different variable names, dimension names, and dimension
orders than EAM for the COSP histogram variables:

- `FISCCP1_COSP(time, cosp_prs, cosp_tau, lat, lon)` -> `isccp_ctptau(time, cosp_tau, cosp_prs, lat, lon)`
- `CLMODIS(time, cosp_prs, cosp_tau_modis, lat, lon)` -> `modis_ctptau(time, cosp_tau, cosp_prs, lat, lon)`
- `CLD_MISR(time, cosp_htmisr, cosp_tau, lat, lon)` -> `misr_cthtau(time, cosp_tau, cosp_cth, lat, lon)`

Note that EAMxx orders the histogram axes as `(tau, prs/cth)`, the reverse of
EAM. `cosp_histogram_standardize()` transposes them back to `(prs/cth, tau)`.

> **Note:** Earlier EAMxx output had these dimensions but **no coordinate
> values**, which required a separate preprocessing step to inject defaults.
> The current EAMxx output already ships coordinate values for `cosp_tau`,
> `cosp_prs`, and `cosp_cth`, so **no preprocessing is needed**. The old
> `preprocess_eamxx_data.py` script has been removed.

## Workflow

### Step 1: Inspect the data (optional)

```bash
python inspect_data.py
```

Confirms the COSP variables and that `cosp_tau`, `cosp_prs`, and `cosp_cth`
carry coordinate values.

### Step 2: Run diagnostics

```bash
python test_eamxx_cosp_histogram.py
```

Runs the `cosp_histogram` set and the `lat_lon` set (COSP cloud-fraction
maps) on the EAMxx output, driven by `cosp_model_vs_obs.cfg`. CALIPSO
variables are excluded because EAMxx output does not provide them.

## Coordinate Values

From the EAMxx specification (bin centers):

```python
cosp_prs = [90000, 74000, 62000, 50000, 37500, 24500, 9000]  # Pa
cosp_tau = [0.15, 0.8, 2.45, 6.5, 16.2, 41.5, 100]  # unitless
cosp_cth = [0, 250, 750, 1250, 1750, 2250, 2750, 3500, 4500, 6000,
            8000, 10000, 12000, 14500, 16000, 18000]  # meters
```

## Files

- `test_eamxx_cosp_histogram.py` - Run full diagnostic suite
- `cosp_model_vs_obs.cfg` - lat_lon + cosp_histogram diagnostics (CALIPSO excluded)
- `inspect_data.py` - Inspect data structure
- `README.md` - This file

## Code Changes

The following changes were made to support EAMxx COSP variables:

1. **e3sm_diags/derivations/formulas_cosp.py**
   - Added `"cosp_cth"` to `CLOUD_HIST_MAP["prs"]["keys"]`
   - Transpose EAMxx `(tau, prs/cth)` axes back to `(prs/cth, tau)`
   - Added `"misr_cthtau"` to MISR variable list
   - Added `"isccp_ctptau"` and `"modis_ctptau"` to ISCCP/MODIS variable list
   - Added `"cosp_cth": 1000` to `PRS_UNIT_ADJ_MAP` (MISR height m -> km) so
     `cosp_bin_sum` cloud-level subsetting works for `lat_lon` MISR variables

2. **e3sm_diags/derivations/derivations.py**
   - Mapped new variables to `cosp_histogram_standardize` (histograms) and
     `cosp_bin_sum` (lat_lon cloud-fraction `CLD*_TAU*` variables)

3. **tests/e3sm_diags/derivations/test_formulas_cosp.py**
   - Added test coverage for new variables (standardize and bin sum)

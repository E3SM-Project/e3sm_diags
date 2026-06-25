# New ENSO diagnostics (ported from a-prime)

New `plot_type`s were added to the existing `enso_diags` set, porting
diagnostics from [a-prime](../../../../a-prime) that e3sm_diags previously
lacked.

| `plot_type` | a-prime source | Driver | Plot |
|---|---|---|---|
| `nino_index_timeseries` | `plot_multiple_index` | `run_diag_nino_index_timeseries` | `plot_nino_index_timeseries` |
| `seasonality` | `plot_multiple_index_seasonality` | `run_diag_seasonality` | `plot_seasonality` |
| `interannual_variability` | `plot_stddev` | `run_diag_interannual_variability` | `plot_map` (reused) |
| `equatorial_soi` | `plot_multiple_index_same_plot` (`SOI_Nino`) | `run_diag_equatorial_soi` | `plot_equatorial_soi` |
| `lead_lag` | `plot_regress_lead_lag_index_field` (`ENSO Evolution`) | `run_diag_lead_lag` | `plot_lead_lag` |

The two index plot types stack the three Niño regions (NINO3, NINO34, NINO4) as
subplot rows and reuse the existing, already-validated Niño index computation
(`calculate_nino_index_model` / `calculate_nino_index_obs` in
`e3sm_diags/driver/enso_diags_driver.py`) — only the plotting transforms are
new.

The Niño index is the area-averaged SST (or TS, if SST is unavailable) anomaly
over a Niño region, with the monthly climatology removed. For model–vs–obs
runs the reference index is the built-in HadISST record in
`e3sm_diags/driver/default_diags/enso_NINO{3,34,4}.long.data`.

## `nino_index_timeseries`

The monthly Niño index anomaly as a time series, with the test case in the left
column and the reference in the right column. Each panel shows the monthly
index, its mean, and a 13-month moving-average smoothing curve; the title
reports the mean and standard deviation. This visualizes the amplitude,
irregularity, and any drift of simulated ENSO variability.

## `seasonality`

The **per-calendar-month standard deviation** of the Niño index (12 values,
Jan–Dec), with test and reference overlaid in each region panel. This diagnoses
**ENSO phase-locking** — the tendency of ENSO events to mature in boreal winter
(Nov–Jan) and decay in boreal spring (Mar–May). Observations show a
characteristic U-shape: a spring minimum (the "spring predictability barrier")
and a late-year maximum. A model can reproduce the right ENSO amplitude but the
wrong seasonal timing; this plot checks whether the std-dev curve peaks in the
right months. A flat curve or a wrong-season peak signals a coupled-feedback or
mean-state bias.

## `interannual_variability`

The six panels of a-prime's `plot_stddev` over a region (e.g. the tropical band
`20S20N`), produced as **two** three-panel map figures that reuse the existing
`plot_map`:

- **Climatological mean** (`...-mean-...`): the day-weighted annual SST
  climatology (test / reference / difference). The mean is computed with the
  shared e3sm_diags `climo(..., "ANN")` routine.
- **Interannual standard deviation** (`...-std-...`): the standard deviation of
  the day-weighted annual means across years (test / reference / difference).
  The annual means come from xcdat's `temporal.group_average(freq="year")`,
  which `climo` cannot provide because it collapses the time axis.

The standard-deviation map is the headline ENSO diagnostic: it shows the
spatial footprint of interannual SST variability, which in observations peaks as
a tongue along the eastern equatorial Pacific cold-tongue / Niño regions. A
model with too-weak or mislocated variability, or a double-ITCZ-style bias,
shows up directly. SST is derived as `TS` (in degC) masked to ocean points; for
model–vs–obs the reference is the gridded HadISST SST.

## `equatorial_soi`

The **Equatorial Southern Oscillation Index (EQSOI)** overlaid with the Niño3.4
SST anomaly index, as a two-row × two-column (test / reference) time series. The
correlation between the two indices is annotated in each panel.

The EQSOI is the difference of the **standardized** monthly sea level pressure
(`PSL`) anomalies between the eastern equatorial Pacific (`EPAC`,
−5–5°N, 230–280°E) and the equatorial Indian Ocean / Indonesia (`INDO`,
−5–5°N, 90–140°E): for each region the area-average monthly anomaly is computed,
then standardized (subtract mean, divide by std); EQSOI = z(EPAC) − z(INDO)
(unitless). This is a pressure-based ENSO index that is independent of SST, so
its strong anticorrelation with Niño3.4 (≈ −0.8) is a check on the simulated
**Walker-circulation / Bjerknes feedback**: a positive EQSOI (high pressure in
the east, low in the west) accompanies La Niña, and vice versa. A weak or
positive correlation signals a broken atmosphere–ocean coupling.

`EPAC`/`INDO` were added to `REGION_SPECS`. The Niño3.4 index reuses
`calculate_nino_index_*` and is plotted as the raw anomaly (not standardized),
matching a-prime's plot-time `--stdize 0`. For model–vs–obs the EQSOI reference
is observed `PSL` (ERA5 time series, resolved via the existing `psl`
derivation) and the Niño3.4 reference is the built-in HadISST record; both must
cover the same reference years for the correlation.

## `lead_lag`

The **lead-lag regression and correlation** of a field anomaly on the Niño3.4
index, as tiled global maps with the lags (`-8, -4, 0, 4, 8` months) as rows and
the test case, reference, and difference as columns. Two figures are produced
per variable: the regression coefficients (with significance hatching on the
test and reference columns) and the correlations. Positive lags indicate the
Niño index leading the field. This is a port of a-prime's
`plot_regress_lead_lag_index_field` ("ENSO Evolution") diagnostic.

At each lag the field anomaly and the Niño index are sliced relative to each
other (`lag >= 0`: field at `t+lag` on index at `t`; `lag < 0`: the reverse),
then `xs.linslope` / `xs.pearson_r` give the regression and correlation, with
`pearson_r_p_value < 0.05` marking significance. The diagnostic shows the ENSO
life cycle — how anomalies develop before, peak at, and decay after the SST
maximum — which is a check on the simulated coupled feedbacks.

**Variables** mirror the `regression_map` set (`TS, PRECT, TAUX, TAUY, LHFLX,
SHFLX, NET_FLUX_SRF`, plus `CLD*` for model–vs–model). Because `lead_lag` at
lag 0 is the same computation as the `regression_map` set, the contour levels are
taken from the cfg (`contour_levels`, same as the corresponding `regression_map` block);
correlation levels are fixed to ±1. The field is **not** land-masked (e.g. `TS`
includes land/sea ice), and the maps are global, matching a-prime. For
model–vs–obs the reference fields are ERA5 (and GPCP_v3.2 for `PRECT`) and the
Niño3.4 reference is the built-in HadISST record.

> **Note:** `lead_lag` is by far the slowest set — a global regression at every
> grid point for five lags and both cases, times two figures, per variable.

## Configuration

Each plot type is a single config block in the default diags
(`e3sm_diags/driver/default_diags/enso_diags_model_vs_obs.cfg` and
`..._model_vs_model.cfg`):

```
[#]
sets = ["enso_diags"]
plot_type = "nino_index_timeseries"   # or "seasonality"
case_id = "NINO-index"           # or "NINO-seasonality"
variables = ["TS"]
ref_name = "HadISST"
reference_name = "HadISST"
seasons = ["ANN"]

[#]
sets = ["enso_diags"]
plot_type = "interannual_variability"
case_id = "SST-variability"
variables = ["SST"]
regions = ["20S20N"]
ref_name = "HadISST"
reference_name = "HadISST"
seasons = ["ANN"]
test_colormap = "viridis"
reference_colormap = "viridis"
diff_colormap = "bwr"

[#]
sets = ["enso_diags"]
plot_type = "equatorial_soi"
case_id = "EQSOI"
variables = ["PSL"]
ref_name = "ERA5"                 # obs PSL; Nino3.4 ref is built-in HadISST
reference_name = "ERA5 (PSL), HadISST (Nino3.4)"
seasons = ["ANN"]

[#]
sets = ["enso_diags"]
plot_type = "lead_lag"
case_id = "TS-lead-lag"           # one block per variable, like the map set
variables = ["TS"]
regions = ["global"]
ref_name = "ERA5"                 # GPCP_v3.2 for PRECT; Nino3.4 ref is HadISST
reference_name = "ERA5 (TS), HadISST (Nino3.4)"
seasons = ["ANN"]
test_colormap = "diverging_bwr.rgb"
contour_levels = [-1, -0.8, -0.6, -0.4, -0.2, 0.2, 0.4, 0.6, 0.8, 1]
figsize = [14, 13]
```

For the index plot types the stacked regions default to
`["NINO3", "NINO34", "NINO4"]` via the `nino_regions` attribute on
`EnsoDiagsParameter`; `interannual_variability` instead uses the standard
`regions` attribute for the map domain.

## Running the local test

`run_enso_diags.py` sets the data paths/years:

```bash
conda activate e3sm_diags_dev_py313
python run_enso_diags.py
```

It runs model–vs–obs on the v3.LR.amip_0101 simulation (1995–2004), with
observations from the standard e3sm_diags obs directory
(`/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/time-series`). The index
plots use the built-in HadISST Niño index as the reference; the
`interannual_variability` maps use the gridded HadISST SST; `equatorial_soi`
uses ERA5 PSL; and `lead_lag` uses ERA5 (and GPCP_v3.2 for `PRECT`) fields. The
viewer is written to the NERSC portal directory configured in `results_dir`.

> **Note:** `lead_lag` dominates the runtime (several minutes per variable). To
> iterate on the other sets quickly, run a `-d` cfg that omits the `lead_lag`
> blocks.

## Verifying faithfulness to a-prime

`verify_vs_aprime.py` reimplements a-prime's algorithms in plain numpy on the
same input files and diffs them against the netCDF each diagnostic saves.

- **Index / seasonality** (region box → cos-lat area average →
  monthly-climatology anomaly): both match to machine precision (~1e-13) for the
  open-ocean regions NINO34 and NINO4. NINO3 differs by <0.1% because its eastern
  edge abuts South America, where xcdat's bounds-aware area weighting (used by
  the diagnostic) handles partial/land cells differently than the simple
  point-wise cos-lat weighting in the verification script — a weighting nuance,
  not an algorithmic difference.
- **Equatorial SOI** (per-region PSL anomaly → standardize → EPAC − INDO):
  matches the numpy reimplementation to ~5e-6, the same bounds-aware vs. cos-lat
  area-weighting nuance as the index above (smaller here since the PSL regions
  are open ocean with no land mask).
- **Interannual variability** (day-weighted annual means → mean & std across
  years): the standard-deviation map matches to machine precision (~1e-15). The
  climatological-mean map agrees to ~1e-6 degC because the diagnostic's `climo`
  accumulates its day-weighted average in the field's native float32 over all
  months, whereas the verification computes a float64 mean of annual means
  (mathematically identical for a noleap calendar) — again a floating-point /
  weighting nuance, not an algorithmic difference.
- **Lead-lag**: not reimplemented in numpy here. The diagnostic reuses the same
  anomaly + `xs.linslope` / `pearson_r` machinery already validated for the
  `regression_map` set, and `lead_lag` at **lag 0 is identical to the `regression_map` set**
  (same field, same Niño index). That equivalence is the intended cross-check
  (and is also the basis of the planned ERA5-vs-ERA-Interim reference comparison).

```bash
python verify_vs_aprime.py
```

# New ENSO diagnostics (ported from a-prime)

New `plot_type`s were added to the existing `enso_diags` set, porting
diagnostics from [a-prime](../../../../a-prime) that e3sm_diags previously
lacked.

| `plot_type` | a-prime source | Driver | Plot |
|---|---|---|---|
| `index_timeseries` | `plot_multiple_index` | `run_diag_index_timeseries` | `plot_index_timeseries` |
| `seasonality` | `plot_multiple_index_seasonality` | `run_diag_seasonality` | `plot_seasonality` |
| `interannual_variability` | `plot_stddev` | `run_diag_interannual_variability` | `plot_map` (reused) |

The two index plot types stack the three Niño regions (NINO3, NINO34, NINO4) as
subplot rows and reuse the existing, already-validated Niño index computation
(`calculate_nino_index_model` / `calculate_nino_index_obs` in
`e3sm_diags/driver/enso_diags_driver.py`) — only the plotting transforms are
new.

The Niño index is the area-averaged SST (or TS, if SST is unavailable) anomaly
over a Niño region, with the monthly climatology removed. For model–vs–obs
runs the reference index is the built-in HadISST record in
`e3sm_diags/driver/default_diags/enso_NINO{3,34,4}.long.data`.

## `index_timeseries`

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

## Configuration

Each plot type is a single config block in the default diags
(`e3sm_diags/driver/default_diags/enso_diags_model_vs_obs.cfg` and
`..._model_vs_model.cfg`):

```
[#]
sets = ["enso_diags"]
plot_type = "index_timeseries"   # or "seasonality"
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
```

For the index plot types the stacked regions default to
`["NINO3", "NINO34", "NINO4"]` via the `nino_regions` attribute on
`EnsoDiagsParameter`; `interannual_variability` instead uses the standard
`regions` attribute for the map domain.

## Running the local test

`run_enso_aprime_diags.py` sets the data paths/years; `run_enso_aprime_diags.cfg`
limits the run to the ported blocks. Run with the dev conda env on `PYTHONPATH`
so the working tree (not the installed package) is used:

```bash
conda activate e3sm_diags_dev_py313
PYTHONPATH=/global/u2/c/chengzhu/e3sm_diags:$PYTHONPATH \
  python run_enso_aprime_diags.py -d run_enso_aprime_diags.cfg
```

It runs model–vs–obs on a 10-year E3SM v2 piControl time series
(years 0051–0060). The index plots use the built-in HadISST Niño index as the
reference; the `interannual_variability` maps use the gridded HadISST SST. The
viewer is written to the NERSC portal directory configured in `results_dir`.

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
- **Interannual variability** (day-weighted annual means → mean & std across
  years): the standard-deviation map matches to machine precision (~1e-15). The
  climatological-mean map agrees to ~1e-6 degC because the diagnostic's `climo`
  accumulates its day-weighted average in the field's native float32 over all
  months, whereas the verification computes a float64 mean of annual means
  (mathematically identical for a noleap calendar) — again a floating-point /
  weighting nuance, not an algorithmic difference.

```bash
PYTHONPATH=/global/u2/c/chengzhu/e3sm_diags:$PYTHONPATH \
  python verify_vs_aprime.py
```

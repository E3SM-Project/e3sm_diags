# ERA5 reference data for e3sm_diags

A single workflow that retrieves ERA5 monthly means from the ECMWF Climate Data
Store (CDS), converts them to the units and conventions e3sm_diags expects, and
builds the time-series and climatology files served from
`observations/Atm/{time-series,climatology}/ERA5`.

Rerun it with a later `--end-year` whenever new ERA5 years become available.

## Why this replaces the original scripts

The 1979-2019 dataset was assembled from two different sources:

| Source | Script | Variables |
| ------ | ------ | --------- |
| obs4MIPs / CREATE-IP (CMORized ERA5) | `../create_ERA5_climo.sh` | the CMIP-named fields (`pr`, `ta`, `rlut`, …) |
| Direct CDS downloads | `../create_ERA5_ext_climo.sh` | fields obs4MIPs does not carry (`si10`, `d2m`, `sp`, `t2m`, `tp`, `cp`, `lsp`, `e`, `ro`, `vimd`) |
| Direct CDS download, regridded to 1 degree | `../create_ERA5_U_climo.sh` | `ua` for the QBO diagnostic |

Those scripts hard-coded CDS request IDs (`adaptor.mars.internal-…nc`), so
extending the record meant redoing the downloads by hand. They are kept for
provenance; new work should use this directory.

Everything now comes from the CDS and is described by one table,
[`era5_variables.yml`](era5_variables.yml), which records for each variable the
CDS fields it is built from, the conversion formula, its output metadata, and —
in the `legacy` key — how the original dataset produced the same variable.

## Prerequisites

1. `cdsapi`, which is in `conda-env/dev.yml`:

   ```bash
   conda env create -f conda-env/dev.yml    # or: conda install -c conda-forge cdsapi
   conda activate e3sm_diags_dev
   ```

2. A CDS account and a personal access token. Register at
   <https://cds.climate.copernicus.eu>, accept the ERA5 licence terms on both
   the [single-level](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels-monthly-means)
   and [pressure-level](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-pressure-levels-monthly-means)
   monthly-means dataset pages, then copy your token from
   <https://cds.climate.copernicus.eu/profile> into `~/.cdsapirc`:

   ```
   url: https://cds.climate.copernicus.eu/api
   key: <your-personal-access-token>
   ```

   The legacy `…/api/v2` endpoint and UID:key credentials were retired in 2024
   and no longer work.

## Workflow

```bash
cd analysis_data_preprocess/ERA5

# 1. Retrieve raw monthly means (resumable; existing files are skipped).
python era5_pipeline.py download --start-year 1979 --end-year 2025

# 2. Convert to per-variable time series in e3sm_diags conventions.
python era5_pipeline.py process --start-year 1979 --end-year 2025

# 3. Build ANN/DJF/MAM/JJA/SON and monthly climatologies.
python era5_pipeline.py climo --start-year 1979 --end-year 2025
```

Output under `--base-dir` (default
`/global/cfs/cdirs/e3sm/zhang40/analysis_data_e3sm_diags/ERA5_v2`):

```
raw/            era5_<short_name>_<years>.nc     as downloaded from the CDS
time_series/    <variable>_197901_202512.nc      one file per variable
climatology/    ERA5_<season>_<yyyymm>_<yyyymm>_climo.nc   all variables merged
```

Useful flags: `-v pr ta` to work on a subset of variables, `--dry-run` to list
the retrievals without submitting them, `--complevel 0` to write uncompressed
files, `--overwrite` to reprocess.

## Validating a change

Before regenerating the full record, reproduce a year that the original dataset
already covers and compare:

```bash
python era5_pipeline.py download --start-year 2019 --end-year 2019
python era5_pipeline.py process  --start-year 2019 --end-year 2019
python era5_pipeline.py compare  --start-year 2019 --end-year 2019
```

`compare` lines up the two datasets on the months they share and reports
`max|diff|`, RMSE and both global means per variable.

Expect exact agreement for the analysis fields (`ta`, `ua`, `va`, `zg`, `hus`,
`hur`, `wap`, `ps`, `psl`, `tas`, `ts`, `uas`, `vas`, `tauu`, `tauv`) and small
differences for the flux and precipitation fields, which the original dataset
derived from daily accumulations (`strd/86400`, `cp/86.4`, …) while this
workflow uses ERA5's mean-rate fields (`msdwlwrf`, `mcpr`, …). A large
difference, a sign flip or a round-number ratio means the formula in
`era5_variables.yml` needs attention.

## Notes on the ERA5 data

- **ERA5T.** The most recent ~3 months are preliminary (ERA5T) and arrive with
  an extra `expver` dimension. `process` prefers `expver=1` (ERA5) and falls
  back to `expver=5` (ERA5T) only where ERA5 is not yet available. Reprocess the
  affected years once they are finalized.
- **ERA5.1.** The CDS serves ERA5.1 for 2000-2006, which corrects a cold bias in
  the lower stratosphere. The original 1979-2019 files predate that, so
  stratospheric temperatures over those years will not match exactly.
- **Accumulated fields.** In the monthly-means product, accumulations (`tp`,
  `cp`, `lsp`, `e`, `ro`) are reported as a mean *daily* total, so `tp` is
  m day-1 despite a units attribute of `m`. They are carried through unchanged
  for backward compatibility; use `pr`, `prc` and `evspsbl` for anything
  quantitative.
- **Time bounds.** Time stamps are placed at mid-month with bounds spanning the
  whole month, which is what xCDAT needs to weight climatologies correctly.
- **Grid.** Native 0.25 degree (721x1440), latitude ordered south-to-north, all
  37 ERA5 pressure levels ordered surface-to-top — the same conventions as the
  original CMORized files, so no diagnostic needs to change.

## Adding a variable

Add an entry to `era5_variables.yml` naming the CDS fields it needs, the formula
over their short names, and its output metadata, then run `download`, `process`
and `climo` with `-v <variable>`. Nothing in the pipeline script needs to change.

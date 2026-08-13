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

Output under `--base-dir` (default `$SCRATCH/analysis_data_e3sm_diags/ERA5_v2`):

```
raw/            era5_<short_name>_<years>.nc     as downloaded from the CDS
time_series/    <variable>_197901_202512.nc      one file per variable
climatology/    ERA5_<season>_<yyyymm>_<yyyymm>_climo.nc   all variables merged
```

Useful flags: `-v pr ta` to work on a subset of variables, `--dry-run` to list
the retrievals without submitting them, `--complevel 0` to write uncompressed
files, `--overwrite` to reprocess.

## Long runs

The full record is ~1700 retrievals and takes many hours, so do not run it in a
foreground shell. Either detach it from the login session:

```bash
setsid nohup python -u era5_pipeline.py download --start-year 1979 --end-year 2025 \
    > $SCRATCH/era5_download.log 2>&1 < /dev/null &
```

or submit it as a batch job. Compute nodes reach the CDS only through the NERSC
proxy:

```bash
#SBATCH --qos=shared --time=24:00:00 --constraint=cpu
export https_proxy=http://proxy.nersc.gov:3128
export http_proxy=http://proxy.nersc.gov:3128
python -u era5_pipeline.py download --start-year 1979 --end-year 2025
```

Either way the step is resumable: finished files are skipped and each retrieval
is written to a `.part` file that is renamed only once the download completes,
so an interrupted run never leaves a truncated file behind.

## Disk space

On the native 0.25 degree grid (721x1440) a 2D month is 4.2 MB and a 3D month
(37 levels) is 154 MB, so the full 1979-2025 record is large:

| | 2019 only | 1979-2025 (564 months) |
| --- | --- | --- |
| `raw/` (35 2D + 8 3D fields, deflated by the CDS) | ~10 GB | ~450 GB |
| `time_series/` (45 variables, `--complevel 1`) | ~10 GB | ~500 GB |
| `climatology/` (17 files, all variables merged) | 24 GB | 24 GB |

That is why the default `--base-dir` is on scratch. Delete `raw/` once the
variables built from it have been processed, then copy `time_series/` and
`climatology/` to their permanent home:

```bash
cp -r $SCRATCH/analysis_data_e3sm_diags/ERA5_v2/{time_series,climatology} \
      /global/cfs/cdirs/e3sm/diagnostics/observations/Atm/...
```

Scratch is purged periodically, so do not leave the only copy of a finished run
there.

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

### Results for 2019

43 of the 45 variables agree with the 1979-2019 dataset to float32 roundoff,
including the flux and precipitation fields (≤0.007% relative difference in the
global mean). The two that do not are described below; in both cases the new
output is not the file at fault, so neither blocks a full run.

#### `sp` is displaced in the original file after 2013-11

`observations/Atm/time-series/ERA5/sp_197901_201912.nc` is bit-identical to the
CMORized `ps_197901_201912.nc` from 1979-01 through 2013-11. Its last 73 months,
2013-12 through 2019-12, hold the field rolled 16 grid cells (4.00 degrees) east:
an FFT cross-correlation against `ps` gives `dlat=0`, `dlon=16`, with the
correlation rising from 0.921 to 0.9999 and the RMSE dropping from 3792 Pa to
147 Pa once the roll is undone.

A longitudinal roll preserves area averages, so global means still agree to
0.006 Pa and the error is invisible to any global-mean check. Locally it is
large: 99175.8 Pa over the Tibetan plateau (lat 31, lon 79.5) where the true
surface pressure is 51437 Pa.

This propagates into the `ERA5_ext` climatologies e3sm_diags ships, where `sp`
differs from `ps` by up to 6939 Pa — the 73/492 dilution of the monthly error —
with 27.7% of cells off by more than 100 Pa. `sp` has one consumer,
`("d2m", "sp"): qsat` in `e3sm_diags/derivations/derivations.py`, which derives
`QREFHT` because ERA5 ships no `huss`; that runs as a shipped default in
`lat_lon_model_vs_obs.cfg` against `ref_name = "ERA5_ext"`. The resulting
`QREFHT` error reaches 1.11 g/kg (11.6%) over high terrain, against difference
contours that start at 0.25 g/kg, so it plots as an apparent model bias. The
global mean is nearly unaffected (7.0642 vs 7.0627 g/kg), which is presumably
why it went unnoticed.

The `sp` this workflow downloads is correct, so regenerating the dataset fixes
`QREFHT`. `compare` will keep reporting a large `sp` difference until the
original file is replaced — that difference is the bug, not a regression.

#### `vimd` disagrees beyond the unit change

The original file is in kg m-2 day-1, matching its `units = "kg m**-2"`
attribute, while this workflow writes kg m-2 s-1 — a factor of 86400. The two
disagree by more than that factor, and the moisture budget favours the
*original*: since `dW/dt + div Q = E - P`, monthly `vimd` should track
`evspsbl - pr`, and against this workflow's own `evspsbl - pr` for 2019 the
original correlates 0.982 with a slope of 85305 (within 1.3% of 86400) while
the new download correlates only 0.853 with a slope of 0.961. Both global means
are ~0 once the units are matched, as a divergence integral requires.

So the CDS field named in `era5_variables.yml` may not be the one the ext script
downloaded — worth resolving before anyone relies on `vimd`. Nothing in
e3sm_diags reads it today.

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

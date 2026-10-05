# ERA5 reference data for e3sm_diags

A single workflow that retrieves ERA5 monthly means from the ECMWF Climate Data
Store (CDS), converts them to the units and conventions e3sm_diags expects, and
builds the time-series and climatology files served under
`observations/Atm/`.

It replaces the scripts that produced the original 1979-2019 dataset
(`../create_ERA5_climo.sh`, `../create_ERA5_ext_climo.sh`,
`../create_ERA5_U_climo.sh`), which mixed obs4MIPs files with hand-run CDS
downloads. Those are kept for provenance. Every variable is now described by
one table, [`era5_variables.yml`](era5_variables.yml): the CDS fields it is
built from, the conversion formula, its output metadata and, under `legacy`,
how the original dataset produced it.

## Prerequisites

1. `cdsapi`, needed by `download` only: `pip install cdsapi` (or
   `conda install -c conda-forge cdsapi`) in the e3sm_diags environment.

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
raw/            era5_<short_name>_<year>-<year>.nc   as downloaded, one file per field-year
time_series/    <variable>_197901_202512.nc          one file per variable
climatology/    ERA5_<season>_<yyyymm>_<yyyymm>_climo.nc   all variables merged
```

Useful flags: `-v pr ta` to work on a subset of variables, `--dry-run` to list
the retrievals without submitting them, `--complevel 0` to write uncompressed
files, `--overwrite` to reprocess.

## Running on Perlmutter

Run the full record as batch jobs, submitted from a log directory:

```bash
conda activate e3sm_diags_dev_py313
export ERA5_DIR=$HOME/e3sm_diags/analysis_data_preprocess/ERA5
mkdir -p $SCRATCH/analysis_data_e3sm_diags/ERA5_v2/slurm
cd $SCRATCH/analysis_data_e3sm_diags/ERA5_v2/slurm

# download: one year per task, free `xfer` QOS (~1-2 days, mostly CDS queueing)
sbatch --array=1979-2025%3 $ERA5_DIR/submit_download.sh

# process: one variable per task, `shared` QOS
sbatch --array=0-43%12 $ERA5_DIR/submit_process.sh

# climo: one job over all variables
sbatch --qos=shared --constraint=cpu --cpus-per-task=16 --mem=96G --time=06:00:00 \
    --job-name=era5_climo --output=%x_%j.log \
    --wrap "python -u $ERA5_DIR/era5_pipeline.py climo --start-year 1979 --end-year 2025"
```

What the 1979-2025 run took:

- **process:** about 3 minutes per surface field and 1-2 hours per pressure-level
  field (`ta ua va wap hus hur zg tro3`). The fields built from two sources
  (`rsus`, `rsut`, `rsutcs`, `clwvi`) ran out of memory at 32 GB and had to be
  resubmitted with more.
- **climo:** 58 minutes, peaking at about 96 GB.

Every stage is resumable: finished outputs are skipped, and downloads are
written to a `.part` file that is renamed only once complete. Resubmit the same
array to retry.

## Why `download` bundles variables

The CDS queues each request for an hour or more regardless of its size, so
`download` groups every field still missing from a year into one retrieval per
CDS dataset (surface and pressure levels), about 94 for the full record instead
of about 2000. It then splits the result into one file per field.

Two things in the CDS response are handled in `split_bundle`:

- A bundle spanning both ECMWF streams arrives as a **zip**, one netCDF per
  stream. Analysis and forecast fields are stamped at different hours, so the
  members must never be merged.
- The CDS **renames** the mean-rate fields (`msdwlwrf` arrives as
  `avg_sdlwrf`). The renames are irregular, so `cds_short_names` in
  `era5_variables.yml` lists them all.

## Notes on the data

- **Agreement with the 1979-2019 dataset.** Over 1979-2019 the new time series
  match the original to float32 roundoff, with the same grid, time axis, units
  and sign conventions. There are two exceptions; see PR #1079 for the analysis.
  - **`sp`**: the original file is rolled 16 grid cells east from 2013-12
    onward. The new `sp` is correct. This also fixes the `QREFHT` that
    e3sm_diags derives from `d2m` and `sp`.
  - **`vimd`** is not carried over. The CDS field delivers `vimdf`, which does
    not match the original, and nothing in e3sm_diags reads it.
- **Climatologies follow `ncclimo -a sdd`**, like the original files and the
  model climatologies they are compared with:
  - a monthly climatology weights every year equally;
  - a season weights its months by the non-leap calendar (February is always
    28 days).
- **ERA5T.** The most recent ~3 months are preliminary (ERA5T). `process` uses
  final ERA5 where it exists and falls back to ERA5T otherwise. Reprocess those
  months once they are finalized.
- **ERA5.1.** The CDS serves ERA5.1 for 2000-2006, which corrects a cold bias in
  the lower stratosphere. The original files predate it, so stratospheric
  temperatures over those years differ.
- **Accumulated fields** (`tp`, `cp`, `lsp`, `e`, `ro`) are mean daily totals,
  so `tp` is m day-1 despite a units attribute of `m`. They are kept for
  backward compatibility; use `pr`, `prc` and `evspsbl` for anything
  quantitative.
- **Grid.** Native 0.25 degree (721x1440). Latitude runs south to north and the
  37 pressure levels run from the surface up, the same conventions as the
  original files.

## Disk space

| | 1979-2025 (564 months) |
| --- | --- |
| `raw/` | 237 GB (1974 files) |
| `time_series/` | ~250 GB (44 files) |
| `climatology/` | 12.4 GB (17 files) |

That is why the default `--base-dir` is on scratch. Copy `time_series/` and
`climatology/` to their permanent home once validated; scratch is purged.

## Extending the record

Process each new year into its own file and add it next to the existing time
series. e3sm_diags opens every `<variable>_<yyyymm>_<yyyymm>.nc` file in the
directory together and joins them along time, so `rlds_197901_202512.nc` plus
`rlds_202601_202612.nc` read as one 1979-2026 record.

```bash
sbatch --array=2026 $ERA5_DIR/submit_download.sh   # new year only (~5 GB)
ERA5_START_YEAR=2026 ERA5_END_YEAR=2026 sbatch --array=0-43%12 $ERA5_DIR/submit_process.sh
```

Only the new year's raw files are needed; earlier years do not have to be on
disk. Then copy the 44 new `time_series/` files into the published time-series
directory.

- **Wait for final ERA5.** Process a year only once its December is final ERA5
  rather than ERA5T, roughly three months after the year ends.
- **Files must not overlap in time.** If a longer range is ever rebuilt into a
  single file, remove the yearly files it covers.
- **Climatology.** `climo` reads a single time-series file per variable, so it
  cannot build a climatology from yearly files. The climatology is a fixed
  reference period and does not need to change with every new year.

## Adding a variable

Add an entry to `era5_variables.yml` naming the CDS fields it needs, the formula
over their short names, and its output metadata, then run `download` and
`process` with `-v <variable>`. Nothing in the pipeline script needs to change.

Run `climo` **without** `-v`. Each climatology file holds every variable, and
`climo -v` rewrites it with only the selected ones.

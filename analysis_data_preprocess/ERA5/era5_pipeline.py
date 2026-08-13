#!/usr/bin/env python3
"""
Unified retrieval and processing workflow for the e3sm_diags ERA5 reference data.

Replaces the three scripts that produced the original 1979-2019 dataset
(`create_ERA5_climo.sh`, `create_ERA5_ext_climo.sh`, `create_ERA5_U_climo.sh`),
which pulled from a mix of obs4MIPs/CREATE-IP and ad-hoc CDS downloads. Every
variable now comes from the ECMWF Climate Data Store, driven by the table in
`era5_variables.yml`, so the dataset can be regenerated or extended by rerunning
the same command with a different year range.

Sub-commands
------------
download : retrieve raw monthly-mean fields from the CDS
process  : convert raw fields to e3sm_diags time series (one file per variable)
climo    : build ANN/DJF/MAM/JJA/SON + monthly climatologies from the time series
compare  : check new time series against the original files over shared years

Examples
--------
Validate the workflow against the existing dataset using a single year::

    python era5_pipeline.py download --start-year 2019 --end-year 2019
    python era5_pipeline.py process  --start-year 2019 --end-year 2019
    python era5_pipeline.py compare  --start-year 2019 --end-year 2019

Produce the full extended dataset::

    python era5_pipeline.py download --start-year 1979 --end-year 2025
    python era5_pipeline.py process  --start-year 1979 --end-year 2025
    python era5_pipeline.py climo    --start-year 1979 --end-year 2025

Requirements
------------
`cdsapi` (in `conda-env/dev.yml`) plus a CDS personal access token in
`~/.cdsapirc`. See the README in this directory.
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import xarray as xr
import yaml

SCRIPT_NAME = Path(__file__).name
CONFIG_PATH = Path(__file__).parent / "era5_variables.yml"

# Default working directory, holding `raw/`, `time_series/` and `climatology/`.
# The full record needs ~1 TB while the raw downloads are still around, so this
# defaults to scratch; copy `time_series/` and `climatology/` to their permanent
# home once a run is validated.
DEFAULT_BASE_DIR = (
    Path(
        os.environ.get(
            "SCRATCH", "/global/cfs/cdirs/e3sm/zhang40/analysis_data_e3sm_diags"
        )
    )
    / "analysis_data_e3sm_diags/ERA5_v2"
)

# The original 1979-2019 dataset, used by `compare`.
DEFAULT_REFERENCE_DIR = Path(
    "/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/time-series/ERA5"
)

# CDS dimension names mapped onto the ones the output uses. Applied both to the
# raw downloads and, in `compare`, to the original files written by the ext
# scripts, which kept the CDS names.
RENAME_DIMS = {
    "valid_time": "time",
    "latitude": "lat",
    "longitude": "lon",
    "pressure_level": "plev",
    "level": "plev",
    "isobaricInhPa": "plev",
}

# ERA5 is served on a 0.25 degree lat/lon grid; the output keeps that grid with
# latitude ordered south-to-north and pressure ordered surface-to-top, matching
# the original CMORized files.
MONTHS = [f"{m:02d}" for m in range(1, 13)]
SEASONS = {
    "DJF": [12, 1, 2],
    "MAM": [3, 4, 5],
    "JJA": [6, 7, 8],
    "SON": [9, 10, 11],
    "ANN": list(range(1, 13)),
}

logger = logging.getLogger("era5_pipeline")


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
def load_config(path: Path = CONFIG_PATH) -> dict[str, Any]:
    """Read the variable table.

    Parameters
    ----------
    path : Path
        Location of ``era5_variables.yml``.

    Returns
    -------
    dict
        Parsed configuration with ``datasets`` and ``variables`` keys.
    """
    with open(path) as f:
        config = yaml.safe_load(f)

    for key in ("datasets", "variables"):
        if key not in config:
            raise ValueError(f"{path} is missing the '{key}' section")

    return config


def select_variables(config: dict[str, Any], names: list[str] | None) -> dict[str, Any]:
    """Return the requested subset of the variable table.

    Raises
    ------
    ValueError
        If a requested variable is not defined in the table.
    """
    variables = config["variables"]

    if not names:
        return variables

    unknown = sorted(set(names) - set(variables))
    if unknown:
        raise ValueError(
            f"Unknown variable(s) {unknown}. Defined variables: {sorted(variables)}"
        )

    return {name: variables[name] for name in names}


def year_chunks(start_year: int, end_year: int, chunk_years: int) -> list[list[int]]:
    """Split a year range into retrieval chunks, aligned on decade boundaries."""
    years = list(range(start_year, end_year + 1))

    if chunk_years <= 1:
        return [[year] for year in years]

    chunks: list[list[int]] = []
    for year in years:
        # Start a new chunk at each decade boundary so reruns with a different
        # end year reuse the files already downloaded.
        if not chunks or year % chunk_years == 0 or len(chunks[-1]) >= chunk_years:
            chunks.append([])
        chunks[-1].append(year)

    return chunks


# -----------------------------------------------------------------------------
# download
# -----------------------------------------------------------------------------
def collect_sources(
    variables: dict[str, Any], config: dict[str, Any]
) -> list[tuple[str, str, str]]:
    """List the unique (dataset, CDS variable, short name) triples to retrieve.

    Several output variables share a source field (``sp`` feeds both ``ps`` and
    ``sp``, ``msdwlwrf`` feeds both ``rlds`` and ``rlus``), so the raw fields are
    de-duplicated before retrieval.
    """
    sources: dict[tuple[str, str], str] = {}

    for name, spec in variables.items():
        dataset = spec["dataset"]
        if dataset not in config["datasets"]:
            raise ValueError(f"Variable '{name}' refers to unknown dataset '{dataset}'")

        for cds_var, short_name in spec["sources"].items():
            sources[(dataset, cds_var)] = short_name

    return [(dataset, cds_var, short) for (dataset, cds_var), short in sources.items()]


def raw_path(raw_dir: Path, short_name: str, years: list[int]) -> Path:
    """Path of a raw download covering ``years``."""
    return raw_dir / f"era5_{short_name}_{years[0]}-{years[-1]}.nc"


def download(args: argparse.Namespace, config: dict[str, Any]) -> None:
    """Retrieve raw monthly-mean fields from the CDS.

    Existing files are left alone, so an interrupted retrieval can be resumed by
    rerunning the same command.
    """
    import cdsapi

    variables = select_variables(config, args.variables)
    sources = collect_sources(variables, config)
    raw_dir = args.base_dir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)

    requests: list[tuple[str, dict[str, Any], Path]] = []
    for dataset, cds_var, short_name in sorted(sources):
        dataset_spec = config["datasets"][dataset]

        for years in year_chunks(
            args.start_year, args.end_year, dataset_spec.get("chunk_years", 1)
        ):
            target = raw_path(raw_dir, short_name, years)
            if target.exists():
                logger.info("Already downloaded, skipping: %s", target.name)
                continue

            request: dict[str, Any] = {
                "product_type": [dataset_spec["product_type"]],
                "variable": [cds_var],
                "year": [str(year) for year in years],
                "month": MONTHS,
                "time": ["00:00"],
                "data_format": "netcdf",
                "download_format": "unarchived",
            }
            if "pressure_levels" in dataset_spec:
                request["pressure_level"] = [
                    str(level) for level in dataset_spec["pressure_levels"]
                ]

            requests.append((dataset_spec["cds_name"], request, target))

    logger.info("%d retrieval(s) queued", len(requests))
    if args.dry_run:
        for cds_name, request, target in requests:
            logger.info("%s <- %s %s", target.name, cds_name, request["variable"])
        return

    client = cdsapi.Client()
    for index, (cds_name, request, target) in enumerate(requests, start=1):
        logger.info("[%d/%d] Retrieving %s", index, len(requests), target.name)
        # Download to a partial file first so an interrupted retrieval is not
        # mistaken for a complete one on the next run.
        partial = target.with_suffix(".nc.part")
        client.retrieve(cds_name, request, str(partial))
        partial.rename(target)


# -----------------------------------------------------------------------------
# process
# -----------------------------------------------------------------------------
def normalize_raw(ds: xr.Dataset, short_name: str) -> xr.DataArray:
    """Put a raw CDS file into e3sm_diags coordinate conventions.

    Renames the CDS coordinate names, drops the ERA5T ``expver`` and ensemble
    ``number`` dimensions, orders latitude south-to-north and pressure
    surface-to-top, and converts pressure from hPa to Pa.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset as downloaded from the CDS. It holds a single data variable.
    short_name : str
        Name to give the data variable, as used in the ``formula`` entries.

    Returns
    -------
    xr.DataArray
        The renamed data variable on ``time``/``plev``/``lat``/``lon``.
    """
    ds = ds.rename({old: new for old, new in RENAME_DIMS.items() if old in ds.dims})

    # Requests spanning the ERA5/ERA5T boundary can return both streams stacked
    # on an `expver` dimension. Prefer ERA5 (expver 1) and fall back to ERA5T
    # (expver 5) only where ERA5 is not yet available. The current CDS instead
    # labels each month with an `expver` coordinate, which needs no merging.
    if "expver" in ds.dims:
        expver = [str(v).lstrip("0") for v in ds["expver"].values]
        era5 = ds.isel(expver=expver.index("1")) if "1" in expver else None
        era5t = ds.isel(expver=expver.index("5")) if "5" in expver else None
        if era5 is None:
            ds = era5t
        elif era5t is None:
            ds = era5
        else:
            ds = era5.fillna(era5t)

    if "number" in ds.dims:
        ds = ds.isel(number=0, drop=True)

    # `number` (ensemble member) and `expver` (ERA5 vs ERA5T) also arrive as
    # scalar or time-varying coordinates, which do not belong in the output.
    ds = ds.drop_vars([name for name in ("number", "expver") if name in ds.coords])

    data_vars = [name for name in ds.data_vars if ds[name].ndim >= 3]
    if len(data_vars) != 1:
        raise ValueError(
            f"Expected exactly one field in the raw file, found {data_vars}"
        )
    da = ds[data_vars[0]].rename(short_name)

    # ERA5 is served north-to-south; the e3sm_diags reference files are
    # south-to-north.
    if da["lat"][0] > da["lat"][-1]:
        da = da.isel(lat=slice(None, None, -1))

    if "plev" in da.dims:
        # The CDS reports levels in hPa. The original files use Pa ordered from
        # the surface upward.
        if float(da["plev"].max()) <= 1100.0:
            da = da.assign_coords(plev=da["plev"].astype("float64") * 100.0)
        da = da.sortby("plev", ascending=False)

    return da.astype("float32")


def open_source(raw_dir: Path, short_name: str, years: list[int]) -> xr.DataArray:
    """Open every raw chunk holding ``short_name`` and concatenate along time."""
    paths = sorted(raw_dir.glob(f"era5_{short_name}_*.nc"))
    if not paths:
        raise FileNotFoundError(
            f"No raw files for '{short_name}' in {raw_dir}. Run the download step first."
        )

    arrays = []
    for path in paths:
        # `chunks={}` keeps the file's own chunking, so a 3D field streams
        # through dask instead of being read into memory whole.
        with xr.open_dataset(path, decode_timedelta=False, chunks={}) as ds:
            arrays.append(normalize_raw(ds, short_name))

    da = xr.concat(arrays, dim="time").sortby("time")
    da = da.sel(time=da["time"].dt.year.isin(years))

    expected = len(years) * 12
    if da.sizes["time"] != expected:
        raise ValueError(
            f"'{short_name}' has {da.sizes['time']} months, expected {expected} "
            f"for {years[0]}-{years[-1]}"
        )

    return da


def evaluate_formula(formula: str, sources: dict[str, xr.DataArray]) -> xr.DataArray:
    """Evaluate a variable's ``formula`` over its source fields.

    The formulas in ``era5_variables.yml`` are simple arithmetic over the source
    short names, e.g. ``msdwlwrf - msnlwrf`` or ``z / 9.80665``.
    """
    return eval(formula, {"__builtins__": {}}, dict(sources))  # noqa: S307


def month_midpoints(years: list[int]) -> tuple[pd.DatetimeIndex, np.ndarray]:
    """Build mid-month time stamps and month-boundary bounds.

    e3sm_diags (via xCDAT) weights climatologies by time bounds, so the bounds
    must span the whole month. The original files use the same convention.
    """
    starts = pd.date_range(
        f"{years[0]}-01-01", f"{years[-1]}-12-01", freq="MS", inclusive="both"
    )
    ends = starts + pd.offsets.MonthBegin(1)
    midpoints = starts + (ends - starts) / 2
    bounds = np.stack([starts.values, ends.values], axis=1)

    return midpoints, bounds


def build_dataset(
    da: xr.DataArray, name: str, spec: dict[str, Any], years: list[int]
) -> xr.Dataset:
    """Attach e3sm_diags metadata, mid-month times and coordinate bounds."""
    midpoints, time_bounds = month_midpoints(years)
    da = da.assign_coords(time=midpoints).rename(name)

    da.attrs = {
        key: spec[key]
        for key in ("standard_name", "long_name", "units", "comment")
        if key in spec
    }
    da.attrs["original_name"] = spec["formula"]
    da.attrs["institution"] = "ECMWF"

    ds = da.to_dataset()
    ds["time_bnds"] = xr.DataArray(time_bounds, dims=["time", "bnds"])
    ds["time"].attrs["bounds"] = "time_bnds"

    for axis, coord in (("X", "lon"), ("Y", "lat")):
        ds[coord].attrs.update(
            {
                "standard_name": "longitude" if coord == "lon" else "latitude",
                "long_name": "longitude" if coord == "lon" else "latitude",
                "units": "degrees_east" if coord == "lon" else "degrees_north",
                "axis": axis,
            }
        )
    if "plev" in ds.dims:
        ds["plev"].attrs.update(
            {
                "standard_name": "air_pressure",
                "long_name": "pressure",
                "units": "Pa",
                "positive": "down",
                "axis": "Z",
            }
        )

    ds.attrs = {
        "title": "ERA5 monthly means processed for e3sm_diags",
        "institution": "European Centre for Medium-Range Weather Forecasts",
        "source": (
            "ECMWF ERA5 reanalysis, monthly averaged reanalysis on a 0.25 degree "
            "lat/lon grid, retrieved from the Copernicus Climate Data Store. "
            "Generated using Copernicus Climate Change Service information; "
            "neither the European Commission nor ECMWF is responsible for any "
            "use that may be made of the Copernicus information or data it "
            "contains."
        ),
        "references": "https://doi.org/10.24381/cds.f17050d7",
        "frequency": "mon",
        "history": (
            f"{datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}: "
            f"{SCRIPT_NAME}: {name} = {spec['formula']} "
            f"from ERA5 monthly means for {years[0]}-{years[-1]}"
        ),
    }

    return ds


def time_series_path(directory: Path, name: str, years: list[int]) -> Path:
    """Path of a processed time series, using the e3sm_diags naming convention."""
    return directory / f"{name}_{years[0]}01_{years[-1]}12.nc"


def write_dataset(ds: xr.Dataset, path: Path, complevel: int) -> None:
    """Write a dataset with the encoding the reference files use."""
    encoding: dict[str, dict[str, Any]] = {
        "time": {
            "units": "days since 1900-01-01",
            "calendar": "gregorian",
            "dtype": "float64",
        },
        "time_bnds": {
            "units": "days since 1900-01-01",
            "calendar": "gregorian",
            "dtype": "float64",
        },
    }
    for name in ds.data_vars:
        if name == "time_bnds":
            continue
        encoding[str(name)] = {
            "dtype": "float32",
            "_FillValue": 1.0e20,
            "zlib": complevel > 0,
            "complevel": complevel,
        }

    path.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(path, unlimited_dims=["time"], encoding=encoding)


def process(args: argparse.Namespace, config: dict[str, Any]) -> None:
    """Convert raw CDS downloads into per-variable e3sm_diags time series."""
    variables = select_variables(config, args.variables)
    raw_dir = args.base_dir / "raw"
    out_dir = args.base_dir / "time_series"
    years = list(range(args.start_year, args.end_year + 1))

    for name, spec in variables.items():
        path = time_series_path(out_dir, name, years)
        if path.exists() and not args.overwrite:
            logger.info("Already processed, skipping: %s", path.name)
            continue

        logger.info("Processing %s = %s", name, spec["formula"])
        sources = {
            short: open_source(raw_dir, short, years)
            for short in spec["sources"].values()
        }
        da = evaluate_formula(spec["formula"], sources)
        ds = build_dataset(da, name, spec, years)
        write_dataset(ds, path, args.complevel)
        logger.info("Wrote %s", path)


# -----------------------------------------------------------------------------
# climo
# -----------------------------------------------------------------------------
def monthly_climatology(da: xr.DataArray) -> xr.DataArray:
    """Average each calendar month over all years.

    Months are weighted by their length so that February in leap years carries
    slightly more weight, matching how `ncclimo` builds monthly climatologies.
    """
    weights = da["time"].dt.days_in_month.astype("float64")
    weighted = (da * weights).groupby("time.month").sum(dim="time", skipna=True)
    norm = weights.groupby("time.month").sum(dim="time")

    return weighted / norm


def season_mean(
    monthly: xr.DataArray, months: list[int], years: list[int]
) -> xr.DataArray:
    """Combine monthly climatologies into a season, weighted by month length.

    Uses `ncclimo`'s "seasonally discontinuous December" convention: DJF pools
    every December in the period with every January and February, rather than
    dropping the unpaired months at the ends of the record.
    """
    # Month lengths averaged over the period, so that leap years are reflected
    # in the February weight.
    lengths = {
        month: float(
            np.mean(
                [
                    pd.Period(f"{year}-{month:02d}", freq="M").days_in_month
                    for year in years
                ]
            )
        )
        for month in months
    }
    total = sum(lengths.values())

    selected = monthly.sel(month=months)
    weights = xr.DataArray(
        [lengths[month] / total for month in months],
        coords={"month": months},
        dims=["month"],
    )

    return (selected * weights).sum(dim="month", skipna=True)


def climo_filename(season: str, years: list[int]) -> str:
    """Climatology file name, following the original ERA5 dataset's convention."""
    if season == "ANN" or season == "DJF":
        first, last = 1, 12
    elif season in SEASONS:
        first, last = SEASONS[season][0], SEASONS[season][-1]
    else:
        first = last = int(season)

    return f"ERA5_{season}_{years[0]}{first:02d}_{years[-1]}{last:02d}_climo.nc"


def climo(args: argparse.Namespace, config: dict[str, Any]) -> None:
    """Build seasonal and monthly climatologies from the processed time series."""
    variables = select_variables(config, args.variables)
    ts_dir = args.base_dir / "time_series"
    out_dir = args.base_dir / "climatology"
    years = list(range(args.start_year, args.end_year + 1))

    periods = list(SEASONS) + MONTHS
    outputs: dict[str, dict[str, xr.DataArray]] = {period: {} for period in periods}

    for name in variables:
        path = time_series_path(ts_dir, name, years)
        if not path.exists():
            raise FileNotFoundError(f"{path} not found. Run the process step first.")

        logger.info("Computing climatology for %s", name)
        with xr.open_dataset(path, chunks={"time": 12}) as ds:
            monthly = monthly_climatology(ds[name]).compute()

        for season, months in SEASONS.items():
            outputs[season][name] = season_mean(monthly, months, years)
        for month in MONTHS:
            outputs[month][name] = monthly.sel(month=int(month), drop=True)

    for period, fields in outputs.items():
        ds = xr.Dataset(fields)
        ds.attrs = {
            "title": f"ERA5 {period} climatology processed for e3sm_diags",
            "yrs_averaged": f"{years[0]}-{years[-1]}",
            "history": (
                f"{datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}: "
                f"{SCRIPT_NAME}: {period} climatology over {years[0]}-{years[-1]}"
            ),
        }

        path = out_dir / climo_filename(period, years)
        path.parent.mkdir(parents=True, exist_ok=True)
        encoding = {
            str(name): {
                "dtype": "float32",
                "_FillValue": 1.0e20,
                "zlib": args.complevel > 0,
                "complevel": args.complevel,
            }
            for name in ds.data_vars
        }
        ds.to_netcdf(path, encoding=encoding)
        logger.info("Wrote %s", path)


# -----------------------------------------------------------------------------
# compare
# -----------------------------------------------------------------------------
def compare(args: argparse.Namespace, config: dict[str, Any]) -> None:
    """Compare new time series against the original files over shared months.

    Reports the maximum absolute difference, RMSE and both global means for every
    variable that exists in both datasets, so unit, sign and grid-convention
    mistakes surface immediately.
    """
    variables = select_variables(config, args.variables)
    ts_dir = args.base_dir / "time_series"
    years = list(range(args.start_year, args.end_year + 1))

    rows: list[tuple[str, str]] = []
    for name in variables:
        new_path = time_series_path(ts_dir, name, years)
        old_paths = sorted(args.reference_dir.glob(f"{name}_??????_??????.nc"))

        if not new_path.exists():
            rows.append((name, "new time series missing"))
            continue
        if not old_paths:
            rows.append((name, "no counterpart in the original dataset"))
            continue

        with (
            xr.open_dataset(new_path, chunks={"time": 1}) as ds_new,
            xr.open_dataset(old_paths[0], chunks={"time": 1}) as ds_old,
        ):
            new = ds_new[name]
            old = ds_old[name]

            # The ext-script files keep the raw CDS dimension names. Without
            # this rename the two arrays share no horizontal dimension and the
            # subtraction below broadcasts to a 5-D array instead of aligning.
            old = old.rename({d: r for d, r in RENAME_DIMS.items() if d in old.dims})

            # The originals cover 1979-2019; line the two up on the months they
            # share, comparing by year and month rather than exact time stamps.
            stamp = lambda da: da["time"].dt.year * 100 + da["time"].dt.month  # noqa: E731
            shared = np.intersect1d(stamp(new).values, stamp(old).values)
            if shared.size == 0:
                rows.append((name, "no overlapping months"))
                continue

            new = new.isel(time=np.isin(stamp(new).values, shared))
            old = old.isel(time=np.isin(stamp(old).values, shared))

            if new.dims != old.dims:
                rows.append((name, f"dims {new.dims} vs {old.dims}"))
                continue

            if new.shape != old.shape:
                rows.append((name, f"shape {new.shape} vs {old.shape}"))
                continue

            diff = (new - old.assign_coords(time=new["time"])).compute()
            max_abs = float(abs(diff).max())
            rmse = float(np.sqrt((diff**2).mean()))
            new_mean = float(new.mean())
            old_mean = float(old.mean())
            scale = max(abs(new_mean), abs(old_mean), 1e-30)

            rows.append(
                (
                    name,
                    f"max|diff|={max_abs:.6g}  rmse={rmse:.6g}  "
                    f"mean new={new_mean:.6g} old={old_mean:.6g}  "
                    f"rel={(new_mean - old_mean) / scale:+.3%}  "
                    f"months={shared.size}",
                )
            )
        logger.info("%-8s %s", *rows[-1])

    print(
        f"\nComparison over {args.start_year}-{args.end_year} vs {args.reference_dir}"
    )
    for name, summary in rows:
        print(f"  {name:<8} {summary}")


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    for command, function in (
        ("download", download),
        ("process", process),
        ("climo", climo),
        ("compare", compare),
    ):
        sub = subparsers.add_parser(command, help=function.__doc__.splitlines()[0])
        sub.set_defaults(function=function)
        sub.add_argument("--start-year", type=int, default=1979)
        sub.add_argument("--end-year", type=int, default=2025)
        sub.add_argument(
            "-v",
            "--variables",
            nargs="+",
            help="Subset of variables to handle (default: all in era5_variables.yml)",
        )
        sub.add_argument("--base-dir", type=Path, default=DEFAULT_BASE_DIR)

        if command == "download":
            sub.add_argument(
                "--dry-run",
                action="store_true",
                help="List the retrievals that would be submitted, then stop",
            )
        if command in ("process", "climo"):
            sub.add_argument(
                "--complevel",
                type=int,
                default=1,
                help="netCDF deflate level; 0 disables compression (default: 1)",
            )
        if command == "process":
            sub.add_argument(
                "--overwrite",
                action="store_true",
                help="Reprocess variables whose output file already exists",
            )
        if command == "compare":
            sub.add_argument(
                "--reference-dir", type=Path, default=DEFAULT_REFERENCE_DIR
            )

    return parser


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )

    args = build_parser().parse_args(argv)
    if args.end_year < args.start_year:
        raise ValueError("--end-year must not precede --start-year")

    args.function(args, load_config())

    return 0


if __name__ == "__main__":
    sys.exit(main())

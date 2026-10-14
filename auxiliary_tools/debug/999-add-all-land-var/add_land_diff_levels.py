"""Add zero-centered diff_levels to the lat_lon_land model_vs_model cfg (#1017).

The land group asked that difference maps use the same 5th/95th percentile
cutoff as the contour maps, with zero centered and distinct colors for
positive and negative values.

Levels are static so plots stay comparable across runs. For each variable the
bound is max(|p5|, |p95|) of the difference field pooled over ANN, DJF, MAM,
JJA and SON, rounded up to a "nice" value. Eleven levels span -bound..bound so
zero is the middle boundary and no color bin straddles zero.

The difference fields are the ``*_diff.nc`` files written by running the
``lat_lon_land`` set with ``save_netcdf = True`` for v3.LR.historical_0051
(1987-1988) vs v3.LR.piControl (0401-0450).

Usage:
    python add_land_diff_levels.py <cfg> <results_dir> <report.json>
"""

import json
import math
import re
import sys
from pathlib import Path

import numpy as np
import xarray as xr

SEASONS = ["ANN", "DJF", "MAM", "JJA", "SON"]
N_PER_SIDE = 5
NICE_MANTISSAS = [1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0]
# plot/utils._get_c_levels_and_norm pads levels with -1e8 and 1e8.
SENTINEL = 1.0e8
MARKER = (
    "# diff_levels: zero-centered 5th/95th percentile of v3.LR.historical_0051 "
    "(1987-1988) minus v3.LR.piControl 0401-0450, all seasons"
)


def nice_ceil(x: float) -> float:
    """Round x up to m * 10**e with m in NICE_MANTISSAS."""
    exponent = math.floor(math.log10(x))
    for mantissa in NICE_MANTISSAS:
        candidate = mantissa * 10**exponent
        if candidate >= x * (1 - 1e-12):
            return float(f"{candidate:.6g}")

    return float(f"{10 ** (exponent + 1):.6g}")


def symmetric_levels(values: np.ndarray) -> list[float] | None:
    """Return 2 * N_PER_SIDE + 1 levels spanning -bound..bound, or None."""
    values = values[np.isfinite(values)]
    if values.size == 0:
        return None

    p5, p95 = np.percentile(values, [5, 95])
    bound = max(abs(p5), abs(p95))
    if not bound > 0:
        return None

    bound = nice_ceil(bound)
    if bound >= SENTINEL:
        return None

    step = bound / N_PER_SIDE
    return [float(f"{k * step:.6g}") for k in range(-N_PER_SIDE, N_PER_SIDE + 1)]


def format_levels(levels: list[float]) -> str:
    """Format levels so the cfg parser (yaml.safe_load) reads them as floats."""
    tokens = []
    for x in levels:
        s = f"{x:.6g}"
        if "e" in s and "." not in s:
            s = s.replace("e", ".0e", 1)
        tokens.append(s)

    return "[" + ",".join(tokens) + "]"


def main(cfg_path: str, results_dir: str, report_path: str) -> None:
    blocks = Path(cfg_path).read_text().split("[#]")
    report = []

    for i, block in enumerate(blocks):
        var_match = re.search(r'^variables = \["([^"]+)"\]', block, re.M)
        if var_match is None:
            continue

        var = var_match.group(1)
        case_id = re.search(r'^case_id = "([^"]+)"', block, re.M).group(1)
        entry = {"var": var, "case_id": case_id}

        if re.search(r"^diff_levels = ", block, re.M):
            entry["status"] = "has diff_levels, unchanged"
            report.append(entry)
            continue

        pooled = []
        for season in SEASONS:
            nc = Path(results_dir, "lat_lon_land", case_id).glob(
                f"*-{var}-{season}-global_diff.nc"
            )
            nc = next(nc, None)
            if nc is None:
                continue
            with xr.open_dataset(nc) as ds:
                da = ds[var] if var in ds else ds[list(ds.data_vars)[0]]
                pooled.append(da.values.ravel())

        if not pooled:
            entry["status"] = "no diff data, auto"
            report.append(entry)
            continue

        levels = symmetric_levels(np.concatenate(pooled))
        if levels is None:
            entry["status"] = "degenerate or out of range, auto"
            report.append(entry)
            continue

        line = f"diff_levels = {format_levels(levels)}"
        anchor = re.search(r"^contour_levels = .*$", block, re.M) or re.search(
            r"^diff_colormap = .*$", block, re.M
        )
        block = block[: anchor.end()] + "\n" + line + block[anchor.end() :]
        block = block.rstrip("\n") + "\n" + MARKER + "\n\n"
        blocks[i] = block

        entry.update(status="added", levels=levels)
        report.append(entry)

    Path(cfg_path).write_text("[#]".join(blocks))
    Path(report_path).write_text(json.dumps(report, indent=1))

    counts: dict[str, int] = {}
    for entry in report:
        counts[entry["status"]] = counts.get(entry["status"], 0) + 1
    print(counts)


if __name__ == "__main__":
    main(*sys.argv[1:4])

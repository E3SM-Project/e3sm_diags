"""Checks that `climo` reproduces `ncclimo`'s averaging conventions.

The original 1979-2019 ERA5 climatology was built with `ncclimo -a sdd`, so the
replacement workflow has to weight months the same way or the new files will
not line up with the ones they replace. Two conventions are easy to get wrong,
and both were wrong at one point:

* a monthly climatology weights every year equally, even though February is
  longer in leap years;
* a season weights its months by the fixed non-leap calendar, so February
  counts as 28 days regardless of how many leap years the period holds.

`tests/` is for the `e3sm_diags` package, and `pyproject.toml` limits pytest to
`tests/e3sm_diags`, so this lives beside the script it covers. Run it with:

    pytest analysis_data_preprocess/ERA5/test_climo.py
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import xarray as xr

from era5_pipeline import MONTH_LENGTHS, SEASONS, monthly_climatology, season_mean

YEARS = [1979, 1980, 1981, 1982]  # 1980 is a leap year


def _series(values_by_month: dict[int, float] | None = None) -> xr.DataArray:
    """A monthly series over `YEARS`, 1.0 everywhere unless told otherwise."""
    time = pd.date_range(f"{YEARS[0]}-01-01", f"{YEARS[-1]}-12-01", freq="MS")
    data = np.ones(len(time), dtype="float64")

    if values_by_month:
        for i, stamp in enumerate(time):
            if stamp.month in values_by_month:
                data[i] = values_by_month[stamp.month]

    return xr.DataArray(data, coords={"time": time}, dims=["time"])


def test_monthly_climatology_weights_every_year_equally():
    # February 1980 is a leap February; give it a value the others do not have
    # so a length weighting would show up as a pull towards it.
    da = _series()
    da[da["time"].dt.month.values == 2] = [0.0, 10.0, 0.0, 0.0]

    february = float(monthly_climatology(da).sel(month=2))

    # Plain mean of the four years. Weighting by 29/28 days would give 2.5885.
    assert february == pytest.approx(10.0 / 4)


def test_season_mean_uses_fixed_non_leap_month_lengths():
    monthly = monthly_climatology(_series({12: 0.0, 1: 0.0, 2: 1.0}))

    djf = float(season_mean(monthly, SEASONS["DJF"]))

    # February's share of DJF is 28 / (31 + 31 + 28), not 28.25 / 90.25.
    assert djf == pytest.approx(28 / 90)


def test_season_weights_do_not_depend_on_the_period():
    """A period with more leap years must not shift the weights."""
    short = monthly_climatology(_series({12: 0.0, 1: 0.0, 2: 1.0}))
    assert float(season_mean(short, SEASONS["DJF"])) == pytest.approx(28 / 90)


def test_annual_mean_weights_the_whole_calendar():
    monthly = monthly_climatology(_series({m: 0.0 for m in range(2, 13)}))

    assert float(season_mean(monthly, SEASONS["ANN"])) == pytest.approx(31 / 365)


def test_month_lengths_are_the_non_leap_calendar():
    assert sum(MONTH_LENGTHS.values()) == 365
    assert MONTH_LENGTHS[2] == 28

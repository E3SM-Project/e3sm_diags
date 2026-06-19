"""Numerically verify the index_timeseries output against an a-prime-style,
hand-coded numpy computation of the nino index.

Reproduces a-prime's algorithm (get_reg_box -> cos-lat area average ->
remove_seasonal_cycle_monthly_data) with raw numpy and compares to the netCDF
saved by the new e3sm_diags index_timeseries diagnostic.
"""

import numpy as np
import xcdat as xc

TS_FILE = (
    "/global/cfs/cdirs/e3sm/e3sm_diags/postprocessed_e3sm_v2_data_for_e3sm_diags/"
    "20210528.v2rc3e.piControl.ne30pg2_EC30to60E2r2.chrysalis/time-series/rgr/"
    "TS_005101_006012.nc"
)
RESULTS = "/global/cfs/cdirs/e3sm/www/chengzhu/tests/enso_index_timeseries_test"
NC_OUT = f"{RESULTS}/enso_diags/NINO-index/nino-index-timeseries_test.nc"
NC_SEASONALITY = f"{RESULTS}/enso_diags/NINO-seasonality/nino-index-seasonality_test.nc"

# a-prime get_reg_box boxes (lat_ll, lat_ul, lon_ll, lon_ul).
BOXES = {
    "NINO3": (-5, 5, 210, 270),
    "NINO34": (-5, 5, 190, 240),
    "NINO4": (-5, 5, 160, 210),
}


def aprime_index(ds, var, box):
    lat_ll, lat_ul, lon_ll, lon_ul = box
    sub = ds.sel(lat=slice(lat_ll, lat_ul), lon=slice(lon_ll, lon_ul))[var]

    lat = sub["lat"].values
    # a-prime area-weights with the model area field; on a rectilinear grid this
    # reduces to cos(lat) weighting (longitude weights are uniform).
    wgts = np.cos(np.deg2rad(lat))[:, None] * np.ones(sub["lon"].size)[None, :]

    field = sub.values  # (time, lat, lon)
    area_avg = np.array(
        [np.sum(field[t] * wgts) / np.sum(wgts) for t in range(field.shape[0])]
    )

    # remove_seasonal_cycle_monthly_data
    nyrs = area_avg.shape[0] // 12
    yrs = np.arange(nyrs)
    anom = area_avg + np.nan
    for m in range(12):
        idx = 12 * yrs + m
        anom[idx] = area_avg[idx] - np.mean(area_avg[idx])
    return anom


ds = xc.open_dataset(TS_FILE)
out = xc.open_dataset(NC_OUT)
out_seasonality = xc.open_dataset(NC_SEASONALITY)

print("Index time series:")
print(f"{'region':8s} {'max|diff|':>12s} {'rms diff':>12s} {'a-prime std':>12s}")
aprime_indices = {}
for region, box in BOXES.items():
    ref = aprime_index(ds, "TS", box)
    aprime_indices[region] = ref
    got = out[region].values
    diff = got - ref
    print(
        f"{region:8s} {np.nanmax(np.abs(diff)):12.3e} "
        f"{np.sqrt(np.nanmean(diff**2)):12.3e} {np.nanstd(ref):12.4f}"
    )

print("\nSeasonality (per-month std dev):")
print(f"{'region':8s} {'max|diff|':>12s} {'rms diff':>12s}")
for region in BOXES:
    index = aprime_indices[region]
    nyrs = index.shape[0] // 12
    # a-prime: std across years for each calendar month (Jan..Dec).
    ref_seasonality = np.nanstd(index[: nyrs * 12].reshape(nyrs, 12), axis=0)
    got_seasonality = out_seasonality[region].values
    diff = got_seasonality - ref_seasonality
    print(
        f"{region:8s} {np.nanmax(np.abs(diff)):12.3e} "
        f"{np.sqrt(np.nanmean(diff**2)):12.3e}"
    )

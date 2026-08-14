"""Numerically verify the ported enso_diags plot types against a-prime-style,
hand-coded numpy computations.

Reproduces a-prime's algorithms with raw numpy and compares to the netCDF saved
by the new e3sm_diags diagnostics:

* nino_index_timeseries / seasonality: get_reg_box -> cos-lat area average ->
  remove_seasonal_cycle_monthly_data.
* interannual_variability: compute_reg_seasonal_climo_and_stddev (day-weighted
  annual means -> mean & std dev across years).
"""

import numpy as np
import xcdat as xc

DATA = (
    "/global/cfs/cdirs/e3sm/e3sm_diags/postprocessed_e3sm_v2_data_for_e3sm_diags/"
    "20210528.v2rc3e.piControl.ne30pg2_EC30to60E2r2.chrysalis/time-series/rgr/"
)
TS_FILE = f"{DATA}/TS_005101_006012.nc"
OCNFRAC_FILE = f"{DATA}/OCNFRAC_005101_006012.nc"
PSL_FILE = f"{DATA}/PSL_005101_006012.nc"

RESULTS = "/global/cfs/cdirs/e3sm/www/chengzhu/tests/enso_aprime_diags_test"
NC_OUT = f"{RESULTS}/enso_diags/NINO-index/nino-index-timeseries_test.nc"
NC_SEASONALITY = f"{RESULTS}/enso_diags/NINO-seasonality/nino-index-seasonality_test.nc"
VAR_DIR = f"{RESULTS}/enso_diags/SST-variability"
# The mean and std fields are saved as SST_mean / SST_std in one _test.nc file.
NC_VAR = f"{VAR_DIR}/interannual-variability-sst-20s20n_test.nc"

# a-prime get_reg_box boxes (lat_ll, lat_ul, lon_ll, lon_ul).
BOXES = {
    "NINO3": (-5, 5, 210, 270),
    "NINO34": (-5, 5, 190, 240),
    "NINO4": (-5, 5, 160, 210),
}

# a-prime EQSOI sea level pressure regions.
EQSOI_BOXES = {
    "EPAC": (-5, 5, 230, 280),
    "INDO": (-5, 5, 90, 140),
}
NC_EQSOI = f"{RESULTS}/enso_diags/EQSOI/equatorial-soi_test.nc"


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


def aprime_climo_and_std(box_lat=(-20, 20)):
    """Reimplement a-prime compute_reg_seasonal_climo_and_stddev for SST.

    SST is derived the same way as e3sm_diags: TS converted to degC, masked to
    ocean points (OCNFRAC >= 0.9). The series is aggregated into day-weighted
    annual means, then reduced across years into the climatological mean and the
    interannual standard deviation.
    """
    ds_ts = xc.open_dataset(TS_FILE)
    ds_ocn = xc.open_dataset(OCNFRAC_FILE)

    ts = ds_ts.sel(lat=slice(*box_lat))["TS"]
    ocnfrac = ds_ocn.sel(lat=slice(*box_lat))["OCNFRAC"]

    # SST in degC, masked over non-ocean (OCNFRAC < 0.9).
    sst = (ts - 273.15).where(ocnfrac >= 0.9)
    field = sst.values  # (time, lat, lon)

    # Day weights from the time bounds (days per month).
    time_bnds = ds_ts.bounds.get_bounds(axis="T").values
    day_wgts = (time_bnds[:, 1] - time_bnds[:, 0]) / np.timedelta64(1, "D")

    nt = field.shape[0]
    nyrs = nt // 12
    annual = np.empty((nyrs,) + field.shape[1:])
    for y in range(nyrs):
        sl = slice(12 * y, 12 * (y + 1))
        w = day_wgts[sl]
        annual[y] = np.sum(field[sl] * w[:, None, None], axis=0) / np.sum(w)

    return np.nanmean(annual, axis=0), np.nanstd(annual, axis=0)


print("\nEquatorial SOI (EQSOI = z(EPAC) - z(INDO) from PSL):")
print(f"{'field':8s} {'max|diff|':>12s} {'rms diff':>12s}")
ds_psl = xc.open_dataset(PSL_FILE)


def aprime_standardize(field):
    # standardize_time_series: (x - mean) / std (population std, ddof=0).
    return (field - np.nanmean(field)) / np.nanstd(field)


# a-prime: per-region anomaly -> standardize -> difference (EPAC minus INDO).
z_epac = aprime_standardize(aprime_index(ds_psl, "PSL", EQSOI_BOXES["EPAC"]))
z_indo = aprime_standardize(aprime_index(ds_psl, "PSL", EQSOI_BOXES["INDO"]))
ref_eqsoi = z_epac - z_indo

ds_eqsoi = xc.open_dataset(NC_EQSOI)
got_eqsoi = ds_eqsoi["EQSOI"].values
diff = got_eqsoi - ref_eqsoi
print(
    f"{'EQSOI':8s} {np.nanmax(np.abs(diff)):12.3e} "
    f"{np.sqrt(np.nanmean(diff**2)):12.3e}"
)

print("\nInterannual variability (SST, 20S20N):")
print(f"{'field':8s} {'max|diff|':>12s} {'rms diff':>12s}")
ref_mean, ref_std = aprime_climo_and_std()
ds_var = xc.open_dataset(NC_VAR)
for label, ref, nc_var in [
    ("mean", ref_mean, "SST_mean"),
    ("std", ref_std, "SST_std"),
]:
    got = ds_var[nc_var].values
    diff = got - ref
    print(
        f"{label:8s} {np.nanmax(np.abs(diff)):12.3e} "
        f"{np.sqrt(np.nanmean(diff**2)):12.3e}"
    )

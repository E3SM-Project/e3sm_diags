# Running E3SM Diags on EAMxx Output

## Prerequisites

0. Secure an interactive compute node and activate the E3SM-Unified environment:

```
salloc --nodes 1 --qos interactive --time 02:00:00 --constraint cpu --account e3sm
source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh
```

(The version of E3SM Diags (v3) that has EAMxx variable support is available in E3SM-Unified v1.11 (Mid Feb 2025 release.)

## EAMxx ne256pg2 Case

1. Use `test_eamxx_ne256.sh` to generate climatology, diurnal cycle climo, and daily time-series from ne256pg2 EAMxx output:

```
bash test_eamxx_ne256.sh
```

This script uses NCO/ncclimo to remap and generate:
- Monthly climatology files
- Diurnal cycle (3-hourly) climatology files
- Per-variable 3-hourly time-series, then averaged to daily

2. Run E3SM Diags:

```
python run_e3sm_diags_eamxx256_climo.py
```

This runs the following diagnostic sets: lat_lon, zonal_mean_xy, zonal_mean_2d, zonal_mean_2d_stratosphere, polar, meridional_mean_2d, annual_cycle_zonal_mean, tropical_subseasonal, and precip_pdf.

## Archived: Earlier EAMxx ne256pg2 Case (rainfrac1)

The original scripts for the `ne256pg2_ne256pg2.F20TR-SCREAMv1.rainfrac1` case are preserved in `case1_rainfrac1/`:

- `case1_rainfrac1/nco.sh`: NCO remapping script
- `case1_rainfrac1/run_e3sm_diags_1996.py`: compare 1996 climatology to available 1990 obs climatology
- `case1_rainfrac1/run_e3sm_diags_climo.py`: compare 1996 climatology to pre-calculated obs climatology

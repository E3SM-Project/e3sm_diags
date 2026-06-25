# Running E3SM Diags with Arbitrary Start Month (S2D)

This example demonstrates using the `start_month` parameter to compute climatologies with a non-January annual cycle start. This is useful for seasonal-to-decadal (S2D) diagnostics where the annual cycle may start in a month other than January.

## Setup

0. Secure an interactive compute node and activate the E3SM-Unified environment:

```
salloc --nodes 1 --qos interactive --time 02:00:00 --constraint cpu --account e3sm

source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh
```

## Prepare Data

1. Use ncclimo to generate per-variable monthly time series from E3SM output. An example script is provided in `run_nco.sh`. The `--mth_srt` and `--mth_end` flags define the annual cycle boundaries (e.g., May through April):

```
bash run_nco.sh
```

## Run Diagnostics

2. Run the diagnostics with `start_month` set to match the ncclimo cycle:

```
python run_e3sm_diags_s2d.py
```

## Key Parameter

`param.start_month = 5` sets the annual cycle to start in May. Each cycle year then covers May of year N through April of year N+1. The default is 1 (January).

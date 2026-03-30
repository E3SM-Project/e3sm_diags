# Example 11: AIRS Spectral OLR Diagnostics

Compare E3SM model output against AIRS satellite spectral OLR observations.

## Overview

This diagnostic set evaluates E3SM spectral outgoing longwave radiation (OLR)
against observations derived from the Atmospheric Infrared Sounder (AIRS).
It covers broadband and spectral band fluxes (bands 02 and 06), longwave cloud
forcing, and fractional contributions of each spectral band.

The AIRS spectral OLR data and diagnostic methodology were contributed by
Professor Xianglei Huang's group at the University of Michigan.

## Variables

| Variable | Description |
|---|---|
| FLUT / FLUTC | Broadband all-sky / clear-sky OLR |
| FLSU02 / FLSU06 | Band 02 (350-500 cm-1) / Band 06 (820-980 cm-1) all-sky OLR |
| FLSUCLR02 / FLSUCLR06 | Band 02 / Band 06 clear-sky OLR |
| LWCF | TOA longwave cloud forcing |
| LWCF02 / LWCF06 | Band 02 / Band 06 cloud forcing |
| FLSU02_FRAC / FLSU06_FRAC | Fractional contribution of band to all-sky OLR |
| FLSUCLR02_FRAC / FLSUCLR06_FRAC | Fractional contribution of band to clear-sky OLR |
| LWCF02_FRAC / LWCF06_FRAC | Fractional contribution of band to LW cloud forcing |

## How to run

1. Update the paths in `ex11.py` (test data, reference data, output directory).
2. Run: `python ex11.py -d diags.cfg`

## User guide

For details on interpreting the diagnostics, see:
`E3SM_spectralOLR_diagnostics.pdf`

"""Module for defining regions used for spatial subsetting.

This module is the refactored, Xarray version of `default_regions.py`, which is
based on `cdutil`.
"""

# "lat": The latitude domain for subsetting a variable.
# "lon": The longitude domain for subsetting a variable.
# "value": The lower limit for masking.
region_specs = {
    "global": {},
    "land": {"value": 0.65},
    "ocean": {"value": 0.65},
    # North American Monsoon
    "NAMM": {"lat": (0.0, 45.0), "lon": (210.0, 310.0)},
}

"""Module for defining regions used for spatial subsetting."""

# A dictionary storing the specifications for each region.
# "lat": The latitude domain for subsetting a variable, (lon_west, lon_east).
# "lon": The longitude domain for subsetting a variable (lat_west, lat_east).
# "value": The lower limit for masking.
REGION_SPECS = {
    "global": {},
    "NHEX": {"lat": (30.0, 90)},
    "SHEX": {"lat": (-90.0, -30)},
    "TROPICS": {"lat": (-30.0, 30)},
    "TRMM_region": {"lat": (-38.0, 38)},
    "90S50S": {"lat": (-90.0, -50)},
    "50S20S": {"lat": (-50.0, -20)},
    "20S20N": {"lat": (-20.0, 20)},
    "50S50N": {"lat": (-50.0, 50)},
    "5S5N": {"lat": (-5.0, 5)},
    "20N50N": {"lat": (20.0, 50)},
    "50N90N": {"lat": (50.0, 90)},
    "60S90N": {"lat": (-60.0, 90)},
    "60S60N": {"lat": (-60.0, 60)},
    "75S75N": {"lat": (-75.0, 75)},
    "ocean": {"value": 0.65},
    "ocean_seaice": {"value": 0.65},
    "land": {"value": 0.65},
    "land_60S90N": {"value": 0.65, "lat": (-60.0, 90)},
    "ocean_TROPICS": {"value": 0.65, "lat": (-30.0, 30)},
    "land_NHEX": {"value": 0.65, "lat": (30.0, 90)},
    "land_SHEX": {"value": 0.65, "lat": (-90.0, -30)},
    "land_TROPICS": {"value": 0.65, "lat": (-30.0, 30)},
    "ocean_NHEX": {"value": 0.65, "lat": (30.0, 90)},
    "ocean_SHEX": {"value": 0.65, "lat": (-90.0, -30)},
    # follow AMWG polar range,more precise selector
    "polar_N": {"lat": (50.0, 90.0)},
    "polar_S": {"lat": (-90.0, -55.0)},
    # To match AMWG results, the bounds is not as precise in this case
    # 'polar_N_AMWG':{'domain': Selector("lat":(50., 90.))},
    # 'polar_S_AMWG':{'domain': Selector("lat":(-90., -55.))},
    # Below is for modes of variability
    "NAM": {"lat": (20.0, 90), "lon": (-180, 180)},
    "NAO": {"lat": (20.0, 80), "lon": (-90, 40)},
    "SAM": {"lat": (-20.0, -90), "lon": (0, 360)},
    "PNA": {"lat": (20.0, 85), "lon": (120, 240)},
    "PDO": {"lat": (20.0, 70), "lon": (110, 260)},
    # Below is for monsoon domains
    # All monsoon domains
    "AllM": {"lat": (-45.0, 45.0), "lon": (0.0, 360.0)},
    # North American Monsoon
    "NAMM": {"lat": (0, 45.0), "lon": (210.0, 310.0)},
    # South American Monsoon
    "SAMM": {"lat": (-45.0, 0.0), "lon": (240.0, 330.0)},
    # North African Monsoon
    "NAFM": {"lat": (0.0, 45.0), "lon": (310.0, 60.0)},
    # South African Monsoon
    "SAFM": {"lat": (-45.0, 0.0), "lon": (0.0, 90.0)},
    # Asian Summer Monsoon
    "ASM": {"lat": (0.0, 45.0), "lon": (60.0, 180.0)},
    # Australian Monsoon
    "AUSM": {"lat": (-45.0, 0.0), "lon": (90.0, 160.0)},
    # Below is for NINO domains.
    "NINO3": {"lat": (-5.0, 5.0), "lon": (210.0, 270.0)},
    "NINO34": {"lat": (-5.0, 5.0), "lon": (190.0, 240.0)},
    "NINO4": {"lat": (-5.0, 5.0), "lon": (160.0, 210.0)},
    # Below is for additional domains for diurnal cycle of precipitation
    "W_Pacific": {"lat": (-20.0, 20.0), "lon": (90.0, 180.0)},
    "CONUS": {"lat": (25.0, 50.0), "lon": (-125.0, -65.0)},
    "Amazon": {"lat": (-20.0, 5.0), "lon": (-80.0, -45.0)},
    # Below is for RRM(regionally refined model) domains.
    # 'CONUS_RRM': {'domain': "lat":(20., 50., 'ccb'), "lon":(-125., -65., 'ccb'))},For RRM dataset, negative value won't work
    "CONUS_RRM": {"lat": (20.0, 50.0), "lon": (235.0, 295.0)},
    # Below is for debugging. A smaller latitude range reduces processing time.
    "DEBUG": {"lat": (-2.0, 2)},
}

# A dictionary storing ARM site specifications with specific coordinates.
# Select nearest grid point to ARM site coordinate.
# "lat": The latitude point.
# "lon": The longitude point.
# "description": The description of the ARM site.
ARM_SITE_SPECS = {
    "sgpc1": {"lat": 36.4, "lon": -97.5, "description": "97.5W 36.4N Oklahoma ARM"},
    "nsac1": {"lat": 71.3, "lon": -156.6, "description": "156.6W 71.3N Barrow ARM"},
    "twpc1": {"lat": -2.1, "lon": 147.4, "description": "147.4E 2.1S Manus ARM"},
    "twpc2": {"lat": -0.5, "lon": 166.9, "description": "166.9E 0.5S Nauru ARM"},
    "twpc3": {"lat": -12.4, "lon": 130.9, "description": "130.9E 12.4S Darwin ARM"},
    "enac1": {
        "lat": 39.1,
        "lon": -28.0,
        "description": "28.0E 39.1N Graciosa Island ARM",
    },
}

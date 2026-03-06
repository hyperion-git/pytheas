"""
Pytheas -- body-tide prediction on a normal-gravity baseline

Named after Pytheas of Massalia (c. 325 BC), the Greek explorer who made the
first recorded systematic observations connecting ocean tides to the Moon.

Computes the normal-gravity background plus lunisolar body tides along an
arbitrary measurement axis, using exact Newtonian point-mass forcing, Meeus
ephemerides, and a single scalar elastic-response factor.

Quoted ~200-1000 nGal accuracy applies to the vertical body-tide channel at
quiet inland sites.  Pytheas is not a complete terrestrial gravity model:
loading, local anomalies, and full IERS constituent-/component-dependent
response corrections are omitted.
Dependency: numpy only.
"""

__version__ = "3.4.0"

from ._core import (
    # Main API
    compute_g,
    compute_timeseries,
    GravityResult,
    TimeSeries,
    # Lab frame
    GravityField,
    LabFrame,
    # Building blocks
    normal_gravity,
    sun_position_ecef,
    moon_position_ecef,
    geodetic_to_ecef,
    enu_basis,
    measurement_axis,
    tidal_acceleration,
    julian_date,
    gmst_rad,
    # Constants
    GM_MOON,
    GM_SUN,
    GM_E,
    A_WGS84,
    F_WGS84,
    B_WGS84,
    E2,
    OMEGA,
    AU,
    H2,
    K2,
    DELTA_GRAV,
    GAMMA_E,
    GAMMA_P,
)

__all__ = [
    "compute_g",
    "compute_timeseries",
    "GravityResult",
    "TimeSeries",
    "GravityField",
    "LabFrame",
    "normal_gravity",
    "sun_position_ecef",
    "moon_position_ecef",
    "geodetic_to_ecef",
    "enu_basis",
    "measurement_axis",
    "tidal_acceleration",
    "julian_date",
    "gmst_rad",
    "GM_MOON",
    "GM_SUN",
    "GM_E",
    "A_WGS84",
    "F_WGS84",
    "B_WGS84",
    "E2",
    "OMEGA",
    "AU",
    "H2",
    "K2",
    "DELTA_GRAV",
    "GAMMA_E",
    "GAMMA_P",
]

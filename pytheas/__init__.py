"""
Pytheas -- gravitational acceleration g(t) at a point on Earth

Named after Pytheas of Massalia (c. 325 BC), the Greek explorer who made the
first recorded systematic observations connecting ocean tides to the Moon.

Computes the total gravitational acceleration projected along an arbitrary
measurement axis, including normal gravity (WGS84), lunar and solar tidal
acceleration (exact Newtonian, Meeus ephemeris), and elastic Earth
amplification (IERS 2010 Love numbers).

Accuracy: ~10-100 nGal (1e-10 to 1e-9 m/s^2) for inland sites.
Dependency: numpy only.
"""

__version__ = "3.1.1"

from ._core import (
    # Main API
    compute_g,
    compute_timeseries,
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

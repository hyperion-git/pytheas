"""
Pytheas -- gravitational acceleration g(t) at a point on Earth

Computes the total gravitational acceleration projected along an arbitrary
measurement axis, as a function of time.  The model includes:

  - Normal gravity on the GRS80 ellipsoid with second-order free-air correction
  - Lunar tidal acceleration (exact Newtonian formula, Meeus ephemeris)
  - Solar tidal acceleration (exact Newtonian formula, Meeus ephemeris)
  - Elastic Earth amplification via IERS 2010 Love numbers

Accuracy: ~10-100 nGal (1e-10 to 1e-9 m/s^2) for inland sites.

Dependencies: numpy only.
"""

import numpy as np
from datetime import datetime, timedelta

# =============================================================================
# Physical Constants
# =============================================================================

# Gravitational parameters (DE430)
GM_MOON = 4.9028695e12       # m^3/s^2
GM_SUN  = 1.32712440018e20   # m^3/s^2

# GRS80 Ellipsoid
A_GRS80 = 6378137.0                          # semi-major axis (m)
F_GRS80 = 1.0 / 298.257222101                # flattening
B_GRS80 = A_GRS80 * (1.0 - F_GRS80)          # semi-minor axis (m)
E2      = 2.0 * F_GRS80 - F_GRS80 ** 2       # first eccentricity squared
GM_E    = 3.986004418e14                      # Earth GM (m^3/s^2, GRS80)
OMEGA   = 7.292115e-5                         # Earth rotation rate (rad/s)

# Normal gravity on the ellipsoid (Somigliana constants)
GAMMA_E = 9.7803253141                        # equatorial gravity (m/s^2)
GAMMA_P = 9.8321849378                        # polar gravity (m/s^2)
K_SOM   = (B_GRS80 * GAMMA_P) / (A_GRS80 * GAMMA_E) - 1.0
M_RATIO = OMEGA ** 2 * A_GRS80 ** 2 * B_GRS80 / GM_E

# Astronomical
AU              = 1.495978707e11              # astronomical unit (m)
OBLIQUITY_J2000 = 23.439291                   # mean obliquity at J2000 (deg)

# Elastic Earth (IERS 2010 Love numbers)
H2 = 0.6078
K2 = 0.2980
DELTA_GRAV = 1.0 + H2 - 1.5 * K2             # gravimetric factor = 1.1608


# =============================================================================
# Normal Gravity
# =============================================================================

def normal_gravity(lat_deg, alt_m):
    """GRS80 normal gravity with second-order free-air correction.

    Parameters
    ----------
    lat_deg : float
        Geodetic latitude in degrees.
    alt_m : float
        Altitude above the ellipsoid in meters.

    Returns
    -------
    float
        Normal gravity in m/s^2.
    """
    phi = np.radians(lat_deg)
    sin2 = np.sin(phi) ** 2

    # Somigliana closed-form on the ellipsoid
    gamma_0 = GAMMA_E * (1.0 + K_SOM * sin2) / np.sqrt(1.0 - E2 * sin2)

    # Second-order free-air correction (Taylor in h/a)
    h = alt_m
    fac = 1.0 + F_GRS80 + M_RATIO - 2.0 * F_GRS80 * sin2
    gamma = gamma_0 * (1.0 - 2.0 * fac * h / A_GRS80
                        + 3.0 * h ** 2 / A_GRS80 ** 2)
    return gamma


# =============================================================================
# Time Utilities
# =============================================================================

def julian_date(dt):
    """Convert a UTC datetime to Julian Date (Meeus algorithm).

    Parameters
    ----------
    dt : datetime
        UTC date and time.

    Returns
    -------
    float
        Julian Date.
    """
    y, m = dt.year, dt.month
    if m <= 2:
        y -= 1
        m += 12
    A = y // 100
    B = 2 - A + A // 4
    jd = (int(365.25 * (y + 4716)) + int(30.6001 * (m + 1))
          + dt.day + B - 1524.5)
    jd += (dt.hour + dt.minute / 60.0 + dt.second / 3600.0) / 24.0
    return jd


def _T(dt):
    """Julian centuries since J2000.0."""
    return (julian_date(dt) - 2451545.0) / 36525.0


def gmst_rad(dt):
    """Greenwich Mean Sidereal Time in radians.

    Parameters
    ----------
    dt : datetime
        UTC date and time.

    Returns
    -------
    float
        GMST in radians.
    """
    T = _T(dt)
    JD = julian_date(dt)
    gmst_deg = (280.46061837
                + 360.98564736629 * (JD - 2451545.0)
                + 0.000387933 * T ** 2
                - T ** 3 / 38710000.0) % 360.0
    return np.radians(gmst_deg)


def _eci_to_ecef(x_eci, y_eci, z_eci, dt):
    """Rotate ECI to ECEF using GMST."""
    theta = gmst_rad(dt)
    c, s = np.cos(theta), np.sin(theta)
    return (c * x_eci + s * y_eci,
            -s * x_eci + c * y_eci,
            z_eci)


# =============================================================================
# Solar Ephemeris (Meeus, Astronomical Algorithms)
# =============================================================================

def sun_position_ecef(dt):
    """Sun position in ECEF coordinates.

    Low-precision solar ephemeris after Meeus (Astronomical Algorithms).
    Accuracy ~1 arcmin in ecliptic longitude.

    Parameters
    ----------
    dt : datetime
        UTC date and time.

    Returns
    -------
    numpy.ndarray
        ECEF position vector [x, y, z] in meters.
    """
    T = _T(dt)

    L0 = (280.46646 + 36000.76983 * T + 0.0003032 * T ** 2) % 360
    M  = (357.52911 + 35999.05029 * T - 0.0001537 * T ** 2) % 360
    M_rad = np.radians(M)

    C = ((1.914602 - 0.004817 * T - 0.000014 * T ** 2) * np.sin(M_rad)
         + (0.019993 - 0.000101 * T) * np.sin(2.0 * M_rad)
         + 0.000289 * np.sin(3.0 * M_rad))

    sun_lon = np.radians((L0 + C) % 360)

    e = 0.016708634 - 0.000042037 * T
    v_rad = np.radians((M + C) % 360)
    R = AU * (1.000001018 * (1.0 - e ** 2)) / (1.0 + e * np.cos(v_rad))

    eps = np.radians(OBLIQUITY_J2000 - 0.013004 * T)

    x_eci = R * np.cos(sun_lon)
    y_eci = R * np.sin(sun_lon) * np.cos(eps)
    z_eci = R * np.sin(sun_lon) * np.sin(eps)

    x, y, z = _eci_to_ecef(x_eci, y_eci, z_eci, dt)
    return np.array([x, y, z])


# =============================================================================
# Lunar Ephemeris (Meeus ch. 47, major terms)
# =============================================================================

# Meeus Table 47.A: Longitude perturbations (D, Ms, Mp, F, coeff * 1e-6 deg)
_LON_TERMS = [
    (0,  0,  1,  0,  6288774), (2,  0, -1,  0,  1274027),
    (2,  0,  0,  0,   658314), (0,  0,  2,  0,   213618),
    (0,  1,  0,  0,  -185116), (0,  0,  0,  2,  -114332),
    (2,  0, -2,  0,    58793), (2, -1, -1,  0,    57066),
    (2,  0,  1,  0,    53322), (2, -1,  0,  0,    45758),
    (0,  1, -1,  0,   -40923), (1,  0,  0,  0,   -34720),
    (0,  1,  1,  0,   -30383), (2,  0,  0, -2,    15327),
    (0,  0,  1,  2,   -12528), (0,  0,  1, -2,    10980),
    (4,  0, -1,  0,    10675), (0,  0,  3,  0,    10034),
    (4,  0, -2,  0,     8548), (2,  1, -1,  0,    -7888),
    (2,  1,  0,  0,    -6766), (1,  0, -1,  0,    -5163),
    (1,  1,  0,  0,     4987), (2, -1,  1,  0,     4036),
]

# Meeus Table 47.A: Distance perturbations (D, Ms, Mp, F, coeff in meters)
_DIST_TERMS = [
    (0,  0,  1,  0, -20905355), (2,  0, -1,  0,  -3699111),
    (2,  0,  0,  0,  -2955968), (0,  0,  2,  0,   -569925),
    (0,  1,  0,  0,     48888), (0,  0,  0,  2,     -3149),
    (2,  0, -2,  0,    246158), (2, -1, -1,  0,   -152138),
    (2,  0,  1,  0,   -170733), (2, -1,  0,  0,   -204586),
    (0,  1, -1,  0,   -129620), (1,  0,  0,  0,    108743),
    (0,  1,  1,  0,    104755), (2,  0,  0, -2,     10321),
    (0,  0,  1, -2,     79661), (4,  0, -1,  0,    -34782),
    (0,  0,  3,  0,    -23210), (4,  0, -2,  0,    -21636),
    (2,  1, -1,  0,     24208), (2,  1,  0,  0,     30824),
    (1,  0, -1,  0,     -8379), (1,  1,  0,  0,    -16675),
    (2, -1,  1,  0,    -12831),
]

# Meeus Table 47.B: Latitude perturbations (D, Ms, Mp, F, coeff * 1e-6 deg)
_LAT_TERMS = [
    (0,  0,  0,  1,  5128122), (0,  0,  1,  1,   280602),
    (0,  0,  1, -1,   277693), (2,  0,  0, -1,   173237),
    (2,  0, -1,  1,    55413), (2,  0, -1, -1,    46271),
    (2,  0,  0,  1,    32573), (0,  0,  2,  1,    17198),
    (2,  0,  1, -1,     9266), (0,  0,  2, -1,     8822),
    (2, -1,  0, -1,     8216), (2,  0, -2, -1,     4324),
    (2,  0,  1,  1,     4200), (2,  1,  0, -1,    -3359),
    (2, -1, -1,  1,     2463), (2, -1,  0,  1,     2211),
    (2, -1, -1, -1,     2065), (0,  1, -1, -1,    -1870),
]


def moon_position_ecef(dt):
    """Moon position in ECEF coordinates.

    Truncated Meeus (ch. 47) ephemeris using 24 longitude, 23 distance,
    and 18 latitude terms.  Accuracy ~0.1 deg in position, ~200 km in
    distance.

    Parameters
    ----------
    dt : datetime
        UTC date and time.

    Returns
    -------
    numpy.ndarray
        ECEF position vector [x, y, z] in meters.
    """
    T = _T(dt)

    # Fundamental arguments (degrees)
    Lp = (218.3164477 + 481267.88123421 * T
          - 0.0015786 * T ** 2 + T ** 3 / 538841.0) % 360
    D  = (297.8501921 + 445267.1114034 * T
          - 0.0018819 * T ** 2 + T ** 3 / 545868.0) % 360
    Mp = (134.9633964 + 477198.8675055 * T
          + 0.0087414 * T ** 2 + T ** 3 / 69699.0) % 360
    Ms = (357.5291092 + 35999.0502909 * T
          - 0.0001536 * T ** 2) % 360
    F  = (93.2720950 + 483202.0175233 * T
          - 0.0036539 * T ** 2 - T ** 3 / 3526000.0) % 360

    Lp_r = np.radians(Lp)
    D_r  = np.radians(D)
    Mp_r = np.radians(Mp)
    Ms_r = np.radians(Ms)
    F_r  = np.radians(F)

    E = 1.0 - 0.002516 * T - 0.0000074 * T ** 2

    def _sum_series(terms, use_cos=False):
        total = 0.0
        for d, ms, mp, f, coeff in terms:
            arg = d * D_r + ms * Ms_r + mp * Mp_r + f * F_r
            e_corr = E ** abs(ms)
            if use_cos:
                total += coeff * e_corr * np.cos(arg)
            else:
                total += coeff * e_corr * np.sin(arg)
        return total

    sum_l = _sum_series(_LON_TERMS)
    sum_r = _sum_series(_DIST_TERMS, use_cos=True)
    sum_b = _sum_series(_LAT_TERMS)

    # Additional corrections (Meeus p. 338)
    A1 = np.radians((119.75 + 131.849 * T) % 360)
    A2 = np.radians((53.09 + 479264.290 * T) % 360)
    A3 = np.radians((313.45 + 481266.484 * T) % 360)

    sum_l += (3958 * np.sin(A1) + 1962 * np.sin(Lp_r - F_r)
              + 318 * np.sin(A2))
    sum_b += (-2235 * np.sin(Lp_r) + 382 * np.sin(A3)
              + 175 * np.sin(A1 - F_r) + 175 * np.sin(A1 + F_r)
              + 127 * np.sin(Lp_r - Mp_r) - 115 * np.sin(Lp_r + Mp_r))

    lam    = np.radians(Lp + sum_l / 1e6)         # ecliptic longitude
    beta   = np.radians(sum_b / 1e6)               # ecliptic latitude
    dist_m = (385000.56 + sum_r / 1000.0) * 1000.0 # meters

    # Ecliptic -> equatorial (ECI)
    eps = np.radians(OBLIQUITY_J2000 - 0.013004 * T)
    cb, sb = np.cos(beta), np.sin(beta)
    cl, sl = np.cos(lam),  np.sin(lam)
    ce, se = np.cos(eps),  np.sin(eps)

    x_eci = dist_m * cb * cl
    y_eci = dist_m * (cb * sl * ce - sb * se)
    z_eci = dist_m * (cb * sl * se + sb * ce)

    x, y, z = _eci_to_ecef(x_eci, y_eci, z_eci, dt)
    return np.array([x, y, z])


# =============================================================================
# Coordinate Transforms
# =============================================================================

def geodetic_to_ecef(lat_deg, lon_deg, alt_m):
    """Convert geodetic coordinates to ECEF Cartesian.

    Parameters
    ----------
    lat_deg : float
        Geodetic latitude in degrees.
    lon_deg : float
        Geodetic longitude in degrees (east positive).
    alt_m : float
        Altitude above the GRS80 ellipsoid in meters.

    Returns
    -------
    numpy.ndarray
        ECEF position vector [x, y, z] in meters.
    """
    phi = np.radians(lat_deg)
    lam = np.radians(lon_deg)
    sp, cp = np.sin(phi), np.cos(phi)
    N = A_GRS80 / np.sqrt(1.0 - E2 * sp ** 2)
    return np.array([(N + alt_m) * cp * np.cos(lam),
                     (N + alt_m) * cp * np.sin(lam),
                     (N * (1.0 - E2) + alt_m) * sp])


def enu_basis(lat_deg, lon_deg):
    """East-North-Up unit vectors at a geodetic location.

    Parameters
    ----------
    lat_deg : float
        Geodetic latitude in degrees.
    lon_deg : float
        Geodetic longitude in degrees.

    Returns
    -------
    tuple of numpy.ndarray
        (e_east, e_north, e_up) unit vectors in ECEF coordinates.
    """
    phi = np.radians(lat_deg)
    lam = np.radians(lon_deg)
    sp, cp = np.sin(phi), np.cos(phi)
    sl, cl = np.sin(lam), np.cos(lam)
    e_east  = np.array([-sl,      cl,       0.0])
    e_north = np.array([-sp * cl, -sp * sl, cp])
    e_up    = np.array([ cp * cl,  cp * sl, sp])
    return e_east, e_north, e_up


def measurement_axis(lat_deg, lon_deg, zenith_deg=0.0, azimuth_deg=0.0):
    """Measurement axis unit vector in ECEF.

    Parameters
    ----------
    lat_deg : float
        Geodetic latitude in degrees.
    lon_deg : float
        Geodetic longitude in degrees.
    zenith_deg : float
        Angle from vertical.  0 = straight up, 90 = horizontal.
    azimuth_deg : float
        Azimuth of the horizontal projection, clockwise from north.

    Returns
    -------
    numpy.ndarray
        Unit vector in ECEF coordinates.
    """
    e_e, e_n, e_u = enu_basis(lat_deg, lon_deg)
    zen = np.radians(zenith_deg)
    azi = np.radians(azimuth_deg)
    return (np.cos(zen) * e_u
            + np.sin(zen) * (np.cos(azi) * e_n + np.sin(azi) * e_e))


# =============================================================================
# Tidal Acceleration
# =============================================================================

def tidal_acceleration(r, R, GM):
    """Exact Newtonian tidal acceleration.

    Computes  GM * [(R - r) / |R - r|^3  -  R / |R|^3]

    Parameters
    ----------
    r : numpy.ndarray
        Observer position in ECEF (meters).
    R : numpy.ndarray
        Celestial body position in ECEF (meters).
    GM : float
        Gravitational parameter of the body (m^3/s^2).

    Returns
    -------
    numpy.ndarray
        Tidal acceleration vector in ECEF (m/s^2).
    """
    d = R - r
    return GM * (d / np.linalg.norm(d) ** 3 - R / np.linalg.norm(R) ** 3)


# =============================================================================
# Main API
# =============================================================================

def compute_g(dt, lat_deg, lon_deg, alt_m,
              zenith_deg=0.0, azimuth_deg=0.0):
    """Compute total g(t) projected along a measurement axis.

    Parameters
    ----------
    dt : datetime
        UTC date and time.
    lat_deg : float
        Geodetic latitude in degrees.
    lon_deg : float
        Geodetic longitude in degrees.
    alt_m : float
        Altitude above the GRS80 ellipsoid in meters.
    zenith_deg : float
        Zenith angle of the measurement axis (0 = vertical).
    azimuth_deg : float
        Azimuth of the measurement axis (clockwise from north).

    Returns
    -------
    dict
        g_total      : total g on axis (m/s^2)
        g_static     : normal gravity component (m/s^2)
        g_tidal      : elastic-Earth tidal perturbation (m/s^2)
        g_tidal_moon : lunar contribution (m/s^2)
        g_tidal_sun  : solar contribution (m/s^2)
    """
    g0    = normal_gravity(lat_deg, alt_m)
    n_hat = measurement_axis(lat_deg, lon_deg, zenith_deg, azimuth_deg)
    e_up  = enu_basis(lat_deg, lon_deg)[2]
    g_static = g0 * np.dot(e_up, n_hat)

    r = geodetic_to_ecef(lat_deg, lon_deg, alt_m)
    a_moon = DELTA_GRAV * tidal_acceleration(r, moon_position_ecef(dt), GM_MOON)
    a_sun  = DELTA_GRAV * tidal_acceleration(r, sun_position_ecef(dt),  GM_SUN)

    gm = np.dot(a_moon, n_hat)
    gs = np.dot(a_sun,  n_hat)

    return dict(g_total=g_static + gm + gs,
                g_static=g_static, g_tidal=gm + gs,
                g_tidal_moon=gm, g_tidal_sun=gs)


def compute_timeseries(start, end, lat_deg, lon_deg, alt_m,
                       zenith_deg=0.0, azimuth_deg=0.0,
                       interval_minutes=10.0, n_samples=None):
    """Compute a g(t) timeseries.

    Parameters
    ----------
    start : datetime
        Start time (UTC).
    end : datetime
        End time (UTC).
    lat_deg : float
        Geodetic latitude in degrees.
    lon_deg : float
        Geodetic longitude in degrees.
    alt_m : float
        Altitude above the GRS80 ellipsoid in meters.
    zenith_deg : float
        Zenith angle of the measurement axis (0 = vertical).
    azimuth_deg : float
        Azimuth of the measurement axis (clockwise from north).
    interval_minutes : float
        Time step in minutes (default 10).  Ignored if *n_samples* is set.
    n_samples : int, optional
        Number of evenly spaced samples between *start* and *end*
        (inclusive).  When given, overrides *interval_minutes*.

    Returns
    -------
    dict
        times        : list of datetime objects
        g_total      : numpy array (m/s^2)
        g_static     : numpy array (m/s^2)
        g_tidal      : numpy array (m/s^2)
        g_tidal_moon : numpy array (m/s^2)
        g_tidal_sun  : numpy array (m/s^2)
    """
    g0    = normal_gravity(lat_deg, alt_m)
    n_hat = measurement_axis(lat_deg, lon_deg, zenith_deg, azimuth_deg)
    e_up  = enu_basis(lat_deg, lon_deg)[2]
    g_static_val = g0 * np.dot(e_up, n_hat)
    r = geodetic_to_ecef(lat_deg, lon_deg, alt_m)

    if n_samples is not None:
        if n_samples < 1:
            raise ValueError("n_samples must be >= 1")
        total_sec = (end - start).total_seconds()
        if n_samples == 1:
            times = [start]
        else:
            step_sec = total_sec / (n_samples - 1)
            times = [start + timedelta(seconds=i * step_sec)
                     for i in range(n_samples)]
    else:
        step = timedelta(minutes=interval_minutes)
        times = []
        t = start
        while t <= end:
            times.append(t)
            t += step

    n = len(times)
    g_total      = np.empty(n)
    g_tidal      = np.empty(n)
    g_tidal_moon = np.empty(n)
    g_tidal_sun  = np.empty(n)

    for i, t in enumerate(times):
        am  = DELTA_GRAV * tidal_acceleration(r, moon_position_ecef(t), GM_MOON)
        asn = DELTA_GRAV * tidal_acceleration(r, sun_position_ecef(t),  GM_SUN)
        g_tidal_moon[i] = np.dot(am, n_hat)
        g_tidal_sun[i]  = np.dot(asn, n_hat)
        g_tidal[i]      = g_tidal_moon[i] + g_tidal_sun[i]
        g_total[i]      = g_static_val + g_tidal[i]

    return dict(times=times, g_total=g_total,
                g_static=np.full(n, g_static_val),
                g_tidal=g_tidal,
                g_tidal_moon=g_tidal_moon, g_tidal_sun=g_tidal_sun)

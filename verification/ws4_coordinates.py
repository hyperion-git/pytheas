"""
WS4: Verification of coordinate transforms and GMST in Pytheas.

Tests:
  1. Geodetic-to-ECEF at known stations (quadrant, radius, order-of-magnitude)
  2. GMST at 12+ epochs vs independent ERA computation
  3. ENU basis orthogonality and round-trip consistency
  4. GMST formula coefficient verification against Meeus eq. 12.4

Run:
  /home/xeal/.local/bin/micromamba run -n sci-base python verification/ws4_coordinates.py
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
from datetime import datetime
from pytheas._core import (
    geodetic_to_ecef, enu_basis, gmst_rad, julian_date, _T,
    A_WGS84, B_WGS84
)

PASS = 0
FAIL = 0

def check(name, condition, detail=""):
    global PASS, FAIL
    if condition:
        PASS += 1
        print(f"  [PASS] {name}")
    else:
        FAIL += 1
        print(f"  [FAIL] {name}  -- {detail}")


# =====================================================================
# 1. ECEF at known geodetic stations
# =====================================================================
print("=" * 70)
print("1. GEODETIC-TO-ECEF AT KNOWN STATIONS")
print("=" * 70)

R_EARTH = 6371e3  # mean Earth radius (m)

def _ref_geodetic_to_ecef(lat_deg, lon_deg, alt_m):
    """Independent WGS84 geodetic-to-ECEF (for cross-checking Pytheas)."""
    a = 6378137.0
    f = 1.0 / 298.257223563
    e2 = 2.0 * f - f * f
    phi = np.radians(lat_deg)
    lam = np.radians(lon_deg)
    sp, cp = np.sin(phi), np.cos(phi)
    N = a / np.sqrt(1.0 - e2 * sp**2)
    x = (N + alt_m) * cp * np.cos(lam)
    y = (N + alt_m) * cp * np.sin(lam)
    z = (N * (1.0 - e2) + alt_m) * sp
    return np.array([x, y, z])


stations = [
    ("Greenwich Observatory", 51.4769, -0.0005, 50.0),
    ("Paris Observatory",     48.8361,  2.3366, 67.0),
    ("Sydney Observatory",   -33.8599, 151.2044, 40.0),
    ("Null Island",            0.0,     0.0,      0.0),
    ("North Pole",            90.0,     0.0,      0.0),
    ("South Pole",           -90.0,     0.0,      0.0),
]

for name, lat, lon, alt in stations:
    ecef = geodetic_to_ecef(lat, lon, alt)
    ecef_ref = _ref_geodetic_to_ecef(lat, lon, alt)
    r = np.linalg.norm(ecef)
    diff = np.linalg.norm(ecef - ecef_ref)

    print(f"\n  {name}: lat={lat}, lon={lon}, alt={alt}m")
    print(f"    Pytheas ECEF:    ({ecef[0]:.1f}, {ecef[1]:.1f}, {ecef[2]:.1f})")
    print(f"    Reference ECEF:  ({ecef_ref[0]:.1f}, {ecef_ref[1]:.1f}, {ecef_ref[2]:.1f})")
    print(f"    Difference:      {diff:.2e} m")
    print(f"    Radius:          {r/1e3:.3f} km")

    # Exact match against independent implementation (sub-mm)
    check(f"{name}: matches independent WGS84 (diff < 1e-6 m)",
          diff < 1e-6,
          f"diff = {diff:.2e} m")

    # Radius should be ~6350-6400 km
    check(f"{name}: radius in [6350,6400] km",
          6350e3 < r < 6400e3,
          f"r = {r/1e3:.1f} km")

    # Quadrant checks
    if lat > 0:
        check(f"{name}: Z > 0 (northern hemisphere)", ecef[2] > 0, f"Z = {ecef[2]:.0f}")
    elif lat < 0:
        check(f"{name}: Z < 0 (southern hemisphere)", ecef[2] < 0, f"Z = {ecef[2]:.0f}")
    else:
        check(f"{name}: Z ~ 0 (equator)", abs(ecef[2]) < 1.0, f"Z = {ecef[2]:.0f}")

    if 0 < lon < 180:
        check(f"{name}: Y > 0 (east of Greenwich)", ecef[1] > 0, f"Y = {ecef[1]:.0f}")
    elif abs(lon) < 1.0:
        check(f"{name}: Y near zero (near prime meridian)", abs(ecef[1]) < 1000,
              f"Y = {ecef[1]:.0f}")


# =====================================================================
# 2. GMST AT MULTIPLE EPOCHS
# =====================================================================
print("\n" + "=" * 70)
print("2. GMST AT MULTIPLE EPOCHS vs EARTH ROTATION ANGLE (ERA)")
print("=" * 70)

def era_rad(dt):
    """Earth Rotation Angle (IAU 2006).
    ERA = 2*pi*(0.7790572732640 + 1.00273781191135448 * Du)
    where Du = JD(UT1) - 2451545.0
    We approximate UT1 ~ UTC.
    """
    JD = julian_date(dt)
    Du = JD - 2451545.0
    theta = 2.0 * np.pi * (0.7790572732640 + 1.00273781191135448 * Du)
    return theta % (2.0 * np.pi)


def gmst_meeus_deg(dt):
    """Independently compute GMST using Meeus eq. 12.4 coefficients."""
    JD = julian_date(dt)
    T = (JD - 2451545.0) / 36525.0
    gmst = (280.46061837
            + 360.98564736629 * (JD - 2451545.0)
            + 0.000387933 * T**2
            - T**3 / 38710000.0) % 360.0
    return gmst


epochs = [
    datetime(2000, 1, 1, 12, 0, 0),   # J2000.0
    datetime(2005, 3, 15, 6, 0, 0),
    datetime(2010, 6, 15, 0, 0, 0),
    datetime(2012, 12, 21, 12, 0, 0),
    datetime(2015, 1, 1, 0, 0, 0),
    datetime(2017, 6, 21, 18, 0, 0),
    datetime(2020, 3, 20, 3, 30, 0),
    datetime(2022, 9, 22, 12, 0, 0),
    datetime(2024, 1, 1, 0, 0, 0),
    datetime(2024, 6, 15, 12, 0, 0),
    datetime(2025, 6, 21, 0, 0, 0),
    datetime(2025, 12, 31, 23, 59, 59),
    datetime(2030, 1, 1, 0, 0, 0),
]

print(f"\n  {'Epoch':>24s}  {'GMST(Pytheas)':>14s}  {'ERA':>14s}  {'Diff(arcsec)':>14s}")
print("  " + "-" * 70)

max_diff_arcsec = 0.0
for dt in epochs:
    gmst_pytheas = np.degrees(gmst_rad(dt))
    era_deg = np.degrees(era_rad(dt))

    # Difference: GMST - ERA = equation of the equinoxes + accumulated precession
    # For modern dates this is dominated by the precession term, typically ~few degrees
    # But the *rate* difference over short intervals should be small
    # More relevant: compare Pytheas GMST against our independent Meeus implementation
    gmst_indep = gmst_meeus_deg(dt)
    diff_arcsec = (gmst_pytheas - gmst_indep) * 3600.0
    # Wrap to [-180, 180] degrees first
    diff_arcsec = ((gmst_pytheas - gmst_indep + 180) % 360 - 180) * 3600.0

    if abs(diff_arcsec) > max_diff_arcsec:
        max_diff_arcsec = abs(diff_arcsec)

    print(f"  {str(dt):>24s}  {gmst_pytheas:14.8f}  {era_deg:14.8f}  {diff_arcsec:14.6f}")

check("GMST matches independent Meeus implementation (max diff < 0.001 arcsec)",
      max_diff_arcsec < 0.001,
      f"max diff = {max_diff_arcsec:.6f} arcsec")

# J2000.0 specific check: GMST should be ~280.46061837 deg
dt_j2000 = datetime(2000, 1, 1, 12, 0, 0)
gmst_j2000 = np.degrees(gmst_rad(dt_j2000))
print(f"\n  J2000.0 GMST: {gmst_j2000:.8f} deg (expected ~280.46061837 deg)")
check("J2000.0 GMST matches constant term",
      abs(gmst_j2000 - 280.46061837) < 0.001,
      f"got {gmst_j2000:.10f}")

# Compare GMST vs ERA: they should agree to within ~1 degree for modern dates
# (difference is equation of equinoxes + precession offset)
print("\n  GMST vs ERA comparison (checking reasonableness):")
for dt in [datetime(2000, 1, 1, 12, 0, 0), datetime(2024, 1, 1, 0, 0, 0)]:
    gmst_val = np.degrees(gmst_rad(dt))
    era_val = np.degrees(era_rad(dt))
    diff_deg = ((gmst_val - era_val + 180) % 360 - 180)
    print(f"  {str(dt):>24s}: GMST-ERA = {diff_deg:.4f} deg = {diff_deg*3600:.1f} arcsec")
    # The GMST-ERA difference should be small (within a degree or so)
    check(f"GMST-ERA at {dt.year}: |diff| < 1.5 deg",
          abs(diff_deg) < 1.5,
          f"diff = {diff_deg:.4f} deg")


# =====================================================================
# 3. ENU ROUND-TRIP TEST
# =====================================================================
print("\n" + "=" * 70)
print("3. ENU BASIS ORTHOGONALITY AND ROUND-TRIP")
print("=" * 70)

test_locs = [
    (51.4769, -0.0005, "Greenwich"),
    (48.8361, 2.3366, "Paris"),
    (-33.8599, 151.2044, "Sydney"),
    (0.0, 0.0, "Null Island"),
    (90.0, 0.0, "North Pole"),
    (-90.0, 0.0, "South Pole"),
    (45.0, 90.0, "Mid-latitude 90E"),
]

for lat, lon, label in test_locs:
    e, n, u = enu_basis(lat, lon)

    # Build rotation matrix R = [e; n; u] (rows are basis vectors)
    R = np.vstack([e, n, u])

    # Orthogonality: R^T R = I
    RtR = R.T @ R
    ortho_err = np.max(np.abs(RtR - np.eye(3)))

    # Determinant = +1 (right-handed)
    det = np.linalg.det(R)

    # Round-trip: generate random ECEF vector, project to ENU, back to ECEF
    rng = np.random.default_rng(42)
    v_ecef = rng.normal(size=3) * 1000.0  # random vector in meters
    v_enu = R @ v_ecef
    v_ecef_back = R.T @ v_enu
    roundtrip_err = np.max(np.abs(v_ecef - v_ecef_back))

    print(f"\n  {label} (lat={lat}, lon={lon}):")
    print(f"    max|R^T R - I| = {ortho_err:.2e}")
    print(f"    det(R)         = {det:.10f}")
    print(f"    round-trip err = {roundtrip_err:.2e} m")

    check(f"{label}: orthogonality (err < 1e-14)", ortho_err < 1e-14,
          f"err = {ortho_err:.2e}")
    check(f"{label}: det(R) = +1 (right-handed)", abs(det - 1.0) < 1e-14,
          f"det = {det:.15f}")
    check(f"{label}: round-trip (err < 1e-10 m)", roundtrip_err < 1e-10,
          f"err = {roundtrip_err:.2e}")

# Verify that e_up at a station points roughly radially outward
print("\n  Checking e_up alignment with radial direction:")
for lat, lon, label in [(51.4769, -0.0005, "Greenwich"), (-33.8599, 151.2044, "Sydney")]:
    _, _, u = enu_basis(lat, lon)
    ecef = geodetic_to_ecef(lat, lon, 0.0)
    r_hat = ecef / np.linalg.norm(ecef)
    # On an ellipsoid, e_up is the normal to the ellipsoid, not exactly radial
    # But the angle should be small (< 0.2 deg)
    cos_angle = np.dot(u, r_hat)
    angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1, 1)))
    print(f"  {label}: angle(e_up, radial) = {angle_deg:.4f} deg")
    check(f"{label}: e_up near radial (< 0.2 deg)", angle_deg < 0.2,
          f"angle = {angle_deg:.4f} deg")


# =====================================================================
# 4. VERIFY GMST FORMULA COEFFICIENTS (Meeus eq. 12.4)
# =====================================================================
print("\n" + "=" * 70)
print("4. GMST FORMULA COEFFICIENTS (Meeus eq. 12.4)")
print("=" * 70)

# The formula in Pytheas (lines 136-139 of _core.py):
#   gmst_deg = (280.46061837
#               + 360.98564736629 * (JD - 2451545.0)
#               + 0.000387933 * T^2
#               - T^3 / 38710000.0) % 360.0
#
# Meeus, "Astronomical Algorithms" (2nd ed), eq. 12.4:
#   theta_0 = 280.46061837 + 360.98564736629*(JD - 2451545.0)
#           + 0.000387933*T^2 - T^3/38710000.0
#
# These are the standard coefficients. Let's verify them numerically.

# Extract coefficients by probing the function at specific times
# At J2000.0 (JD=2451545.0, T=0): GMST = 280.46061837
dt0 = datetime(2000, 1, 1, 12, 0, 0)
JD0 = julian_date(dt0)
T0 = _T(dt0)
gmst0 = np.degrees(gmst_rad(dt0))

print(f"\n  J2000.0: JD = {JD0:.1f}, T = {T0:.10f}")
print(f"  GMST at J2000.0 = {gmst0:.10f} deg")

# Reference coefficients from Meeus eq. 12.4
c0 = 280.46061837
c1 = 360.98564736629
c2 = 0.000387933
c3 = 1.0 / 38710000.0

print(f"\n  Meeus eq. 12.4 coefficients:")
print(f"    c0 = {c0}")
print(f"    c1 = {c1}")
print(f"    c2 = {c2}")
print(f"    c3 = 1/{38710000.0:.1f} = {c3:.15e}")

# Verify: at J2000.0, only c0 contributes (T=0, JD-J2000=0)
check("Constant term c0 = 280.46061837",
      abs(gmst0 - c0) < 1e-6,
      f"GMST(J2000) = {gmst0}")

# Verify daily rate: compute GMST 1 sidereal day apart
# One mean solar day later: JD increases by 1, so the linear contribution
# should add 360.98564736629 deg (i.e., one full rotation + precession)
dt1 = datetime(2000, 1, 2, 12, 0, 0)
gmst1 = np.degrees(gmst_rad(dt1))
daily_advance = (gmst1 - gmst0) % 360.0 + 360.0  # should be ~360.986
# Account for wrapping
raw_diff = gmst1 - gmst0
# The actual advance includes multiple full rotations; extract just the expected value
T1 = _T(dt1)
JD1 = julian_date(dt1)
expected_advance = c1 * (JD1 - JD0) + c2 * (T1**2 - T0**2) - c3 * (T1**3 - T0**3)
actual_advance = gmst1 - gmst0
# Both may have had modulo applied; compare predicted vs actual
pred_gmst1 = (c0 + c1 * (JD1 - 2451545.0) + c2 * T1**2 - c3 * T1**3) % 360.0
print(f"\n  GMST at J2000+1d: Pytheas={gmst1:.8f}, predicted={pred_gmst1:.8f}")
check("Daily rate c1: GMST(J2000+1d) matches prediction",
      abs(gmst1 - pred_gmst1) < 1e-6,
      f"diff = {abs(gmst1 - pred_gmst1):.2e} deg")

# Verify T^2 term: compare at T=1 century (year 2100)
dt_2100 = datetime(2100, 1, 1, 12, 0, 0)
JD_2100 = julian_date(dt_2100)
T_2100 = _T(dt_2100)
gmst_2100 = np.degrees(gmst_rad(dt_2100))
pred_2100 = (c0 + c1 * (JD_2100 - 2451545.0) + c2 * T_2100**2 - c3 * T_2100**3) % 360.0
print(f"\n  GMST at 2100-01-01 12:00: Pytheas={gmst_2100:.8f}, predicted={pred_2100:.8f}")
print(f"  T = {T_2100:.6f} centuries, T^2 term = {c2 * T_2100**2:.6f} deg")
print(f"  T^3 term = {c3 * T_2100**3:.6f} deg")
check("T^2/T^3 terms: GMST at 2100 matches full formula",
      abs(gmst_2100 - pred_2100) < 1e-6,
      f"diff = {abs(gmst_2100 - pred_2100):.2e} deg")

# Cross-check: the daily rate 360.98564736629 should equal
# 360 * 1.00273781191135448 (ratio of sidereal to solar day)
# 360 * 1.00273781191135448 = 360.98561228808...
# Actually Meeus uses a slightly different rate that absorbs precession.
# Let's just verify the number is consistent.
sidereal_rate = 360.0 * 1.00273781191135448  # from ERA definition
print(f"\n  ERA-based daily rate:   {sidereal_rate:.11f} deg/day")
print(f"  Meeus GMST daily rate: {c1:.11f} deg/day")
print(f"  Difference:            {(c1 - sidereal_rate)*3600:.4f} arcsec/day")
print(f"  (This difference accumulates the precession rate ~50 arcsec/year)")
precession_arcsec_per_year = (c1 - sidereal_rate) * 3600 * 365.25
print(f"  Implied precession:    {precession_arcsec_per_year:.1f} arcsec/year (expect ~50)")
check("Precession rate from c1 is ~50 arcsec/year",
      45 < precession_arcsec_per_year < 55,
      f"got {precession_arcsec_per_year:.1f}")


# =====================================================================
# SUMMARY
# =====================================================================
print("\n" + "=" * 70)
TOTAL = PASS + FAIL
print(f"SUMMARY: {PASS}/{TOTAL} passed, {FAIL}/{TOTAL} failed")
print("=" * 70)

if FAIL > 0:
    sys.exit(1)

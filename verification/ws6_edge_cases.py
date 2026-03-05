"""
WS6: Numerical Stability and Edge Case Tests
=============================================
Tests extreme inputs for NaN/Inf, symmetry, continuity, and wrapping.

Sections:
  1. Near-polar latitudes
  2. Extreme altitudes
  3. Date range (Meeus ephemeris limits)
  4. Longitude wrapping
  5. Latitude symmetry
  6. Timeseries continuity
  7. Comprehensive NaN/Inf sweep
"""

import sys
import numpy as np
from datetime import datetime, timedelta

sys.path.insert(0, "/home/xeal/dev/pytheas")
from pytheas import (
    compute_g, compute_timeseries, normal_gravity,
    geodetic_to_ecef, enu_basis, moon_position_ecef, sun_position_ecef,
    julian_date, gmst_rad, LabFrame,
)

PASS = "PASS"
FAIL = "FAIL"
WARN = "WARN"

n_pass = 0
n_fail = 0
n_warn = 0


def section(title):
    print(f"\n{'=' * 72}")
    print(f"  {title}")
    print(f"{'=' * 72}\n")


def check(name, condition, detail=""):
    global n_pass, n_fail
    status = PASS if condition else FAIL
    if condition:
        n_pass += 1
    else:
        n_fail += 1
    tag = f"  [{status}] {name}"
    if detail:
        tag += f"  -- {detail}"
    print(tag)
    return condition


def warn(name, detail=""):
    global n_warn
    n_warn += 1
    tag = f"  [{WARN}] {name}"
    if detail:
        tag += f"  -- {detail}"
    print(tag)


def is_finite(x):
    """Check finiteness for scalars, arrays, or dataclass results."""
    if isinstance(x, np.ndarray):
        return np.all(np.isfinite(x))
    try:
        return np.isfinite(float(x))
    except (TypeError, ValueError):
        return False


# =========================================================================
# 1. Near-polar latitudes
# =========================================================================
section("1. Near-Polar Latitudes")

ref_dt = datetime(2025, 3, 20, 12, 0, 0)
polar_lats = [89.999, 89.9999, 90.0, 90.0001, -89.999, -90.0]

for lat in polar_lats:
    label = f"lat={lat:+.4f}"

    # normal_gravity
    ng = normal_gravity(lat, 0.0)
    check(f"{label}: normal_gravity finite", is_finite(ng), f"g={ng:.8f}")
    check(f"{label}: normal_gravity in [9.78, 9.84]", 9.78 <= ng <= 9.84,
          f"g={ng:.8f}")

    # geodetic_to_ecef
    pos = geodetic_to_ecef(lat, 0.0, 0.0)
    check(f"{label}: geodetic_to_ecef finite", is_finite(pos),
          f"|r|={np.linalg.norm(pos):.1f} m")

    # enu_basis orthonormality
    e, n, u = enu_basis(lat, 0.0)
    ortho_ok = True
    rh_ok = True
    for v in (e, n, u):
        if abs(np.linalg.norm(v) - 1.0) > 1e-12:
            ortho_ok = False
    if abs(np.dot(e, n)) > 1e-12 or abs(np.dot(e, u)) > 1e-12 or abs(np.dot(n, u)) > 1e-12:
        ortho_ok = False
    cross = np.cross(e, n)
    if np.linalg.norm(cross - u) > 1e-10:
        rh_ok = False
    check(f"{label}: enu_basis orthonormal", ortho_ok)
    check(f"{label}: enu_basis right-handed", rh_ok,
          f"|ExN - U| = {np.linalg.norm(cross - u):.2e}")

    # compute_g should not crash or produce NaN
    try:
        res = compute_g(ref_dt, lat, 0.0, 0.0)
        check(f"{label}: compute_g finite", is_finite(res.g_total),
              f"g_total={res.g_total:.8f}")
    except Exception as exc:
        check(f"{label}: compute_g no exception", False, str(exc))


# =========================================================================
# 2. Extreme altitudes
# =========================================================================
section("2. Extreme Altitudes")

alt_cases = [
    (-500, "underground"),
    (0, "sea level"),
    (100, "100 m"),
    (1000, "1 km"),
    (8849, "Everest"),
    (35_786_000, "GEO orbit"),
    (384_400_000, "Moon distance"),
]

prev_g = None
for alt, desc in alt_cases:
    ng = normal_gravity(45.0, alt)
    check(f"alt={desc}: normal_gravity finite", is_finite(ng), f"g={ng:.6f}")

    # Monotonic decrease for moderate altitudes
    if prev_g is not None and alt <= 8849:
        check(f"alt={desc}: g decreases with altitude", ng < prev_g,
              f"g={ng:.8f} vs prev={prev_g:.8f}")

    if alt <= 8849:
        prev_g = ng

    # compute_g should not crash
    try:
        res = compute_g(ref_dt, 45.0, 0.0, alt)
        check(f"alt={desc}: compute_g finite", is_finite(res.g_total),
              f"g_total={res.g_total:.6f}")
    except Exception as exc:
        check(f"alt={desc}: compute_g no exception", False, str(exc))

# Check extreme altitude behavior
for alt, desc in alt_cases:
    if alt > 8849:
        ng = normal_gravity(45.0, alt)
        if ng < 0:
            warn(f"alt={desc}: negative normal_gravity (free-air breakdown)",
                 f"g={ng:.6f}")
        elif ng > 20:
            warn(f"alt={desc}: unreasonably large normal_gravity",
                 f"g={ng:.6f}")


# =========================================================================
# 3. Date range
# =========================================================================
section("3. Date Range (Meeus Ephemeris)")

date_cases = [
    datetime(1900, 1, 1, 0, 0, 0),
    datetime(1950, 6, 15, 12, 0, 0),
    datetime(2000, 1, 1, 12, 0, 0),
    datetime(2025, 3, 20, 0, 0, 0),
    datetime(2050, 1, 1, 0, 0, 0),
    datetime(2100, 1, 1, 0, 0, 0),
    datetime(2200, 1, 1, 0, 0, 0),
]

for dt in date_cases:
    label = dt.strftime("%Y-%m-%d")

    # julian_date
    jd = julian_date(dt)
    check(f"{label}: julian_date finite", is_finite(jd), f"JD={jd:.4f}")
    # J2000 = 2451545.0, range sanity
    check(f"{label}: julian_date reasonable",
          2_000_000 < jd < 2_600_000, f"JD={jd:.4f}")

    # moon_position_ecef
    try:
        moon = moon_position_ecef(dt)
        moon_dist = np.linalg.norm(moon)
        check(f"{label}: moon position finite", is_finite(moon),
              f"|r|={moon_dist/1e6:.0f} km")
        # Moon distance should be roughly 356,000-407,000 km
        check(f"{label}: moon distance plausible",
              300e6 < moon_dist < 450e6,
              f"{moon_dist/1e6:.0f} km")
    except Exception as exc:
        check(f"{label}: moon_position_ecef no exception", False, str(exc))

    # sun_position_ecef
    try:
        sun = sun_position_ecef(dt)
        sun_dist = np.linalg.norm(sun)
        check(f"{label}: sun position finite", is_finite(sun),
              f"|r|={sun_dist/1e9:.1f} Gm")
        # Sun distance ~147-152 Gm
        check(f"{label}: sun distance plausible",
              140e9 < sun_dist < 160e9,
              f"{sun_dist/1e9:.1f} Gm")
    except Exception as exc:
        check(f"{label}: sun_position_ecef no exception", False, str(exc))

    # compute_g
    try:
        res = compute_g(dt, 48.0, 11.0, 500.0)
        check(f"{label}: compute_g finite", is_finite(res.g_total),
              f"g_total={res.g_total:.8f}")
    except Exception as exc:
        check(f"{label}: compute_g no exception", False, str(exc))


# =========================================================================
# 4. Longitude wrapping
# =========================================================================
section("4. Longitude Wrapping")

lon_pairs = [
    (0.0, 360.0, "lon=0 vs lon=360"),
    (-90.0, 270.0, "lon=-90 vs lon=270"),
    (-180.0, 180.0, "lon=-180 vs lon=180"),
]

for lon_a, lon_b, desc in lon_pairs:
    # geodetic_to_ecef
    pos_a = geodetic_to_ecef(45.0, lon_a, 0.0)
    pos_b = geodetic_to_ecef(45.0, lon_b, 0.0)
    diff = np.linalg.norm(pos_a - pos_b)
    check(f"{desc}: geodetic_to_ecef match", diff < 0.01,
          f"diff={diff:.2e} m")

    # compute_g
    res_a = compute_g(ref_dt, 45.0, lon_a, 0.0)
    res_b = compute_g(ref_dt, 45.0, lon_b, 0.0)
    g_diff = abs(res_a.g_total - res_b.g_total)
    check(f"{desc}: compute_g match", g_diff < 1e-12,
          f"diff={g_diff:.2e} m/s^2")

# All listed longitudes should be finite
all_lons = [-180, -90, 0, 90, 180, 270, 360]
for lon in all_lons:
    res = compute_g(ref_dt, 45.0, lon, 0.0)
    check(f"lon={lon}: compute_g finite", is_finite(res.g_total),
          f"g_total={res.g_total:.8f}")


# =========================================================================
# 5. Latitude symmetry
# =========================================================================
section("5. Latitude N/S Symmetry")

sym_lats = [0, 15, 30, 45, 60, 75, 90]

for lat in sym_lats:
    gn = normal_gravity(lat, 0.0)
    gs = normal_gravity(-lat, 0.0)
    check(f"lat=+/-{lat}: normal_gravity symmetric",
          abs(gn - gs) < 1e-15,
          f"N={gn:.10f}, S={gs:.10f}, diff={abs(gn-gs):.2e}")


# =========================================================================
# 6. Timeseries continuity
# =========================================================================
section("6. Timeseries Continuity (3 days, 1-min steps)")

ts_start = datetime(2025, 3, 20, 0, 0, 0)
ts_end = datetime(2025, 3, 23, 0, 0, 0)

ts = compute_timeseries(ts_start, ts_end, 48.0, 11.0, 500.0,
                        interval_minutes=1.0)

g_arr = ts.g_total
n_samples = len(g_arr)
print(f"  Samples: {n_samples}")
print(f"  g range: [{np.min(g_arr):.8f}, {np.max(g_arr):.8f}] m/s^2")

# Check for NaN/Inf
nan_count = np.sum(~np.isfinite(g_arr))
check("timeseries: no NaN/Inf", nan_count == 0,
      f"{nan_count} non-finite values out of {n_samples}")

# Check continuity: max step-to-step change
dg = np.abs(np.diff(g_arr))
max_dg = np.max(dg)

max_dg_nGal = max_dg / 1e-11
print(f"  Max |dg/dt| = {max_dg_nGal:.2f} nGal/min")

# Physical bound: peak-to-peak ~150,000 nGal with ~12.4-hour dominant period
# gives max rate ~ pi * 150000 / (12.4*60) ~ 633 nGal/min.
# Use 700 nGal/min as a generous threshold for genuine discontinuities.
threshold = 700e-11  # 700 nGal in m/s^2

check("timeseries: continuity (max step < 700 nGal/min)",
      max_dg < threshold,
      f"max_dg = {max_dg_nGal:.2f} nGal/min")

# Also check tidal components
for name, arr in [("g_tidal", ts.g_tidal),
                  ("g_tidal_moon", ts.g_tidal_moon),
                  ("g_tidal_sun", ts.g_tidal_sun)]:
    nan_count = np.sum(~np.isfinite(arr))
    check(f"timeseries {name}: no NaN/Inf", nan_count == 0,
          f"{nan_count} non-finite")


# =========================================================================
# 7. Comprehensive NaN/Inf sweep
# =========================================================================
section("7. Comprehensive NaN/Inf Sweep")

sweep_lats = [-90, -89.999, -45, 0, 45, 89.999, 90]
sweep_lons = [-180, 0, 180, 360]
sweep_alts = [-500, 0, 8849, 35_786_000]

n_combo = 0
n_finite = 0
failures = []

for lat in sweep_lats:
    for lon in sweep_lons:
        for alt in sweep_alts:
            n_combo += 1
            label = f"lat={lat}, lon={lon}, alt={alt}"

            # normal_gravity
            ng = normal_gravity(lat, alt)
            if not is_finite(ng):
                failures.append(f"normal_gravity({label}): {ng}")

            # geodetic_to_ecef
            pos = geodetic_to_ecef(lat, lon, alt)
            if not is_finite(pos):
                failures.append(f"geodetic_to_ecef({label}): {pos}")

            # enu_basis
            e, n_vec, u = enu_basis(lat, lon)
            for vname, v in [("e", e), ("n", n_vec), ("u", u)]:
                if not is_finite(v):
                    failures.append(f"enu_basis.{vname}({label}): {v}")

            # compute_g (skip extreme alts to save time)
            if alt <= 8849:
                try:
                    res = compute_g(ref_dt, lat, lon, alt)
                    if not is_finite(res.g_total):
                        failures.append(f"compute_g({label}): g_total={res.g_total}")
                    else:
                        n_finite += 1
                except Exception as exc:
                    failures.append(f"compute_g({label}): EXCEPTION {exc}")
            else:
                # For extreme alts, just test that it doesn't crash
                try:
                    res = compute_g(ref_dt, lat, lon, alt)
                    if is_finite(res.g_total):
                        n_finite += 1
                    else:
                        warn(f"compute_g({label}): non-finite g_total={res.g_total}")
                except Exception as exc:
                    failures.append(f"compute_g({label}): EXCEPTION {exc}")

print(f"  Tested {n_combo} combinations")
print(f"  Finite compute_g results: {n_finite}")

if failures:
    check("NaN/Inf sweep: all outputs finite", False,
          f"{len(failures)} failures")
    for f in failures:
        print(f"    - {f}")
else:
    check("NaN/Inf sweep: all outputs finite", True,
          f"{n_combo} combos, zero failures")


# =========================================================================
# 8. LabFrame at edge cases
# =========================================================================
section("8. LabFrame Edge Cases")

edge_labs = [
    (90.0, 0.0, 0.0, "North Pole"),
    (-90.0, 0.0, 0.0, "South Pole"),
    (0.0, 0.0, 0.0, "Equator 0 lon"),
    (0.0, 180.0, 0.0, "Equator dateline"),
    (45.0, 0.0, 8849.0, "Mid-lat Everest alt"),
]

for lat, lon, alt, desc in edge_labs:
    try:
        lab = LabFrame(lat, lon, alt)
        fld = lab.field(ref_dt)
        check(f"LabFrame {desc}: g finite", is_finite(fld.g),
              f"|g|={np.linalg.norm(fld.g):.6f}")
        check(f"LabFrame {desc}: T finite", is_finite(fld.T),
              f"diag(T)={np.diag(fld.T)}")
        check(f"LabFrame {desc}: omega finite", is_finite(fld.omega))

        # eom should not crash
        a = fld.eom([0, 0, 0], [0, 0, 0])
        check(f"LabFrame {desc}: eom finite", is_finite(a))

    except Exception as exc:
        check(f"LabFrame {desc}: no exception", False, str(exc))


# =========================================================================
# Summary
# =========================================================================
section("SUMMARY")
total = n_pass + n_fail
print(f"  {n_pass}/{total} checks passed")
print(f"  {n_fail}/{total} checks FAILED")
print(f"  {n_warn} warnings")

if n_fail > 0:
    print(f"\n  *** {n_fail} FAILURES detected ***")
    sys.exit(1)
else:
    print("\n  All edge-case checks passed.")
    sys.exit(0)

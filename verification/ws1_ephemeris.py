#!/usr/bin/env python
"""
WS1: Pytheas ephemeris verification against JPL Horizons DE441.

Compares moon_position_ecef() and sun_position_ecef() against JPL Horizons
geocentric ICRF vector positions.  Since ECEF is a z-rotation of ECI,
frame-independent quantities are:
  - Geocentric distance |r|
  - Z-component (declination proxy, preserved by z-rotation)
  - Declination angle: arcsin(z/|r|)

Attempts live fetch from JPL Horizons API; falls back to hardcoded DE441
reference values on failure.
"""

import sys
import json
import ssl
import urllib.request
import urllib.parse
import numpy as np
from datetime import datetime, timedelta
from pathlib import Path

# -- Add pytheas to path --
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from pytheas._core import moon_position_ecef, sun_position_ecef

# For plots
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter


# ============================================================================
# JPL Horizons API fetch
# ============================================================================

def _horizons_ssl_ctx():
    """SSL context that tolerates corporate/WSL certificate issues."""
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    return ctx


def fetch_horizons(command, start="2024-01-01", stop="2025-12-31",
                   step="1 MONTHS"):
    """Fetch geocentric ICRF vector positions from JPL Horizons.

    Returns list of (datetime, x_km, y_km, z_km) or None on failure.
    """
    params = {
        "format": "json",
        "COMMAND": f"'{command}'",
        "CENTER": "'500@399'",
        "EPHEM_TYPE": "VECTORS",
        "OBJ_DATA": "NO",
        "VEC_TABLE": "2",
        "REF_SYSTEM": "ICRF",
        "REF_PLANE": "FRAME",
        "START_TIME": f"'{start}'",
        "STOP_TIME": f"'{stop}'",
        "STEP_SIZE": f"'{step}'",
        "VEC_LABELS": "NO",
    }
    url = ("https://ssd.jpl.nasa.gov/api/horizons.api?"
           + urllib.parse.urlencode(params, quote_via=urllib.parse.quote))
    try:
        ctx = _horizons_ssl_ctx()
        resp = urllib.request.urlopen(url, timeout=30, context=ctx)
        data = json.loads(resp.read())
    except Exception as e:
        print(f"  [WARN] Horizons API error: {e}")
        return None

    if "result" not in data:
        print(f"  [WARN] No 'result' in response: {list(data.keys())}")
        return None

    result = data["result"]
    soe = result.find("$$SOE")
    eoe = result.find("$$EOE")
    if soe < 0 or eoe < 0:
        print("  [WARN] No $$SOE/$$EOE markers in response")
        return None

    block = result[soe + 5:eoe].strip()
    entries = []
    lines = block.split("\n")
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        if "A.D." in line:
            # Parse: "2460310.500000000 = A.D. 2024-Jan-01 00:00:00.0000 TDB"
            after_ad = line.split("A.D.")[1].strip()
            parts = after_ad.split()
            date_str = parts[0]       # "2024-Jan-01"
            time_str = parts[1][:8]   # "00:00:00"
            try:
                dt = datetime.strptime(date_str + " " + time_str,
                                       "%Y-%b-%d %H:%M:%S")
            except ValueError:
                i += 1
                continue
            # Next line: X Y Z (km)
            i += 1
            if i < len(lines):
                vals = lines[i].strip().split()
                if len(vals) >= 3:
                    x = float(vals[0])
                    y = float(vals[1])
                    z = float(vals[2])
                    entries.append((dt, x, y, z))
        i += 1

    if not entries:
        print("  [WARN] Parsed 0 entries from Horizons response")
        return None

    print(f"  [OK] Fetched {len(entries)} epochs from Horizons")
    return entries


# ============================================================================
# Hardcoded DE441 reference data (from verified JPL Horizons queries)
# ============================================================================
# Geocentric ICRF (J2000 equatorial), units: km
# Format: (datetime, X, Y, Z)

MOON_HARDCODED = [
    (datetime(2024, 1, 1, 0, 0, 0),   -3.679525311419403E+05,  1.427749731429682E+05,  8.934228118085049E+04),
    (datetime(2024, 2, 1, 0, 0, 0),   -3.794441905598769E+05, -1.217969281133812E+05, -5.404168732891611E+04),
    (datetime(2024, 3, 1, 0, 0, 0),   -3.047794061123109E+05, -2.297847597615251E+05, -1.160346733941420E+05),
    (datetime(2024, 4, 1, 0, 0, 0),   -1.729808728858708E+04, -3.380369459332119E+05, -1.833989312917340E+05),
    (datetime(2024, 5, 1, 0, 0, 0),    2.124949070670388E+05, -2.695873807499340E+05, -1.520821685021791E+05),
    (datetime(2024, 6, 1, 0, 0, 0),    3.683053100782324E+05, -1.057009909977951E+04, -1.456014827837195E+04),
    (datetime(2024, 7, 1, 0, 0, 0),    2.987660891339004E+05,  1.982135583152877E+05,  1.013051731959705E+05),
    (datetime(2024, 8, 1, 0, 0, 0),    1.454325221888769E+04,  3.385090214490337E+05,  1.832956760255061E+05),
    (datetime(2024, 9, 1, 0, 0, 0),   -2.842620959102936E+05,  2.445828706699075E+05,  1.368053910670627E+05),
    (datetime(2024, 10, 1, 0, 0, 0),  -3.973452891771129E+05,  6.860925730142037E+04,  4.215473811211392E+04),
    (datetime(2024, 11, 1, 0, 0, 0),  -3.395751117227855E+05, -1.947459810487748E+05, -1.025087892396738E+05),
    (datetime(2024, 12, 1, 0, 0, 0),  -1.617063866868304E+05, -3.193096051662483E+05, -1.719635790319248E+05),
    (datetime(2025, 1, 1, 0, 0, 0),    1.520523605713538E+05, -3.078236314310748E+05, -1.668798858747519E+05),
    (datetime(2025, 2, 1, 0, 0, 0),    3.547881170375207E+05, -8.666465192598013E+04, -4.583670657885905E+04),
    (datetime(2025, 3, 1, 0, 0, 0),    3.602941433419737E+05, -3.651416391512818E+04, -1.805044980283984E+04),
    (datetime(2025, 4, 1, 0, 0, 0),    2.454482148126173E+05,  2.315919964997484E+05,  1.281359663481913E+05),
    (datetime(2025, 5, 1, 0, 0, 0),    3.054550482385845E+04,  3.210665129518766E+05,  1.750773889393343E+05),
    (datetime(2025, 6, 1, 0, 0, 0),   -2.715799234455151E+05,  2.396883531917965E+05,  1.268727794582563E+05),
    (datetime(2025, 7, 1, 0, 0, 0),   -3.869569259184093E+05,  7.007467527785245E+04,  3.167086475374909E+04),
    (datetime(2025, 8, 1, 0, 0, 0),   -3.383055393434293E+05, -1.906013419981420E+05, -1.103757034534148E+05),
    (datetime(2025, 9, 1, 0, 0, 0),   -9.285930384322161E+04, -3.416982409293890E+05, -1.883065658881358E+05),
    (datetime(2025, 10, 1, 0, 0, 0),   1.286548406946826E+05, -3.276488648101639E+05, -1.756072374609563E+05),
    (datetime(2025, 11, 1, 0, 0, 0),   3.454182341515119E+05, -1.354335235830779E+05, -6.521909375703422E+04),
    (datetime(2025, 12, 1, 0, 0, 0),   3.575453219118375E+05,  6.998882397408022E+04,  4.746851033844583E+04),
]

SUN_HARDCODED = [
    (datetime(2024, 1, 1, 0, 0, 0),    2.481099325965390E+07, -1.330334521311550E+08, -5.766810624017932E+07),
    (datetime(2024, 2, 1, 0, 0, 0),    9.721451787005229E+07, -1.016422084200529E+08, -4.406122075247613E+07),
    (datetime(2024, 3, 1, 0, 0, 0),    1.397742662283037E+08, -4.526662373785767E+07, -1.962358943210541E+07),
    (datetime(2024, 4, 1, 0, 0, 0),    1.465255382642640E+08,  2.715671643902075E+07,  1.177068710573139E+07),
    (datetime(2024, 5, 1, 0, 0, 0),    1.141340593681229E+08,  9.032649341896772E+07,  3.915416940941525E+07),
    (datetime(2024, 6, 1, 0, 0, 0),    5.022212022184425E+07,  1.313326759801844E+08,  5.693049608345818E+07),
    (datetime(2024, 7, 1, 0, 0, 0),   -2.464397118193815E+07,  1.377008675035678E+08,  5.969178077550165E+07),
    (datetime(2024, 8, 1, 0, 0, 0),   -9.536966553407273E+07,  1.084065060025865E+08,  4.699376311594316E+07),
    (datetime(2024, 9, 1, 0, 0, 0),   -1.406846254024134E+08,  5.028350738382646E+07,  2.179844838325282E+07),
    (datetime(2024, 10, 1, 0, 0, 0),  -1.483359738689474E+08, -1.904939774848031E+07, -8.256583874035988E+06),
    (datetime(2024, 11, 1, 0, 0, 0),  -1.158586058346020E+08, -8.520880312977119E+07, -3.693633446405344E+07),
    (datetime(2024, 12, 1, 0, 0, 0),  -5.301698654340477E+07, -1.263075472631757E+08, -5.475264045732966E+07),
    (datetime(2025, 1, 1, 0, 0, 0),    2.673066229892559E+07, -1.327246809685359E+08, -5.753486058267201E+07),
    (datetime(2025, 2, 1, 0, 0, 0),    9.868437168826628E+07, -1.004703354482892E+08, -4.355317480190057E+07),
    (datetime(2025, 3, 1, 0, 0, 0),    1.395523282499756E+08, -4.584127077907413E+07, -1.987242172913667E+07),
    (datetime(2025, 4, 1, 0, 0, 0),    1.466517885910457E+08,  2.656135948804406E+07,  1.151331514036667E+07),
    (datetime(2025, 5, 1, 0, 0, 0),    1.145577272135454E+08,  8.986568053010082E+07,  3.895493986807102E+07),
    (datetime(2025, 6, 1, 0, 0, 0),    5.083222443804654E+07,  1.311267381302641E+08,  5.684110688594246E+07),
    (datetime(2025, 7, 1, 0, 0, 0),   -2.400825742134599E+07,  1.377890812256956E+08,  5.972930041537907E+07),
    (datetime(2025, 8, 1, 0, 0, 0),   -9.487022255787645E+07,  1.087696726283600E+08,  4.714996401921754E+07),
    (datetime(2025, 9, 1, 0, 0, 0),   -1.404523085048895E+08,  5.082643798664138E+07,  2.203270185244682E+07),
    (datetime(2025, 10, 1, 0, 0, 0),  -1.484277473129742E+08, -1.846969715139007E+07, -8.005826458529575E+06),
    (datetime(2025, 11, 1, 0, 0, 0),  -1.162674360126390E+08, -8.474778723513289E+07, -3.673606153111023E+07),
    (datetime(2025, 12, 1, 0, 0, 0),  -5.363216352801231E+07, -1.260926872018423E+08, -5.465841350520901E+07),
]


# ============================================================================
# Comparison engine
# ============================================================================

def compute_ref_invariants(ref_data):
    """Extract frame-independent quantities from ICRF XYZ data.

    Returns list of (datetime, distance_km, z_km).
    """
    out = []
    for entry in ref_data:
        dt, x, y, z = entry
        dist = np.sqrt(x**2 + y**2 + z**2)
        out.append((dt, dist, z))
    return out


def compare_body(name, ref_invariants, position_fn, dist_tol_km,
                 angle_tol_deg):
    """Compare Pytheas positions against reference invariants.

    Parameters
    ----------
    ref_invariants : list of (datetime, distance_km, z_km)
    """
    results = []
    for dt, dist_ref, z_ref in ref_invariants:
        pos = position_fn(dt)
        dist_py = np.linalg.norm(pos) / 1e3   # km
        z_py = pos[2] / 1e3                    # km

        dist_err = dist_py - dist_ref
        z_err = z_py - z_ref

        # Declination: arcsin(z/r) -- frame-independent
        dec_ref = np.degrees(np.arcsin(np.clip(z_ref / dist_ref, -1, 1)))
        dec_py = np.degrees(np.arcsin(np.clip(z_py / dist_py, -1, 1)))
        dec_err = dec_py - dec_ref

        results.append({
            "dt": dt,
            "dist_ref": dist_ref,
            "dist_py": dist_py,
            "dist_err": dist_err,
            "z_ref": z_ref,
            "z_py": z_py,
            "z_err": z_err,
            "dec_ref": dec_ref,
            "dec_py": dec_py,
            "dec_err": dec_err,
        })

    dist_errs = np.array([r["dist_err"] for r in results])
    z_errs = np.array([r["z_err"] for r in results])
    dec_errs = np.array([r["dec_err"] for r in results])

    stats = {
        "name": name,
        "n_epochs": len(results),
        "dist_rms_km": float(np.sqrt(np.mean(dist_errs**2))),
        "dist_max_km": float(np.max(np.abs(dist_errs))),
        "dist_mean_km": float(np.mean(dist_errs)),
        "z_rms_km": float(np.sqrt(np.mean(z_errs**2))),
        "z_max_km": float(np.max(np.abs(z_errs))),
        "z_mean_km": float(np.mean(z_errs)),
        "dec_rms_deg": float(np.sqrt(np.mean(dec_errs**2))),
        "dec_max_deg": float(np.max(np.abs(dec_errs))),
        "dec_mean_deg": float(np.mean(dec_errs)),
        "dist_tol_km": dist_tol_km,
        "angle_tol_deg": angle_tol_deg,
    }

    stats["dist_pass"] = stats["dist_max_km"] < dist_tol_km
    stats["angle_pass"] = stats["dec_max_deg"] < angle_tol_deg
    stats["overall_pass"] = stats["dist_pass"] and stats["angle_pass"]

    return results, stats


def print_report(results, stats):
    """Print formatted comparison report."""
    name = stats["name"]
    print(f"\n{'='*76}")
    print(f"  {name} Ephemeris Verification ({stats['n_epochs']} epochs)")
    print(f"{'='*76}")

    print(f"\n  {'Date':20s} {'Dist Ref (km)':>15s} {'Dist Err (km)':>14s} "
          f"{'Z Err (km)':>12s} {'Dec Err (deg)':>14s}")
    print(f"  {'-'*20} {'-'*15} {'-'*14} {'-'*12} {'-'*14}")
    for r in results:
        dt_str = r["dt"].strftime("%Y-%m-%d %H:%M")
        print(f"  {dt_str:20s} {r['dist_ref']:15.1f} {r['dist_err']:+14.2f} "
              f"{r['z_err']:+12.2f} {r['dec_err']:+14.4f}")

    print(f"\n  Summary Statistics:")
    print(f"    Distance RMS:  {stats['dist_rms_km']:10.2f} km")
    print(f"    Distance Max:  {stats['dist_max_km']:10.2f} km  "
          f"(tol: {stats['dist_tol_km']:.0f} km)  "
          f"{'PASS' if stats['dist_pass'] else 'FAIL'}")
    print(f"    Distance Bias: {stats['dist_mean_km']:+10.2f} km")
    print(f"    Z-comp RMS:    {stats['z_rms_km']:10.2f} km")
    print(f"    Z-comp Max:    {stats['z_max_km']:10.2f} km")
    print(f"    Dec RMS:       {stats['dec_rms_deg']:10.4f} deg "
          f"({stats['dec_rms_deg']*60:.2f} arcmin)")
    print(f"    Dec Max:       {stats['dec_max_deg']:10.4f} deg "
          f"({stats['dec_max_deg']*60:.2f} arcmin)  "
          f"(tol: {stats['angle_tol_deg']:.4f} deg)  "
          f"{'PASS' if stats['angle_pass'] else 'FAIL'}")
    print(f"    Dec Bias:      {stats['dec_mean_deg']:+10.4f} deg "
          f"({stats['dec_mean_deg']*60:+.2f} arcmin)")
    print(f"\n    Overall: {'PASS' if stats['overall_pass'] else 'FAIL'}")


def plot_errors(moon_results, sun_results, moon_stats, sun_stats):
    """Generate verification error plots."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("Pytheas Ephemeris vs JPL Horizons DE441", fontsize=14,
                 fontweight="bold")

    # -- Moon distance error --
    ax = axes[0, 0]
    dates = [r["dt"] for r in moon_results]
    errs = [r["dist_err"] for r in moon_results]
    ax.plot(dates, errs, "o-", color="steelblue", markersize=4)
    ax.axhline(0, color="gray", linewidth=0.5)
    ax.axhline(+200, color="red", linewidth=0.5, linestyle="--", alpha=0.5,
               label="200 km tolerance")
    ax.axhline(-200, color="red", linewidth=0.5, linestyle="--", alpha=0.5)
    ax.set_ylabel("Distance Error (km)")
    ax.set_title(f"Moon Distance Error (max={moon_stats['dist_max_km']:.0f} km)")
    ax.legend(fontsize=8)
    ax.xaxis.set_major_formatter(DateFormatter("%Y-%m"))
    ax.tick_params(axis="x", rotation=45)

    # -- Moon declination error --
    ax = axes[0, 1]
    dec_errs = [r["dec_err"] * 60 for r in moon_results]  # arcmin
    ax.plot(dates, dec_errs, "o-", color="steelblue", markersize=4)
    ax.axhline(0, color="gray", linewidth=0.5)
    ax.axhline(+6, color="red", linewidth=0.5, linestyle="--", alpha=0.5,
               label="0.1 deg tolerance")
    ax.axhline(-6, color="red", linewidth=0.5, linestyle="--", alpha=0.5)
    ax.set_ylabel("Declination Error (arcmin)")
    ax.set_title(f"Moon Declination Error "
                 f"(max={moon_stats['dec_max_deg']*60:.1f} arcmin)")
    ax.legend(fontsize=8)
    ax.xaxis.set_major_formatter(DateFormatter("%Y-%m"))
    ax.tick_params(axis="x", rotation=45)

    # -- Sun fractional distance error --
    ax = axes[1, 0]
    dates_s = [r["dt"] for r in sun_results]
    errs_pct = [r["dist_err"] / r["dist_ref"] * 100 for r in sun_results]
    ax.plot(dates_s, errs_pct, "o-", color="darkorange", markersize=4)
    ax.axhline(0, color="gray", linewidth=0.5)
    ax.axhline(+0.01, color="red", linewidth=0.5, linestyle="--", alpha=0.5,
               label="0.01% tolerance")
    ax.axhline(-0.01, color="red", linewidth=0.5, linestyle="--", alpha=0.5)
    ax.set_ylabel("Distance Error (%)")
    ax.set_title(f"Sun Fractional Distance Error "
                 f"(max={max(abs(e) for e in errs_pct):.5f}%)")
    ax.legend(fontsize=8)
    ax.xaxis.set_major_formatter(DateFormatter("%Y-%m"))
    ax.tick_params(axis="x", rotation=45)

    # -- Sun declination error --
    ax = axes[1, 1]
    dec_errs_s = [r["dec_err"] * 60 for r in sun_results]  # arcmin
    ax.plot(dates_s, dec_errs_s, "o-", color="darkorange", markersize=4)
    ax.axhline(0, color="gray", linewidth=0.5)
    ax.axhline(+1, color="red", linewidth=0.5, linestyle="--", alpha=0.5,
               label="1 arcmin tolerance")
    ax.axhline(-1, color="red", linewidth=0.5, linestyle="--", alpha=0.5)
    ax.set_ylabel("Declination Error (arcmin)")
    ax.set_title(f"Sun Declination Error "
                 f"(max={sun_stats['dec_max_deg']*60:.2f} arcmin)")
    ax.legend(fontsize=8)
    ax.xaxis.set_major_formatter(DateFormatter("%Y-%m"))
    ax.tick_params(axis="x", rotation=45)

    plt.tight_layout()
    out_path = Path(__file__).resolve().parent / "ws1_ephemeris.png"
    plt.savefig(out_path, dpi=150)
    print(f"\n  Plot saved: {out_path}")
    plt.close()


# ============================================================================
# Main
# ============================================================================

def main():
    print("Pytheas Ephemeris Verification (WS1)")
    print("=" * 76)

    # -- Try live Horizons fetch --
    print("\nAttempting JPL Horizons API fetch...")
    moon_live = fetch_horizons("301")
    sun_live = fetch_horizons("10")

    use_live = moon_live is not None and sun_live is not None

    if use_live:
        print("\n  Using LIVE Horizons data")
        moon_ref = moon_live
        sun_ref = sun_live
    else:
        print("\n  Falling back to HARDCODED DE441 reference data")
        moon_ref = MOON_HARDCODED
        sun_ref = SUN_HARDCODED

    # Compute frame-independent invariants from XYZ reference data
    moon_inv = compute_ref_invariants(moon_ref)
    sun_inv = compute_ref_invariants(sun_ref)

    # -- Moon comparison --
    # Pass criteria: distance < 200 km, direction < 0.1 deg
    moon_results, moon_stats = compare_body(
        "Moon", moon_inv, moon_position_ecef,
        dist_tol_km=200.0, angle_tol_deg=0.1)
    print_report(moon_results, moon_stats)

    # -- Sun comparison --
    # Pass criteria: distance < 0.01% (~15000 km), direction < 1 arcmin
    sun_results, sun_stats = compare_body(
        "Sun", sun_inv, sun_position_ecef,
        dist_tol_km=15000.0, angle_tol_deg=1.0 / 60.0)
    print_report(sun_results, sun_stats)

    # -- Sun fractional distance --
    frac_errs = [abs(r["dist_err"] / r["dist_ref"] * 100)
                 for r in sun_results]
    sun_frac_pass = max(frac_errs) < 0.01
    print(f"\n  Sun Fractional Distance Error:")
    print(f"    Max: {max(frac_errs):.6f}%  (tol: 0.01%)  "
          f"{'PASS' if sun_frac_pass else 'FAIL'}")

    # -- Systematic bias vs T --
    print(f"\n  Systematic Bias Analysis:")
    t0 = datetime(2024, 1, 1)
    for label, results in [("Moon", moon_results), ("Sun", sun_results)]:
        t_sec = np.array([(r["dt"] - t0).total_seconds() for r in results])
        dist_e = np.array([r["dist_err"] for r in results])
        dec_e = np.array([r["dec_err"] for r in results])
        if len(t_sec) > 2:
            pd = np.polyfit(t_sec, dist_e, 1)
            slope_d = pd[0] * 365.25 * 86400
            print(f"    {label} distance trend: {slope_d:+.2f} km/year")
            pp = np.polyfit(t_sec, dec_e, 1)
            slope_p = pp[0] * 365.25 * 86400
            print(f"    {label} dec trend:      {slope_p:+.6f} deg/year "
                  f"({slope_p*60:+.4f} arcmin/year)")

    # -- Generate plots --
    try:
        plot_errors(moon_results, sun_results, moon_stats, sun_stats)
    except Exception as e:
        print(f"\n  [WARN] Plot generation failed: {e}")

    # -- Final verdict --
    print(f"\n{'='*76}")
    print(f"  FINAL VERDICT")
    print(f"{'='*76}")

    checks = [
        ("Moon distance < 200 km", moon_stats["dist_pass"]),
        ("Moon direction < 0.1 deg (6 arcmin)", moon_stats["angle_pass"]),
        ("Sun distance < 0.01%", sun_frac_pass),
        ("Sun direction < 1 arcmin", sun_stats["angle_pass"]),
    ]
    all_pass = all(p for _, p in checks)
    for label, passed in checks:
        print(f"    {'PASS' if passed else 'FAIL'}  {label}")

    print(f"\n    {'ALL CHECKS PASSED' if all_pass else 'SOME CHECKS FAILED'}")
    print(f"{'='*76}")

    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())

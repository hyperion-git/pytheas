# Pytheas

**Gravitational acceleration g(t) at a point on Earth**

Named after [Pytheas of Massalia](https://en.wikipedia.org/wiki/Pytheas)
(c. 325 BC), the Greek explorer who sailed from the Mediterranean to
Britain and made the first recorded systematic observations connecting
ocean tides to the phases of the Moon.

Pytheas computes the total gravitational acceleration projected along an
arbitrary measurement axis at any location on Earth, as a function of
time.  It is a single-file, dependency-light tool designed for quick
tidal predictions in precision accelerometry, gravimetry, and sensor
characterization.

**Dependency:** numpy only (matplotlib optional, for `--plot`).


## Installation

### Quick install (pip)

```bash
pip install pytheas
```

### From source with micromamba (recommended)

```bash
git clone https://github.com/hyperion-git/pytheas.git
cd pytheas

# Option A: use the included environment file
micromamba create -f environment.yml -y
micromamba activate pytheas

# Option B: create manually
micromamba create -n pytheas python=3.12 numpy matplotlib pytest -c conda-forge -y
micromamba activate pytheas

# Install in development mode
pip install -e .

# Verify
pytheas --lat 48.14 --lon 11.58 --alt 500
python -m pytest tests/ -v
```

### From source with pip only

```bash
git clone https://github.com/hyperion-git/pytheas.git
cd pytheas
pip install -e ".[test]"
```


## Quick Start

### Python API

```python
from datetime import datetime
from pytheas import compute_g, compute_timeseries

# Single epoch -- vertical sensor at Munich, 500 m altitude
result = compute_g(
    datetime(2025, 3, 20, 12, 0),
    lat_deg=48.14, lon_deg=11.58, alt_m=500.0,
)
print(f"g       = {result['g_total']:.10f} m/s^2")
print(f"g_tidal = {result['g_tidal'] * 1e6:.4f} um/s^2")

# 48-hour timeseries
data = compute_timeseries(
    start=datetime(2025, 3, 20),
    end=datetime(2025, 3, 22),
    lat_deg=48.14, lon_deg=11.58, alt_m=500.0,
    interval_minutes=10.0,
)
# data['times'], data['g_total'], data['g_tidal'], ...
```

### Tilted sensor

```python
# Sensor tilted 15 deg from vertical toward east at Ulm, Eselsberg
result = compute_g(
    datetime(2025, 3, 20, 12, 0),
    lat_deg=48.42, lon_deg=9.96, alt_m=620.0,
    zenith_deg=15.0, azimuth_deg=90.0,
)
print(f"g       = {result['g_total']:.10f} m/s^2")
print(f"g_tidal = {result['g_tidal'] * 1e6:.4f} um/s^2")
```

The `zenith_deg` parameter is the angle from vertical (0 = straight up,
90 = horizontal) and `azimuth_deg` is clockwise from north.  For a
perfectly vertical sensor both default to 0.  The static component
scales as `cos(zenith)` while the tidal projection depends on the full
3D geometry.

### Horizontal sensor

```python
# North-pointing horizontal accelerometer at Ulm, Eselsberg
result = compute_g(
    datetime(2025, 3, 20, 12, 0),
    lat_deg=48.42, lon_deg=9.96, alt_m=620.0,
    zenith_deg=90.0, azimuth_deg=0.0,
)
print(f"horizontal tidal = {result['g_tidal'] * 1e6:.4f} um/s^2")
```

### N-sample timeseries

```python
# 500 evenly spaced samples over a 48-hour window at Ulm, Eselsberg
data = compute_timeseries(
    start=datetime(2025, 3, 20),
    end=datetime(2025, 3, 22),
    lat_deg=48.42, lon_deg=9.96, alt_m=620.0,
    n_samples=500,
)
print(f"{len(data['times'])} points, dt ~ "
      f"{(data['times'][1] - data['times'][0]).total_seconds():.1f} s")
```

When `n_samples` is given it overrides the default `interval_minutes`
cadence.  The samples are inclusive of both endpoints.

### Command line

```bash
# Default: 48-hour timeseries from now, vertical axis
pytheas --lat 48.14 --lon 11.58 --alt 500

# Specify time window
pytheas --lat 48.14 --lon 11.58 --alt 500 --start 2025-03-20 --hours 72

# Horizontal axis, CSV export, plot
pytheas --lat 48.14 --lon 11.58 --alt 500 --zenith 90 --azimuth 0 \
        --csv output.csv --plot
```

Or run as a module: `python -m pytheas --lat 48.14 --lon 11.58`.


## What It Computes

The output `g_total` along the measurement axis is:

```
g_total(t) = g_static + delta * [a_tidal_moon(t) + a_tidal_sun(t)] . n_hat
```

| Component | Source | Typical magnitude |
|-----------|--------|-------------------|
| `g_static` | GRS80 normal gravity + free-air correction | 9.78 -- 9.83 m/s^2 |
| `g_tidal_moon` | Exact Newtonian lunar tidal acceleration | ~1.1 um/s^2 peak |
| `g_tidal_sun` | Exact Newtonian solar tidal acceleration | ~0.5 um/s^2 peak |
| `delta` | Elastic Earth amplification (gravimetric factor) | 1.1608 |


## Physics Model

### Normal gravity

The static gravitational acceleration on the GRS80 reference ellipsoid is
computed using the Somigliana closed-form formula:

```
gamma_0 = gamma_e * (1 + k * sin^2(phi)) / sqrt(1 - e^2 * sin^2(phi))
```

where `gamma_e = 9.7803253141 m/s^2` (equatorial) and
`gamma_p = 9.8321849378 m/s^2` (polar).  The altitude correction uses
the second-order free-air gradient:

```
gamma(h) = gamma_0 * [1 - 2*(1+f+m-2f*sin^2(phi))*h/a + 3*h^2/a^2]
```

This is accurate to <1 nGal for altitudes below 5 km.

### Tidal acceleration

The tidal acceleration from a celestial body at position **R** acting on
an observer at **r** is computed using the exact Newtonian formula:

```
a_tidal = GM * [(R - r) / |R - r|^3  -  R / |R|^3]
```

This is the difference between the gravitational pull at the observer's
position and at the Earth's center.  No tidal-tensor (gradient)
approximation is used, avoiding the ~2.5% truncation error inherent in
the first-order expansion.

### Ephemerides

**Sun:** Low-precision Meeus (Astronomical Algorithms) ephemeris.
Ecliptic longitude accuracy ~1 arcmin, distance via true anomaly.
Results in <1 nGal tidal error.

**Moon:** Truncated Meeus (ch. 47) ephemeris with:
- 24 longitude terms (Table 47.A)
- 23 distance terms (Table 47.A)
- 18 latitude terms (Table 47.B)
- A1, A2, A3 additive corrections
- Eccentricity correction factor E

Position accuracy ~0.1 deg, distance ~200 km.  The tidal acceleration
error is ~0.1% of the lunar signal, or ~1 nGal.

### Elastic Earth correction

The solid Earth is not rigid: it deforms under the tidal potential.  This
amplifies the surface gravity change by the gravimetric factor:

```
delta = 1 + h2 - (3/2) * k2 = 1.1608
```

using IERS 2010 Love numbers h2 = 0.6078, k2 = 0.2980.  A single
frequency-independent value is used; the true factor has a ~1% resonance
near the K1 frequency due to the free core nutation (FCN), contributing
up to ~10 nGal of error near that constituent.

### Coordinate system

All computation is done in ECEF (Earth-Centered, Earth-Fixed) Cartesian
coordinates.  The ECI-to-ECEF rotation uses GMST with the IAU polynomial.
The measurement axis is specified via zenith angle (0 = vertical up,
90 = horizontal) and azimuth (clockwise from north), decomposed into the
local East-North-Up frame.


## Accuracy Budget

For a vertical sensor at a continental (inland) site:

| Error source | Magnitude | Notes |
|-------------|-----------|-------|
| Lunar ephemeris (position) | ~1 nGal | 0.1 deg = 0.1% of tidal signal |
| Lunar ephemeris (distance) | ~1 nGal | 200 km / 385000 km ~ 0.05% |
| Solar ephemeris | <1 nGal | 1 arcmin, negligible |
| Elastic Earth (no FCN) | ~10 nGal | 1% near K1 only |
| Normal gravity | <1 nGal | GRS80 + 2nd-order free-air |
| Ocean loading (not modeled) | 1--50 nGal | negligible inland, significant at coast |
| Atmospheric pressure (not modeled) | ~3 nGal/hPa | needs barometer data |
| Polar motion (not modeled) | <1 nGal | |
| Planetary tides (not modeled) | <0.1 nGal | Jupiter, Venus combined |

**Total for inland sites: ~10 nGal (1e-8 m/s^2).**

Coastal sites may see up to 50 nGal additional error from unmodeled ocean
loading.


## API Reference

### `compute_g(dt, lat_deg, lon_deg, alt_m, zenith_deg=0, azimuth_deg=0)`

Single-epoch computation.  Returns a dict with keys:

- `g_total` -- total g projected on axis (m/s^2)
- `g_static` -- normal gravity component (m/s^2)
- `g_tidal` -- total tidal perturbation (m/s^2)
- `g_tidal_moon` -- lunar tidal contribution (m/s^2)
- `g_tidal_sun` -- solar tidal contribution (m/s^2)

### `compute_timeseries(start, end, lat_deg, lon_deg, alt_m, ..., n_samples=None)`

Multi-epoch computation over `[start, end]`.  By default uses
`interval_minutes=10.0` cadence; pass `n_samples=N` to get exactly N
evenly spaced samples instead.  Returns the same keys as `compute_g`,
each as a numpy array, plus `times` (list of datetime objects).

### Building blocks

| Function | Description |
|----------|-------------|
| `normal_gravity(lat_deg, alt_m)` | GRS80 gravity with free-air correction |
| `sun_position_ecef(dt)` | Sun ECEF position (m) |
| `moon_position_ecef(dt)` | Moon ECEF position (m) |
| `geodetic_to_ecef(lat, lon, alt)` | Geodetic to Cartesian |
| `enu_basis(lat, lon)` | Local East-North-Up unit vectors |
| `measurement_axis(lat, lon, zenith, azimuth)` | Sensor axis vector |
| `tidal_acceleration(r, R, GM)` | Exact tidal acceleration |
| `julian_date(dt)` | UTC datetime to Julian Date |
| `gmst_rad(dt)` | Greenwich Mean Sidereal Time |

### Constants

All physical constants are accessible as module attributes:

```python
import pytheas
pytheas.GM_MOON   # 4.9028695e12  m^3/s^2
pytheas.GM_SUN    # 1.32712440018e20 m^3/s^2
pytheas.A_GRS80   # 6378137.0 m
pytheas.DELTA_GRAV  # 1.1608
pytheas.H2, pytheas.K2  # Love numbers
```


## Limitations

Effects **not** included that would be needed for sub-10-nGal accuracy:

- **Ocean tidal loading:** Requires site-specific loading coefficients from
  models like FES2014 or GOT4.10c.  Dominant at coastal sites (up to 50 nGal).
- **Atmospheric pressure loading:** ~3 nGal per hPa of barometric pressure
  anomaly.  Requires real-time barometer data.
- **FCN Love number correction:** The gravimetric factor has a ~1% resonance
  near the K1 tidal frequency.  A frequency-dependent delta(f) table would
  reduce this error.
- **JPL planetary ephemeris:** Using DE440/DE441 instead of Meeus would give
  sub-arcsecond positions, reducing ephemeris error to <0.01 nGal.
- **Polar motion and LOD:** Earth orientation parameters from IERS bulletins.


## License

MIT

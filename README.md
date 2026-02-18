# Pytheas

**Gravitational acceleration g(t) at a point on Earth**

Named after [Pytheas of Massalia](https://en.wikipedia.org/wiki/Pytheas)
(c. 325 BC), the Greek explorer who first systematically connected
ocean tides to the phases of the Moon.

Pytheas computes time-dependent gravitational acceleration projected
along an arbitrary measurement axis at any point on Earth.  A
single-file, dependency-light tool, it targets tidal predictions for
precision accelerometry, gravimetry, and sensor characterization.

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
pytheas --lat 48.42 --lon 9.96 --alt 620
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

# Single epoch -- vertical sensor at Ulm, Eselsberg (620 m)
result = compute_g(
    datetime(2025, 3, 20, 12, 0),
    lat_deg=48.42, lon_deg=9.96, alt_m=620.0,
)
print(f"g       = {result['g_total']:.10f} m/s^2")  # ≈ 9.8073743815 m/s²
print(f"g_tidal = {result['g_tidal'] * 1e6:.4f} µm/s^2")

# 48-hour timeseries
data = compute_timeseries(
    start=datetime(2025, 3, 20),
    end=datetime(2025, 3, 22),
    lat_deg=48.42, lon_deg=9.96, alt_m=620.0,
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
print(f"g_tidal = {result['g_tidal'] * 1e6:.4f} µm/s^2")
```

The `zenith_deg` parameter specifies the angle from vertical (0 = up,
90 = horizontal); `azimuth_deg` specifies the bearing clockwise from
north.  Both default to 0 (vertical sensor).  The static component
scales as `cos(zenith)`, while the tidal projection depends on the full
3D geometry.

### Horizontal sensor

```python
# North-pointing horizontal accelerometer at Ulm, Eselsberg
result = compute_g(
    datetime(2025, 3, 20, 12, 0),
    lat_deg=48.42, lon_deg=9.96, alt_m=620.0,
    zenith_deg=90.0, azimuth_deg=0.0,
)
print(f"horizontal tidal = {result['g_tidal'] * 1e6:.4f} µm/s^2")
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

When `n_samples` is given, it overrides `interval_minutes`.
Both endpoints are included.

### Command line

```bash
# Default: 48-hour timeseries from now, vertical axis
pytheas --lat 48.42 --lon 9.96 --alt 620

# Specify time window
pytheas --lat 48.42 --lon 9.96 --alt 620 --start 2025-03-20 --hours 72

# Horizontal axis, CSV export, plot
pytheas --lat 48.42 --lon 9.96 --alt 620 --zenith 90 --azimuth 0 \
        --csv output.csv --plot
```

Or run as a module: `python -m pytheas --lat 48.42 --lon 9.96`.


## What It Computes

The output `g_total` along the measurement axis is:

```
g_total(t) = g_static + delta * [a_tidal_moon(t) + a_tidal_sun(t)] . n_hat
```

| Component | Source | Typical magnitude |
|-----------|--------|-------------------|
| `g_static` | GRS80 normal gravity + free-air correction | 9.78 -- 9.83 m/s² |
| `g_tidal_moon` | Exact Newtonian lunar tidal acceleration | ~1.1 µm/s² peak |
| `g_tidal_sun` | Exact Newtonian solar tidal acceleration | ~0.5 µm/s² peak |
| `delta` | Elastic Earth amplification (gravimetric factor) | 1.1608 |


## Physics Model

### Normal gravity

The static gravitational acceleration on the GRS80 reference ellipsoid
follows the Somigliana closed-form formula:

```
gamma_0 = gamma_e * (1 + k * sin^2(phi)) / sqrt(1 - e^2 * sin^2(phi))
```

where `gamma_e = 9.7803253141 m/s²` (equatorial) and
`gamma_p = 9.8321849378 m/s²` (polar).  Altitude dependence uses the
second-order free-air gradient:

```
gamma(h) = gamma_0 * [1 - 2*(1+f+m-2f*sin^2(phi))*h/a + 3*h^2/a^2]
```

Accuracy is <1 nGal below ~400 m altitude and <20 nGal below 1 km.

### Tidal acceleration

The tidal acceleration from a celestial body at position **R** on an
observer at **r** follows the exact Newtonian expression:

```
a_tidal = GM * [(R - r) / |R - r|^3  -  R / |R|^3]
```

This computes the difference between the gravitational acceleration at
the observer and at the Earth's center, without the tidal-tensor
(gradient) approximation.  Using the exact formula avoids the ~2.5%
truncation error of the first-order expansion.

### Ephemerides

**Sun:** Low-precision Meeus (*Astronomical Algorithms*) ephemeris.
Ecliptic longitude accurate to ~1 arcmin; distance derived from true
anomaly.  Tidal error contribution: <1 nGal.

**Moon:** Truncated Meeus (ch. 47) ephemeris including:
- 24 longitude terms (Table 47.A)
- 23 distance terms (Table 47.A)
- 18 latitude terms (Table 47.B)
- A1, A2, A3 additive corrections
- Eccentricity correction factor E

Position accurate to ~0.1 deg; distance accurate to ~200 km.  The
resulting tidal error is ~0.1% of the lunar signal (~1 nGal).

### Elastic Earth correction

The solid Earth deforms under the tidal potential, amplifying the
surface gravity change by the gravimetric factor:

```
delta = 1 + h2 - (3/2) * k2 = 1.1608
```

using IERS 2010 Love numbers h2 = 0.6078, k2 = 0.2980.  This single
frequency-independent value neglects the ~1% resonance near the K1
frequency caused by the free core nutation (FCN), introducing up to
~10 nGal of error at that constituent.

### Coordinate system

All calculations use ECEF (Earth-Centered, Earth-Fixed) Cartesian
coordinates.  The ECI-to-ECEF rotation applies GMST from the IAU
polynomial.  The measurement axis, specified by zenith angle
(0 = up, 90 = horizontal) and azimuth (clockwise from north), is
decomposed in the local East-North-Up frame.


## Accuracy Budget

For a vertical sensor at a continental (inland) site:

| Error source | Magnitude | Notes |
|-------------|-----------|-------|
| Lunar ephemeris (position) | ~1 nGal | 0.1 deg = 0.1% of tidal signal |
| Lunar ephemeris (distance) | ~1 nGal | 200 km / 385000 km ~ 0.05% |
| Solar ephemeris | <1 nGal | 1 arcmin, negligible |
| Elastic Earth (no FCN) | ~10 nGal | 1% near K1 only |
| Normal gravity | <1 nGal | GRS80 + 2nd-order free-air |
| Ocean loading (not modeled) | 1--3000 nGal | negligible inland, up to ~3 µGal at coast |
| Atmospheric pressure (not modeled) | ~300 nGal/hPa | needs barometer data |
| Polar motion (not modeled) | ~10--20 nGal | 48-hour variation; static offset ~5 µGal |
| Planetary tides (not modeled) | ~1--6 nGal | Venus dominates at inferior conjunction |

**Total for inland sites: ~10 nGal (1e-10 m/s²).**

At coastal sites, unmodeled ocean loading can add up to ~3 µGal of error.


## API Reference

### `compute_g(dt, lat_deg, lon_deg, alt_m, zenith_deg=0, azimuth_deg=0)`

Single-epoch computation.  Returns a dict with keys:

- `g_total` -- total g projected on axis (m/s²)
- `g_static` -- normal gravity component (m/s²)
- `g_tidal` -- total tidal perturbation (m/s²)
- `g_tidal_moon` -- lunar tidal contribution (m/s²)
- `g_tidal_sun` -- solar tidal contribution (m/s²)

### `compute_timeseries(start, end, lat_deg, lon_deg, alt_m, ..., n_samples=None)`

Multi-epoch computation over `[start, end]`.  Default cadence is
`interval_minutes=10.0`; pass `n_samples=N` for exactly N evenly spaced
samples instead.  Returns the same keys as `compute_g` (each as a numpy
array), plus `times` (list of datetime objects).

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

The following effects, omitted here, would be needed for sub-10-nGal accuracy:

- **Ocean tidal loading:** Requires site-specific coefficients from models
  such as FES2014 or GOT4.10c.  Dominant at coastal sites (up to 50 nGal).
- **Atmospheric pressure loading:** ~3 nGal per hPa of barometric anomaly.
  Requires real-time barometer data.
- **FCN Love number correction:** A frequency-dependent delta(f) table would
  capture the ~1% gravimetric resonance near K1.
- **JPL planetary ephemeris:** DE440/DE441 would provide sub-arcsecond
  positions, reducing ephemeris error to <0.01 nGal.
- **Polar motion and LOD:** Earth orientation parameters from IERS bulletins.


## License

MIT

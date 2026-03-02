# Pytheas Roadmap

## Current: v3.2 — Lightweight Reference Tool

Single-file, NumPy-only library. ~10 nGal accuracy for inland sites.
Physics: WGS84 normal gravity, Meeus ephemerides, exact Newtonian tides, elastic Earth (IERS 2010).

---

## Phase 1 — Streamline & Harden (next)

Goal: Better API ergonomics, comprehensive testing, no new physics.

### 1a. Dataclass API
- [ ] `GravityResult` dataclass for `compute_g` return type
- [ ] `TimeSeries` dataclass (array fields) for `compute_timeseries`
- [ ] `dataclasses.asdict()` for backward compatibility
- [ ] Update CLI and README examples

### 1b. External Validation
- [ ] Validate Moon position against JPL Horizons (query via API or cached reference)
- [ ] Validate Sun position against JPL Horizons
- [ ] Validate tidal acceleration against IERS/ETERNA reference predictions
- [ ] Validate normal gravity against published benchmark values (IGSN71 stations)

### 1c. Property-Based Tests
- [ ] Symmetry: tidal acceleration antisymmetric under body reflection
- [ ] Limiting cases: tidal → 0 as R → ∞, tidal → gradient approx as r/R → 0
- [ ] Conservation: g_tidal = g_tidal_moon + g_tidal_sun (decomposition identity)
- [ ] Coordinate invariants: ENU basis orthonormality, right-handedness for all lat/lon
- [ ] Normal gravity: monotonic increase with latitude, equator < pole

---

## Phase 2 — Modular Accuracy Upgrades (~1 nGal target)

Goal: Optional modules that improve accuracy when external data/deps are available.
New dependencies allowed: astropy, jplephem.

### 2a. JPL Ephemerides (`pytheas.ephem`)
- Replace Meeus with DE440 via `jplephem`
- Fallback to built-in Meeus when jplephem unavailable
- Expected improvement: Moon position 200 km → <1 km; tidal error ~1 nGal → <0.01 nGal

### 2b. Atmospheric Pressure Loading (`pytheas.atmosphere`)
- Simple admittance model: δg ≈ -0.3 µGal/hPa × (P - P_ref)
- User supplies local pressure or a time series
- P_ref from standard atmosphere at station altitude
- Expected correction: up to ~3 µGal (≈ 300 nGal) during weather systems

### 2c. Ocean Tide Loading (`pytheas.ocean`)
- Precomputed loading coefficients per station (from OLFG/SPOTL or bundled table)
- 11 major tidal constituents (M2, S2, N2, K2, K1, O1, P1, Q1, Mf, Mm, Ssa)
- Hardisp-style harmonic synthesis
- Expected correction: 0–50 nGal inland, up to ~5 µGal at coastal sites

### 2d. Earth Orientation (`pytheas.orientation`)
- Polar motion (x_p, y_p) from IERS Bulletin A
- UT1-UTC correction
- Requires periodic data file updates or network fetch
- Expected correction: <1 nGal for tides, matters more for absolute positioning

### 2e. Free Core Nutation (`pytheas.fcn`)
- Analytic resonance model (no external data)
- Modifies body tide Love numbers near diurnal band
- Expected correction: ~1–2 nGal at diurnal frequencies

### Module Priority Order
1. Atmospheric pressure (biggest correction, simplest to implement)
2. Ocean loading (biggest correction at coastal sites)
3. JPL ephemerides (cleanest accuracy gain, well-supported packages)
4. FCN (pure math, small correction)
5. Earth orientation (smallest correction, most operational burden)

---

## Phase 3 — Future Considerations (not planned)

- Solid Earth pole tide
- Non-tidal ocean/atmosphere loading (de-aliasing)
- Instrument response modeling (spring, superconducting gravimeter)
- Real-time mode (streaming pressure + tidal prediction)
- Comparison/benchmarking CLI against ETERNA, TSoft, Tsoft

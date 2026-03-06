# Implementation Notes

What the code computes, how it computes it, and how accurate the result is.

For mathematical derivations, see [Theory](theory.md).  For usage examples, see [Tutorial](tutorial.md).


## 1. Pipeline Overview

A call to `compute_g(dt, lat_deg, lon_deg, alt_m, zenith_deg, azimuth_deg)` executes the following chain:

```
(lat, lon, alt) ──► geodetic_to_ecef ──► r_ecef
                                          │
dt ──► moon_position_ecef ──► R_moon ─────┤
   └─► sun_position_ecef  ──► R_sun  ─────┤
                                           ▼
                              tidal_acceleration(r, R, GM)
                              × DELTA_GRAV  (elastic amplification)
                                           │
(lat, lon) ──► enu_basis ──► (E, N, U)    │
(zenith, azimuth) ──► measurement_axis ──► n_hat
                                           │
                              dot(a_tidal, n_hat)
                                           │
(lat, alt) ──► normal_gravity ──► γ        │
                              γ · cos(θ) + tidal projection
                                           ▼
                                       GravityResult
```

`compute_g()` and `compute_timeseries()` report the scalar reading on the chosen measurement axis, so a vertical sensor returns a positive static value $\gamma$.  `LabFrame.field()`, by contrast, returns the signed ENU gravity vector itself, with $g_U < 0$.


## 2. Normal Gravity

### Somigliana Formula

Normal gravity on the WGS84 ellipsoid:

$$\gamma_0(\varphi) = \gamma_e \frac{1 + k \sin^2\varphi}{\sqrt{1 - e^2 \sin^2\varphi}}$$

### WGS84 Constants

| Symbol | Value | Description |
|--------|-------|-------------|
| $a$ | 6 378 137.0 m | Semi-major axis |
| $f$ | 1/298.257 223 563 | Flattening |
| $b$ | 6 356 752.314 245 m | Semi-minor axis |
| $e^2$ | 0.006 694 379 9901 | First eccentricity squared |
| $GM$ | 3.986 004 418 $\times 10^{14}$ m$^3$/s$^2$ | Geocentric gravitational parameter |
| $\Omega$ | 7.292 115 $\times 10^{-5}$ rad/s | Rotation rate |
| $\gamma_e$ | 9.780 325 3141 m/s$^2$ | Equatorial normal gravity |
| $\gamma_p$ | 9.832 184 9378 m/s$^2$ | Polar normal gravity |

### Free-Air Correction

Altitude dependence uses the second-order Taylor expansion:

$$\gamma(\varphi, h) = \gamma_0 \left[1 - \frac{2(1+f+m-2f\sin^2\varphi)}{a}\,h + \frac{3}{a^2}\,h^2\right]$$

where $m = \Omega^2 a^2 b / GM$.

The first-order free-air gradient is about $-308.6\ \mu\text{Gal/m}$ (equivalently $-3.086\ \mu\text{m}\,\text{s}^{-2}\,\text{m}^{-1}$).  At 500 m altitude, the quadratic term changes $\gamma$ by about $1.8\times10^{-7}$ m/s$^2$ ($\approx 1.8\times10^4$ nGal): tiny compared with total gravity, but still worth retaining because the code already has the closed-form second-order expression.


## 3. Ephemerides

### Moon: Meeus ch. 47

The lunar ephemeris uses the truncated series from Meeus, *Astronomical Algorithms*:

- **24** longitude terms (Table 47.A) → ecliptic longitude
- **23** distance terms (Table 47.A) → geocentric distance
- **18** latitude terms (Table 47.B) → ecliptic latitude
- Additive corrections: **3** longitude terms and **6** latitude terms (using A1/A2/A3 and combinations with $L'$, $F$, $M'$)
- **Eccentricity correction** factor $E = 1 - 0.002516\,T - 0.0000074\,T^2$

The five fundamental arguments ($L'$, $D$, $M'$, $M_s$, $F$) are polynomials in Julian centuries $T$ since J2000.0.

**Position accuracy:** ~0.1$\degree$ (~700 km at mean distance).
**Distance accuracy:** ~200 km out of 385 000 km (~0.05%).

The tidal acceleration scales as $R^{-3}$, so a 0.05% distance error produces ~0.15% tidal error: for a 100 $\mu$Gal lunar tide this is ~0.15 $\mu$Gal (~150 nGal).  Position error contributes a comparable order.

### Sun: Low-Precision Meeus

The solar ephemeris uses:

- **Mean longitude** $L_0$ and **mean anomaly** $M$ as polynomials in $T$
- **Equation of center** $C$: three terms in $\sin(M)$, $\sin(2M)$, $\sin(3M)$
- **Distance** from the vis-viva relation: $R = a(1-e^2)/(1 + e\cos\nu)$
- **Ecliptic-to-equatorial** rotation using mean obliquity $\varepsilon$

**Accuracy:** ~1 arcmin in ecliptic longitude, ~0.01% in distance.

Solar tides are ~2.2$\times$ smaller than lunar tides, and the ephemeris is proportionally more accurate.  **Tidal error: typically O(10 nGal).**

### Coordinate Chain

Both ephemerides produce positions in **ecliptic coordinates** (geocentric), which are then transformed:

1. Ecliptic → Equatorial (ECI) via obliquity rotation
2. ECI → ECEF via GMST rotation

GMST uses the IAU polynomial in Julian centuries.


## 4. Tidal Acceleration

### Exact Newtonian Formula

$$\mathbf{a}_\text{tidal} = GM\left[\frac{\mathbf{R} - \mathbf{r}}{|\mathbf{R} - \mathbf{r}|^3} - \frac{\mathbf{R}}{R^3}\right]$$

This computes the difference between gravitational pull at the observer ($\mathbf{r}$) and at the Earth's center (origin), without expanding in $r/R$.

### Why Not the Gradient Approximation?

The standard tidal tensor (gradient) approximation:

$$a_i^\text{grad} \approx -\frac{GM}{R^3}\left(\delta_{ij} - 3\hat{R}_i\hat{R}_j\right) r_j$$

truncates at order $(r/R)^1$.  For the Moon with $r/R \approx 1/60$, the next-order term contributes:

$$\frac{a_\text{exact} - a_\text{grad}}{a_\text{exact}} \sim \frac{3r}{2R} \approx 2.5\%$$

At a peak tidal signal of ~1 $\mu$m/s$^2$ (~100 $\mu$Gal), this is **~2.5 $\mu$Gal (~2500 nGal)**.  The exact formula costs negligibly more to compute, so there is no reason to approximate.

### Elastic Amplification

The tidal acceleration is multiplied by the gravimetric factor before use:

$$\mathbf{a}_\text{measured} = \delta \cdot \mathbf{a}_\text{tidal}$$

with $\delta = 1 + h_2 - \tfrac{3}{2}k_2 = 1.1608$ (see Section 7).


## 5. Coordinate Transforms

### Geodetic → ECEF

$$\begin{pmatrix} x \\ y \\ z \end{pmatrix} = \begin{pmatrix} (N + h)\cos\varphi\cos\lambda \\ (N + h)\cos\varphi\sin\lambda \\ [N(1 - e^2) + h]\sin\varphi \end{pmatrix}$$

where $N = a / \sqrt{1 - e^2\sin^2\varphi}$ is the prime vertical radius of curvature.

### ENU Basis

The local East-North-Up unit vectors in ECEF:

$$\hat{\mathbf{e}}_E = (-\sin\lambda,\ \cos\lambda,\ 0)$$
$$\hat{\mathbf{e}}_N = (-\sin\varphi\cos\lambda,\ -\sin\varphi\sin\lambda,\ \cos\varphi)$$
$$\hat{\mathbf{e}}_U = (\cos\varphi\cos\lambda,\ \cos\varphi\sin\lambda,\ \sin\varphi)$$

These form an orthonormal right-handed triad.

### Measurement Axis

Given zenith angle $\theta$ (from vertical) and azimuth $\alpha$ (from north, clockwise):

$$\hat{\mathbf{n}} = \cos\theta\,\hat{\mathbf{e}}_U + \sin\theta\left(\cos\alpha\,\hat{\mathbf{e}}_N + \sin\alpha\,\hat{\mathbf{e}}_E\right)$$

### ECEF ↔ ENU Rotation

The rotation matrix $R = [\hat{\mathbf{e}}_E;\ \hat{\mathbf{e}}_N;\ \hat{\mathbf{e}}_U]$ converts ECEF vectors to ENU:

$$\mathbf{v}_\text{ENU} = R\,\mathbf{v}_\text{ECEF}$$

Rank-2 tensors transform as $T_\text{ENU} = R\,T_\text{ECEF}\,R^T$.


## 6. LabFrame Internals

`LabFrame.__init__` precomputes:

| Quantity | Expression | Time dependence |
|----------|------------|-----------------|
| $\mathbf{r}$ | `geodetic_to_ecef(lat, lon, alt)` | Static |
| ENU basis | `enu_basis(lat, lon)` | Static |
| $\gamma$ | `normal_gravity(lat, alt)` | Static |
| $\boldsymbol{\Omega}_\text{ENU}$ | $(0,\ \Omega\cos\varphi,\ \Omega\sin\varphi)$ | Static |
| $T_\text{Earth}$ | `_earth_gradient_tensor(lat, alt)` | Static |

`LabFrame.field(dt)` adds time-dependent tidal contributions:

1. Compute $\mathbf{R}_\text{Moon}(t)$, $\mathbf{R}_\text{Sun}(t)$ in ECEF
2. Compute tidal acceleration vectors in ECEF, multiply by $\delta$
3. Rotate tidal vectors to ENU
4. Assemble $\mathbf{g} = (0, 0, -\gamma) + \mathbf{a}_\text{Moon}^\text{ENU} + \mathbf{a}_\text{Sun}^\text{ENU}$
5. If `order >= 1`: compute tidal gradient tensors in ECEF, rotate to ENU, add to $T_\text{Earth}$
6. Return `GravityField(g, omega, T, ...)`

### Earth Gradient Tensor

The static Earth gradient tensor is diagonal in ENU:

$$T_\text{Earth} = \operatorname{diag}(T_H,\ T_H,\ T_{UU})$$

where $T_{UU} = \partial g_U/\partial h = -\partial\gamma/\partial h > 0$ near sea level (gravity weakens going up), consistent with the acceleration-gradient convention used by the tidal tensor and [Theory](theory.md) Eq. (5.3).  The horizontal components satisfy the vacuum Poisson trace condition:

$$T_{EE} + T_{NN} + T_{UU} = 2\Omega^2$$

Setting $T_{EE} = T_{NN} = T_H$:

$$T_H = \frac{2\Omega^2 - T_{UU}}{2}$$

Off-diagonal terms of order $O(f)$ (flattening) are deferred.

### Tidal Gradient Tensor

The tidal gradient tensor for a body at $\mathbf{R}$ with $\mathbf{d} = \mathbf{R} - \mathbf{r}$:

$$T_{ij} = GM\left[-\frac{\delta_{ij}}{|\mathbf{d}|^3} + \frac{3\,d_i\,d_j}{|\mathbf{d}|^5}\right]$$

This tensor is symmetric and traceless (vacuum Laplace equation $\nabla^2\Phi = 0$).


## 7. Accuracy Budget

### Modeled Error Sources

| Source | Error | Notes |
|--------|-------|-------|
| Lunar ephemeris (position) | ~100–300 nGal | 0.1$\degree$ angular error on a ~100 $\mu$Gal signal |
| Lunar ephemeris (distance) | ~150 nGal | 200 km / 385 000 km = 0.05% distance $\Rightarrow$ 0.15% tidal |
| Solar ephemeris | ~10–100 nGal | Smaller signal but comparable angular/distance uncertainties |
| Elastic Earth (no FCN) | ~100–1000 nGal | 1% Love-number mismatch on major constituents |
| Normal gravity | < 1 nGal | WGS84 Somigliana + 2nd-order free-air |
| Tidal formula (exact, not gradient) | < 0.1 nGal | Full expression, no truncation |
| **Total modeled** | **~200–1000 nGal** | Usually dominated by ephemeris + FCN treatment |

### Unmodeled Effects

| Effect | Magnitude | Notes |
|--------|-----------|-------|
| Ocean tidal loading | 1–3000 nGal | < 1 nGal deep inland; up to ~3 $\mu$Gal at coast |
| Atmospheric pressure | ~300 nGal/hPa | Requires barometer data |
| Polar motion | ~10–20 nGal | 48-hr variation; static offset ~5 $\mu$Gal |
| Planetary tides (Venus) | 1–6 nGal | At inferior conjunction |
| FCN resonance (K1) | ~100–1000 nGal | Frequency-dependent Love number |
| LOD variations | < 1 nGal | Earth rotation rate fluctuations |
| Non-tidal ocean/atmosphere | Variable | Weather, currents |

### Bottom Line

For a **vertical sensor** at an **inland continental site** below ~1 km altitude:

- **Total modeled error: typically sub-$\mu$Gal** (~10$^2$–10$^3$ nGal)
- Error is usually dominated by ephemeris accuracy and frequency-independent Love numbers near K1
- At coastal sites, unmodeled ocean loading can add up to ~3 $\mu$Gal

For horizontal or tilted sensors, the same accuracy applies to the projected tidal component.  The static projection $\gamma\cos\theta$ is exact to the precision of the Somigliana formula.


## 8. Validation

### External References

**JPL Horizons (DE441).**  Moon and Sun positions from JPL's DE441 planetary ephemeris are used as truth values in the regression tests (2024--2025 epochs):

- Moon distance agreement: < 200 km
- Moon declination agreement: < 0.15$\degree$
- Sun distance agreement: < 10 000 km (~0.007%)
- Sun declination agreement: < 0.2$\degree$

**GRS80 normal gravity.**  The WGS84 Somigliana formula (which uses GRS80-compatible constants) is validated against published GRS80 values at the equator and pole.

### Test Suite

The test suite (`tests/test_pytheas.py`) contains **79 test methods** across **13 test classes**, with additional parameterized cases in the external-validation and property sections:

| Category | Tests | What is checked |
|----------|-------|-----------------|
| Constants | 6 | GM values, $\delta$, Somigliana self-consistency |
| Normal gravity | 5 | Equator/pole values, latitude trend, altitude dependence, free-air gradient |
| Time utilities | 2 | Julian date, GMST at known epochs |
| Coordinates | 7 | ECEF conversion, ENU orthonormality, right-handedness, axis vectors |
| Ephemeris | 4 | Sun/Moon distances, position ranges, external validation |
| Tidal acceleration | 3 | Order of magnitude, Moon/Sun ratio, center-of-mass limit |
| Full pipeline | 6 | Decomposition sums, gravimetric amplification, zenith projection |
| Timeseries | 4 | Length, variation range, spectral periods, static constancy |
| Edge cases | 4 | South pole, date line, high altitude, negative longitude |
| LabFrame | 21 | Rotation vector, tensor symmetry/trace, EOM, Coriolis, timeseries |
| Dataclass API | 6 | Frozen behavior, `asdict` conversion |
| External validation | 5 | JPL Horizons Moon/Sun, GRS80 gravity |
| Physics properties | 6 | Scaling limits, ENU handedness, gravity bounds, Coriolis $\perp$ velocity |

### Property-Based Tests

Several tests use **property-based** (randomized) assertions rather than fixed reference values:

- **Tidal scaling:** $a_\text{tidal}$ decreases with distance (1x, 10x, 100x checks)
- **ENU orthonormality/right-handedness:** Verified at 7 global locations spanning quadrants
- **Gravity bounds:** $g \in [9.76, 9.84]$ m/s$^2$ for any latitude/altitude in range
- **Coriolis perpendicularity:** $(\boldsymbol{\Omega} \times \mathbf{v}) \cdot \mathbf{v} = 0$ for 20 random velocities

### Running Tests

```bash
pip install -e ".[test]"
pytest tests/ -v
```

The full parameterized suite should pass.

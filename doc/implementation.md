# Implementation Notes

What the code computes, how it computes it, and how accurate each output is.

For mathematical derivations, see [Theory](theory.md).  For usage examples, see [Tutorial](tutorial.md).

Pytheas is a body-tide predictor on top of a normal-gravity baseline.  The current implementation combines WGS84/Somigliana normal gravity with exact Newtonian Moon/Sun point-mass tides and a single scalar elastic-response factor.  It does not attempt to be a complete terrestrial gravity model.

**Included in the code path.**

- Normal gravity on the WGS84 ellipsoid with second-order free-air correction
- Exact Newtonian lunisolar body tides from point masses
- Rotation of vectors and tensors into the local ENU frame
- A single scalar gravimetric factor $\delta = 1 + h_2 - \tfrac{3}{2}k_2$

**Not included.**

- Static gravity anomalies beyond normal gravity (geoid, terrain, local mass distribution): affects the absolute background field and static tensor at a real site
- Ocean, atmospheric, and hydrological loading: can add site-dependent signals from tens of nGal to $\mu$Gal
- Full IERS constituent-dependent Earth response, including FCN-sensitive diurnal corrections and horizontal/tensor transfer factors: this is the main reason the quoted sub-$\mu$Gal accuracy applies only to the vertical body-tide channel
- A full non-diagonal Earth gradient model: the returned tensor is useful for local linearization, not precision gradiometry


## Notation

| Symbol | Code | Description |
|--------|------|-------------|
| $\varphi$ | `lat_deg` | Geodetic latitude |
| $\lambda$ | `lon_deg` | Geodetic longitude |
| $h$ | `alt_m` | Altitude above ellipsoid |
| $\gamma$ | `g_normal` | Normal gravity magnitude |
| $\mathbf{g}$ | `field.g` | Gravity vector (ENU, $g_U < 0$) |
| $\mathbf{T}$ | `field.T` | Gravity gradient tensor (ENU) |
| $\boldsymbol{\Omega}$ | `field.omega` | Earth rotation vector (ENU) |
| $\delta$ | `DELTA_GRAV` | Scalar gravimetric factor $1 + h_2 - \frac{3}{2}k_2$ |
| $h_2$, $k_2$ | `H2`, `K2` | Degree-2 Love numbers |
| $\hat{\mathbf{n}}$ | `n_hat` | Measurement axis unit vector |
| $\theta$ | `zenith_deg` | Zenith angle (0 = vertical) |
| $\alpha$ | `azimuth_deg` | Azimuth (clockwise from north) |


## 1. Pipeline Overview

A call to `compute_g(dt, lat_deg, lon_deg, alt_m, zenith_deg, azimuth_deg)` executes the following chain:

```
(lat, lon, alt) ──► geodetic_to_ecef ──► r_ecef
                                          │
dt ──► moon_position_ecef ──► R_moon ─────┤
   └─► sun_position_ecef  ──► R_sun  ─────┤
                                           ▼
                              tidal_acceleration(r, R, GM)
                              × DELTA_GRAV  (scalar elastic-response approximation)
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

**Position accuracy:** ~0.1° (~700 km at mean distance).
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

### Time Scale

The code uses UTC throughout, but the Meeus ephemerides are parameterized by Terrestrial Time (TT).  The difference TT $-$ UTC $\approx$ 69 s in 2025 (and growing with each leap second) is not corrected.  At the Moon's angular rate of $\sim$0.5°/h, 69 s produces $\sim$0.01° position error, which is small compared with the Meeus truncation error ($\sim$0.1°) but should be noted for future high-precision work.

### Measured Ephemeris Accuracy

Against JPL Horizons DE441 reference positions at four test epochs (2024--2025):

| Body | Distance error | Declination error |
|------|---------------|-------------------|
| Moon | 1.6--36.7 km (of 385 000 km) | $\le$ 0.146° |
| Sun  | $\le$ 5 830 km (of 150 000 000 km) | $\le$ 0.144° |

These are tighter than the Meeus textbook estimates, confirming that the truncated series is adequate for the $10^2$--$10^3$ nGal accuracy target.


## 4. Tidal Acceleration

### Exact Newtonian Formula

$$\mathbf{a}_\text{tidal} = GM\left[\frac{\mathbf{R} - \mathbf{r}}{|\mathbf{R} - \mathbf{r}|^3} - \frac{\mathbf{R}}{R^3}\right]$$

This computes the difference between gravitational pull at the observer ($\mathbf{r}$) and at the Earth's center (origin), without expanding in $r/R$.

**Numerical note.**  The subtraction $(\mathbf{R}-\mathbf{r})/|\mathbf{R}-\mathbf{r}|^3 - \mathbf{R}/R^3$ involves two nearly equal vectors, but catastrophic cancellation is not a concern: for the Moon, $r/R \approx 1/60$, so the subtraction loses at most $\sim$1.7 decimal digits, well within the $\sim$15.7 digits available in IEEE 754 double precision.

### Why Not the Gradient Approximation?

The standard tidal tensor (gradient) approximation:

$$a_i^\text{grad} \approx -\frac{GM}{R^3}\left(\delta_{ij} - 3\hat{R}_i\hat{R}_j\right) r_j$$

truncates at order $(r/R)^1$.  For the Moon with $r/R \approx 1/60$, the next-order term contributes:

$$\frac{a_\text{exact} - a_\text{grad}}{a_\text{exact}} \sim \frac{3r}{2R} \approx 2.5\%$$

At a peak tidal signal of ~1 $\mu$m/s$^2$ (~100 $\mu$Gal), this is **~2.5 $\mu$Gal (~2500 nGal)**.  The exact formula costs negligibly more to compute, so there is no reason to approximate.

### Elastic Amplification

Pytheas multiplies the rigid-Earth tidal acceleration by a single gravimetric factor before use:

$$\mathbf{a}_\text{measured} = \delta \cdot \mathbf{a}_\text{tidal}$$

with $\delta = 1 + h_2 - \tfrac{3}{2}k_2 = 1.1608$ (see [Theory](theory.md), Section 7).

This scalar factor is a simplification.  In the full IERS 2010 treatment (Chapters 6--7), the Earth response depends on constituent, latitude, and observable, and horizontal components use different transfer factors than vertical gravity.  Applying one $\delta$ to the full vector and to every tensor element is therefore only a first-order approximation; it is weakest in the K1/FCN band, in horizontal components, and in gravity-gradient tensors.


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

Off-diagonal terms of order $O(f)$ (flattening) are deferred.  This diagonal closure is adequate for the linearized equation of motion but sets an order-10 Eotvos floor on the total tensor.

### Tidal Gradient Tensor

The tidal gradient tensor for a body at $\mathbf{R}$ with $\mathbf{d} = \mathbf{R} - \mathbf{r}$:

$$T_{ij} = GM\left[-\frac{\delta_{ij}}{|\mathbf{d}|^3} + \frac{3\,d_i\,d_j}{|\mathbf{d}|^5}\right]$$

This tensor is symmetric and traceless (vacuum Laplace equation $\nabla^2\Phi = 0$).


## 7. Accuracy Budget

### Vertical body tide (best-supported channel)

| Source | Error | Notes |
|--------|-------|-------|
| Lunar ephemeris (position) | ~100–300 nGal | 0.1° angular error on a ~100 $\mu$Gal signal |
| Lunar ephemeris (distance) | ~150 nGal | 200 km / 385 000 km = 0.05% distance $\Rightarrow$ 0.15% tidal |
| Solar ephemeris | ~10–100 nGal | Smaller signal but comparable angular/distance uncertainties |
| Scalar elastic response ($\delta$, no constituent table) | ~100–1000 nGal | Dominant near diurnal constituents; quoted for the vertical body tide |
| Normal gravity baseline | < 1 nGal | WGS84 Somigliana + 2nd-order free-air |
| Tidal forcing formula (exact point mass) | < 0.1 nGal | Full expression, no $r/R$ truncation |
| **Total modeled vertical body tide** | **~200–1000 nGal** | Quiet inland site; excludes loading and weather |

### Effects not modeled in the vertical budget

| Effect | Magnitude | Notes |
|--------|-----------|-------|
| Ocean tidal loading | 1–3000 nGal | < 1 nGal deep inland; up to ~3 $\mu$Gal at coast |
| Atmospheric pressure | ~300 nGal/hPa | Requires barometer data |
| Polar motion | ~10–20 nGal | 48-hr variation; static offset ~5 $\mu$Gal |
| Planetary tides (Venus) | 1–6 nGal | At inferior conjunction |
| LOD variations | < 1 nGal | Earth rotation rate fluctuations |
| Non-tidal ocean/atmosphere | Variable | Weather, currents |

### Horizontal components and tilted sensors

The rigid Newtonian forcing geometry is still computed exactly, but the Earth response is not: Pytheas uses the same scalar $\delta$ that is derived for the vertical gravimetric response.  The full IERS treatment uses constituent- and latitude-dependent Love/Shida numbers for horizontal motion and tilt.  Horizontal components should therefore be treated as qualitative or low-order predictions, not as sub-$\mu$Gal products.

### Gradient tensor

The external-body tensor is exact for a rigid point mass, but the total tensor is limited by two implementation choices: the same scalar $\delta$ is applied to every tensor element, and the static Earth tensor is represented by a diagonal ENU approximation constrained only by the trace.  That is adequate for local linearization of the equation of motion, but not for precision gradiometry.  Expect order-10 Eotvos fidelity for the total tensor, with the diagonal static-Earth approximation setting the floor.

### Bottom Line

- **Vertical body tide:** ~200–1000 nGal is a defensible modeled-error range for a quiet inland site
- **Horizontal components:** uncertainty is substantially larger because the elastic response is not modeled with the required component/frequency dependence
- **Gradient tensor:** useful at the ~10 Eotvos level, not as a calibrated gravity-gradient product


## 8. Validation

### External References

**JPL Horizons (DE441).**  Moon and Sun positions from JPL's DE441 planetary ephemeris are used as truth values in the regression tests (2024--2025 epochs):

- Moon distance agreement: < 200 km
- Moon declination agreement: < 0.15°
- Sun distance agreement: < 10 000 km (~0.007%)
- Sun declination agreement: < 0.2°

**WGS84 normal gravity.**  The Somigliana formula is validated against the WGS84 defining parameters: $\gamma_e = 9.780\,325\,3141$ m/s$^2$ (equatorial) and $\gamma_p = 9.832\,184\,9378$ m/s$^2$ (polar), with agreement to better than 1 nGal at the boundary latitudes.

### Test Suite

The test suite (`tests/test_pytheas.py`) contains **113 tests** across **13 test classes**, including parameterized cases in the external-validation and property sections:

| Category | Tests | What is checked |
|----------|-------|-----------------|
| Constants | 6 | GM values, $\delta$, Somigliana self-consistency |
| Normal gravity | 5 | Equator/pole values, latitude trend, altitude dependence, free-air gradient |
| Time utilities | 4 | Julian date, GMST, timezone normalization |
| Coordinates | 8 | ECEF conversion, ENU orthonormality, right-handedness, axis vectors, ECI→ECEF rotation |
| Ephemeris | 4 | Sun/Moon distances, position ranges, external validation |
| Tidal acceleration | 3 | Order of magnitude, Moon/Sun ratio, center-of-mass limit |
| Full pipeline | 6 | Decomposition sums, gravimetric amplification, zenith projection |
| Timeseries | 6 | Length, variation, spectral periods, static constancy, input validation |
| Edge cases | 4 | South pole, date line, high altitude, negative longitude |
| LabFrame | 23 | Rotation vector, tensor symmetry/trace, EOM, Coriolis, timeseries, input validation |
| Dataclass API | 8 | Frozen behavior, `asdict` conversion, array immutability |
| External validation | 19 | JPL Horizons Moon/Sun (parameterized), WGS84 gravity |
| Physics properties | 17 | Scaling limits, ENU handedness (7 locations), gravity bounds, Coriolis $\perp$ velocity |

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

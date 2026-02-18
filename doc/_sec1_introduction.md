# 1. Introduction

## 1.1 The Question

Bolt a gravimeter to bedrock in Munich ($\varphi = 48.14^\circ$N, $\lambda = 11.58^\circ$E, $h = 500$ m above WGS84). The display settles near 9.807 m/s$^2$, but the last decimal places drift by $\sim$1 $\mu$m/s$^2$ over twelve hours, tracing a nearly sinusoidal oscillation that repeats -- not quite -- twice per day. Three contributions account for this signal:

**Static field.** Earth's gravitational attraction plus the centrifugal acceleration from rotation, determined by $\varphi$ and $h$. Time-independent; magnitude $\sim$9.8 m/s$^2$.
**Lunar tide.** Differential gravitational pull of the Moon between observer and Earth's center. The quadrupolar (degree-2) structure produces a semi-diurnal period of $\sim$12.4 h; peak amplitude $\sim$1.1 $\mu$m/s$^2$.
**Solar tide.** Same mechanism, but weaker by the $R^{-3}$ tidal scaling: peak $\sim$0.5 $\mu$m/s$^2$.

The measured acceleration is their sum, with tidal terms amplified by the solid Earth's elastic response ($\delta = 1.160808$) and projected onto the sensor axis. A superconducting gravimeter resolves parts-per-billion accelerations, so the tidal signal exceeds its noise floor by three orders of magnitude. The Pytheas library -- named after Pytheas of Massalia (c. 325 BC), the first to record a systematic Moon--tide connection -- implements the formulas derived in this document.

---

## 1.2 The Formula at a Glance

$$\boxed{g_{\text{total}}(t) = \gamma(\varphi, h)\,(\hat{\mathbf{e}}_U \cdot \hat{\mathbf{n}}) + \delta\left[\mathbf{a}_{\text{Moon}}(t) + \mathbf{a}_{\text{Sun}}(t)\right] \cdot \hat{\mathbf{n}}} \tag{1.1}$$

**$\gamma(\varphi, h)$: Normal gravity.** The static baseline set by Earth's mass and rotation. Ranges from 9.780 m/s$^2$ at the equator to 9.832 m/s$^2$ at the poles, decreasing $\sim$0.3 mGal/m with elevation. Computed via the Somigliana closed-form expression on WGS84 with a second-order free-air correction (Section 3).

**$\mathbf{a}_{\text{Moon}}(t)$, $\mathbf{a}_{\text{Sun}}(t)$: Tidal acceleration vectors.** The difference between each body's gravitational pull at the observer and at Earth's center:

$$\mathbf{a}_{\text{tidal}} = GM\left[\frac{\mathbf{R} - \mathbf{r}}{|\mathbf{R} - \mathbf{r}|^3} - \frac{\mathbf{R}}{R^3}\right] \tag{1.2}$$

with $\mathbf{R}$ the body's position, $\mathbf{r}$ the observer's, $GM$ the gravitational parameter (Section 6).

**$\delta$: Gravimetric factor.** The solid Earth deforms elastically under tidal stress, amplifying the observed signal. Direct gravity change, free-air displacement, and mass redistribution combine to give $\delta = 1.160808$ (Section 7).

**$\hat{\mathbf{n}}$: Measurement axis.** Unit vector along the sensor direction, expressed in the local ENU frame (Section 2). For a vertical gravimeter ($\hat{\mathbf{n}} = \hat{\mathbf{e}}_U$), equation (1.1) reduces to:

$$g_{\text{total}}(t) = \gamma(\varphi, h) + \delta\left[\mathbf{a}_{\text{Moon}}(t) + \mathbf{a}_{\text{Sun}}(t)\right] \cdot \hat{\mathbf{e}}_U \tag{1.3}$$

The decomposition $g = g_{\text{static}} + g_{\text{tidal}}$ has perturbation-theoretic structure: a time-independent baseline plus corrections suppressed by $\sim$$10^{-7}$, encoding the instantaneous lunar and solar positions.

![Figure 12: Anatomy of g(t) showing the decomposition into static normal gravity baseline at 9.807402 m/s², lunar tidal oscillation with 1.5 μm/s² peak-to-peak amplitude, and smaller solar tidal contribution, both amplified by the elastic factor δ = 1.16.](figures/fig12_anatomy.png)

Figure 12: (a) $g_{\text{static}} = 9.807402$ m/s$^2$ (Somigliana + free-air, Munich, 500 m). (b) Lunar tide, p-p $\sim$1.5 $\mu$m/s$^2$. (c) Solar tide, $\sim$$2\times$ smaller. Both include $\delta = 1.160808$.

---

## 1.3 A 48-Hour Portrait

![Figure 1: 48-hour time series of gravitational acceleration g(t) at Munich starting 2025-03-20 00:00 UTC, showing the total acceleration dominated by the static baseline (top panel) and the tidal residual revealing the M2 semi-diurnal lunar tide with 12.4-hour period plus solar contribution (bottom panel).](figures/fig01_g_timeseries.png)

**Figure 1.** 48 hours of $g(t)$ at Munich beginning 2025-03-20 00:00 UTC (equinox). Panel (a): total $g$ — tidal modulation is invisible at full scale. Panel (b): after removing the static baseline, the M2 semi-diurnal tide emerges with $\sim$12.4 h period; the lunar component (dashed) is $\sim$$2\times$ the solar (dotted); their weekly interference produces spring–neap modulation.

The 12.4 h period arises because the Moon advances $\sim 13^\circ$/day along the ecliptic, stretching the semi-diurnal period to half a lunar day ($\approx$12 h 25 min). Two daily maxima reflect the degree-2 tidal potential $\propto \cos^2\theta$, which peaks at both $\theta = 0$ (sub-lunar point) and $\theta = \pi$ (anti-lunar point). Atmospheric pressure loading, instrumental drift, and microseismic noise are not modeled (Section 1.5).

---

## 1.4 Roadmap

Evaluating equation (1.1) requires six ingredients, each developed in its own section:

**Section 2 -- Coordinates.** WGS84 ellipsoid-to-ECEF conversion ($N(\varphi)$), GMST rotation between Earth-fixed and inertial frames, and the ENU basis that defines $\hat{\mathbf{n}}$.
**Section 3 -- Normal Gravity.** $\gamma(\varphi, h)$ via the Somigliana formula with a second-order free-air expansion in $h/a$. Equator-to-pole variation $\sim$0.05 m/s$^2$.
**Section 4 -- Lunar Position.** Meeus ephemeris: $\sim$65 perturbation terms built from five fundamental arguments. Dominant corrections: equation of center ($6.29^\circ$), evection ($1.27^\circ$), variation ($0.66^\circ$). Accuracy $\sim$10 arcsec.
**Section 5 -- Solar Position.** Near-Keplerian orbit ($e \approx 0.017$) with a three-term equation of center from the $O(e^3)$ Kepler expansion. Accuracy $\sim$1 arcmin ($<$1 nGal).
**Section 6 -- Tidal Acceleration.** Exact formula (1.2) and $R^{-3}$ scaling law. The gradient approximation truncates at $\sim$2,800 nGal; Pytheas uses the exact form.
**Section 7 -- Elastic Response.** $\delta = 1 + h_2 - \frac{3}{2}k_2 = 1.1608$: displacement ($h_2 \approx 0.61$) amplifies the signal, while mass redistribution ($\frac{3}{2}k_2 \approx 0.45$) partially offsets it; the net amplification is 16%.

Each section includes worked examples and Pytheas code.

---

## 1.5 The Error Budget at a Glance

Overall model accuracy is $\sim$10 nGal inland and $\sim$50 nGal near coastlines (1 nGal = $10^{-11}$ m/s$^2$). Two scaling relations govern the error analysis:

$$a_{\text{tidal}} \sim \frac{2\,GM\,r}{R^3} \tag{1.4}$$

$$\delta = 1 + h_2 - \frac{3}{2}\,k_2 \tag{1.5}$$

with IERS 2010 Love numbers $h_2 = 0.6078$, $k_2 = 0.2980$, giving $\delta = 1.160808$ (Section 7). The cubic distance dependence in equation (1.4) explains why the nearby Moon dominates the far more massive Sun (Section 6).

**Dominant error:** treating $\delta$ as frequency-independent. Near the K1 constituent ($\sim$23.93 h), Free Core Nutation resonance alters $\delta$ by $\sim$1%, contributing $\sim$10 nGal. **Secondary errors:** lunar ephemeris uncertainty ($\sim$10 arcsec) introduces $\sim$1 nGal RMS; the solar ephemeris ($\sim$1 arcmin) contributes $<$1 nGal; normal gravity is exact on the ellipsoid to machine precision.

![Figure 13: Error budget ranking nine sources by magnitude on a logarithmic axis, showing modeled effects (ephemeris, normal gravity) below 1 nGal, frequency-independent Love number treatment at ~10 nGal, and unmodeled effects (ocean loading, atmospheric pressure) setting the accuracy floor at 10–50 nGal.](figures/fig13_error_budget.png)

Figure 13 ranks all error sources: modeled effects contribute $<$1 nGal; the frequency-independent $\delta$ approximation contributes $\sim$10 nGal; unmodeled ocean loading (1--50 nGal) and atmospheric pressure ($\sim$3 nGal/hPa) set a separate accuracy floor. The following effects are deliberately omitted, as each requires external data incompatible with the NumPy-only architecture:

- **Ocean tidal loading** (1--50 nGal): requires site-specific coefficients from FES2014.
- **Atmospheric pressure loading** ($\sim$3 nGal/hPa): requires real-time barometric data.
- **Polar motion** ($<$1 nGal): requires IERS orientation parameters.
- **Planetary tides** ($<$0.1 nGal): negligible magnitude.

Setting $\delta = 1$ (rigid Earth) introduces $\sim$18,000 nGal of error, so the elastic correction is essential. Superconducting gravimeters resolve $\sim$0.1--1 nGal; the 10 nGal floor of Pytheas represents a deliberate balance between accuracy and first-principles simplicity.

---

## 1.6 Notation and Conventions

**Units.** SI throughout. Angles in radians unless otherwise stated; time in UTC. Tidal signals are quoted in $\mu$Gal ($10^{-8}$ m/s$^2$) or nGal ($10^{-11}$ m/s$^2$). Typical signal amplitude is $O(100\;\mu\text{Gal})$; accuracy floor is 10 nGal.

**Notation.**

| Symbol | Meaning | Defined in |
|--------|---------|------------|
| $\varphi$, $\lambda$, $h$ | Geodetic latitude, longitude, altitude | Section 2 |
| $\gamma(\varphi, h)$ | Normal gravity (Somigliana + free-air) | Section 3 |
| $\mathbf{R}_{\text{Moon}}(t)$, $\mathbf{R}_{\text{Sun}}(t)$ | ECEF position vectors of Moon, Sun | Sections 4, 5 |
| $\mathbf{r}$ | Observer's ECEF position vector | Section 2 |
| $\mathbf{a}_{\text{tidal}}$ | Tidal acceleration vector | Section 6 |
| $\delta$ | Gravimetric factor (= 1.1608) | Section 7 |
| $\hat{\mathbf{n}}$ | Measurement axis unit vector | Section 2 |
| $\hat{\mathbf{e}}_U$, $\hat{\mathbf{e}}_N$, $\hat{\mathbf{e}}_E$ | Local Up, North, East unit vectors | Section 2 |
| GMST | Greenwich Mean Sidereal Time | Section 2 |

**Code.** All implementations use Python/NumPy; variable names mirror the mathematical notation (`H2`, `K2`, `DELTA_GRAV`). Top-level API:

```python
from datetime import datetime
from pytheas import compute_g

result = compute_g(
    dt=datetime(2025, 3, 20, 12, 0, 0),  # UTC time
    lat_deg=48.14,                         # Munich latitude
    lon_deg=11.58,                         # Munich longitude
    alt_m=500.0,                           # altitude above WGS84
)

print(f"g_total  = {result['g_total']:.6f} m/s²")
print(f"g_static = {result['g_static']:.6f} m/s²")
print(f"g_tidal  = {result['g_tidal'] * 1e6:.3f} μm/s²")
```

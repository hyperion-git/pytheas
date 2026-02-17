# 1. Introduction: What Does a Gravimeter Read?

## 1.1 The Question

Bolt a gravimeter to bedrock in Munich -- geodetic latitude 48.14 degrees north, longitude 11.58 degrees east, altitude 500 meters above the GRS80 ellipsoid -- and read off the gravitational acceleration. The display settles near 9.807 m/s$^2$, exactly as every introductory physics course promises. But watch the last few decimal places. They will not hold still. Over the next twelve hours, the reading drifts by roughly a micron per second squared, tracing a smooth, nearly sinusoidal oscillation that repeats -- not quite -- twice per day. What is happening?

The answer involves three distinct contributors, each operating on a different timescale.

**Static field.** The gravitational attraction of the Earth itself, modified by the centrifugal acceleration due to Earth's rotation. This depends on the observer's latitude $\varphi$ (through the centrifugal force and the ellipsoidal shape) and altitude $h$ (through the inverse-square falloff). For a fixed instrument, this contribution is constant in time. Its magnitude is approximately 9.8 m/s$^2$.

**Lunar tide.** The Moon's gravitational pull is not uniform across the Earth. The pull at the observer's location differs from the pull at Earth's center, producing a residual acceleration -- the tidal acceleration -- that varies as the Moon moves through its orbit. The Moon completes an orbit in approximately 27.3 days, but because the tidal pattern has a quadrupolar (degree-2) structure, the dominant period is approximately 12.4 hours (semi-diurnal). The peak magnitude is approximately 1.1 $\mu$m/s$^2$.

**Solar tide.** The same mechanism applies to the Sun. Although the Sun is vastly more massive than the Moon, it is also much farther away. The tidal acceleration depends on the gradient of the gravitational field (scaling as $1/R^3$, not $1/R^2$), so the Sun's contribution is smaller: approximately 0.5 $\mu$m/s$^2$ at peak.

The total measured acceleration is the sum of these three contributions, with the tidal terms amplified by the elastic response of the solid Earth and projected onto the sensor's measurement axis.

To put the scale in perspective: imagine standing on a bathroom scale while the Moon passes overhead. The Moon's tidal pull changes your apparent weight by roughly 0.1 microgram out of 70 kilograms. You will never feel it. But a superconducting gravimeter -- a niobium sphere levitated in a persistent-current magnetic trap, cooled to 4 kelvin -- resolves accelerations to parts per billion. The tidal signal is a thousand times larger than its noise floor. What we are about to compute is not a theoretical curiosity; it is a signal that is routinely measured, catalogued, and used to calibrate some of the most sensitive instruments in experimental physics.

The library that implements the formulas derived in this document is named after Pytheas of Massalia (circa 325 BC), the Greek navigator and astronomer who sailed from the Mediterranean to Britain and, in doing so, became the first person known to have recorded a systematic connection between the Moon and the ocean tides. Twenty-three centuries later, we understand the mechanism in complete quantitative detail. This document derives it from first principles.

---

## 1.2 The Formula at a Glance

The complete gravitational acceleration measured by a surface gravimeter, as a function of time, is:

$$\boxed{g_{\text{total}}(t) = \gamma(\varphi, h)\,(\hat{\mathbf{e}}_U \cdot \hat{\mathbf{n}}) + \delta\left[\mathbf{a}_{\text{moon}}(t) + \mathbf{a}_{\text{sun}}(t)\right] \cdot \hat{\mathbf{n}}} \tag{1.1}$$

Every symbol earns its place. Let us walk through them one at a time.

**$\gamma(\varphi, h)$: Normal gravity.** This is what gravity would be if the Moon and Sun did not exist -- the static baseline produced by the Earth's own mass distribution and rotation. It depends on geodetic latitude $\varphi$ (ranging from 9.780 m/s$^2$ at the equator to 9.832 m/s$^2$ at the poles) and altitude $h$ above the reference ellipsoid (decreasing by approximately 0.3 mGal per meter of elevation). The formula is given by Somigliana's closed-form expression on the GRS80 ellipsoid, extended to altitude by a second-order free-air correction. We derive it in Section 3.

**$\mathbf{a}_{\text{moon}}(t)$ and $\mathbf{a}_{\text{sun}}(t)$: Tidal acceleration vectors.** These are the residual gravitational accelerations produced by the Moon and Sun -- the difference between each body's pull at the observer's location and its pull at Earth's center. The Moon contributes approximately 1.1 $\mu$m/s$^2$ at peak; the Sun approximately 0.5 $\mu$m/s$^2$. Both are three-dimensional vectors in Earth-Centered Earth-Fixed (ECEF) coordinates, computed from the exact Newtonian tidal formula:

$$\mathbf{a}_{\text{tidal}} = GM\left[\frac{\mathbf{R} - \mathbf{r}}{|\mathbf{R} - \mathbf{r}|^3} - \frac{\mathbf{R}}{R^3}\right] \tag{1.2}$$

where $\mathbf{R}$ is the position of the external body (Moon or Sun), $\mathbf{r}$ is the observer's position, and $GM$ is the body's gravitational parameter. We derive this formula from first principles in Section 6.

**$\delta$: The gravimetric factor.** The Earth is not rigid. The same tidal forces that tug on the gravimeter also deform the entire planet, raising and lowering the crust by tens of centimeters. This elastic response amplifies the measured tidal signal by approximately 16%. The gravimetric factor $\delta \approx 1.16$ is a dimensionless number that encodes three competing effects -- direct tidal gravity, free-air displacement, and mass redistribution -- into a single multiplicative correction. We derive it in Section 7.

**$\hat{\mathbf{n}}$: The measurement axis.** A gravimeter measures a scalar, not a vector: it registers the component of acceleration along its sensitive axis. The unit vector $\hat{\mathbf{n}}$ specifies the direction in which the sensor is pointed, parametrized by a zenith angle and an azimuth in the local East-North-Up frame. We define it in Section 2.

For the common case of a vertical gravimeter ($\hat{\mathbf{n}} = \hat{\mathbf{e}}_U$, the local "up" direction), the formula simplifies to:

$$g_{\text{total}}(t) = \gamma(\varphi, h) + \delta\left[\mathbf{a}_{\text{moon}}(t) + \mathbf{a}_{\text{sun}}(t)\right] \cdot \hat{\mathbf{e}}_U \tag{1.3}$$

A quick dimensional check confirms the structure: $\gamma$ carries units of m/s$^2$; each $\mathbf{a}_{\text{tidal}}$ is also m/s$^2$; $\delta$ is dimensionless; $\hat{\mathbf{n}}$ is a dimensionless unit vector; and the dot product yields a scalar. Everything checks out.

The decomposition $g(t) = g_{\text{static}} + g_{\text{tidal}}$ has the same additive structure as perturbation theory in general: a large, time-independent baseline plus small, time-varying corrections. The tidal terms are the "perturbation," suppressed by a factor of roughly $10^{-7}$ relative to the static part. Yet it is precisely these tiny corrections that encode the positions and masses of the Moon and Sun, and that make gravimetry a window into celestial mechanics.

![Figure 12: Anatomy of g(t) showing the decomposition into static normal gravity baseline at 9.807402 m/s², lunar tidal oscillation with 1.5 μm/s² peak-to-peak amplitude, and smaller solar tidal contribution, both amplified by the elastic factor δ = 1.16.](figures/fig12_anatomy.png)

Figure 12 makes this additive structure visually explicit. Panel (a) displays $g_{\text{static}}$ as a horizontal line at 9.807402 m/s$^2$ -- the normal gravity baseline computed from Somigliana's formula plus the free-air correction for Munich's 500-meter altitude. Panel (b) shows the lunar tidal contribution oscillating in $\mu$m/s$^2$, with a peak-to-peak amplitude of approximately 1.5 $\mu$m/s$^2$. Panel (c) adds the solar tidal contribution, smaller by roughly a factor of two. Both tidal panels include the elastic amplification factor $\delta = 1.16$ already applied.

---

## 1.3 A 48-Hour Portrait

What does the formula predict in practice? Figure 1 answers this question with a 48-hour snapshot of $g(t)$ at Munich, starting at 2025-03-20 00:00 UTC (near the March equinox, chosen so that the Sun crosses the celestial equator and the solar tide takes a particularly symmetric form).

Panel (a) displays the total gravitational acceleration in m/s$^2$. At this scale -- the full value of $g$ -- the tidal modulation is invisible. The reading hovers near 9.8074 m/s$^2$, and even a careful eye cannot distinguish the tidal oscillation from the thickness of the plotted line. The static baseline dominates by seven orders of magnitude.

Panel (b) tells a different story. Here we subtract the static baseline and zoom in by a factor of $10^7$, revealing the tidal residual in $\mu$m/s$^2$. The dominant rhythm is unmistakable: a nearly sinusoidal oscillation with a period of approximately 12.4 hours. This is the M2 lunar semi-diurnal tide, the single largest constituent of the tidal signal. The plot separates the total tidal signal (solid line) from its lunar (dashed) and solar (dotted) components; the lunar contribution is roughly twice the solar one, and the two interfere constructively and destructively on a timescale of weeks, producing the familiar spring-tide and neap-tide modulation.

Why 12.4 hours and not 12? Because the Moon is not stationary. It advances approximately 13 degrees per day along its orbit, so the tidal pattern -- which is locked to the Moon's position, not the Sun's -- slips relative to the solar day. The semi-diurnal period is half a lunar day, or about 12 hours and 25 minutes. Two tidal bulges pass under the observer each lunar day: one on the side of the Earth nearest the Moon, one on the far side.

And why two bulges, not one? The tidal field is quadrupolar: it has spherical harmonic degree 2. The near side of the Earth is pulled harder toward the Moon than the center; the far side is pulled less hard. Both extremes produce a positive vertical tidal acceleration relative to the freely falling center. The angular dependence of the degree-2 tidal potential goes as $\cos^2\theta$ (where $\theta$ is the angle between the observer and the sub-lunar point), which has two maxima -- at $\theta = 0$ (directly beneath the Moon) and $\theta = \pi$ (directly opposite). The "two high tides per day" that oceanographers have documented since antiquity follow directly from this quadrupolar symmetry.

In a real gravimeter record, the 48-hour window would also show atmospheric pressure effects (a few nGal per hectopascal of pressure change), instrumental drift, and microseismic noise. We are deliberately not modeling any of these -- our formula captures only the gravitational signal from the static Earth and the two dominant tidal sources. The effects we omit, and why, are catalogued in Subsection 1.5.

![Figure 1: 48-hour timeseries of gravitational acceleration g(t) at Munich starting 2025-03-20 00:00 UTC, showing the total acceleration dominated by the static baseline (top panel) and the tidal residual revealing the M2 semi-diurnal lunar tide with 12.4-hour period plus solar contribution (bottom panel).](figures/fig01_g_timeseries.png)

---

## 1.4 Roadmap

To evaluate $g_{\text{total}}(t)$ from equation (1.1), we need six ingredients: a coordinate system, the static baseline, the Moon's position, the Sun's position, the tidal acceleration formula, and the elastic correction factor. Each of these occupies one section of this document. Here is what each section derives, and why.

**Section 2 -- Setting Up Coordinates.** Before we can compute anything, we need to convert a GPS position (latitude, longitude, altitude) into the three-dimensional vectors that the tidal formula demands. This section derives the GRS80 ellipsoid-to-Cartesian conversion -- including the prime vertical radius of curvature $N(\varphi)$ that accounts for Earth's oblate shape -- the Greenwich Mean Sidereal Time (GMST) rotation that connects the Earth-fixed frame (where the observer sits) to the inertial frame (where celestial bodies are defined), and the East-North-Up basis vectors that define the local measurement axis $\hat{\mathbf{n}}$. The coordinate pipeline is the plumbing of the entire calculation: invisible when it works, catastrophic when it doesn't.

**Section 3 -- Gravity on a Quiet Earth.** What would gravity be at our location if the Moon and Sun did not exist? This section derives the static baseline $\gamma(\varphi, h)$, starting from Newton's shell theorem for a spherical Earth, adding the centrifugal correction from rotation (approximately 0.034 m/s$^2$ at the equator), and incorporating the oblate ellipsoidal shape to arrive at Somigliana's closed-form formula. The free-air correction extends the result to altitude $h$ via a second-order Taylor expansion in $h/a$, where $a$ is Earth's equatorial radius. The equator-to-pole gravity variation of approximately 0.05 m/s$^2$ -- half a percent of $g$ -- emerges naturally from the interplay of geometry and rotation.

**Section 4 -- Where is the Moon?** The Moon's position at any given instant is the single most important input to the tidal calculation, and it is also the hardest to compute. The three-body problem (Earth, Moon, Sun) has no closed-form solution, so the Moon's trajectory is represented as a sum of approximately 65 sinusoidal perturbation terms, each with an amplitude, frequency, and phase built from five fundamental angles that track the Moon's mean longitude, elongation from the Sun, solar and lunar anomalies, and argument of latitude. The largest terms -- the equation of the center (6.29 degrees), evection (1.27 degrees), and variation (0.66 degrees) -- capture the dominant physics of the elliptical orbit and solar perturbations. This section derives the Meeus analytical ephemeris that produces the Moon's ECEF position to approximately 10 arcsecond accuracy.

**Section 5 -- Where is the Sun?** The Sun is mercifully easier than the Moon. Earth's orbit around the Sun is nearly Keplerian -- the eccentricity is only $e \approx 0.017$, and planetary perturbations are negligible at our accuracy level -- so the Sun's ecliptic longitude is captured by just three terms in the equation of center: $1.915^\circ \sin M + 0.020^\circ \sin 2M + 0.0003^\circ \sin 3M$, derived from a perturbative expansion of Kepler's equation to third order in eccentricity. The angular accuracy of approximately 1 arcminute translates to less than 1 nGal error in tidal gravity. This section derives the Meeus low-precision solar ephemeris and explains why the Sun's contribution, despite being only 46% the strength of the Moon's, is worth including.

**Section 6 -- The Tidal Acceleration.** This is the heart of the model. Starting from Newton's law of gravitation, this section derives the exact tidal acceleration formula (equation 1.2) as the difference between the external body's gravitational pull at the observer and at Earth's center. The $R^{-3}$ scaling law -- which explains why the closer Moon dominates the more massive Sun -- emerges from a Taylor expansion of the exact formula. The gradient approximation, widely used in textbooks, is shown to introduce approximately 2,800 nGal truncation error for the Moon -- unacceptable for precision work -- which is why Pytheas uses the exact formula at negligible additional computational cost.

**Section 7 -- The Earth Fights Back.** The solid Earth deforms under tidal stress, raising and lowering the crust by up to 30 centimeters. This section derives the gravimetric factor $\delta = 1 + h_2 - \frac{3}{2}k_2 = 1.16$ from the three competing effects: the direct tidal gravity on a rigid Earth (baseline), the free-air displacement as the crust rides upward on the tidal bulge (amplifying the signal by $h_2 \approx 0.61$), and the mass redistribution beneath the bulge (reducing the signal by $\frac{3}{2}k_2 \approx 0.45$). The net 16% amplification is the result of displacement winning over mass redistribution.

Each section includes worked numerical examples and code snippets from the Pytheas library, so the reader can verify every derived formula computationally. The error budget (Figure 13), presented in summary in the next subsection, quantifies the accuracy of each ingredient and identifies where the model's limitations lie.

---

## 1.5 The Error Budget at a Glance

How good is this model? The headline number: approximately 10 nGal accuracy for inland sites, approximately 50 nGal for coastal sites. One nanoGal is $10^{-11}$ m/s$^2$ -- roughly one part in $10^{12}$ of the gravitational acceleration at Earth's surface. An accuracy of 10 nGal means we are predicting $g(t)$ to approximately one part in $10^{10}$ of $g$.

To appreciate what limits the accuracy, we need two formulas that the later sections will derive in full. The first is the tidal magnitude scaling:

$$a_{\text{tidal}} \sim \frac{2\,GM\,r}{R^3} \tag{1.4}$$

where $GM$ is the gravitational parameter of the external body, $r$ is the Earth's radius, and $R$ is the distance to the body. This expression -- derived from a first-order Taylor expansion of equation (1.2) -- reveals the critical $R^{-3}$ dependence that makes the closer Moon dominate the more massive Sun (Section 6).

The second is the gravimetric factor:

$$\delta = 1 + h_2 - \frac{3}{2}\,k_2 \tag{1.5}$$

where $h_2 = 0.6078$ and $k_2 = 0.2980$ are the IERS 2010 Love numbers for the degree-2 body tide (Section 7). Its numerical value is $\delta = 1.1608$.

The single largest error source in the model is the frequency-independent treatment of $\delta$. The Love numbers $h_2$ and $k_2$ are, strictly speaking, functions of tidal frequency. Near the K1 tidal constituent (period approximately 23.93 hours), the Free Core Nutation -- a normal mode involving differential precession of the fluid outer core and solid mantle -- produces a resonance that alters $\delta$ by up to approximately 1%. Applied to the total tidal signal of approximately 1,000 nGal at the K1 frequency, this produces an error of approximately 10 nGal. This is the dominant limitation of the single-$\delta$ approach.

The ephemeris errors are secondary. The Meeus lunar ephemeris achieves approximately 10 arcsecond accuracy in ecliptic longitude and sub-arcsecond accuracy in parallax; when propagated through the tidal formula and averaged over the semi-diurnal cycle, the RMS error is approximately 1 nGal. The solar ephemeris, at approximately 1 arcminute accuracy, contributes less than 1 nGal. The normal gravity formula (Somigliana with second-order free-air correction) is exact on the reference ellipsoid to machine precision and accurate to better than 1 nGal for altitudes up to several kilometers.

![Figure 13: Error budget ranking nine sources by magnitude on a logarithmic axis, showing modeled effects (ephemeris, normal gravity) below 1 nGal, frequency-independent Love number treatment at ~10 nGal, and unmodeled effects (ocean loading, atmospheric pressure) setting the accuracy floor at 10–50 nGal.](figures/fig13_error_budget.png)

Figure 13 ranks all error sources visually. The horizontal bar chart displays nine contributions on a logarithmic axis, color-coded by modeling status: the modeled sources (ephemeris, normal gravity) cluster below 1 nGal; the frequency-independent Love number treatment sits at approximately 10 nGal; and the unmodeled effects -- ocean tidal loading (1--50 nGal depending on proximity to the coast) and atmospheric pressure loading (approximately 3 nGal per hectopascal) -- set a separate accuracy floor for sites near the ocean or exposed to large barometric swings.

Several physical effects are deliberately omitted from the formula:

- **Ocean tidal loading** (1--50 nGal at coastal sites): moving water masses gravitationally attract nearby land gravimeters; requires site-specific loading coefficients from ocean tide models such as FES2014.
- **Atmospheric pressure loading** (approximately 3 nGal/hPa): pressure changes redistribute air mass and flex the crust; requires real-time barometric data.
- **Polar motion** ($<$1 nGal): wandering of Earth's rotation axis (Chandler wobble) alters the local centrifugal acceleration; requires Earth orientation parameters from IERS bulletins.
- **Planetary tides** ($<$0.1 nGal): Jupiter and Venus produce tidal accelerations at Earth, but at levels completely negligible for this model.

These omissions are deliberate design choices, not oversights. Each could be added as a correction layer, but each requires external data (ocean models, barometric readings, orientation parameters) that would compromise the library's self-contained, numpy-only architecture.

As a limiting case, consider what happens if we set $\delta = 1$ -- treating the Earth as perfectly rigid. The error from ignoring the elastic response would be approximately 16% of the tidal signal, or roughly 18,000 nGal. That is three orders of magnitude larger than the 10 nGal accuracy floor. The elastic correction is not a refinement; it is essential.

For comparison, a superconducting gravimeter has a noise floor of approximately 0.1--1 nGal. Our model is accurate enough to serve as the dominant theoretical prediction for tidal gravity at inland sites, but not accurate enough for the highest-precision tidal analysis, which requires frequency-dependent Love numbers, ocean loading corrections, and atmospheric admittance models. The 10 nGal floor is a deliberate sweet spot: accurate enough for most applications, simple enough to derive from first principles in a self-contained document.

---

## 1.6 Prerequisites, Notation, and Conventions

**Prerequisites.** This document assumes familiarity with vector calculus (dot products, cross products, coordinate transformations between Cartesian frames), Newtonian mechanics (gravitational force law, centrifugal acceleration in rotating frames), and Taylor series expansions. No quantum mechanics, no general relativity, no signal processing. Everything that follows is classical Newtonian gravity in Euclidean geometry. A reader comfortable with upper-division undergraduate physics or first-year graduate coursework will find the derivations self-contained.

**Units.** SI units throughout. Angles are in radians unless stated otherwise; time is in UTC. Gravity perturbations are frequently quoted in the conventional unit of the Gal, named after Galileo Galilei: 1 Gal = 1 cm/s$^2$. In practice, tidal signals are measured in microGal ($1\;\mu\text{Gal} = 10^{-8}$ m/s$^2$) or nanoGal ($1\;\text{nGal} = 10^{-11}$ m/s$^2$ = 0.01 nm/s$^2$). The tidal signal we are computing is of order 100 $\mu$Gal = 100,000 nGal. The model's accuracy floor of 10 nGal corresponds to $10^{-11}$ m/s$^2$ -- roughly the weight of a single grain of pollen.

**Notation.** The following table collects the principal symbols used throughout the document, with the section where each is defined:

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

**Conventions.** Code snippets use Python with NumPy; variable names mirror the mathematical notation (e.g., `H2`, `K2`, `DELTA_GRAV`, `normal_gravity`). The Pytheas library is the reference implementation of every formula in this document. Each section includes code showing the corresponding Pytheas function, so the reader can verify every derivation numerically. As a preview of the interface, here is the top-level API call that evaluates equation (1.1):

```python
from datetime import datetime
from pytheas import compute_g

result = compute_g(
    dt=datetime(2025, 3, 20, 12, 0, 0),  # UTC time
    lat_deg=48.14,                         # Munich latitude
    lon_deg=11.58,                         # Munich longitude
    alt_m=500.0,                           # altitude above GRS80
)

print(f"g_total  = {result['g_total']:.6f} m/s²")
print(f"g_static = {result['g_static']:.6f} m/s²")
print(f"g_tidal  = {result['g_tidal'] * 1e6:.3f} μm/s²")
```

Every formula derived in Sections 2--7 has a corresponding function in the library. The derivations explain the physics; the code confirms the numbers.

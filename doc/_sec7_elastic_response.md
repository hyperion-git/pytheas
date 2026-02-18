# 7. Elastic Response of the Solid Earth

The preceding sections treated the Earth as rigid. In reality, the Earth is an elastic solid that deforms under tidal stress, amplifying the gravimeter's tidal signal by about 16%. This section derives the correction factor.


## 7.1 Body Tides and Elastic Deformation

### The Body Tide

The tidal force acts on every mass element of the Earth, not just a test mass at the surface. Its quadrupolar pattern stretches the Earth along the Earth--Moon line and compresses it transversely, raising a solid-body bulge -- the **body tide**. The crust rises and falls by ~30 cm (lunar) and ~15 cm (solar) twice daily.

### Why the Deformation Is Imperceptible

The body tide has a wavelength comparable to the Earth's circumference, so over human-scale distances the displacement is uniform and imperceptible. A gravimeter, however, measures gravitational acceleration at a fixed coordinate point (latitude, longitude, ellipsoidal height). When the solid Earth deforms, two quantities change at that point, both affecting the reading.


## 7.2 Two Competing Effects on Measured Gravity

A gravimeter bolted to the crust rides with the deforming Earth. Two competing effects alter the measured acceleration:

**Effect 1: Free-air effect.**
The tidal bulge lifts the gravimeter radially outward by $\Delta r$, reducing $g$ by the same mechanism as the free-air decrease of gravity with altitude.

**Effect 2: Mass redistribution.**
Deformation redistributes mass, altering the gravitational potential -- and hence the acceleration -- at the measurement point.

These two effects partially oppose each other. Quantifying them requires the **Love numbers**.


## 7.3 Love Numbers

The Love numbers (A. Love, 1911) are dimensionless parameters that characterize the Earth's deformational response to an applied tidal potential.

### The Tidal Potential

The tidal potential is dominated by its degree-2 (quadrupolar) spherical harmonic component $V_2$, whose relevant properties are:

- $V_2$ is a degree-2 spherical harmonic in the angular coordinates
- At the Earth's surface ($r = R$), $V_2$ depends on the zenith angle of the tide-raising body and its distance
- The radial tidal acceleration on a rigid Earth equals $-\partial V_2 / \partial r$

All results below are expressed as ratios relative to $V_2$; its explicit form is not needed.

### The Displacement Love Number $h_2$

$h_2$ quantifies the radial surface displacement in response to $V_2$:

$$\Delta r = h_2 \, \frac{V_2}{g}$$

A perfectly fluid body in hydrostatic equilibrium would have $\Delta r = V_2/g$, so $h_2$ measures the fraction of this fluid limit attained:
- $h_2 = 1$: perfectly fluid body (full hydrostatic adjustment)
- $h_2 = 0$: perfectly rigid body (no deformation)
- Real Earth: $h_2 \approx 0.61$ (61% of the fluid limit)

### The Potential Love Number $k_2$

The additional gravitational potential at the surface due to the deformation is:

$$V_{\text{deformation}} = k_2 \, V_2$$

- $k_2 = 0$: rigid body
- Real Earth: $k_2 \approx 0.30$

### Why Only Degree 2?

Higher-degree terms scale as $(R/d)^l$ with $R/d \sim 1/60$ (Moon) and $R/d \sim 4 \times 10^{-5}$ (Sun), so they are strongly suppressed. The degree-3 lunar contribution amounts to ~1.7% of degree 2 -- negligible at our target accuracy.


## 7.4 Deriving the Gravimetric Factor $\delta$

The measured tidal gravity on an elastic Earth combines three contributions: (1) rigid-Earth tidal gravity, (2) the free-air correction from surface displacement, and (3) the potential change from mass redistribution.

### Step 1: Rigid-Earth Tidal Gravity

The radial gravity perturbation is $\Delta g_{\text{rigid}} = -\partial V_2/\partial r$. The tidal potential from a distant body takes the form $V_2 = -\frac{GM}{r_M}\left(\frac{r}{r_M}\right)^2 P_2(\cos\theta)$ at the observer's position $r$. Because $V_2$ is proportional to $r^2$ (from the degree-2 Legendre polynomial), its radial derivative satisfies $\partial V_2/\partial r = 2V_2/r$. Evaluating at the Earth's surface $r = R$:

$$\boxed{\Delta g_{\text{rigid}} = -\frac{2 \, V_2}{R}}$$

### Step 2: Gravity Change from Radial Displacement (Free-Air Effect)

The gravimeter is displaced outward by $\Delta r = h_2 V_2/g$. For a spherical Earth, the free-air gradient is $\partial g/\partial r = -2g/R$ (from $g = GM/r^2$). The resulting gravity change is:

$$\Delta g_{\text{disp}} = \frac{\partial g}{\partial r} \, \Delta r = -\frac{2g}{R} \cdot \frac{h_2 V_2}{g}$$

$$\boxed{\Delta g_{\text{disp}} = -\frac{2 \, h_2 \, V_2}{R}}$$

The sign matches $\Delta g_{\text{rigid}}$: as the gravimeter rides the bulge upward, it moves away from Earth's center, so $g$ decreases further. Displacement therefore *amplifies* $|\Delta g|$.

### Step 3: Gravity Change from Mass Redistribution (Potential Effect)

The deformation adds a potential $k_2 V_2$ at the surface. Since the rigid-Earth tidal potential was already included in Step 1, only the deformation-generated part $k_2 V_2$ remains. This potential is generated by internal mass redistribution and decays as $r^{-(l+1)} = r^{-3}$ outside the surface:

$$\frac{\partial (k_2 V_2)}{\partial r}\bigg|_{\text{ext}} = -3 \, \frac{k_2 V_2}{R}$$

The gravity perturbation $\Delta g_{\text{pot}} = -\partial(k_2 V_2)/\partial r$:

$$\boxed{\Delta g_{\text{pot}} = +\frac{3 \, k_2 \, V_2}{R}}$$

This term is *positive* -- opposite in sign to $\Delta g_{\text{rigid}}$. The extra mass concentrated beneath the gravimeter by the bulge partially counteracts the tidal stretching, *reducing* $|\Delta g|$.

### Step 4: Combining All Three Contributions

$$\Delta g_{\text{elastic}} = \Delta g_{\text{rigid}} + \Delta g_{\text{disp}} + \Delta g_{\text{pot}}$$

$$= -\frac{2V_2}{R} - \frac{2 h_2 V_2}{R} + \frac{3 k_2 V_2}{R} = -\frac{2V_2}{R}\left(1 + h_2 - \frac{3}{2} k_2\right)$$

$$\Delta g_{\text{elastic}} = \Delta g_{\text{rigid}} \cdot \underbrace{\left(1 + h_2 - \frac{3}{2}\, k_2\right)}_{\delta}$$

The gravimetric factor is:

$$\boxed{\delta = 1 + h_2 - \frac{3}{2}\, k_2}$$

### Summary of the Three Contributions

| Term | Origin | Expression | Effect on $\lvert\Delta g\rvert$ |
|------|--------|------------|----------------------|
| 1 | Rigid-Earth tidal gravity | $-2V_2/R$ | Baseline |
| $+h_2$ | Displacement (free-air) | $-2h_2 V_2/R$ | Amplifies |
| $-\frac{3}{2}k_2$ | Mass redistribution | $+3k_2 V_2/R$ | Reduces |

Displacement dominates because the gravimeter rides the tidal bulge upward. Mass redistribution partially cancels this effect, since the deformed Earth concentrates mass beneath the bulge.


## 7.5 Numerical Values

### IERS 2010 Love Numbers

The IERS 2010 Conventions give nominal degree-2 Love numbers for a spherically symmetric, non-rotating, elastic, isotropic (SNREI) Earth model:

$$h_2 = 0.6078, \qquad k_2 = 0.2980$$

### The Gravimetric Factor

$$\delta = 1 + 0.6078 - 1.5 \times 0.2980$$

$$\delta = 1 + 0.6078 - 0.4470$$

$$\boxed{\delta = 1.160808}$$

### Physical Interpretation

The tidal signal is **16% larger** than on a rigid Earth:

- The rigid-Earth baseline: **1.0000** (100%)
- The displacement (free-air) amplification: **+0.6078** (+60.8%)
- The mass redistribution reduction: **-0.4470** (-44.7%)
- **Net amplification: +0.1608 (+16.1%)**

Although displacement dominates, mass redistribution substantially counteracts it, leaving a net amplification far smaller than $h_2$ alone would suggest. Figure 11 overlays the rigid and elastic predictions; the shaded region highlights the 16% amplification.

![Figure 11: Rigid Earth vs. elastic Earth tidal gravity signal. The shaded region between the curves shows the 16% amplification from the elastic response (gravimetric factor δ = 1.1608).](figures/fig11_rigid_vs_elastic.png)

### How Pytheas Applies This

Pytheas applies the gravimetric factor as a multiplicative correction to the tidal acceleration:

```python
H2 = 0.6078
K2 = 0.2980
DELTA_GRAV = 1.0 + H2 - 1.5 * K2   # = 1.1608

a_moon = DELTA_GRAV * tidal_acceleration(r, moon_position_ecef(dt), GM_MOON)
a_sun  = DELTA_GRAV * tidal_acceleration(r, sun_position_ecef(dt),  GM_SUN)
```

The `tidal_acceleration` function computes the rigid-Earth value; multiplying by $\delta$ gives the measured tidal gravity change on the elastic Earth.


## 7.6 Omitted Effects and Their Magnitudes

The derivation above assumes a single, frequency-independent gravimetric factor $\delta$. Several effects have been omitted.

### Frequency Dependence of Love Numbers

The Love numbers depend on the forcing frequency. The values 0.6078 and 0.2980 are appropriate for semi-diurnal tides (~12 h period). At longer periods the mantle flows viscously, increasing $h_2$. At diurnal periods, a resonance significantly alters both Love numbers.

### The Free Core Nutation Resonance

Near the K1 tidal frequency (~23.93 h period), tidal forcing resonates with the **Free Core Nutation (FCN)** -- a normal mode in which the fluid outer core and solid mantle precess around slightly different axes. Close to this resonance, the Love numbers vary rapidly and $\delta$ can deviate from 1.1608 by up to ~1% (~10 nGal). This constitutes the **single largest error source** in Pytheas. A proper correction would require frequency-dependent Love numbers for each tidal constituent.

### Higher-Degree Love Numbers

The degree-3 lunar tidal component (~1.7% of degree 2) involves Love numbers $h_3$ and $k_3$. Its effect on gravity is negligible at the 10 nGal level.

### Latitude Dependence

The Love numbers assume a spherically symmetric Earth. Lateral variations in elastic properties cause the effective Love numbers to vary with location (typically $<$1%), an effect usually subsumed into ocean loading corrections at coastal sites.

### What These Omissions Cost

| Omission | Typical error | Notes |
|----------|--------------|-------|
| Frequency-independent $\delta$ (no FCN) | $\sim$10 nGal | Largest near K1 frequency |
| No degree-3 Love numbers | $<$1 nGal | Suppressed by $(R/d)$ |
| No lateral variation | $<$1 nGal | Subsumed into ocean loading |

For Pytheas's target accuracy of ~10 nGal at inland sites, a single gravimetric factor $\delta = 1.160808$ is adequate. The FCN resonance remains the dominant limitation of this approach.

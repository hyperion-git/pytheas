# Section 7: The Earth Fights Back (Elastic Response)

So far, we have computed the tidal acceleration that the Moon and Sun exert at a point on Earth's surface, treating the Earth as a perfectly rigid body. But the Earth is not rigid. It is an elastic solid that deforms under stress -- and the tidal forces from the Moon and Sun are exactly the kind of stress that makes it flex. This deformation changes what a gravimeter actually measures, amplifying the tidal signal by about 16%. This section derives the correction factor.


## 7.1 The Earth is not rigid

### The body tide

Picture the Moon pulling on the Earth. In Sections 5--6 we focused on the differential gravitational acceleration (the tidal acceleration) acting on a test mass at the surface. But the same tidal force field acts on every piece of the Earth itself. The result: the entire solid Earth deforms.

The tidal force field has a quadrupolar pattern -- it stretches the Earth along the Earth--Moon line and squeezes it perpendicular to that line. The solid Earth responds by bulging toward and away from the Moon, just as the ocean does, but with a much smaller amplitude. The crust rises and falls by roughly 30 cm twice per day in response to the lunar tide (and roughly 15 cm for the solar tide). This is the **body tide** or **solid Earth tide**.

### Why you don't notice

You might wonder: if the ground moves 30 cm, why don't we feel it? The answer is that everything around you -- the floor, the walls, the table, the nearby landscape -- all move together. The body tide has a wavelength comparable to the Earth's circumference; over any human-scale distance, the displacement is essentially uniform. You cannot see a 30 cm rise when the entire city rises by 30 cm.

But a gravimeter can detect it. A gravimeter measures the gravitational acceleration at a fixed coordinate point (a specific latitude, longitude, and height above the ellipsoid). When the solid Earth deforms, two things change at that coordinate point, and both affect the reading.


## 7.2 How deformation changes what you measure

A gravimeter bolted to the crust rides with the deforming Earth. The deformation produces two competing effects on the measured gravity:

**Effect 1: The ground moves up (free-air effect).**
The tidal bulge lifts the gravimeter radially outward by some displacement $\Delta r$. Being farther from the center of the Earth, the gravimeter experiences a weaker gravitational pull. This is the same free-air effect that causes gravity to decrease as you climb a mountain: moving away from Earth's center reduces $g$.

**Effect 2: Mass redistributes (potential effect).**
When the Earth deforms, mass shifts around. The bulge that lifts you also means there is more mass beneath some places and less beneath others. This redistribution of mass changes the gravitational potential -- and hence the gravitational acceleration -- at your location. This is a separate effect from the displacement.

These two effects partially oppose each other: the free-air effect reduces your measured tidal gravity change, while the mass redistribution can either increase or decrease it depending on the details. To quantify them, we need the **Love numbers**.


## 7.3 Love numbers: parameterizing deformation

The Love numbers are dimensionless parameters that describe how much the Earth deforms in response to an applied tidal potential, compared to a reference case. They were introduced by Augustus Love in 1911.

### The tidal potential

The tidal potential from an external body (Moon or Sun) can be expanded in spherical harmonics. The dominant term is degree $l = 2$ (the quadrupolar term), which we call $V_2$. This is a function of position on the Earth's surface and of time (because the Moon and Sun move). For what follows, we only need to know:

- $V_2$ is a degree-2 spherical harmonic in the angular coordinates
- At the Earth's surface ($r = R$), $V_2$ depends on the zenith angle of the tide-raising body and the body's distance
- The tidal acceleration on a rigid Earth equals $-\partial V_2 / \partial r$ in the radial direction

We will not need the explicit form of $V_2$. Everything is expressed as ratios relative to $V_2$.

### The displacement Love number $h_2$

The Love number $h_2$ quantifies the radial displacement of the Earth's surface in response to $V_2$:

$$\Delta r = h_2 \, \frac{V_2}{g}$$

where $g$ is the surface gravity. This is the *definition* of $h_2$.

**Physical meaning:** On a perfectly fluid body with no rigidity (which deforms until it reaches hydrostatic equilibrium with the tidal potential), the radial displacement would be exactly $V_2 / g$. The Love number $h_2$ is the ratio of the actual displacement to this fluid limit:
- $h_2 = 1$: perfectly fluid body (deforms completely)
- $h_2 = 0$: perfectly rigid body (no deformation)
- For the real Earth: $h_2 \approx 0.61$, meaning the Earth deforms to about 61% of the fluid limit

### The potential Love number $k_2$

The deformation of the Earth redistributes mass, which produces an additional gravitational potential at the surface. The Love number $k_2$ quantifies this:

$$V_{\text{deformation}} = k_2 \, V_2$$

That is, the additional potential caused by the Earth's deformation is $k_2$ times the applied tidal potential. Again, this is a definition.

**Physical meaning:**
- $k_2 = 0$: rigid body (no deformation, no additional potential)
- For the real Earth: $k_2 \approx 0.30$, meaning the deformation-induced potential is about 30% of the applied tidal potential

### Why only degree 2?

The tidal potential is dominated by the degree-2 ($l=2$) term. Higher-degree terms ($l=3, 4, \ldots$) exist but are much smaller because the tidal potential involves the ratio $(R/d)^l$ where $d$ is the distance to the Moon or Sun. Since $R/d \sim 1/60$ for the Moon and $R/d \sim 4 \times 10^{-5}$ for the Sun, higher powers are strongly suppressed. The degree-3 contribution to lunar tidal gravity is roughly $(R/d) \sim 1.7\%$ of degree 2 -- small enough to ignore at our accuracy level.


## 7.4 Deriving the gravimetric factor $\delta$

We now derive the factor $\delta$ by which the elastic response modifies the measured tidal gravity change. We do this in three steps: (1) compute the rigid-Earth tidal gravity, (2) add the free-air correction from displacement, and (3) add the effect of mass redistribution.

### Step 1: Tidal gravity on a rigid Earth

On a rigid Earth, the only effect of the tidal potential $V_2$ is the direct gravitational perturbation. The change in radial gravity is the negative radial gradient of the potential:

$$\Delta g_{\text{rigid}} = -\frac{\partial V_2}{\partial r}$$

Since $V_2$ is a degree-2 spherical harmonic, it scales as $r^2$ inside the Earth (for the part generated by external sources). Therefore:

$$\frac{\partial V_2}{\partial r} = \frac{2 \, V_2}{r}$$

Evaluating at the surface $r = R$:

$$\boxed{\Delta g_{\text{rigid}} = -\frac{2 \, V_2}{R}}$$

This is our baseline -- the tidal gravity change if the Earth were perfectly rigid.

### Step 2: Gravity change from radial displacement (free-air effect)

On an elastic Earth, the gravimeter is displaced radially outward by $\Delta r = h_2 V_2 / g$. Moving upward in the Earth's gravity field decreases the measured gravity. The rate of decrease is the free-air gradient:

$$\frac{\partial g}{\partial r} \approx -\frac{2g}{R}$$

This is the familiar result from Newtonian gravity: for a spherical Earth, $g = GM/r^2$, so $\partial g / \partial r = -2GM/r^3 = -2g/r$, evaluated at $r = R$.

The gravity change due to the displacement is:

$$\Delta g_{\text{disp}} = -\frac{\partial g}{\partial r} \, \Delta r$$

Wait -- we need the sign carefully. If the gravimeter moves *up* (increasing $r$), gravity *decreases*. Since $\partial g/\partial r$ is negative (gravity weakens with altitude), and $\Delta r$ is positive (upward displacement), the change in measured gravity is:

$$\Delta g_{\text{disp}} = \frac{\partial g}{\partial r} \, \Delta r = -\frac{2g}{R} \cdot \frac{h_2 V_2}{g}$$

$$\boxed{\Delta g_{\text{disp}} = -\frac{2 \, h_2 \, V_2}{R}}$$

Note the sign: negative. The displacement *reduces* the measured gravity perturbation (the gravimeter moves away from the center, weakening the pull). But recall that $\Delta g_{\text{rigid}}$ is also negative (for a tide-raising potential with $V_2 > 0$ at the sub-lunar point, the tidal acceleration is directed outward, reducing the inward component of $g$). So the displacement effect has the *same sign* as the rigid-Earth tidal gravity -- it makes $|\Delta g|$ *larger*. The ground rising toward the Moon carries the gravimeter up, making the measured $g$ drop further.

### Step 3: Gravity change from mass redistribution (potential effect)

The deformation of the Earth produces an additional gravitational potential $k_2 V_2$ at the surface. The total potential perturbation at the surface of the deformed Earth is therefore:

$$V_{\text{total}} = V_2 + k_2 V_2 = (1 + k_2) \, V_2$$

But wait -- we already accounted for the $V_2$ part in Step 1 (the rigid-Earth tidal gravity). The *additional* gravity perturbation from the mass redistribution alone comes from the $k_2 V_2$ part.

The additional potential $k_2 V_2$ is generated by the deformed Earth itself (an internal source). For an internal degree-$l$ potential, the radial dependence at and above the surface goes as $r^{-(l+1)}$. So for $l = 2$:

$$\frac{\partial (k_2 V_2)}{\partial r} = -(l+1) \, \frac{k_2 V_2}{r} = -\frac{3 \, k_2 V_2}{R}$$

The gravity perturbation from mass redistribution is:

$$\Delta g_{\text{pot}} = -\frac{\partial (k_2 V_2)}{\partial r} = \frac{3 \, k_2 V_2}{R}$$

Wait -- let us be careful with signs and the direction of the gradient.

The radial component of gravitational acceleration is $g_r = -\partial V/\partial r$ (pointing inward, toward the center). A perturbation $\delta V = k_2 V_2$ that decays as $r^{-3}$ outside the surface has:

$$\frac{\partial (k_2 V_2)}{\partial r}\bigg|_{\text{ext}} = -3 \, \frac{k_2 V_2}{R}$$

Therefore the associated gravity perturbation (change in radial acceleration) is:

$$\Delta g_{\text{pot}} = -\frac{\partial (k_2 V_2)}{\partial r} = +\frac{3 \, k_2 V_2}{R}$$

$$\boxed{\Delta g_{\text{pot}} = +\frac{3 \, k_2 \, V_2}{R}}$$

The sign is *positive* -- opposite to $\Delta g_{\text{rigid}}$. The mass redistribution *opposes* the tidal gravity signal. Physically, the bulge of extra mass beneath you partially pulls you back down, counteracting the tidal stretching.

### Step 4: Combine all three contributions

The total gravity change on the elastic Earth is:

$$\Delta g_{\text{elastic}} = \Delta g_{\text{rigid}} + \Delta g_{\text{disp}} + \Delta g_{\text{pot}}$$

Substituting:

$$\Delta g_{\text{elastic}} = -\frac{2V_2}{R} - \frac{2 h_2 V_2}{R} + \frac{3 k_2 V_2}{R}$$

Factor out $-2V_2/R = \Delta g_{\text{rigid}}$:

$$\Delta g_{\text{elastic}} = -\frac{2V_2}{R}\left(1 + h_2 - \frac{3}{2} k_2\right)$$

$$\Delta g_{\text{elastic}} = \Delta g_{\text{rigid}} \cdot \underbrace{\left(1 + h_2 - \frac{3}{2}\, k_2\right)}_{\delta}$$

The gravimetric factor is:

$$\boxed{\delta = 1 + h_2 - \frac{3}{2}\, k_2}$$

This is the key result. It tells us that the measured tidal gravity change on the elastic Earth equals $\delta$ times the rigid-Earth value.

### Summary of the three terms

| Term | Origin | Expression | Effect on $|\Delta g|$ |
|------|--------|------------|----------------------|
| 1 | Rigid-Earth tidal gravity | $-2V_2/R$ | Baseline |
| $+h_2$ | Displacement (free-air) | $-2h_2 V_2/R$ | Amplifies |
| $-\frac{3}{2}k_2$ | Mass redistribution | $+3k_2 V_2/R$ | Reduces |

The displacement amplifies the signal because the gravimeter rides the tidal bulge upward, moving away from Earth's center and experiencing weaker gravity. The mass redistribution partially cancels this because the deformed Earth piles up mass beneath the bulge, pulling the gravimeter back.


## 7.5 Numerical values

### IERS 2010 Love numbers

The International Earth Rotation and Reference Systems Service (IERS) 2010 Conventions give the nominal degree-2 Love numbers for a non-rotating, elastic, isotropic (SNREI) Earth model:

$$h_2 = 0.6078, \qquad k_2 = 0.2980$$

### The gravimetric factor

Substituting into the formula:

$$\delta = 1 + 0.6078 - 1.5 \times 0.2980$$

$$\delta = 1 + 0.6078 - 0.4470$$

$$\boxed{\delta = 1.1608}$$

### Physical interpretation

The tidal signal is **16% larger** than it would be on a rigid Earth. Breaking down the contributions:

- The rigid-Earth baseline: **1.0000** (100%)
- The displacement (free-air) amplification: **+0.6078** (+60.8%)
- The mass redistribution reduction: **-0.4470** (-44.7%)
- **Net amplification: +0.1608 (+16.1%)**

The displacement effect dominates: the ground rises so much that the free-air effect (losing gravity by moving up) is large. But the mass redistribution -- the tidal bulge concentrating mass beneath you -- partially counteracts this. The net result is that the elastic Earth amplifies the tidal signal, but by much less than $h_2$ alone would suggest. Figure 11 overlays the tidal signal for a rigid Earth against the elastic prediction, with the shaded region between the curves making the 16% amplification visually striking.

### How Pytheas applies this

In Pytheas, the gravimetric factor is applied as a simple multiplicative correction to the tidal acceleration:

```python
H2 = 0.6078
K2 = 0.2980
DELTA_GRAV = 1.0 + H2 - 1.5 * K2   # = 1.1608

a_moon = DELTA_GRAV * tidal_acceleration(r, moon_position_ecef(dt), GM_MOON)
a_sun  = DELTA_GRAV * tidal_acceleration(r, sun_position_ecef(dt),  GM_SUN)
```

The `tidal_acceleration` function computes the rigid-Earth value (the exact Newtonian tidal acceleration). Multiplying by $\delta$ gives the measured tidal gravity change on the elastic Earth.


## 7.6 What we are ignoring

The derivation above assumes a single, frequency-independent gravimetric factor $\delta$. The real Earth is more complicated.

### Frequency dependence of Love numbers

The Love numbers $h_2$ and $k_2$ are not constants -- they depend on the frequency of the tidal forcing. The values 0.6078 and 0.2980 are appropriate for the semi-diurnal tides (period ~12 hours). At different frequencies, the Earth's response is slightly different because:

- At very long periods (months to years), the mantle can flow viscously, increasing the effective $h_2$
- At diurnal periods (~24 hours), there is a resonance that significantly alters the Love numbers

### The Free Core Nutation resonance

The most important frequency dependence occurs near the K1 tidal frequency (period $\approx$ 23.93 hours). At this frequency, the tidal forcing resonates with the **Free Core Nutation (FCN)** -- a normal mode of the Earth in which the fluid outer core and the solid mantle precess around slightly different axes.

Near the FCN resonance frequency, the Love numbers change rapidly. The gravimetric factor $\delta$ can deviate from 1.1608 by up to roughly 1%, or about 10 nGal in the gravity signal. This is the **single largest error source** in Pytheas. A proper correction would require using frequency-dependent Love numbers for each tidal constituent -- a table lookup rather than a single number.

### Higher-degree Love numbers

Our derivation used only the degree-2 ($l = 2$) Love numbers, because the tidal potential is dominated by $l = 2$. There exist degree-3 Love numbers ($h_3$, $k_3$), and the lunar tide does have a small degree-3 component ($\sim 1.7\%$ of degree 2). The effect on gravity is much smaller still and is negligible at the 10 nGal level.

### Latitude dependence

The Love numbers we used are for a spherically symmetric Earth model. The real Earth has lateral variations in elastic properties (thinner crust under oceans, thicker under continents, etc.). These cause the effective Love numbers to vary slightly with location, but the effect is small -- typically below 1% -- and is subsumed into ocean loading corrections for coastal sites.

### What these omissions cost

| Omission | Typical error | Notes |
|----------|--------------|-------|
| Frequency-independent $\delta$ (no FCN) | ~10 nGal | Largest near K1 frequency |
| No degree-3 Love numbers | <1 nGal | Suppressed by $(R/d)$ |
| No lateral variation | <1 nGal | Subsumed into ocean loading |

For Pytheas's target accuracy of ~10 nGal at inland sites, using a single gravimetric factor $\delta = 1.1608$ is adequate. The FCN resonance is the dominant limitation of this approach.

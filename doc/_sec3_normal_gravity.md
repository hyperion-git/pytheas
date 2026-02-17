# Section 3: Gravity on a Quiet Earth (Normal Gravity)

Before we bring in the Moon and Sun, we need to know what gravity *would* be at our location if tides did not exist — the "quiet Earth" baseline. This is called **normal gravity**, and computing it properly requires us to account for three things: the mass of the Earth, its rotation, and its shape.

This section builds that baseline from scratch, starting with the simplest possible model (a non-rotating sphere) and adding one physical effect at a time until we arrive at the formula Pytheas actually uses.

---

## 3.1 The non-rotating spherical case

The simplest model of Earth's gravity starts with Newton's law of universal gravitation. Two masses $M$ and $m$ separated by a distance $r$ attract each other with a force

$$F = \frac{G M m}{r^2}$$

where $G = 6.674 \times 10^{-11}$ m$^3$ kg$^{-1}$ s$^{-2}$ is the gravitational constant.

Newton proved a remarkable result known as the **shell theorem**: a uniform spherical shell of mass attracts an external particle exactly as if all the shell's mass were concentrated at its center. Because a solid sphere can be decomposed into nested shells, a uniform (or radially symmetric) sphere of total mass $M_E$ attracts a surface object of mass $m$ with

$$F = \frac{G M_E m}{R^2}$$

where $R$ is the sphere's radius. The gravitational acceleration — force per unit mass — is therefore

$$g = \frac{G M_E}{R^2}$$

This is purely a property of the Earth; the test mass $m$ cancels out.

**Numerical check.** Rather than using $G$ and $M_E$ separately (each of which is known to only about 4 significant figures), geodesists use their product $GM_E = 3.986\,004\,418 \times 10^{14}$ m$^3$ s$^{-2}$, which is known to 10 significant figures from satellite tracking. With a mean Earth radius $R = 6{,}371$ km $= 6.371 \times 10^6$ m:

$$g = \frac{3.986 \times 10^{14}}{(6.371 \times 10^6)^2} = \frac{3.986 \times 10^{14}}{4.059 \times 10^{13}} \approx 9.82 \text{ m/s}^2$$

This is in the right ballpark — about 0.2% from the standard textbook value of $9.81$ m/s$^2$. But the real Earth is neither perfectly spherical nor non-rotating, so the actual value depends on where you stand. To understand why, we need to account for rotation and shape.

---

## 3.2 Why rotation matters

We live on the surface of a spinning planet. An observer standing on the ground is not in an inertial frame — they are traveling in a circle around Earth's rotation axis, completing one revolution per sidereal day. In this rotating frame, a fictitious outward force appears: the **centrifugal force**.

### The centrifugal acceleration

Consider an observer at geodetic latitude $\varphi$ on a spherical Earth of radius $R$. Their distance from the rotation axis is

$$R_\perp = R \cos\varphi$$

In the rotating frame, the outward centrifugal acceleration has magnitude

$$a_\text{cen} = \Omega^2 R_\perp = \Omega^2 R \cos\varphi$$

where $\Omega = 7.292\,115 \times 10^{-5}$ rad/s is Earth's angular velocity.

**How big is this?** At the equator ($\varphi = 0$), where the effect is largest, we use the equatorial radius $R_\text{eq} = 6{,}378{,}137$ m:

$$a_\text{cen}^{\max} = \Omega^2 R_\text{eq} = (7.292 \times 10^{-5})^2 \times 6.378 \times 10^6$$

$$= 5.317 \times 10^{-9} \times 6.378 \times 10^6 \approx 0.0339 \text{ m/s}^2$$

That is about $0.34\%$ of $g \approx 9.82$ m/s$^2$. Small, but easily measurable — it would show up in the fourth significant figure of any gravimeter reading. At the poles ($\varphi = \pm 90°$), $\cos\varphi = 0$ and the centrifugal acceleration vanishes: a person standing on the pole is on the rotation axis and is not traveling in a circle.

### Effective gravity

The gravitational acceleration $\mathbf{g}_\text{grav}$ points toward the center of the Earth (radially inward), while the centrifugal acceleration $\mathbf{a}_\text{cen}$ points outward from the rotation axis (perpendicular to it). These two directions are not exactly opposite except at the equator and poles. What an observer measures as "gravity" — the acceleration of a dropped object — is the vector sum:

$$\mathbf{g}_\text{eff} = \mathbf{g}_\text{grav} + \mathbf{a}_\text{cen}$$

To find the component that reduces the magnitude of gravity (the part that acts "against" gravity), we project the centrifugal acceleration onto the radial direction. At latitude $\varphi$, the centrifugal vector $\Omega^2 R \cos\varphi$ (directed away from the axis) has a radial outward component of $\Omega^2 R \cos^2\!\varphi$. So the effective radial gravity is approximately

$$g_\text{eff}(\varphi) \approx \frac{GM_E}{R^2} - \Omega^2 R \cos^2\!\varphi$$

There is also a small tangential (equatorward) component $\frac{1}{2}\Omega^2 R \sin 2\varphi$ that deflects the local vertical away from the geometric radial direction, but this is a second-order geometric effect that will be absorbed when we move to the correct ellipsoidal shape.

**Summary of rotation effects:**
- Equator: gravity is reduced by $\sim 0.034$ m/s$^2$ ($\sim 0.3\%$)
- Poles: no reduction
- Mid-latitudes: intermediate values, scaling as $\cos^2\!\varphi$

---

## 3.3 Why shape matters

If Earth were non-rotating, its own self-gravity would pull it into a perfect sphere. But it rotates, and the centrifugal effect is strongest at the equator. Over geological time, the Earth has deformed into an **oblate ellipsoid** — it bulges at the equator and is slightly flattened at the poles. This shape is quantified by the **flattening**:

$$f = \frac{a - b}{a}$$

where $a$ is the equatorial (semi-major) radius and $b$ is the polar (semi-minor) radius. For the standard GRS80 reference ellipsoid:

$$f = \frac{1}{298.257} \qquad \Longrightarrow \qquad a - b \approx 21.4 \text{ km}$$

The equatorial radius ($a = 6{,}378{,}137$ m) is about 21 km larger than the polar radius ($b \approx 6{,}356{,}752$ m). That is only a 0.3% difference in radius, but it has a measurable effect on gravity through two competing mechanisms:

1. **Distance effect.** A person at the equator is $\sim$21 km farther from Earth's center than a person at the pole. Since gravity falls off as $1/r^2$, this alone would make equatorial gravity weaker by about $2 \times 21/6371 \approx 0.7\%$.

2. **Mass distribution effect.** The equatorial bulge redistributes mass, placing more of it near the equatorial plane. This partially compensates the distance effect for equatorial observers (more nearby mass pulling sideways and downward) and slightly reduces gravity at the poles (some mass has been moved away from the axis).

The net result of both shape effects combined with the centrifugal reduction is that gravity varies from about $9.78$ m/s$^2$ at the equator to about $9.83$ m/s$^2$ at the poles — a total variation of roughly $0.5\%$, or about $0.05$ m/s$^2$. This is a thousand times larger than the tidal signals we want to model in later sections, so getting the baseline right is essential.

---

## 3.4 Somigliana's formula

Computing gravity on a rotating, oblate Earth from first principles requires solving for the gravitational potential of an ellipsoid — a problem first addressed by Clairaut in the 18th century and refined over the next two centuries. The full derivation involves expanding the gravitational potential in spherical harmonics on an equipotential surface and is well beyond our scope. Instead, we will state the closed-form result and carefully explain what each piece means.

### The level ellipsoid

Geodesists define a **reference ellipsoid** — a mathematically perfect ellipsoidal surface that approximates mean sea level. This surface is chosen to be an **equipotential surface**: the total potential (gravitational + centrifugal) is the same everywhere on it. A ball placed on this surface would not roll anywhere, because there is no potential gradient along the surface.

The standard reference is the **GRS80 ellipsoid** (Geodetic Reference System 1980), defined by four parameters:
- $a = 6{,}378{,}137$ m (equatorial semi-major axis)
- $GM_E = 3.986\,004\,418 \times 10^{14}$ m$^3$ s$^{-2}$ (gravitational parameter)
- $J_2 = 1.082\,63 \times 10^{-3}$ (dynamic form factor — a measure of oblateness)
- $\Omega = 7.292\,115 \times 10^{-5}$ rad/s (angular velocity)

From these, all other quantities (including $b$, $f$, and the gravity values at the equator and poles) can be derived.

### The formula

**Somigliana's formula** (1929) gives the exact gravity on the surface of the level ellipsoid as a function of geodetic latitude $\varphi$:

$$\gamma_0(\varphi) = \frac{\gamma_e \left(1 + k \sin^2\!\varphi\right)}{\sqrt{1 - e^2 \sin^2\!\varphi}}$$

where:

| Symbol | Meaning | GRS80 Value |
|--------|---------|-------------|
| $\gamma_e$ | Gravity at the equator | $9.780\,325\,3141\;\mathrm{m/s^2}$ |
| $\gamma_p$ | Gravity at the pole | $9.832\,184\,9378\;\mathrm{m/s^2}$ |
| $k$ | Somigliana constant | $(b\,\gamma_p)/(a\,\gamma_e) - 1$ |
| $e^2$ | First eccentricity squared | $2f - f^2$ |

Let us unpack each piece to build physical intuition.

**The denominator** $\sqrt{1 - e^2 \sin^2\!\varphi}$: The eccentricity $e$ measures how "squashed" the ellipsoid is. For GRS80, $e^2 = 2f - f^2 \approx 0.00669$. At the equator ($\sin^2\!\varphi = 0$), the denominator is 1. At the pole ($\sin^2\!\varphi = 1$), it becomes $\sqrt{1 - e^2} = b/a \approx 0.99665$. This geometric factor accounts for the fact that the ellipsoidal surface curves differently at different latitudes — the surface normal direction changes, and with it the component of gravity perpendicular to the surface. Dividing by this factor slightly increases gravity at the poles relative to the equator.

**The numerator factor** $1 + k \sin^2\!\varphi$: The constant $k$ encodes the fractional difference between polar and equatorial gravity, corrected by the geometric factor. Numerically:

$$k = \frac{b \cdot \gamma_p}{a \cdot \gamma_e} - 1$$

We can compute this from the GRS80 values. With $b = a(1-f) = 6{,}356{,}752.314$ m:

$$k = \frac{6{,}356{,}752.314 \times 9.832\,184\,9378}{6{,}378{,}137.0 \times 9.780\,325\,3141} - 1 \approx 0.001\,932$$

At the equator, $\sin^2\!\varphi = 0$ and the numerator factor is simply 1, so $\gamma_0(0°) = \gamma_e$. At the pole, $\sin^2\!\varphi = 1$ and the numerator becomes $1 + k$, while the denominator becomes $b/a$, so

$$\gamma_0(90°) = \frac{\gamma_e (1 + k)}{b/a} = \frac{\gamma_e \cdot a (1 + k)}{b} = \frac{a \gamma_e + a\gamma_e k}{b}$$

By the definition of $k$, $a \gamma_e k = b\gamma_p - a\gamma_e$, so the numerator becomes $b\gamma_p$ and we get $\gamma_0(90°) = \gamma_p$. The formula correctly reproduces the known values at both the equator and pole, and smoothly interpolates between them.

**Quick check at 45° latitude:**

$$\sin^2(45°) = 0.5$$

$$\gamma_0(45°) = \frac{9.780\,325 \times (1 + 0.001\,932 \times 0.5)}{\sqrt{1 - 0.006\,694 \times 0.5}} = \frac{9.780\,325 \times 1.000\,966}{\sqrt{0.996\,653}}$$

$$= \frac{9.789\,776}{0.998\,325} \approx 9.806\,2 \text{ m/s}^2$$

This lies between the equatorial and polar values, as expected, and matches the commonly quoted value for mid-latitudes.

**Figure 4** shows how normal gravity varies from equator to pole, comparing three models: the simple spherical approximation (constant at all latitudes), the sphere-plus-rotation model (which accounts for centrifugal acceleration), and the full GRS80 Somigliana formula (which also incorporates the oblate shape of the Earth). The figure demonstrates that rotation accounts for most of the latitude variation, while the ellipsoidal shape provides an additional refinement.

![Figure 4: Normal gravity vs. geodetic latitude comparing three models — spherical Earth (constant), rotation-corrected spherical model, and full GRS80 Somigliana formula on the oblate ellipsoid.](figures/fig04_gravity_vs_latitude.png)

---

## 3.5 Free-air correction (going up)

Somigliana's formula gives gravity *on* the ellipsoid surface — at zero altitude. But real gravimeters, and real atom interferometers, operate at some height $h$ above the ellipsoid. How does gravity change with altitude?

The intuition is straightforward: moving away from the Earth weakens gravity. But the exact rate of decrease depends on latitude (because the Earth is not a sphere), and at heights of a few hundred meters or more, the nonlinear $1/r^2$ falloff means we need more than just a first-order correction.

### Setting up the Taylor expansion

We expand gravity at height $h$ as a Taylor series in $h/a$ (height divided by Earth's semi-major axis — a small dimensionless number for any terrestrial measurement):

$$\gamma(\varphi, h) = \gamma_0(\varphi) \left[1 + c_1 \frac{h}{a} + c_2 \left(\frac{h}{a}\right)^2 + \cdots\right]$$

We need to determine the coefficients $c_1$ and $c_2$.

### First-order coefficient: the gradient of gravity

**Start with the sphere.** For a non-rotating spherical Earth, $g(r) = GM_E/r^2$. The radial gradient is

$$\frac{\partial g}{\partial r} = -\frac{2\,GM_E}{r^3} = -\frac{2g}{r}$$

So gravity decreases at a rate of $2g/r$ per meter of altitude. At the surface, $r = a$, and the first-order correction would be

$$g(h) \approx g_0 \left(1 - \frac{2h}{a}\right)$$

This is the classic "free-air" correction, and for a sphere the first-order coefficient is simply $c_1 = -2$.

**But Earth is not a sphere.** On a rotating oblate ellipsoid, the gradient of gravity involves additional terms from the flattening $f$ and from rotation. The exact first-order coefficient, derived from the theory of the level ellipsoid, is

$$c_1 = -2\left(1 + f + m - 2f\sin^2\!\varphi\right)$$

where $m$ is the **geodetic ratio** of centrifugal to gravitational acceleration at the equator:

$$m = \frac{\Omega^2 a^2 b}{GM_E}$$

Let us compute $m$:

$$m = \frac{(7.292\,115 \times 10^{-5})^2 \times (6{,}378{,}137)^2 \times 6{,}356{,}752.314}{3.986\,004\,418 \times 10^{14}}$$

$$= \frac{5.317 \times 10^{-9} \times 4.068 \times 10^{13} \times 6.357 \times 10^6}{3.986 \times 10^{14}}$$

$$\approx \frac{1.375 \times 10^{12}}{3.986 \times 10^{14}} \approx 0.003\,449$$

The physical meaning of each correction in $c_1$:

- The leading $-2$ comes from the $1/r^2$ falloff, as for a sphere.
- The $f$ term (flattening, $\approx 0.003\,353$) corrects for the ellipsoidal geometry — the surface curves away from a sphere, so the effective "radius" that enters the gradient is modified.
- The $m$ term ($\approx 0.003\,449$) accounts for the decrease of centrifugal acceleration with height. As you go up, you are farther from the rotation axis relative to the Earth's surface, and the centrifugal contribution to effective gravity changes.
- The $-2f\sin^2\!\varphi$ term captures the latitude dependence: at the poles, flattening effects are different than at the equator because the local radius of curvature differs.

**Numerical magnitude of the first-order correction.** At the equator ($\sin^2\!\varphi = 0$):

$$c_1 = -2(1 + 0.003\,353 + 0.003\,449) = -2 \times 1.006\,802 = -2.013\,6$$

Compared to the spherical value of $-2$, the ellipsoidal correction changes the first-order coefficient by about $0.7\%$. The gravity gradient per meter of altitude is

$$\frac{\partial \gamma}{\partial h}\bigg|_{h=0} = -\frac{2\gamma_0 \times \text{fac}}{a} \approx -\frac{2 \times 9.780 \times 1.007}{6{,}378{,}137} \approx -3.086 \times 10^{-6} \text{ m/s}^2\text{/m}$$

In the traditional unit, that is $-0.3086$ mGal/m (where 1 milliGal $= 10^{-5}$ m/s$^2$) — the well-known **free-air gradient**. Over 100 m of altitude, this gives a reduction of about $30.9$ mGal, or $3.09 \times 10^{-4}$ m/s$^2$.

### Second-order coefficient

The second-order coefficient comes from continuing the Taylor expansion of $GM_E/r^2$ to the next term. For the spherical case:

$$g(r) = \frac{GM_E}{r^2} = \frac{GM_E}{(a+h)^2} = \frac{GM_E}{a^2}\frac{1}{(1+h/a)^2}$$

Using the binomial expansion $(1+x)^{-2} = 1 - 2x + 3x^2 - \cdots$:

$$g(h) = g_0\left(1 - 2\frac{h}{a} + 3\frac{h^2}{a^2} - \cdots\right)$$

The second-order coefficient is $c_2 = +3$. On the ellipsoid, the full theory introduces further small corrections to this coefficient as well, but the dominant term remains $+3h^2/a^2$, and the corrections are at the sub-percent level of an already small term. Pytheas uses $c_2 = 3$ without additional corrections, which is adequate for heights up to several kilometers.

### The complete formula

Combining everything, Pytheas computes gravity at height $h$ above the ellipsoid as:

$$\boxed{\gamma(\varphi, h) = \gamma_0(\varphi) \left[1 - \frac{2\left(1 + f + m - 2f\sin^2\!\varphi\right)}{a}\,h + \frac{3\,h^2}{a^2}\right]}$$

where $\gamma_0(\varphi)$ is Somigliana's formula from Section 3.4.

This is a standard result in physical geodesy, found (in various equivalent forms) in textbooks such as Torge & Muller (*Geodesy*, 4th ed.) and Hofmann-Wellenhof & Moritz (*Physical Geodesy*). The formula is exact to second order in $h/a$ on the level ellipsoid.

**Figure 5** illustrates the free-air correction as a function of altitude, showing both the magnitude of the gravity reduction and the difference between the first-order and second-order approximations. Panel (b) highlights how the second-order term, though small at low altitudes, becomes significant above a few hundred meters.

![Figure 5: Free-air correction to gravity as a function of altitude. Panel (a) shows the total gravity reduction; panel (b) shows the difference between first-order and second-order approximations, demonstrating when the quadratic term becomes significant.](figures/fig05_free_air.png)

### Why the second-order term matters

At first glance, $h^2/a^2$ seems negligibly small. Let us check.

**At $h = 100$ m:**

$$\frac{h^2}{a^2} = \frac{10^4}{4.068 \times 10^{13}} \approx 2.5 \times 10^{-10}$$

$$3\gamma_0 \frac{h^2}{a^2} \approx 3 \times 9.8 \times 2.5 \times 10^{-10} \approx 7.3 \times 10^{-9} \text{ m/s}^2 = 0.73 \text{ }\mu\text{Gal}$$

This is below the $\sim 1\,\mu$Gal precision target, so at 100 m the second-order term is barely relevant.

**At $h = 1{,}000$ m:**

$$\frac{h^2}{a^2} = \frac{10^6}{4.068 \times 10^{13}} \approx 2.5 \times 10^{-8}$$

$$3\gamma_0 \frac{h^2}{a^2} \approx 3 \times 9.8 \times 2.5 \times 10^{-8} \approx 7.3 \times 10^{-7} \text{ m/s}^2 = 73 \text{ }\mu\text{Gal}$$

At a kilometer of altitude, omitting the second-order term would introduce an error of $\sim 73\,\mu$Gal — far too large for precision gravimetry. Even at 500 m, the error would be $\sim 18\,\mu$Gal. The second-order term is therefore essential for any site significantly above sea level.

### Numerical example: gravity in Munich

Let us compute normal gravity for a site in Munich, Germany, at latitude $48.1°$ N and altitude 520 m above the ellipsoid.

**Step 1: Somigliana on the ellipsoid.**

$$\sin^2(48.1°) \approx 0.5540$$

$$\gamma_0 = \frac{9.780\,325 \times (1 + 0.001\,932 \times 0.5540)}{\sqrt{1 - 0.006\,694 \times 0.5540}}$$

$$= \frac{9.780\,325 \times 1.001\,070}{\sqrt{0.996\,292}} = \frac{9.790\,792}{0.998\,145} \approx 9.809\,0 \text{ m/s}^2$$

**Step 2: Free-air correction.**

$$\text{fac} = 1 + 0.003\,353 + 0.003\,449 - 2 \times 0.003\,353 \times 0.5540 = 1.003\,087$$

First-order term:
$$-\frac{2 \times 1.003\,087 \times 520}{6{,}378{,}137} = -1.633 \times 10^{-4}$$

Second-order term:
$$\frac{3 \times 520^2}{(6{,}378{,}137)^2} = \frac{811{,}200}{4.068 \times 10^{13}} = 1.994 \times 10^{-8}$$

$$\gamma(48.1°, 520) \approx 9.809\,0 \times (1 - 1.633 \times 10^{-4} + 1.994 \times 10^{-8})$$

$$\approx 9.809\,0 \times 0.999\,837 \approx 9.807\,4 \text{ m/s}^2$$

The free-air correction reduced gravity by about $1.60 \times 10^{-3}$ m/s$^2$ (160 mGal), and the second-order term contributed about $0.2\,\mu$Gal — small at this altitude, but it would grow to tens of $\mu$Gal at higher elevations.

---

### Summary

| Model | Formula | Accuracy |
|-------|---------|----------|
| Spherical, non-rotating | $g = GM_E / R^2$ | ${\sim}\,0.5\%$ |
| + rotation | $g_\text{eff} = GM_E/R^2 - \Omega^2 R \cos^2\!\varphi$ | ${\sim}\,0.2\%$ |
| + shape (Somigliana) | $\gamma_0(\varphi)$ on ellipsoid | exact on surface |
| + altitude (free-air) | $\gamma(\varphi, h)$ with $h/a$ and $h^2/a^2$ terms | ${<}\,1\;\mu\mathrm{Gal}$ to ${\sim}\,1\;\mathrm{km}$ |

Starting from a single number ($9.82$ m/s$^2$), we now have a formula that gives gravity to sub-microGal accuracy anywhere on or near the Earth's surface — a precision of better than one part in $10^{10}$. This is the static baseline onto which tidal perturbations (Section 4) will be superimposed.

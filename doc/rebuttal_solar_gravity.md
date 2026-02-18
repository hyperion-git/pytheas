# Why a Tilted Sensor Cannot Measure Direct Gravitational Pull

## The Claim

> A sensor tilted relative to the local vertical can detect the direct
> gravitational pull of the Sun (~6 × 10⁻³ m/s²) or the Moon
> (~3.3 × 10⁻⁵ m/s²) as a time-varying signal, because the projection
> of the gravitational acceleration onto the sensor axis changes as
> these bodies move across the sky.

This is incorrect. The derivation below shows, in a point-mass model with
no approximations beyond Newtonian mechanics, that the direct pull of
**every** external body cancels exactly against the corresponding
acceleration of the laboratory, leaving only the tidal residual —
regardless of sensor orientation.


## Setup

Four point masses in an inertial frame:

| Body | Mass | Position |
|------|------|----------|
| Sun | $M_S$ | $\mathbf{R}_S(t)$ |
| Moon | $M_M$ | $\mathbf{R}_M(t)$ |
| Earth | $M_E$ | $\mathbf{R}_E(t)$ |
| Test mass | $m$ | $\mathbf{x}(t)$ |

The test mass sits on Earth's surface and is connected to the sensor by a
spring (or equivalent restoring mechanism). The sensor frame is rigidly
attached to the Earth.

The sensor reads the **spring force** $\mathbf{F}_{\text{spring}}$ — the
non-gravitational force required to keep the test mass co-moving with the
lab. This is the only quantity accessible to experiment; gravity itself is
not directly measurable (equivalence principle).


## Step 1: Equation of Motion in the Inertial Frame

Newton's second law for the test mass:

$$
m\,\ddot{\mathbf{x}}
  = \underbrace{-\frac{G M_E\, m}{|\mathbf{x} - \mathbf{R}_E|^3}
    (\mathbf{x} - \mathbf{R}_E)}_{\text{Earth's gravity}}
  \;\underbrace{-\frac{G M_S\, m}{|\mathbf{x} - \mathbf{R}_S|^3}
    (\mathbf{x} - \mathbf{R}_S)}_{\text{Sun's gravity}}
  \;\underbrace{-\frac{G M_M\, m}{|\mathbf{x} - \mathbf{R}_M|^3}
    (\mathbf{x} - \mathbf{R}_M)}_{\text{Moon's gravity}}
  \;+\; \mathbf{F}_{\text{spring}}
$$


## Step 2: Equation of Motion for Earth's Center

Earth's center of mass is accelerated by both external bodies:

$$
M_E\,\ddot{\mathbf{R}}_E
  = -\frac{G M_S\, M_E}{|\mathbf{R}_E - \mathbf{R}_S|^3}
    (\mathbf{R}_E - \mathbf{R}_S)
  \;-\; \frac{G M_M\, M_E}{|\mathbf{R}_E - \mathbf{R}_M|^3}
    (\mathbf{R}_E - \mathbf{R}_M)
$$

Dividing by $M_E$:

$$
\ddot{\mathbf{R}}_E
  = -\frac{G M_S}{|\mathbf{R}_E - \mathbf{R}_S|^3}
    (\mathbf{R}_E - \mathbf{R}_S)
  \;-\; \frac{G M_M}{|\mathbf{R}_E - \mathbf{R}_M|^3}
    (\mathbf{R}_E - \mathbf{R}_M)
$$

This is Earth's acceleration toward the Sun and Moon combined. Every object
attached to the Earth — including the sensor housing, the mounting bracket,
and the reference frame of the measurement — shares this acceleration.


## Step 3: Transform to the Earth-Centered Frame

Define the position of the test mass relative to Earth's center:

$$\mathbf{r} \equiv \mathbf{x} - \mathbf{R}_E$$

and the geocentric positions of the external bodies:

$$\mathbf{R} \equiv \mathbf{R}_S - \mathbf{R}_E, \qquad
  \mathbf{D} \equiv \mathbf{R}_M - \mathbf{R}_E$$

The relative acceleration is $\ddot{\mathbf{r}} = \ddot{\mathbf{x}} - \ddot{\mathbf{R}}_E$.

Substituting from Steps 1 and 2, each external body produces a term of the
form

$$
-\frac{G M_B\,m}{|\mathbf{x} - \mathbf{R}_B|^3}(\mathbf{x} - \mathbf{R}_B)
\;+\; \frac{G M_B\,m}{|\mathbf{R}_E - \mathbf{R}_B|^3}(\mathbf{R}_E - \mathbf{R}_B)
$$

where $B \in \{S, M\}$. Rewriting in geocentric variables
($\mathbf{x} - \mathbf{R}_B = \mathbf{r} - \mathbf{R}_B'$ with
$\mathbf{R}_B' = \mathbf{R}_B - \mathbf{R}_E$):

$$
\boxed{
m\,\ddot{\mathbf{r}}
  = -\frac{G M_E\, m}{|\mathbf{r}|^3}\,\mathbf{r}
  \;+\; G M_S\, m \left[
    \frac{\mathbf{R} - \mathbf{r}}{|\mathbf{R} - \mathbf{r}|^3}
    - \frac{\mathbf{R}}{|\mathbf{R}|^3}
  \right]
  \;+\; G M_M\, m \left[
    \frac{\mathbf{D} - \mathbf{r}}{|\mathbf{D} - \mathbf{r}|^3}
    - \frac{\mathbf{D}}{|\mathbf{D}|^3}
  \right]
  \;+\; \mathbf{F}_{\text{spring}}
}
$$


### Exact cancellation of the direct pull

Consider the Sun's contribution to the boxed equation. It arose from
subtracting Earth's acceleration (Step 2) from the test mass acceleration
(Step 1). Tracing the Sun's terms through the subtraction:

**From the test mass** (Step 1), the Sun contributes an acceleration:

$$-\frac{G M_S}{|\mathbf{x} - \mathbf{R}_S|^3}(\mathbf{x} - \mathbf{R}_S)$$

**From Earth's center** (Step 2), subtracted via
$\ddot{\mathbf{r}} = \ddot{\mathbf{x}} - \ddot{\mathbf{R}}_E$:

$$+\frac{G M_S}{|\mathbf{R}_E - \mathbf{R}_S|^3}(\mathbf{R}_E - \mathbf{R}_S)$$

Substituting $\mathbf{x} - \mathbf{R}_S = \mathbf{r} - \mathbf{R}$ and
$\mathbf{R}_E - \mathbf{R}_S = -\mathbf{R}$, their sum becomes:

$$G M_S \underbrace{\left[\frac{\mathbf{R} - \mathbf{r}}{|\mathbf{R} - \mathbf{r}|^3} - \frac{\mathbf{R}}{R^3}\right]}_{\equiv\;\mathbf{f}(\mathbf{r})}$$

This expression $\mathbf{f}(\mathbf{r})$ is the difference of the Sun's
gravitational field evaluated at two points: the test mass location
($\mathbf{r}$ from Earth's center) and Earth's center itself ($\mathbf{r}=0$).
By construction:

$$\mathbf{f}(\mathbf{0}) = \frac{\mathbf{R}}{R^3} - \frac{\mathbf{R}}{R^3} = \mathbf{0} \qquad\text{(exactly)}$$

This is the exact cancellation. The "direct pull" —
$G M_S\,\mathbf{R}/R^3 = (G M_S / R^2)\,\mathbf{e}_{\mathbf{R}}$ — is the
value of the Sun's field at Earth's center. It appears with opposite signs
in the test mass and Earth equations and cancels identically, without
approximation. What remains,
$\mathbf{f}(\mathbf{r}) = \mathbf{f}(\mathbf{r}) - \mathbf{f}(\mathbf{0})$,
is purely the *variation* of the field across the distance $\mathbf{r}$ —
the tidal acceleration.

The identical cancellation occurs for the Moon, with $\mathbf{D}$ replacing
$\mathbf{R}$. The result generalises to any number of external bodies: for
each, the direct pull is common to test mass and Earth, and subtracts out
exactly.


### Taylor expansion: identifying the surviving term

To determine the magnitude and structure of what survives, expand
$\mathbf{f}(\mathbf{r})$ for a generic body at geocentric position
$\mathbf{R}$ (with $R = |\mathbf{R}|$ and
$\mathbf{e}_{\mathbf{R}} = \mathbf{R}/R$) in powers of the small
parameter $|\mathbf{r}|/R$.

Start with the inverse cube:

$$|\mathbf{R} - \mathbf{r}|^2 = R^2\left(1 - 2\,\frac{\mathbf{e}_{\mathbf{R}} \cdot \mathbf{r}}{R} + \frac{r^2}{R^2}\right)$$

so

$$\frac{1}{|\mathbf{R} - \mathbf{r}|^3} = \frac{1}{R^3}\left(1 + \frac{3\,\mathbf{e}_{\mathbf{R}}\cdot\mathbf{r}}{R} + O\!\left(\frac{r^2}{R^2}\right)\right)$$

Multiply by $(\mathbf{R} - \mathbf{r})$ and keep first-order terms:

$$\frac{\mathbf{R} - \mathbf{r}}{|\mathbf{R} - \mathbf{r}|^3} = \frac{1}{R^3}\Big[\mathbf{R} + 3(\mathbf{e}_{\mathbf{R}}\cdot\mathbf{r})\,\mathbf{e}_{\mathbf{R}} - \mathbf{r}\Big] + O\!\left(\frac{r^2}{R^4}\right)$$

Now subtract $\mathbf{R}/R^3$:

$$\frac{\mathbf{R} - \mathbf{r}}{|\mathbf{R} - \mathbf{r}|^3} - \frac{\mathbf{R}}{R^3} = \frac{1}{R^3}\Big[\cancel{\mathbf{R}} + 3(\mathbf{e}_{\mathbf{R}}\cdot\mathbf{r})\,\mathbf{e}_{\mathbf{R}} - \mathbf{r} - \cancel{\mathbf{R}}\Big] + O\!\left(\frac{r^2}{R^4}\right)$$

The $\mathbf{R}/R^3$ terms — which represent the **direct gravitational
pull** $GM/R^2$ in the direction $\mathbf{e}_{\mathbf{R}}$ — cancel identically.
What survives is the quadrupole (tidal) term:

$$\boxed{\mathbf{a}_{\text{tidal}} = \frac{GM}{R^3}\Big[3\,(\mathbf{e}_{\mathbf{R}}\cdot\mathbf{r})\,\mathbf{e}_{\mathbf{R}} - \mathbf{r}\Big] + O\!\left(\frac{GM\,r^2}{R^4}\right)}$$

This is suppressed by $r/R$ relative to the direct pull. For the Sun,
$r/R \approx 4.3 \times 10^{-5}$; for the Moon, $r/D \approx 1.7 \times 10^{-2}$.

Applying to each body:

$$\mathbf{a}_{\text{tidal}}^{\text{Sun}} = G M_S \left[\frac{\mathbf{R} - \mathbf{r}}{|\mathbf{R} - \mathbf{r}|^3} - \frac{\mathbf{R}}{|\mathbf{R}|^3}\right] \approx \frac{G M_S}{R^3}\Big[3\,(\mathbf{e}_{\mathbf{R}}\cdot\mathbf{r})\,\mathbf{e}_{\mathbf{R}} - \mathbf{r}\Big]$$

$$\mathbf{a}_{\text{tidal}}^{\text{Moon}} = G M_M \left[\frac{\mathbf{D} - \mathbf{r}}{|\mathbf{D} - \mathbf{r}|^3} - \frac{\mathbf{D}}{|\mathbf{D}|^3}\right] \approx \frac{G M_M}{D^3}\Big[3\,(\mathbf{e}_{\mathbf{D}}\cdot\mathbf{r})\,\mathbf{e}_{\mathbf{D}} - \mathbf{r}\Big]$$

To summarise:

- The **exact cancellation** ($\mathbf{f}(\mathbf{0}) = \mathbf{0}$) is
  algebraic: the direct pull drops out of $\mathbf{f}(\mathbf{r})$ because
  $\mathbf{f}(\mathbf{0}) = \mathbf{0}$ by construction. No expansion is needed.
- The **zeroth-order** term in the Taylor expansion ($\mathbf{R}/R^3$)
  confirms this: it cancels identically against the subtracted
  $\mathbf{R}/R^3$.
- The **first-order** term ($\propto r/R^3$) is the leading surviving
  contribution — the tidal acceleration with quadrupolar structure:
  stretching along the Earth–body axis, compression perpendicular,
  trace-free.


## Step 4: What the Sensor Reads

The sensor reads $\mathbf{F}_{\text{spring}}$, which already appears in
the boxed equation from Step 3. To isolate it, we solve that equation
for $\mathbf{F}_{\text{spring}}$ by specifying the acceleration
$\ddot{\mathbf{r}}$ of the test mass.

The test mass is stationary in the lab frame, which co-rotates with
Earth at angular velocity $\boldsymbol{\omega}$. In the non-rotating
geocentric frame used in Step 3, this corresponds to circular motion
with centripetal acceleration:

$$\ddot{\mathbf{r}} = \boldsymbol{\omega} \times (\boldsymbol{\omega} \times \mathbf{r}) = -\omega^2\,\mathbf{r}_\perp$$

where $\mathbf{r}_\perp$ is the component of $\mathbf{r}$ perpendicular
to the rotation axis (pointing outward from that axis).

Substituting into the boxed equation and solving for
$\mathbf{F}_{\text{spring}}$:

$$\boldsymbol{\omega} \times (\boldsymbol{\omega} \times \mathbf{r}) = -\frac{G M_E}{|\mathbf{r}|^3}\,\mathbf{r} + \mathbf{a}_{\text{tidal}}^{\text{Sun}} + \mathbf{a}_{\text{tidal}}^{\text{Moon}} + \frac{\mathbf{F}_{\text{spring}}}{m}$$

$$\mathbf{F}_{\text{spring}} = m\left[\frac{G M_E}{|\mathbf{r}|^3}\,\mathbf{r} + \boldsymbol{\omega} \times (\boldsymbol{\omega} \times \mathbf{r}) - \mathbf{a}_{\text{tidal}}^{\text{Sun}} - \mathbf{a}_{\text{tidal}}^{\text{Moon}}\right]$$

The three contributions have clear physical meaning:

- $GM_E\,\mathbf{r}/|\mathbf{r}|^3$: Earth's gravitational acceleration
  (the spring supports the test mass against this).
- $\boldsymbol{\omega}\times(\boldsymbol{\omega}\times\mathbf{r}) = -\omega^2\mathbf{r}_\perp$:
  centripetal acceleration from Earth's rotation, which *reduces* the
  required spring force (objects weigh less at the equator).
- $-\mathbf{a}_{\text{tidal}}$: the tidal perturbations, acting as small
  corrections.

The sensor projects this onto its measurement axis $\mathbf{e}_{\mathbf{n}}$:

$$g_{\text{measured}} = \frac{\mathbf{F}_{\text{spring}}}{m} \cdot \mathbf{e}_{\mathbf{n}}$$

**Neither the direct solar acceleration** $G M_S / R^2$ **nor the direct
lunar acceleration** $G M_M / D^2$ **appears.** These terms canceled in
Step 3 before we ever solved for $\mathbf{F}_{\text{spring}}$ — the
spring force inherits only the tidal residuals. This holds for any
$\mathbf{e}_{\mathbf{n}}$.


## Step 5: Why Tilting Doesn't Help

The claim implicitly assumes that tilting the sensor changes the
projection of a body's **direct** pull onto $\mathbf{e}_{\mathbf{n}}$:

$$g_{\text{claimed}}^{\text{Sun}} = \frac{G M_S}{R^2}\,(\mathbf{e}_{\mathbf{R}} \cdot \mathbf{e}_{\mathbf{n}}),
\qquad
g_{\text{claimed}}^{\text{Moon}} = \frac{G M_M}{D^2}\,(\mathbf{e}_{\mathbf{D}} \cdot \mathbf{e}_{\mathbf{n}})
\quad\longleftarrow\quad \text{these terms do not exist}$$

But each was canceled by the identical projection of Earth's acceleration
toward that body onto the same axis. Tilting $\mathbf{e}_{\mathbf{n}}$ rotates
all projections equally, because the cancellation is **vectorial** — it
holds component by component, in every direction simultaneously.

To make the argument concrete: the sensor housing accelerates at
$\ddot{\mathbf{R}}_E$ (toward both the Sun and Moon). The test mass, if
released, would also accelerate at $\ddot{\mathbf{R}}_E$ (plus the tidal
corrections). The spring connecting them sees only the difference — which
is the sum of the two tidal terms. This is true regardless of the spring's
orientation.


## Step 6: Numerical Comparison

### Sun

| Quantity | Formula | Value |
|----------|---------|-------|
| Direct pull | $G M_S / R^2$ | 5.93 × 10⁻³ m/s² |
| Earth's acceleration toward Sun (cancels) | same | 5.93 × 10⁻³ m/s² |
| **Tidal residual** | $\sim 2\,G M_S\, r_\oplus / R^3$ | **5.1 × 10⁻⁷ m/s²** |

$$
\frac{a_{\text{tidal}}^{\text{Sun}}}{a_{\text{direct}}^{\text{Sun}}}
= \frac{2\,r_\oplus}{R_{\text{Sun}}}
= \frac{2 \times 6.371 \times 10^6}{1.496 \times 10^{11}}
\approx 8.5 \times 10^{-5}
$$

The direct pull is **11,700×** larger than the measurable tidal signal.

### Moon

| Quantity | Formula | Value |
|----------|---------|-------|
| Direct pull | $G M_M / D^2$ | 3.32 × 10⁻⁵ m/s² |
| Earth's acceleration toward Moon (cancels) | same | 3.32 × 10⁻⁵ m/s² |
| **Tidal residual** | $\sim 2\,G M_M\, r_\oplus / D^3$ | **1.1 × 10⁻⁶ m/s²** |

$$
\frac{a_{\text{tidal}}^{\text{Moon}}}{a_{\text{direct}}^{\text{Moon}}}
= \frac{2\,r_\oplus}{D_{\text{Moon}}}
= \frac{2 \times 6.371 \times 10^6}{3.844 \times 10^{8}}
\approx 3.3 \times 10^{-2}
$$

The direct pull is **30×** larger than the measurable tidal signal. The
Moon's ratio is much less extreme than the Sun's because the Moon is
only ~60 Earth radii away, so the tidal approximation is coarser — but
the cancellation is still exact.

### Summary

| Body | Direct pull | Tidal residual | Ratio | Suppression by |
|------|-------------|----------------|-------|----------------|
| Sun  | 5.93 × 10⁻³ m/s² | 5.1 × 10⁻⁷ m/s² | 11,700× | $R_{\text{Sun}} / 2r_\oplus$ |
| Moon | 3.32 × 10⁻⁵ m/s² | 1.1 × 10⁻⁶ m/s² | 30× | $D_{\text{Moon}} / 2r_\oplus$ |

Note the inversion: the Sun's direct pull is 180× stronger than the
Moon's, but its tidal effect is 2.2× *weaker*. This is the $1/R^3$ vs
$1/R^2$ scaling — the hallmark of a tidal interaction, and a direct
consequence of the cancellation derived above.


## The Point-Mass Model Proves the Cancellation

The crucial point: the cancellation above used **only** Newton's laws
applied to point masses. No general relativity, no equivalence principle,
no appeal to "free fall" as a concept — just subtracting the equation of
motion of Earth's center from the equation of motion of the test mass.

For each external body $B$, the term $G M_B / |\mathbf{R}_B'|^2$ appears
with the same sign, same magnitude, and same direction in both equations,
so it drops out identically.

Anyone who computes $G M_S / R^2 \approx 6 \times 10^{-3}$ m/s² or
$G M_M / D^2 \approx 3.3 \times 10^{-5}$ m/s² and projects it onto a
sensor axis is computing a quantity that is real (those bodies do pull the
test mass that hard) but **unobservable** — because the sensor, the lab,
and the entire Earth are being pulled equally hard in the same direction.
The spring between the test mass and the sensor frame cannot detect a
uniform acceleration field.


## Why the Cancellation Works: Universality of Free Fall

One might object: the test mass $m$ and the Earth $M_E$ have very
different masses, so why should the Sun's effect cancel between them?

The answer is that gravitational **force** is proportional to the
accelerated mass ($F = G M_S\, m / R^2$), so the **acceleration** is
mass-independent:

$$
a = \frac{F}{m} = \frac{G M_S}{R^2}
$$

This is the equality of inertial and gravitational mass
($m_{\text{inert}} = m_{\text{grav}}$), built into Newton's law of
gravitation. The test mass and Earth's center experience the same
gravitational acceleration from the Sun (to leading order in $r/R$),
despite their mass ratio of $\sim 10^{-25}$.

This is **not** Newton's third law (action = reaction between two
interacting bodies). It is a deeper property: two *different* bodies in
the *same* external field accelerate identically. The spring connecting
them therefore reads zero from the external field — only the spatial
*variation* of the field (the tidal gradient) survives.

If any material violated this equality — if some substance fell faster
than others in the Sun's field — the direct pull *would* be measurable as
a composition-dependent residual. This is precisely what Eötvös-type
experiments test, achieving constraints at the $10^{-15}$ level. Their
null results confirm that the cancellation is exact to extraordinary
precision, and that the tidal residual is all that remains.


## Why Converting to Frequency Doesn't Help

One might attempt to circumvent the cancellation by converting the
acceleration into a frequency: use an optomechanical oscillator (a mirror
on a spring, coupled to a cavity) whose resonance frequency depends on
$g$, and beat it against a reference cavity whose frequency depends only
on its length. If the Sun's pull shifts $\omega_{\text{mech}}$ but not
$f_{\text{cav}}$, the beat note should reveal the direct pull.

This fails because the cancellation occurs *before* the readout — at the
level of what forces exist between the test mass and its mount.

### The optomechanical oscillator

A mirror (mass $m$) on a spring (constant $k$), anchored to the lab
ceiling. In equilibrium:

$$k \cdot x_{\text{eq}} = m \cdot (g_{\text{mass}} - a_{\text{anchor}})$$

where $g_{\text{mass}}$ is the total gravitational acceleration at the
mirror, and $a_{\text{anchor}}$ is the acceleration of the anchor point.
From Step 3:

$$g_{\text{mass}} - a_{\text{anchor}} = g_{\text{Earth}} - \omega^2 r_\perp + a_{\text{tidal}}^{\text{Sun}} + a_{\text{tidal}}^{\text{Moon}}$$

The Sun's direct pull appears in both $g_{\text{mass}}$ and
$a_{\text{anchor}}$ and cancels. The equilibrium position — and therefore
$\omega_{\text{mech}} \propto \sqrt{g_{\text{eff}}/L}$ — contains only
the tidal part.

### The reference cavity

A Fabry-Perot cavity has $f = mc/(2L)$. The spacer length $L$ is set by
electromagnetic bond forces. The Sun's uniform field accelerates every
atom equally — no deformation. Only the tidal gradient across $L$ would
stretch it: $\delta L/L \sim (GM_S/R^3)\,L^2/(E/\rho) \sim 10^{-23}$ — unmeasurable.

### The beat note

$\Delta f = f_{\text{mech}} - f_{\text{cav}}$ contains tidal effects from
the oscillator minus negligible tidal deformation of the cavity. The
direct pull $GM_S/R^2$ appears in neither channel. Converting to
frequency does not change *what physical quantity* is being measured.


## What Frequency Measurements *Can* Detect

While no local mechanical measurement can access the direct pull, the
gravitational **potential** (as opposed to force) can be detected through
a fundamentally different channel: the gravitational redshift.

An atomic clock's tick rate depends on the spacetime metric:

$$\frac{d\tau}{dt} = \sqrt{-g_{00}} \approx 1 + \frac{\Phi}{c^2}$$

This is not a force between co-accelerating parts — it is a property of
spacetime itself. However, the equivalence principle guarantees that *all*
local physics is shifted equally. A single clock has nothing local to
compare against. To detect the potential, one needs a **non-local**
reference.

### Clock vs. distant reference (pulsar)

The Sun's potential varies annually due to orbital eccentricity:

$$\frac{\Delta f}{f} = \frac{GM_S}{c^2}\left(\frac{1}{R_{\text{peri}}} - \frac{1}{R_{\text{aph}}}\right) \approx 3.3 \times 10^{-10}$$

This is the **full, direct** solar potential — not tidal, not suppressed
by $r_\oplus/R$. Pulsar timing arrays detect exactly this as the
"Einstein delay" (~1.66 ms annual amplitude).

### Two terrestrial clocks

Two clocks at different locations see different solar potentials, but the
difference is only the tidal potential: $\delta f/f \sim 10^{-17}$.
Modern optical lattice clocks are just reaching this regime.

### Hierarchy of measurements

| Method | Measures | Sun contribution | Magnitude |
|--------|----------|-----------------|-----------|
| Accelerometer | $-\nabla\Phi$ (force) | Tidal only | $5 \times 10^{-7}$ m/s² |
| Optomech. vs. ref. cavity | $-\nabla\Phi$ (force as frequency) | Tidal only | $5 \times 10^{-7}$ m/s² |
| Two local clocks | $\Delta\Phi$ (potential difference) | Tidal potential | $\Delta f/f \sim 10^{-17}$ |
| Clock vs. pulsar | $\Phi$ (absolute potential) | **Full direct** | $\Delta f/f \sim 3 \times 10^{-10}$ |

The dividing line: any instrument where two parts are co-accelerating sees
the cancellation. To escape it, one end of the measurement must be outside
the freely-falling frame.


## What IS Observable by Tilting

Tilting the sensor does reveal genuinely interesting physics, but at the
tidal scale (~10⁻⁷ m/s²), not at the direct-pull scale (~10⁻³ m/s²):

- **Horizontal tidal acceleration.** A vertical sensor sees only the
  vertical component of $\mathbf{a}_{\text{tidal}}$. A tilted sensor
  picks up horizontal components, which have comparable magnitude but
  different time dependence.

- **Angular separation of static and tidal signals.** Normal gravity is
  purely vertical (projects as $g_0 \cos\theta$), while the tidal vector
  has East/North/Up components with independent angular signatures. Scanning
  the sensor azimuth at $\theta = 90°$ gives a pure tidal signal with zero
  static background.

- **Reconstruction of the full tidal vector.** Three non-coplanar
  measurements at one instant determine all four unknowns ($g_0$,
  $a_E$, $a_N$, $a_U$) without a time series.

These effects are real, measurable, and interesting — but they are four
orders of magnitude smaller than the claimed signal.

# Why a Tilted Sensor Cannot Measure Direct Gravitational Pull

## The Claim

> A sensor tilted relative to the local vertical can detect the direct
> gravitational pull of the Sun (~6 × 10⁻³ m/s²) or the Moon
> (~3.3 × 10⁻⁵ m/s²) as a time-varying signal, because the projection
> of the gravitational acceleration onto the sensor axis changes as
> these bodies move across the sky.

This is incorrect. The direct pull cancels exactly against the
acceleration of the laboratory, leaving only the tidal residual —
regardless of sensor orientation. The proof requires nothing beyond
Newton's laws applied to point masses.


## Inertial and Non-Inertial Frames

An **inertial frame** is one in which Newton's second law
$\mathbf{F} = m\,\mathbf{a}$ holds without correction. A body subject
to no forces moves in a straight line at constant velocity.

A **non-inertial frame** is one that accelerates. In such a frame,
Newton's laws require fictitious forces — centrifugal, Coriolis, etc. —
to compensate for the acceleration of the frame itself.

### The laboratory is non-inertial

A laboratory on Earth's surface is non-inertial for two reasons:

1. **Earth rotates**, producing centrifugal and Coriolis
   pseudo-accelerations (~0.03 m/s² at the equator).

2. **Earth freely falls toward the Sun and Moon.** The entire Earth —
   lab, sensor, test mass, and all — accelerates toward the Sun at
   $GM_S/R^2 \approx 5.93 \times 10^{-3}$ m/s² and toward the Moon at
   $GM_M/D^2 \approx 3.3 \times 10^{-5}$ m/s².

The second point is the origin of the cancellation. When you write
equations in the lab frame, you have already subtracted these
accelerations from every term. A sensor physically implements this
subtraction: both ends of the spring accelerate identically in the
Sun's field, so the spring reads zero from that field.

To see this without hidden subtractions, we work in a true inertial
frame and carry every term explicitly.


## The Proof

### Setup

Four point masses in an inertial frame:

| Body | Mass | Position |
|------|------|----------|
| Sun | $M_S$ | $\mathbf{R}_S(t)$ |
| Moon | $M_M$ | $\mathbf{R}_M(t)$ |
| Earth | $M_E$ | $\mathbf{R}_E(t)$ |
| Test mass | $m$ | $\mathbf{x}(t)$ |

The test mass is connected to the sensor by a spring (or equivalent
restoring mechanism). The sensor housing is rigidly attached to the
Earth. The sensor reads the spring force $\mathbf{F}_{\text{spring}}$ —
the non-gravitational force required to keep the test mass co-moving
with the lab.

### Step 1: Test mass equation of motion

In the inertial frame, Newton's second law for the test mass:

$$m\,\ddot{\mathbf{x}}
  = -\frac{G M_E\, m}{|\mathbf{x} - \mathbf{R}_E|^3}
    (\mathbf{x} - \mathbf{R}_E)
  -\frac{G M_S\, m}{|\mathbf{x} - \mathbf{R}_S|^3}
    (\mathbf{x} - \mathbf{R}_S)
  -\frac{G M_M\, m}{|\mathbf{x} - \mathbf{R}_M|^3}
    (\mathbf{x} - \mathbf{R}_M)
  + \mathbf{F}_{\text{spring}}$$

### Step 2: Earth's equation of motion

Earth's center accelerates toward the Sun and Moon:

$$\ddot{\mathbf{R}}_E
  = -\frac{G M_S}{|\mathbf{R}_E - \mathbf{R}_S|^3}
    (\mathbf{R}_E - \mathbf{R}_S)
  - \frac{G M_M}{|\mathbf{R}_E - \mathbf{R}_M|^3}
    (\mathbf{R}_E - \mathbf{R}_M)$$

Every object rigidly attached to the Earth — including the sensor
housing — shares this acceleration.

### Step 3: Subtract to get the relative motion

Define $\mathbf{r} \equiv \mathbf{x} - \mathbf{R}_E$ and the
geocentric positions
$\mathbf{R} \equiv \mathbf{R}_S - \mathbf{R}_E$,
$\mathbf{D} \equiv \mathbf{R}_M - \mathbf{R}_E$. Then
$\ddot{\mathbf{r}} = \ddot{\mathbf{x}} - \ddot{\mathbf{R}}_E$ gives:

$$\boxed{
m\,\ddot{\mathbf{r}}
  = -\frac{G M_E\, m}{|\mathbf{r}|^3}\,\mathbf{r}
  + G M_S\, m \left[
    \frac{\mathbf{R} - \mathbf{r}}{|\mathbf{R} - \mathbf{r}|^3}
    - \frac{\mathbf{R}}{R^3}
  \right]
  + G M_M\, m \left[
    \frac{\mathbf{D} - \mathbf{r}}{|\mathbf{D} - \mathbf{r}|^3}
    - \frac{\mathbf{D}}{D^3}
  \right]
  + \mathbf{F}_{\text{spring}}
}$$

### The cancellation

Each bracketed term has the form

$$\mathbf{f}(\mathbf{r})
  = \frac{\mathbf{R} - \mathbf{r}}{|\mathbf{R} - \mathbf{r}|^3}
  - \frac{\mathbf{R}}{R^3}$$

This is the Sun's gravitational field at the test mass *minus* the
field at Earth's center. At $\mathbf{r} = \mathbf{0}$:

$$\mathbf{f}(\mathbf{0})
  = \frac{\mathbf{R}}{R^3} - \frac{\mathbf{R}}{R^3}
  = \mathbf{0}
  \qquad\text{(exactly)}$$

The direct pull — $GM_S\,\mathbf{R}/R^3$, the uniform field that
accelerates the entire Earth at $5.93 \times 10^{-3}$ m/s² — appears
with equal magnitude and opposite sign in the test mass and Earth
equations, and cancels identically. No approximation is involved. The
same holds for the Moon with $\mathbf{D}$ replacing $\mathbf{R}$.

What survives is $\mathbf{f}(\mathbf{r}) - \mathbf{f}(\mathbf{0})$:
the *variation* of the gravitational field across the baseline
$\mathbf{r}$ — the **tidal acceleration**. To leading order in $r/R$:

$$\boxed{
\mathbf{a}_{\text{tidal}}
  = \frac{GM}{R^3}\Big[
    3\,(\mathbf{e}_{\mathbf{R}}\cdot\mathbf{r})\,\mathbf{e}_{\mathbf{R}} - \mathbf{r}
  \Big]
  + O\!\left(\frac{GM\,r^2}{R^4}\right)
}$$

This is suppressed by a factor $r/R$ relative to the direct pull.

### Derivation of the tidal formula

Expand $|\mathbf{R} - \mathbf{r}|^{-3}$ for $r \ll R$. Write

$$|\mathbf{R} - \mathbf{r}|^2
  = R^2\left(1
    - 2\,\frac{\mathbf{e}_{\mathbf{R}}\cdot\mathbf{r}}{R}
    + \frac{r^2}{R^2}\right)
  \equiv R^2(1 - \epsilon)$$

with $\epsilon = 2\,\mathbf{e}_{\mathbf{R}}\cdot\mathbf{r}/R - r^2/R^2 = O(r/R)$. Then

$$\frac{1}{|\mathbf{R} - \mathbf{r}|^3}
  = \frac{1}{R^3}(1-\epsilon)^{-3/2}
  = \frac{1}{R^3}\left(1
    + \frac{3\,\mathbf{e}_{\mathbf{R}}\cdot\mathbf{r}}{R}
    + O\!\left(\frac{r^2}{R^2}\right)\right)$$

Multiply by $(\mathbf{R} - \mathbf{r})$ and keep terms through first
order:

$$\frac{\mathbf{R} - \mathbf{r}}{|\mathbf{R} - \mathbf{r}|^3}
  = \frac{1}{R^3}\Big[
    \mathbf{R} + 3(\mathbf{e}_{\mathbf{R}}\cdot\mathbf{r})\,\mathbf{e}_{\mathbf{R}} - \mathbf{r}
  \Big]
  + O\!\left(\frac{r^2}{R^4}\right)$$

Subtract $\mathbf{R}/R^3$:

$$\mathbf{f}(\mathbf{r})
  = \frac{\mathbf{R} - \mathbf{r}}{|\mathbf{R} - \mathbf{r}|^3}
  - \frac{\mathbf{R}}{R^3}
  = \frac{1}{R^3}\Big[
    3(\mathbf{e}_{\mathbf{R}}\cdot\mathbf{r})\,\mathbf{e}_{\mathbf{R}} - \mathbf{r}
  \Big]
  + O\!\left(\frac{r^2}{R^4}\right)$$

The $\mathbf{R}/R^3$ terms cancel identically — confirming that the
direct pull drops out — and the leading surviving term is the tidal
quadrupole. $\square$


## Transformation to the Laboratory Frame

The geocentric equation is written in a non-rotating frame. The
laboratory rotates with the Earth at angular velocity
$\boldsymbol{\omega}$. This is a second source of non-inertiality
(see Section 1).

### Rotating-frame equation of motion

Let $\mathbf{r}'$ be the test mass position in the lab frame, related
to its geocentric position by
$\mathbf{r} = \mathcal{R}(t)\,\mathbf{r}'$
where $\mathcal{R}(t)$ is the time-dependent rotation matrix.
Differentiating twice:

$$\ddot{\mathbf{r}}
  = \ddot{\mathbf{r}}'
  + 2\,\boldsymbol{\omega} \times \dot{\mathbf{r}}'
  + \boldsymbol{\omega} \times (\boldsymbol{\omega} \times \mathbf{r}')$$

(the Euler term $\dot{\boldsymbol{\omega}} \times \mathbf{r}'$ vanishes
for uniform rotation). Substituting into the geocentric equation and
expressing all vectors in the lab basis:

$$\boxed{
\ddot{\mathbf{r}}'
  = \underbrace{-\frac{G M_E}{|\mathbf{r}'|^3}\,\mathbf{r}'
    }_{\text{Earth's gravity}}
  + \underbrace{\mathbf{a}_{\text{tidal}}'
    }_{\text{tidal}}
  \underbrace{- \boldsymbol{\omega}
    \times (\boldsymbol{\omega} \times \mathbf{r}')
    }_{\text{centrifugal}}
  \underbrace{- 2\,\boldsymbol{\omega}
    \times \dot{\mathbf{r}}'
    }_{\text{Coriolis}}
  + \frac{\mathbf{F}_{\text{spring}}'}{m}
}$$

The direct solar and lunar pulls do **not** reappear — they canceled
in the geocentric equation before the frame change, and rotating the
coordinate basis cannot restore a term that is already zero.

### What the sensor reads

The test mass sits at rest in the lab:
$\dot{\mathbf{r}}' = \ddot{\mathbf{r}}' = \mathbf{0}$. The Coriolis
term vanishes. Setting the left-hand side to zero:

$$\frac{\mathbf{F}_{\text{spring}}'}{m}
  = \frac{G M_E}{|\mathbf{r}'|^3}\,\mathbf{r}'
  + \boldsymbol{\omega} \times (\boldsymbol{\omega} \times \mathbf{r}')
  - \mathbf{a}_{\text{tidal}}'$$

The sensor projects this onto its measurement axis
$\mathbf{e}_{\mathbf{n}}$:

$$g_{\text{measured}}
  = \frac{\mathbf{F}_{\text{spring}}'}{m} \cdot \mathbf{e}_{\mathbf{n}}$$

The three contributions are:

- $GM_E\,\mathbf{r}'/|\mathbf{r}'|^3$: Earth's gravity (the dominant,
  static term).
- $\boldsymbol{\omega}\times(\boldsymbol{\omega}\times\mathbf{r}')
  = -\omega^2\mathbf{r}'_\perp$: the centrifugal reduction (objects
  weigh less at the equator by ~0.3%).
- $-\mathbf{a}_{\text{tidal}}'$: the tidal perturbation from Sun and
  Moon, of order $10^{-7}$ m/s².

The direct pull $GM_S/R^2$ and $GM_M/D^2$ appear nowhere. This holds
for any $\mathbf{e}_{\mathbf{n}}$.

### Size of each term

| Term | Expression | Magnitude |
|------|-----------|-----------|
| Earth's gravity | $GM_E/r_\oplus^2$ | 9.81 m/s² |
| Centrifugal (equator) | $\omega^2 r_\oplus$ | 3.4 × 10⁻² m/s² |
| Lunar tide | $GM_M\, r_\oplus / D^3$ | 1.1 × 10⁻⁶ m/s² |
| Solar tide | $GM_S\, r_\oplus / R^3$ | 5.1 × 10⁻⁷ m/s² |
| Solar direct pull (canceled) | $GM_S / R^2$ | 5.93 × 10⁻³ m/s² |

The sensor reads a superposition of the first four terms. The fifth —
the direct pull — has been subtracted out by the physics of the
measurement: the lab accelerates with the Earth, and the spring cannot
see it.


## Why Sensor Orientation Is Irrelevant

The cancellation is **vectorial**: the direct pull vanishes as a
three-component vector, not as a particular scalar projection.
Projecting onto any measurement axis $\mathbf{e}_{\mathbf{n}}$:

$$\mathbf{e}_{\mathbf{n}} \cdot \mathbf{f}(\mathbf{0}) = \mathbf{e}_{\mathbf{n}} \cdot \mathbf{0} = 0
  \qquad\text{for all } \mathbf{e}_{\mathbf{n}}$$

Tilting the sensor changes $\mathbf{e}_{\mathbf{n}}$ but cannot make
zero non-zero. The spring force along any axis contains only the tidal
terms.

Physically: the sensor housing and the test mass both accelerate at
$GM_S/R^2$ toward the Sun. The spring connecting them measures only the
*difference* between the accelerations of its two endpoints. A uniform
field produces no difference, regardless of the spring's orientation.

### Sweeping the sensor orientation

Parameterise the measurement axis by zenith angle $\theta$ and
azimuth $\varphi$:

$$\mathbf{e}_{\mathbf{n}}(\theta,\varphi)
  = \sin\theta\cos\varphi\;\mathbf{e}_E
  + \sin\theta\sin\varphi\;\mathbf{e}_N
  + \cos\theta\;\mathbf{e}_U$$

where $\mathbf{e}_E$, $\mathbf{e}_N$, $\mathbf{e}_U$ are the local
East-North-Up unit vectors.

The measured acceleration as a function of orientation is

$$g(\theta,\varphi)
  = g_0\cos\theta
  - a_E'\sin\theta\cos\varphi
  - a_N'\sin\theta\sin\varphi
  - a_U'\cos\theta$$

where $g_0 \equiv GM_E/r_\oplus^2 - \omega^2 r_\perp$ is the static
effective gravity, and $a_E'$, $a_N'$, $a_U'$ are the East, North,
and Up components of the tidal acceleration.

Sweeping $\theta$ and $\varphi$ traces out the full angular
dependence. If the direct pull $GM_S/R^2$ were present, it would
contribute a term $\propto \sin\theta\cos(\varphi - \varphi_\odot)$
with amplitude ~6 × 10⁻³ m/s². Instead, the horizontal terms carry
amplitudes $a_E', a_N' \sim 10^{-7}$ m/s² — the tidal values — and
the vertical tidal correction $a_U'$ is of the same order.

An angular scan therefore confirms the cancellation directly: one
measures $g_0\cos\theta$ plus tidal corrections of order 10⁻⁷ m/s²,
with no 10⁻³ m/s² sinusoidal component at any orientation.


## Numerical Scale

| Body | Direct pull (cancels) | Tidal residual (observable) | Ratio |
|------|----------------------|---------------------------|-------|
| Sun  | 5.93 × 10⁻³ m/s² | 5.1 × 10⁻⁷ m/s² | 11,700× |
| Moon | 3.32 × 10⁻⁵ m/s² | 1.1 × 10⁻⁶ m/s² | 30× |

The Sun's direct pull is 180× stronger than the Moon's, but its tidal
effect is 2.2× *weaker* ($1/R^3$ vs. $1/R^2$ scaling). This inversion
is the hallmark of tidal physics, and a direct consequence of the
cancellation.


## What IS Observable

Tilting the sensor accesses different components of the *tidal* field,
not the direct pull:

- A vertical sensor sees the vertical tidal component
  (~5 × 10⁻⁷ m/s², semi-diurnal).
- A horizontal sensor picks up horizontal tidal components (diurnal).
- Three non-coplanar sensors reconstruct the full tidal tensor.

These effects are real and measurable, but four orders of magnitude
smaller than the claimed signal.

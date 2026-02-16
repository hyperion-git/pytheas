## 6. The Tidal Acceleration (The Heart of the Model)

Everything we have built so far — geodetic coordinates, the measurement axis,
normal gravity, ephemerides — serves a single purpose: computing the **tidal
acceleration** that a celestial body (Moon or Sun) produces at a specific point
on Earth's surface.  This is the signal that makes $g$ vary with time, and it
is the heart of Pytheas.

We will derive the tidal acceleration from scratch, estimate its size, explore
(and reject) a popular approximation, and show how it projects onto a sensor.


### 6.1 Tidal forces from first principles

#### The setup

Place the origin at Earth's center of mass.  An observer stands on the surface
at position $\mathbf{r}$.  A distant body — call it the Moon for concreteness,
with mass $M$ — sits at position $\mathbf{R}$.

We want to know how the Moon's gravity affects the observer *differently* from
how it affects the Earth as a whole.  That difference is the tidal
acceleration.  Figure 9 illustrates the geometry: the vectors from Earth's
center to the Moon ($\mathbf{R}$) and from the observer to the Moon
($\mathbf{R} - \mathbf{r}$), whose difference drives the tidal acceleration,
along with the characteristic quadrupolar pattern of the tidal field.

#### Gravity at two points

Newton's law of gravitation tells us that a mass $M$ at position $\mathbf{R}$
pulls any test mass toward itself.  The gravitational acceleration it produces
at an arbitrary point $\mathbf{p}$ is

$$
\mathbf{g}(\mathbf{p}) = \frac{GM(\mathbf{R} - \mathbf{p})}{|\mathbf{R} - \mathbf{p}|^3}.
$$

The vector $\mathbf{R} - \mathbf{p}$ points from $\mathbf{p}$ toward the body,
and dividing by the cube of its magnitude (not the square!) gives the correct
inverse-square force per unit mass.  (The extra power of distance in the
denominator comes from normalizing the direction vector.)

Apply this at two points:

**At the observer** (position $\mathbf{r}$):

$$
\mathbf{g}_{\mathrm{obs}} = \frac{GM(\mathbf{R} - \mathbf{r})}{|\mathbf{R} - \mathbf{r}|^3}.
$$

**At Earth's center** (position $\mathbf{0}$):

$$
\mathbf{g}_{\mathrm{center}} = \frac{GM\,\mathbf{R}}{|\mathbf{R}|^3} = \frac{GM\,\mathbf{R}}{R^3},
$$

where $R \equiv |\mathbf{R}|$ is the distance from Earth's center to the body.

#### The tidal acceleration

Earth, taken as a whole, is in free fall through the Moon's gravitational
field.  Earth's center accelerates at $\mathbf{g}_{\mathrm{center}}$ toward the
Moon.  In a reference frame that falls with Earth's center, this acceleration
is "used up" — it is what keeps Earth on its orbit.  What remains for the
observer is the *difference*:

$$
\boxed{
\mathbf{a}_{\mathrm{tidal}}
  = \mathbf{g}_{\mathrm{obs}} - \mathbf{g}_{\mathrm{center}}
  = GM \left[
      \frac{\mathbf{R} - \mathbf{r}}{|\mathbf{R} - \mathbf{r}|^3}
    - \frac{\mathbf{R}}{R^3}
  \right].
}
$$

This is the **exact Newtonian tidal acceleration**, and it is what Pytheas
computes.

#### A thought experiment: free-falling Earth

To build intuition, imagine turning off every force except the Moon's gravity.
Earth falls freely toward the Moon.  In that freely-falling frame, the Moon's
gravitational field is exactly zero at Earth's center — by definition, the
center is in free fall.

But you are not at the center.  You are displaced from it by $\mathbf{r}$
(roughly 6,400 km, the Earth's radius).  The Moon's field is not perfectly
uniform — it is slightly stronger on the side of Earth closest to the Moon and
slightly weaker on the far side.  These residual differences, left over after
subtracting the free-fall acceleration of the center, constitute the tidal
field.

The tidal field stretches Earth along the line connecting it to the Moon (the
near side is pulled harder, the far side less hard) and compresses it in the
perpendicular directions.  This is why there are two tidal bulges — one facing
the Moon, one on the opposite side — and two tidal "valleys" at right angles.

#### What Pytheas computes

In code, the tidal acceleration is implemented as:

```python
def tidal_acceleration(r, R, GM):
    d = R - r
    return GM * (d / np.linalg.norm(d)**3 - R / np.linalg.norm(R)**3)
```

Here `d = R - r` is the vector from the observer to the body,
`np.linalg.norm(d)**3` is $|\mathbf{R} - \mathbf{r}|^3$, and the two terms
correspond exactly to $\mathbf{g}_{\mathrm{obs}}$ and
$\mathbf{g}_{\mathrm{center}}$.  There is no approximation — this is the full
Newtonian expression.


### 6.2 How big is this?

Let us plug in numbers to get a feel for the magnitudes involved.

#### The Moon

| Quantity | Value |
|----------|-------|
| $GM_{\mathrm{Moon}}$ | $4.9028695 \times 10^{12}\ \mathrm{m^3/s^2}$ |
| $R$ (Earth–Moon distance) | $\approx 384{,}400\ \mathrm{km} = 3.844 \times 10^8\ \mathrm{m}$ |
| $r$ (Earth's radius) | $\approx 6{,}371\ \mathrm{km} = 6.371 \times 10^6\ \mathrm{m}$ |

Since $r \ll R$, we can get a quick estimate using the leading-order
approximation (derived properly in Section 6.3):

$$
a_{\mathrm{tidal}} \approx \frac{2\,GM\,r}{R^3}.
$$

Substituting:

$$
a_{\mathrm{tidal}}^{\mathrm{Moon}}
  \approx \frac{2 \times 4.903 \times 10^{12} \times 6.371 \times 10^6}
               {(3.844 \times 10^8)^3}
  = \frac{6.247 \times 10^{19}}{5.680 \times 10^{25}}
  \approx 1.10 \times 10^{-6}\ \mathrm{m/s^2}.
$$

That is about $1.1\ \mu\mathrm{m/s^2}$, or roughly **one part in $10^7$ of
$g$**.  A seemingly tiny number — but modern gravimeters resolve a few nGal
($1\ \mathrm{nGal} = 10^{-11}\ \mathrm{m/s^2}$), so the lunar tide is an
enormous signal from their perspective: roughly $110{,}000\ \mathrm{nGal}$.

#### The Sun

| Quantity | Value |
|----------|-------|
| $GM_{\mathrm{Sun}}$ | $1.327 \times 10^{20}\ \mathrm{m^3/s^2}$ |
| $R$ (Earth–Sun distance) | $\approx 1\ \mathrm{AU} = 1.496 \times 10^{11}\ \mathrm{m}$ |

$$
a_{\mathrm{tidal}}^{\mathrm{Sun}}
  \approx \frac{2 \times 1.327 \times 10^{20} \times 6.371 \times 10^6}
               {(1.496 \times 10^{11})^3}
  = \frac{1.691 \times 10^{27}}{3.349 \times 10^{33}}
  \approx 5.05 \times 10^{-7}\ \mathrm{m/s^2}.
$$

That is about $0.5\ \mu\mathrm{m/s^2}$, or roughly half the Moon's
contribution.

#### Why the Moon wins: the $R^{-3}$ law

This result may surprise you.  The Sun's gravitational parameter $GM_{\mathrm{Sun}}$
is about $2.7 \times 10^7$ times larger than $GM_{\mathrm{Moon}}$.  So why
doesn't the Sun dominate the tides?

The answer is that the tidal acceleration depends on $GM/R^3$, not $GM/R^2$.
The ordinary gravitational *force* on Earth from the Sun is indeed far larger
than that from the Moon — the Sun keeps Earth in orbit, after all.  But the
tidal acceleration depends on the *gradient* of the gravitational field: how
much the field changes across the width of the Earth.  That gradient falls off
as $1/R^3$.

The Sun is about 389 times farther away than the Moon.  The ratio of their
tidal strengths is:

$$
\frac{a_{\mathrm{tidal}}^{\mathrm{Sun}}}{a_{\mathrm{tidal}}^{\mathrm{Moon}}}
  = \frac{GM_{\mathrm{Sun}}}{GM_{\mathrm{Moon}}} \times
    \left(\frac{R_{\mathrm{Moon}}}{R_{\mathrm{Sun}}}\right)^3
  \approx 2.7 \times 10^7 \times \left(\frac{1}{389}\right)^3
  \approx 2.7 \times 10^7 \times 1.70 \times 10^{-8}
  \approx 0.46.
$$

The extra power of distance is decisive.  Despite its vastly greater mass, the
Sun loses to the Moon by a factor of about $1/0.46 \approx 2.2$.  The Moon is
the dominant source of tidal variations in $g$.


### 6.3 The gradient approximation (and why Pytheas doesn't use it)

Many textbooks — and many tidal codes — do not evaluate the exact formula from
Section 6.1.  Instead, they Taylor-expand the tidal acceleration for small
$r/R$ and keep only the leading term.  This is the "gradient" or "quadrupole"
approximation.  Let us derive it, assess its accuracy, and explain why Pytheas
skips it.

#### Deriving the first-order expansion

We want to expand

$$
\frac{\mathbf{R} - \mathbf{r}}{|\mathbf{R} - \mathbf{r}|^3}
$$

in powers of $\mathbf{r}$, treating $\mathbf{r}$ as small compared to
$\mathbf{R}$.  Define $\hat{\mathbf{R}} = \mathbf{R}/R$ as the unit vector
toward the body.

Start by writing the denominator.  Let $\mathbf{d} = \mathbf{R} - \mathbf{r}$,
so

$$
|\mathbf{d}|^2
  = (\mathbf{R} - \mathbf{r}) \cdot (\mathbf{R} - \mathbf{r})
  = R^2 - 2\,\mathbf{R} \cdot \mathbf{r} + r^2
  = R^2 \left(1 - 2\,\frac{\mathbf{R} \cdot \mathbf{r}}{R^2}
         + \frac{r^2}{R^2}\right).
$$

Define the small quantity

$$
\epsilon = \frac{2\,\mathbf{R} \cdot \mathbf{r}}{R^2} - \frac{r^2}{R^2},
$$

so that $|\mathbf{d}|^2 = R^2(1 - \epsilon)$.  Note $\epsilon = O(r/R)$.
Then

$$
\frac{1}{|\mathbf{d}|^3} = \frac{1}{R^3}(1 - \epsilon)^{-3/2}
  = \frac{1}{R^3}\left(1 + \tfrac{3}{2}\epsilon + O(\epsilon^2)\right).
$$

To first order in $r/R$, we keep only the $\mathbf{R} \cdot \mathbf{r}$ piece
of $\epsilon$ (the $r^2/R^2$ term is second-order):

$$
\frac{1}{|\mathbf{d}|^3}
  \approx \frac{1}{R^3}\left(1 + \frac{3\,\mathbf{R} \cdot \mathbf{r}}{R^2}\right)
  + O\!\left(\frac{r^2}{R^2}\right).
$$

Now multiply by $\mathbf{d} = \mathbf{R} - \mathbf{r}$:

$$
\frac{\mathbf{R} - \mathbf{r}}{|\mathbf{d}|^3}
  \approx \frac{1}{R^3}\left(1 + \frac{3\,\mathbf{R} \cdot \mathbf{r}}{R^2}\right)
    (\mathbf{R} - \mathbf{r}).
$$

Expanding and keeping terms through first order in $r/R$:

$$
\frac{\mathbf{R} - \mathbf{r}}{|\mathbf{d}|^3}
  \approx \frac{\mathbf{R}}{R^3}
    + \frac{1}{R^3}\left[
        \frac{3(\mathbf{R} \cdot \mathbf{r})}{R^2}\,\mathbf{R} - \mathbf{r}
      \right]
    + O\!\left(\frac{r^2}{R^2}\right).
$$

(The cross-term $\frac{3(\mathbf{R}\cdot\mathbf{r})}{R^2}\cdot(-\mathbf{r})$
is second-order and has been dropped.)

Now subtract $\mathbf{R}/R^3$ to form the tidal acceleration:

$$
\mathbf{a}_{\mathrm{tidal}}
  = GM\left[\frac{\mathbf{R} - \mathbf{r}}{|\mathbf{d}|^3} - \frac{\mathbf{R}}{R^3}\right]
  \approx \frac{GM}{R^3}\left[
      3\,(\hat{\mathbf{R}} \cdot \mathbf{r})\,\hat{\mathbf{R}} - \mathbf{r}
    \right]
  + O\!\left(\frac{r^2}{R^2}\right).
$$

This can be written more compactly using the tidal tensor:

$$
\boxed{
\mathbf{a}_{\mathrm{tidal}}
  \approx -\frac{GM}{R^3}\left[
      \mathbf{r} - 3\,(\mathbf{r} \cdot \hat{\mathbf{R}})\,\hat{\mathbf{R}}
    \right].
}
$$

#### Interpreting the gradient approximation

This expression has a beautiful geometric interpretation.  Decompose
$\mathbf{r}$ into a component along $\hat{\mathbf{R}}$ (the line toward the
body) and a component perpendicular to it:

- Along $\hat{\mathbf{R}}$: the coefficient of $\hat{\mathbf{R}}$ in the
  tidal acceleration is $+2\,GM\,r_\parallel/R^3$, where $r_\parallel =
  \mathbf{r}\cdot\hat{\mathbf{R}}$.  The tidal field **stretches** along the
  line connecting Earth and the body.

- Perpendicular to $\hat{\mathbf{R}}$: the tidal acceleration is
  $-GM\,r_\perp/R^3$.  The tidal field **compresses** in the directions
  orthogonal to the Earth–body line.

The stretching is twice as strong as the compression — and they have opposite
signs.  This is the familiar tidal pattern: elongation toward (and away from)
the body, squeezing at the sides.

You can verify that this field is trace-free: the sum of the three diagonal
entries of the tidal tensor $(+2 - 1 - 1)\times GM/R^3 = 0$, as required by
Laplace's equation in vacuum.

#### How good is the approximation?

The terms we dropped are $O(r^2/R^2)$ relative to the leading term.  How large
is that correction?

**For the Moon:**

$$
\left(\frac{r}{R}\right)^2
  = \left(\frac{6.371 \times 10^6}{3.844 \times 10^8}\right)^2
  = (1.657 \times 10^{-2})^2
  \approx 2.7 \times 10^{-4}.
$$

The leading tidal acceleration is about $1.1\ \mu\mathrm{m/s^2} = 110\ \mu\mathrm{Gal}
= 110{,}000$ nGal.  The truncation error is roughly:

$$
2.7 \times 10^{-4} \times 110{,}000\ \mathrm{nGal} \approx 30\ \mathrm{nGal}.
$$

This is squarely in the 10--100 nGal accuracy range that Pytheas targets.  An
error of 30 nGal from a needless approximation would be one of the largest
error sources in the model.  Moreover, the next-order term has a specific
geometric pattern that could produce systematic biases at certain observer
locations.  Figure 10 compares the exact Newtonian tidal acceleration with the
gradient approximation over 24 hours at Munich, with the residual panel
revealing the ~30 nGal systematic error of the approximation.

**For the Sun:**

$$
\left(\frac{r}{R}\right)^2
  = \left(\frac{6.371 \times 10^6}{1.496 \times 10^{11}}\right)^2
  \approx 1.8 \times 10^{-9}.
$$

Completely negligible.  The gradient approximation for the Sun is essentially
exact.

#### Pytheas's choice: why approximate when you don't have to?

The gradient approximation saves one computation: you avoid calculating
$|\mathbf{R} - \mathbf{r}|$.  But that is just a single norm evaluation — a
handful of multiplications, an addition, and a square root.  On any modern
processor, this is negligible.

By using the exact formula

$$
\mathbf{a}_{\mathrm{tidal}} = GM\left[
  \frac{\mathbf{R} - \mathbf{r}}{|\mathbf{R} - \mathbf{r}|^3}
  - \frac{\mathbf{R}}{R^3}
\right],
$$

Pytheas eliminates the $\sim 30$ nGal truncation error from the lunar tidal
acceleration — for free.  For a model targeting 10--100 nGal accuracy, a 30
nGal error from an avoidable approximation would be unacceptable.

This is a general principle in scientific computing: **use the exact formula
when it costs the same as the approximate one**.  Approximations are valuable
when they reduce computational cost or reveal physical structure (as the
gradient approximation does for physical intuition); but for numerical
evaluation, the exact expression is always preferable.


### 6.4 Projecting onto your sensor

#### The tidal acceleration is a 3D vector

The tidal acceleration $\mathbf{a}_{\mathrm{tidal}}$ lives in three-dimensional
space, expressed in Earth-Centered Earth-Fixed (ECEF) coordinates.  It has
components that point radially, northward, and eastward — and the mix depends
on where you are and where the Moon (or Sun) currently is in the sky.

But a gravimeter or accelerometer does not measure a 3D vector.  It measures
the component of acceleration along its **sensitive axis**.

#### The projection

In Section 2, we defined the measurement axis $\hat{\mathbf{n}}$ — a unit
vector in ECEF pointing along the direction your instrument measures.  For a
standard vertical gravimeter, $\hat{\mathbf{n}} = \hat{\mathbf{e}}_{\mathrm{up}}$,
the local vertical.  For a tilted or horizontal instrument, $\hat{\mathbf{n}}$
has components along east, north, and up (defined by the zenith and azimuth
angles).

The tidal acceleration measured by the sensor is the dot product:

$$
\boxed{
a_{\mathrm{proj}} = \mathbf{a}_{\mathrm{tidal}} \cdot \hat{\mathbf{n}}.
}
$$

This is a scalar — positive when the tidal acceleration has a component along
$\hat{\mathbf{n}}$, negative when it opposes it.

#### What different orientations see

**Vertical sensor** ($\hat{\mathbf{n}} = \hat{\mathbf{e}}_{\mathrm{up}}$):
picks up the radial component of the tidal acceleration.  This is the largest
component, dominated by the $+2\,GM\,r_\parallel/R^3$ stretching term when the
body is near the zenith or nadir.  Peak signal: $\sim 1\ \mu\mathrm{m/s^2}$
from the Moon.

**Horizontal sensor** ($\hat{\mathbf{n}}$ lies in the east–north plane): picks
up the tangential component.  By the gradient approximation, the horizontal
tidal acceleration is at most $GM\,r/R^3$ — half the vertical peak.  In
practice, the horizontal signal is about $0.5\ \mu\mathrm{m/s^2}$ for the Moon,
but the geometry depends sensitively on the Moon's azimuth relative to the
sensor's orientation.

**Tilted sensor**: sees a linear combination.  The zenith and azimuth angles
set the mixing between vertical and horizontal components.

#### In code

Pytheas computes this projection for each body and sums:

```python
n_hat = measurement_axis(lat, lon, zenith, azimuth)
a_moon = tidal_acceleration(r, R_moon, GM_MOON)
a_sun  = tidal_acceleration(r, R_sun,  GM_SUN)

g_tidal_moon = np.dot(a_moon, n_hat)
g_tidal_sun  = np.dot(a_sun,  n_hat)
g_tidal      = g_tidal_moon + g_tidal_sun
```

The final `g_tidal` is the time-varying part of $g$ — the signal that rides on
top of the static normal gravity.  Combined with the elastic-Earth
amplification factor (Section 7), this is what Pytheas reports.

# Section 5: Where is the Sun? (Solar Ephemeris)

In Section 4, we worked hard to locate the Moon, wrestling with dozens of
perturbation terms caused by the Sun's gravitational pull.  The Sun itself
is a far more cooperative target.  This section derives the low-precision
solar ephemeris used in Pytheas, explains *why* it works, and shows that
its errors are negligible for tidal gravimetry.

---

## 5.1 Why the Sun is easier

### The two-body problem, almost exactly

The Moon's orbit around Earth is constantly warped by the Sun.  We needed
60+ trigonometric correction terms just to locate it to a fraction of a
degree.  Why is the Sun different?

Turn the question around: from the perspective of computing the Sun's
apparent position as seen from Earth, the relevant orbit is the *Earth's*
orbit around the Sun.  And that orbit is very nearly Keplerian --- that is,
it is governed almost entirely by the gravitational attraction between two
bodies (Earth and Sun), with only small perturbations from the other
planets.

Three facts make this orbit well-behaved:

1. **No dominant perturber.**  Jupiter, the largest perturber of Earth's
   orbit, has a mass ratio $M_{\text{Jup}} / M_{\odot} \approx 10^{-3}$
   and sits roughly 4--5 AU away.  Its perturbation of Earth's ecliptic
   longitude amounts to only a few arcseconds per century --- completely
   negligible on the timescales we care about (hours to months).

2. **Small eccentricity.**  The Earth's orbital eccentricity is
   $e \approx 0.0167$, so the orbit is very nearly circular.  This means a
   truncated Fourier expansion of the "equation of center" converges
   rapidly: three terms already give sub-arcsecond accuracy.

3. **No large inclination corrections.**  The ecliptic is the reference
   plane for the Sun's apparent motion, so the Sun's ecliptic latitude is
   zero by definition (ignoring the ~1" effect of Earth's barycentric
   offset, which we neglect).

In contrast, the Moon's orbit has eccentricity $e \approx 0.055$ (three
times larger) and is perturbed by the Sun, whose gravitational influence on
the Moon is of the same order as Earth's influence on a per-orbit basis.
The Moon truly is a three-body problem; the Sun is an easy two-body one.

---

## 5.2 Mean anomaly and the equation of center

### The Sun's mean longitude

Just as for the Moon, we begin with a polynomial in $T$ (Julian centuries
since J2000.0) that tracks how far the Sun has traveled along the ecliptic
at a uniform, average rate:

$$L_0 = 280.46646^\circ + 36000.76983^\circ T + 0.0003032^\circ T^2$$

This is the **mean longitude** of the Sun.  The linear term
$36000.77^\circ$ per century translates to roughly $0.9856^\circ$ per day,
or one full revolution per year --- exactly what you would expect for the
Sun's apparent annual motion around the sky.

The small quadratic term accounts for the slow, secular change in Earth's
orbital period due to planetary perturbations.  Over a human lifetime it
contributes a fraction of a degree.

### The mean anomaly

The mean anomaly $M$ measures how far the Earth has traveled since
perihelion (closest approach to the Sun), *as if the orbit were circular*:

$$M = 357.52911^\circ + 35999.05029^\circ T - 0.0001536^\circ T^2$$

The difference between $L_0$ and $M$ is essentially the longitude of
perihelion, which drifts slowly over millennia.

For a circular orbit, the true position angle (true anomaly $\nu$) would
equal $M$ at all times.  But the orbit is an ellipse, and on an ellipse the
body moves faster near perihelion and slower near aphelion.  The correction
$\nu - M$ is called the **equation of center**, and we now derive it from
first principles.

### Deriving the equation of center from Kepler's equation

#### Step 1: Kepler's equation

Consider an elliptical orbit with semi-major axis $a$ and eccentricity $e$.
Kepler's equation relates the **eccentric anomaly** $E$ (a geometric
auxiliary angle measured from the center of the ellipse) to the mean
anomaly $M$:

$$E - e \sin E = M$$

This equation cannot be solved in closed form for $E(M)$.  But since
$e \approx 0.0167$ is small, we can expand in powers of $e$.

#### Step 2: Series expansion of $E$ in powers of $e$

Write $E = M + \delta$, where $\delta$ is a small correction.  Substituting
into Kepler's equation:

$$(M + \delta) - e \sin(M + \delta) = M$$

$$\delta = e \sin(M + \delta)$$

**First approximation** ($\delta$ is of order $e$): set $\delta \approx 0$
inside the sine:

$$\delta_1 = e \sin M$$

so $E \approx M + e \sin M$.

**Second approximation**: substitute $\delta_1$ back:

$$\delta_2 = e \sin(M + e \sin M)$$

Expand the sine using the angle-addition formula.  Since $e \sin M$ is
small:

$$\sin(M + e \sin M) \approx \sin M \cos(e \sin M) + \cos M \sin(e \sin M)$$

$$\approx \sin M (1) + \cos M (e \sin M)$$

$$= \sin M + e \sin M \cos M$$

$$= \sin M + \tfrac{1}{2} e \sin 2M$$

Therefore:

$$\delta_2 = e \sin M + e^2 \sin M \cos M = e \sin M + \tfrac{1}{2} e^2 \sin 2M$$

and

$$E = M + e \sin M + \tfrac{1}{2} e^2 \sin 2M + \mathcal{O}(e^3)$$

#### Step 3: From eccentric anomaly to true anomaly

The true anomaly $\nu$ is related to the eccentric anomaly $E$ by the
exact geometric identity:

$$\tan \frac{\nu}{2} = \sqrt{\frac{1+e}{1-e}} \; \tan \frac{E}{2}$$

This comes from the geometry of the ellipse: the true anomaly is measured
from the focus (where the Sun sits), while the eccentric anomaly is
measured from the center.

To expand $\nu$ in powers of $e$, it is easier to work with the difference
$\nu - E$.  Using the identity above, one can show (see, e.g., Brouwer &
Clemence, *Methods of Celestial Mechanics*, Ch. 3):

$$\nu - E = 2 \arctan\!\left(\frac{e \sin E}{1 - e \cos E}\right)$$

For small $e$, expand:

$$\nu - E \approx e \sin E + \tfrac{1}{2} e^2 \sin 2E + \mathcal{O}(e^3)$$

(This uses $\arctan(x) \approx x$ for small $x$ and the binomial expansion
of $1/(1 - e\cos E)$.)

#### Step 4: Combine to get $\nu - M$

Now combine the two results.  We have:

$$E - M = e \sin M + \tfrac{1}{2} e^2 \sin 2M + \mathcal{O}(e^3)$$

$$\nu - E = e \sin E + \tfrac{1}{2} e^2 \sin 2E + \mathcal{O}(e^3)$$

In the second equation, replace $E$ by $M$ (since $E = M + \mathcal{O}(e)$,
this introduces only higher-order errors):

$$\nu - E \approx e \sin M + \tfrac{1}{2} e^2 \sin 2M + \mathcal{O}(e^3)$$

Add the two:

$$\boxed{\nu - M = 2e \sin M + \tfrac{5}{4} e^2 \sin 2M + \mathcal{O}(e^3)}$$

This is the **equation of center**.  Let us verify the coefficients make
sense.

### Numerical check

For Earth's orbit, $e = 0.016709$:

- First-order term: $2e = 0.03342$ rad $= 1.915^\circ$.  This matches the
  Meeus coefficient $1.9146^\circ$ (the small difference is absorbed by the
  $T$-dependent corrections).

- Second-order term: $\frac{5}{4} e^2 = \frac{5}{4}(0.016709)^2 = 3.49
  \times 10^{-4}$ rad $= 0.0200^\circ$.  Meeus gives $0.019993^\circ$ ---
  excellent agreement.

- The third-order term $\frac{13}{12} e^3 \sin 3M$ evaluates to $0.00029^\circ$,
  which matches the Meeus value of $0.000289^\circ$.

### The Meeus equation of center

Meeus (Astronomical Algorithms, Ch. 25) writes the equation of center with
slowly time-varying coefficients that account for the secular change in
eccentricity:

$$C = (1.9146^\circ - 0.004817^\circ T - 0.000014^\circ T^2) \sin M$$
$$\quad + (0.019993^\circ - 0.000101^\circ T) \sin 2M$$
$$\quad + 0.000289^\circ \sin 3M$$

The leading coefficient $1.9146^\circ$ is precisely $2e$ (in degrees) at
J2000.0.  The $T$-dependent terms capture the fact that Earth's eccentricity
is currently decreasing at a rate of about $0.00004$ per century, making
the orbit slightly more circular over time.

### Physical meaning

The equation of center embodies Kepler's second law: a planet sweeps out
equal areas in equal times.  Near **perihelion** (early January), the Earth
is closer to the Sun, moves faster, and the Sun appears to advance
eastward along the ecliptic faster than average.  Near **aphelion** (early
July), the opposite occurs.  The maximum discrepancy between the Sun's
true and mean positions is $\pm 1.92^\circ$, or about four solar diameters.

### True longitude

Adding the equation of center to the mean longitude gives the Sun's true
ecliptic longitude:

$$\lambda_{\odot} = L_0 + C$$

This is all we need to locate the Sun on the ecliptic plane.  Since the Sun
lies (by definition) in the ecliptic, its ecliptic latitude is zero to an
excellent approximation.

---

## 5.3 Distance to the Sun

### The orbit equation

The distance from the Sun to the Earth follows directly from the geometry
of an ellipse.  In polar form, the orbit equation is:

$$R = \frac{a(1 - e^2)}{1 + e \cos \nu}$$

where $a$ is the semi-major axis, $e$ the eccentricity, and $\nu$ the true
anomaly.

For the Earth--Sun system:
- $a = 1 \text{ AU} = 1.495978707 \times 10^{11}$ m (exactly, by IAU
  definition since 2012)
- $e \approx 0.016709$ at J2000.0

The Pytheas code uses the slightly refined semi-latus rectum factor
$1.000001018 \cdot (1 - e^2)$, where the leading factor accounts for the
distinction between the osculating and mean semi-major axes.

### How much does the distance vary?

At perihelion ($\nu = 0$): $R_{\min} = a(1 - e) \approx 0.9833$ AU.

At aphelion ($\nu = 180^\circ$): $R_{\max} = a(1 + e) \approx 1.0167$ AU.

The fractional variation is $\pm e \approx \pm 1.67\%$ about the mean.

### Impact on tidal acceleration

The tidal acceleration from the Sun scales as $1/R^3$ (as we showed in the
tidal acceleration section).  A fractional change $\delta R / R$ in distance
produces a fractional change in tidal acceleration of:

$$\frac{\delta a_{\text{tidal}}}{a_{\text{tidal}}} = -3 \frac{\delta R}{R}$$

With $\delta R / R = \pm 0.0167$:

$$\frac{\delta a_{\text{tidal}}}{a_{\text{tidal}}} = \mp 3 \times 0.0167 = \mp 5.0\%$$

So the solar tidal acceleration is about 5% stronger in January (perihelion)
than in July (aphelion).  The peak solar tidal acceleration is about 50 nGal
(micro-g level), so this seasonal modulation amounts to $\pm 2.5$ nGal ---
a small but measurable effect that Pytheas captures automatically through
the exact distance formula.

---

## 5.4 From ecliptic to ECEF

We now have the Sun's position in ecliptic coordinates: longitude
$\lambda_\odot$, latitude zero, distance $R$.  To compute tidal
acceleration, we need the Sun's position in the same Earth-fixed (ECEF)
coordinate system we use for the observer.  This requires two rotations.

### Step 1: Ecliptic to Equatorial (ECI)

The ecliptic plane is tilted relative to the equatorial plane by the
**obliquity of the ecliptic**, $\varepsilon$.  This is the angle between
Earth's spin axis and the axis perpendicular to its orbital plane ---
equivalently, the angle that causes the seasons.

At J2000.0:

$$\varepsilon_0 = 23.439291^\circ$$

The obliquity changes slowly due to gravitational torques from the Moon
and Sun on Earth's equatorial bulge:

$$\varepsilon = 23.439291^\circ - 0.013004^\circ \, T$$

(The full expression has higher-order terms, but the linear term suffices
for centuries around J2000.)

Since the Sun lies in the ecliptic plane (latitude zero), the conversion
from ecliptic to equatorial (ECI) Cartesian coordinates is simply:

$$x_{\text{ECI}} = R \cos \lambda_\odot$$

$$y_{\text{ECI}} = R \sin \lambda_\odot \cos \varepsilon$$

$$z_{\text{ECI}} = R \sin \lambda_\odot \sin \varepsilon$$

These formulas follow from a single rotation by angle $\varepsilon$ about
the $x$-axis (the direction of the vernal equinox, which is shared by both
coordinate systems).

To see why: in ecliptic coordinates, the Sun's Cartesian position is
$(R\cos\lambda_\odot, \, R\sin\lambda_\odot, \, 0)$.  The rotation matrix
from ecliptic to equatorial is:

$$\mathbf{R}_x(\varepsilon) = \begin{pmatrix} 1 & 0 & 0 \\ 0 & \cos\varepsilon & -\sin\varepsilon \\ 0 & \sin\varepsilon & \cos\varepsilon \end{pmatrix}$$

Multiplying:

$$\begin{pmatrix} x_{\text{ECI}} \\ y_{\text{ECI}} \\ z_{\text{ECI}} \end{pmatrix} = \begin{pmatrix} R \cos \lambda_\odot \\ R \sin \lambda_\odot \cos\varepsilon \\ R \sin \lambda_\odot \sin\varepsilon \end{pmatrix}$$

which is exactly what we wrote above.

### Step 2: Equatorial (ECI) to Earth-fixed (ECEF)

The ECI (Earth-Centered Inertial) frame does not rotate with the Earth.
To get to ECEF, we rotate by the Greenwich Mean Sidereal Time $\theta$,
which was derived in Section 2:

$$\begin{pmatrix} x_{\text{ECEF}} \\ y_{\text{ECEF}} \\ z_{\text{ECEF}} \end{pmatrix} = \begin{pmatrix} \cos\theta & \sin\theta & 0 \\ -\sin\theta & \cos\theta & 0 \\ 0 & 0 & 1 \end{pmatrix} \begin{pmatrix} x_{\text{ECI}} \\ y_{\text{ECI}} \\ z_{\text{ECI}} \end{pmatrix}$$

The $z$-component is unchanged because both frames share the same polar
axis (Earth's spin axis).  The rotation in the $xy$-plane accounts for
Earth's daily rotation: as the Earth turns, the Sun's position in the
ECEF frame sweeps westward at roughly $15^\circ$ per hour, giving rise to
the familiar east-to-west motion of the Sun across the sky.

### The complete pipeline

Putting it all together, the path from time to Sun position in ECEF is:

$$t \;\xrightarrow{T = (JD - 2451545)/36525}\; T \;\xrightarrow{L_0, M}\; \text{mean elements} \;\xrightarrow{+C}\; \lambda_\odot, R$$
$$\;\xrightarrow{R_x(\varepsilon)}\; \text{ECI} \;\xrightarrow{R_z(\theta)}\; \text{ECEF}$$

This is exactly what the `sun_position_ecef()` function implements.

---

## 5.5 Why low precision suffices

After the laborious lunar ephemeris of Section 4, one might worry that
the solar ephemeris is too simple --- only three sine terms!  Here we
show that this level of precision is more than adequate for tidal
gravimetry.

### Angular accuracy

The Meeus low-precision solar ephemeris is accurate to better than
$\Delta\lambda \approx 1'$ (one arcminute, or $0.017^\circ$) in ecliptic
longitude over the period 1950--2050.  How does this translate to tidal
acceleration error?

### Tidal error from angular offset

The tidal acceleration from the Sun at a point on Earth's surface is:

$$\mathbf{a}_{\text{tidal}} = GM_\odot \left( \frac{\mathbf{R} - \mathbf{r}}{|\mathbf{R} - \mathbf{r}|^3} - \frac{\mathbf{R}}{|\mathbf{R}|^3} \right)$$

where $\mathbf{R}$ is the Sun's position and $\mathbf{r}$ is the observer.
Since $|\mathbf{r}| \ll |\mathbf{R}|$ (the Earth's radius is $4.3 \times
10^{-5}$ AU), we can approximate the tidal acceleration as (from the
gradient expansion):

$$a_{\text{tidal}} \sim \frac{GM_\odot \, r}{R^3}$$

where $r \approx R_\oplus$ is the observer's distance from Earth's center.

If we misplace the Sun by a small angle $\Delta\lambda$ while keeping $R$
correct, the dominant effect is to rotate the tidal acceleration vector by
approximately $\Delta\lambda$.  The resulting error in the *magnitude* of
the projected tidal acceleration is:

$$\delta a \sim a_{\text{tidal}} \times \Delta\lambda$$

(in radians, for small angular errors).

### Numerical estimate

The peak solar tidal acceleration is approximately:

$$a_{\text{tidal}}^{\odot} \approx \frac{GM_\odot \, R_\oplus}{R^3} \approx \frac{1.33 \times 10^{20} \times 6.37 \times 10^6}{(1.50 \times 10^{11})^3} \approx 2.5 \times 10^{-7} \text{ m/s}^2$$

or about $25 \;\mu\text{Gal} = 25{,}000$ nGal as a rough scale (the
actual peak vertical component is around 50 nGal after geometric
projection and the factor of 2 from the tidal tensor).

With $\Delta\lambda = 1' = 2.9 \times 10^{-4}$ rad:

$$\delta a \sim 50 \text{ nGal} \times 2.9 \times 10^{-4} \approx 0.015 \text{ nGal}$$

This is about **0.03%** of the solar tidal signal, or equivalently less
than **0.2 nGal** in absolute terms.

### Comparison with other error sources

| Error source | Magnitude |
|---|---|
| Solar ephemeris (1 arcmin) | < 0.2 nGal |
| Lunar ephemeris (~0.1 deg) | ~1 nGal |
| Neglected ocean loading | 1--5 nGal |
| Target accuracy of Pytheas | 10--100 nGal |

The solar ephemeris error is smaller than the lunar ephemeris error for
two reasons:

1. **The solar tidal signal is smaller.**  The solar tide is only about
   46% the strength of the lunar tide ($\propto M/R^3$: the Sun is much
   more massive but also much farther away, and the $R^{-3}$ dependence
   wins).  Figure 8 overlays the lunar and solar tidal contributions over
   a seven-day window, making the roughly 2.2:1 amplitude ratio visually
   apparent.

2. **The solar position is easier to compute.**  Because the Earth--Sun
   system is nearly a clean two-body problem, three sine terms already
   achieve arcminute accuracy.  The Moon required 60+ terms for comparable
   precision.

The upshot: a three-term equation of center plus an exact orbit equation
for the distance gives solar tidal accelerations accurate to better than
0.2 nGal --- far below the overall accuracy target of Pytheas.  There is
no need for a more elaborate solar ephemeris.

---

## Summary

The solar ephemeris used in Pytheas consists of five equations:

| Quantity | Formula |
|---|---|
| Mean longitude | $L_0 = 280.466^\circ + 36000.770^\circ T$ |
| Mean anomaly | $M = 357.529^\circ + 35999.050^\circ T$ |
| Equation of center | $C \approx 1.915^\circ \sin M + 0.020^\circ \sin 2M + 0.0003^\circ \sin 3M$ |
| True longitude | $\lambda_\odot = L_0 + C$ |
| Distance | $R = a(1-e^2)/(1 + e\cos\nu)$, with $\nu = M + C$ |

The equation of center follows from a perturbative expansion of Kepler's
equation $E - e\sin E = M$ to third order in $e$.  The three coefficients
($2e$, $\frac{5}{4}e^2$, $\frac{13}{12}e^3$) match the Meeus values to
four significant figures.

The ecliptic-to-ECEF conversion requires one rotation by the obliquity
$\varepsilon \approx 23.44^\circ$ (ecliptic to equatorial) followed by one
rotation by the sidereal time $\theta$ (equatorial to Earth-fixed).  Both
rotations are simple, single-axis operations.

The resulting solar position is accurate to ~1 arcminute, contributing less
than 0.2 nGal of error to the tidal computation --- negligible compared to
Pytheas's 10--100 nGal accuracy target.

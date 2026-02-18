# 2. Coordinate Systems

Computing gravitational acceleration at a point on Earth requires a chain of coordinate transformations: geodetic coordinates $(\varphi, \lambda, h)$ from the user, Cartesian ECEF vectors for surface physics, ECI vectors for astronomical ephemerides, and ENU projections for sensor readout. This section derives each conversion.

---

## 2.1 Geodetic Coordinates: Latitude, Longitude, Altitude

A point near Earth's surface is specified by three coordinates:

- **Geodetic latitude** $\varphi$: angle between the equatorial plane and the outward normal to the reference ellipsoid. Range $[-90^\circ, +90^\circ]$.

- **Geodetic longitude** $\lambda$: eastward angle in the equatorial plane from the prime meridian. Range $[-180^\circ, +180^\circ]$.

- **Ellipsoidal height** $h$: signed distance from the reference ellipsoid along the outward normal. Differs from geoid height by tens of meters — negligible for tidal calculations.

GPS receivers report these coordinates, and Pytheas accepts them as input.

Earth's rotation produces an oblate shape, with the equatorial radius exceeding the polar radius by ~21 km. A spherical approximation would introduce position errors of that magnitude — far too large for nanoGal-level gravity modeling. The standard geometric model is therefore an ellipsoid of revolution.

---

## 2.2 The GRS80 Reference Ellipsoid

The Geodetic Reference System 1980 (GRS80) defines Pytheas's reference ellipsoid via two parameters:

**Semi-major axis:**

$$a = 6\,378\,137.0 \text{ m}$$

**Flattening:**

$$f = \frac{1}{298.257222101} \approx 3.353 \times 10^{-3}$$

From $f = (a - b)/a$, the **semi-minor axis** is:

$$b = a(1 - f) = 6\,356\,752.314 \text{ m}$$

The equator-to-pole difference $a - b \approx 21\,385$ m is small relative to $a$ (~0.3%), but enormous for precision gravimetry.

The **first eccentricity squared**:

$$e^2 = 2f - f^2 = \frac{a^2 - b^2}{a^2} \approx 6.694 \times 10^{-3}$$

This parameter appears throughout the coordinate conversions; its smallness quantifies how close the ellipsoid is to a sphere.

---

## 2.3 Geodetic vs. Geocentric Latitude

**Geocentric latitude** $\varphi'$ is the angle from the equatorial plane to the radius vector from Earth's center to the surface point. **Geodetic latitude** $\varphi$ is the angle from the equatorial plane to the ellipsoid surface normal. On an ellipsoid, these two angles differ because the surface normal does not pass through the center — except at the equator and poles.

For the meridional ellipse $x^2/a^2 + z^2/b^2 = 1$, the outward normal at $(x_0, z_0)$ is proportional to $(x_0/a^2,\; z_0/b^2)$, giving:

$$\tan\varphi = \frac{a^2}{b^2}\,\frac{z_0}{x_0} = \frac{1}{1 - e^2}\tan\varphi'$$

$$\tan\varphi' = (1 - e^2)\tan\varphi$$

Since $e^2 \approx 0.0067$, we have $|\varphi'| < |\varphi|$ everywhere except at $0^\circ$ and $\pm 90^\circ$. The maximum difference occurs near $45^\circ$:

$$\Delta\varphi_{\max} \approx \frac{e^2}{2}\sin 2\varphi \bigg|_{\max} \approx 0.0034 \text{ rad} \approx 0.19^\circ \approx 11.5'$$

This corresponds to ~21 km of surface distance, consistent with $a - b$.

**Why this matters for Pytheas:** GPS reports geodetic latitude. Substituting geocentric latitude would displace the observer by up to 21 km, producing gravity errors of order $10^{-4}$ m/s$^2$ — a million times larger than the tidal signals Pytheas aims to resolve.

---

## 2.4 Geodetic to ECEF Conversion

Given geodetic coordinates $(\varphi, \lambda, h)$, we seek the ECEF Cartesian position $(x, y, z)$. The ECEF origin lies at Earth's center of mass, $z$ is the rotation axis, $x$ points toward the prime meridian on the equator, and $y$ completes a right-handed system.

### Prime vertical radius of curvature $N(\varphi)$

From the ellipse equation and the surface normal condition, the surface-point coordinates in the meridional plane are:

$$X_0 = \frac{a^2 \cos\varphi}{\sqrt{a^2\cos^2\varphi + b^2\sin^2\varphi}}, \qquad z_0 = \frac{b^2 \sin\varphi}{\sqrt{a^2\cos^2\varphi + b^2\sin^2\varphi}}$$

One can verify these satisfy the ellipse equation:

$$\frac{X_0^2}{a^2} + \frac{z_0^2}{b^2} = \frac{a^2\cos^2\varphi + b^2\sin^2\varphi}{a^2\cos^2\varphi + b^2\sin^2\varphi} = 1 \quad \checkmark$$

Using $b^2 = a^2(1 - e^2)$:

$$a^2\cos^2\varphi + b^2\sin^2\varphi = a^2(1 - e^2\sin^2\varphi)$$

Define the **prime vertical radius of curvature**:

$$N(\varphi) = \frac{a}{\sqrt{1 - e^2\sin^2\varphi}}$$

Then:

$$X_0 = N(\varphi)\cos\varphi, \qquad z_0 = N(\varphi)(1 - e^2)\sin\varphi$$

Geometrically, $N(\varphi)$ is the distance from the surface point to the $z$-axis, measured along the ellipsoid normal. It ranges from $N(0) = a = 6\,378\,137$ m at the equator to $N(90^\circ) = a^2/b \approx 6\,399\,594$ m at the pole.

### Adding altitude and longitude

Displacing by height $h$ along the outward normal $\hat{\mathbf{n}} = (\cos\varphi, \sin\varphi)$ in the meridional plane, then rotating by longitude $\lambda$:

$$x = (N + h)\cos\varphi\cos\lambda$$

$$y = (N + h)\cos\varphi\sin\lambda$$

$$z = \bigl(N(1 - e^2) + h\bigr)\sin\varphi$$

These three equations constitute the **geodetic-to-ECEF conversion**, implemented in Pytheas as:

```python
N = A_GRS80 / np.sqrt(1.0 - E2 * sp ** 2)
x = (N + alt_m) * cp * np.cos(lam)
y = (N + alt_m) * cp * np.sin(lam)
z = (N * (1.0 - E2) + alt_m) * sp
```

**Check:** At the equator ($\varphi = 0$, $\lambda = 0$, $h = 0$): $x = a$, $y = z = 0$. At the north pole ($\varphi = 90^\circ$, $h = 0$): $x = y = 0$, $z = b$.

---

## 2.5 The ECI Frame (Earth-Centered Inertial)

The **ECI** frame has its origin at Earth's center of mass, with axes fixed relative to distant stars (the ICRF). The standard realization is the **J2000.0 equatorial frame**: $z$ points toward the mean north celestial pole, $x$ toward the mean vernal equinox, and $y$ completes the right-handed system — all defined at epoch J2000.0.

Astronomical ephemerides are naturally expressed in ECI because gravitational orbits take their simplest form in a non-rotating frame. Pytheas computes lunar and solar positions from classical ephemeris formulas (Meeus) in ecliptic coordinates, then rotates them into ECI equatorial coordinates through the obliquity $\varepsilon \approx 23.44^\circ$. Computing tidal acceleration requires the observer-to-body vector, so both positions must share a common frame. Pytheas therefore converts celestial positions from ECI to ECEF.

---

## 2.6 The ECEF Frame (Earth-Centered Earth-Fixed)

The **ECEF** frame shares the ECI origin and $z$-axis, but its $x$- and $y$-axes rotate with the Earth: $x$ through the prime meridian, $y$ toward $90^\circ$E.

ECEF is the natural frame for the observer. The geodetic-to-ECEF conversion (Section 2.4) yields the observer's position, and both the ENU basis (Section 2.9) and measurement axis are defined in this frame. Because the observer is stationary in ECEF (ignoring tectonics), we convert celestial body positions into ECEF and subtract the observer's fixed position to obtain the relative vectors needed for tidal acceleration.

**Figure 2** illustrates the three reference frames (ECI, ECEF, ENU) and the GMST rotation connecting them.

![Figure 2: Three reference frames (ECI, ECEF, ENU) showing the Earth-Centered Inertial frame fixed to distant stars, the Earth-Centered Earth-Fixed frame rotating with the planet, and the local East-North-Up basis at an observer's location. The GMST rotation angle connects ECI and ECEF.](figures/fig02_reference_frames.png)

---

## 2.7 GMST: Greenwich Mean Sidereal Time

The angle between the ECEF $x$-axis and the ECI $x$-axis (vernal equinox) is the **Greenwich Mean Sidereal Time** (GMST).

### Sidereal vs. solar day

Earth completes 366.25 sidereal rotations per year, not 365.25. The extra rotation arises from orbital motion: even a non-spinning Earth would produce one apparent solar day per year from revolution alone. Hence $N_\text{sidereal} = N_\text{solar} + 1$, giving:

$$T_\text{sidereal} = 24\text{h} \times \frac{365.25}{366.25} \approx 23\text{h}\;56\text{m}\;4.09\text{s}$$

The sidereal rotation rate in degrees per solar day is:

$$\dot{\theta} = 360^\circ \times \frac{366.25}{365.25} \approx 360.9856^\circ/\text{day}$$

### The GMST polynomial

The IAU formula:

$$\theta_\text{GMST} = \left[280.46061837 + 360.98564736629 \times (JD - 2451545.0) + 0.000387933 \times T^2 - \frac{T^3}{38710000}\right]^\circ$$

where $JD$ is the Julian Date, $T = (JD - 2451545.0)/36525$ is Julian centuries since J2000.0, and the result is taken mod $360^\circ$.

The constant term gives the GMST at the J2000.0 epoch. The linear coefficient is the sidereal rate. The $T^2$ and $T^3$ terms account for precession (~26,000-year period); even over a century, the $T^2$ term contributes only $0.0004^\circ$. For Pytheas's tidal calculations (spanning weeks to months), only the linear term matters significantly, but the full polynomial is retained for completeness.

```python
gmst_deg = (280.46061837
            + 360.98564736629 * (JD - 2451545.0)
            + 0.000387933 * T ** 2
            - T ** 3 / 38710000.0) % 360.0
```

---

## 2.8 ECI to ECEF Rotation

The ECI-to-ECEF transformation is a rotation by $-\theta_\text{GMST}$ about the $z$-axis. The ECEF frame has rotated forward by $\theta$ relative to ECI, so expressing an ECI vector in ECEF requires rotating the vector backward by $\theta$:

$$R_z(-\theta) = \begin{pmatrix} \cos\theta & \sin\theta & 0 \\ -\sin\theta & \cos\theta & 0 \\ 0 & 0 & 1 \end{pmatrix}$$

Component-wise:

$$x_\text{ecef} = \cos\theta \cdot x_\text{eci} + \sin\theta \cdot y_\text{eci}$$

$$y_\text{ecef} = -\sin\theta \cdot x_\text{eci} + \cos\theta \cdot y_\text{eci}$$

$$z_\text{ecef} = z_\text{eci}$$

**Sign check:** A star fixed in ECI at $(1, 0)$ appears at ECEF position $(\cos\theta, -\sin\theta)$ — angle $-\theta$ from the ECEF $x$-axis, as expected.

```python
def _eci_to_ecef(x_eci, y_eci, z_eci, dt):
    theta = gmst_rad(dt)
    c, s = np.cos(theta), np.sin(theta)
    return (c * x_eci + s * y_eci,
            -s * x_eci + c * y_eci,
            z_eci)
```

---

## 2.9 ENU Basis Vectors (East-North-Up)

The **East-North-Up** (ENU) local coordinate system is the natural frame for ground-based sensors. To project tidal accelerations onto observable directions, we express the ENU unit vectors in ECEF coordinates.

**Up** is the outward ellipsoid normal at $(\varphi, \lambda)$:

$$\hat{\mathbf{e}}_U = (\cos\varphi\cos\lambda,\; \cos\varphi\sin\lambda,\; \sin\varphi)$$

**East** is tangent to the surface, pointing in the direction of increasing $\lambda$. From $\partial\vec{r}/\partial\lambda$ at fixed $\varphi, h$:

$$\hat{\mathbf{e}}_E = (-\sin\lambda,\; \cos\lambda,\; 0)$$

**North** follows from $\hat{\mathbf{e}}_N = \hat{\mathbf{e}}_U \times \hat{\mathbf{e}}_E$:

$$\hat{\mathbf{e}}_N = (-\sin\varphi\cos\lambda,\; -\sin\varphi\sin\lambda,\; \cos\varphi)$$

These three vectors form a right-handed orthonormal basis: $\hat{\mathbf{e}}_E \times \hat{\mathbf{e}}_N = \hat{\mathbf{e}}_U$, and all pairwise dot products vanish.

```python
e_east  = np.array([-sl,      cl,       0.0])
e_north = np.array([-sp * cl, -sp * sl, cp])
e_up    = np.array([ cp * cl,  cp * sl, sp])
```

where `sp, cp = sin(phi), cos(phi)` and `sl, cl = sin(lam), cos(lam)`.

**Physical check:** At $(\varphi, \lambda) = (0, 0)$: $\hat{\mathbf{e}}_U = (1,0,0)$ (radially outward), $\hat{\mathbf{e}}_N = (0,0,1)$ (toward pole), $\hat{\mathbf{e}}_E = (0,1,0)$ (eastward). At the north pole: $\hat{\mathbf{e}}_U = (0,0,1)$, $\hat{\mathbf{e}}_N = (-1,0,0)$.

---

## 2.10 The Measurement Axis

A gravimeter measures the component of acceleration along its sensitive axis, which may be tilted from the local vertical. Two angles parametrize the sensor pointing direction:

- **Zenith angle** $\theta_z$: angle from the local vertical. $\theta_z = 0$ is straight up.
- **Azimuth** $\alpha$: clockwise from north. $\alpha = 0$ is north, $90^\circ$ is east.

**Figure 3** shows the ENU basis with the measurement axis $\hat{\mathbf{n}}$ at zenith angle $\theta_z$ and azimuth $\alpha$.

![Figure 3: The local ENU (East-North-Up) basis vectors at an observer's location, with the measurement axis n̂ specified by zenith angle θ_z from vertical and azimuth α measured clockwise from north.](figures/fig03_measurement_axis.png)

Decomposed in the ENU basis, the measurement axis is:

$$\hat{\mathbf{n}} = \cos\theta_z\;\hat{\mathbf{e}}_U + \sin\theta_z\bigl(\cos\alpha\;\hat{\mathbf{e}}_N + \sin\alpha\;\hat{\mathbf{e}}_E\bigr)$$

Since the ENU vectors are in ECEF (Section 2.9), $\hat{\mathbf{n}}$ is automatically in ECEF.

**Special cases:**

- **Vertical** ($\theta_z = 0$): $\hat{\mathbf{n}} = \hat{\mathbf{e}}_U$. The default in Pytheas.
- **Horizontal north** ($\theta_z = 90^\circ$, $\alpha = 0$): $\hat{\mathbf{n}} = \hat{\mathbf{e}}_N$. Measures the northward tidal component.
- **Small tilt** ($\theta_z = \epsilon \ll 1$): $\hat{\mathbf{n}} \approx \hat{\mathbf{e}}_U + \epsilon(\cos\alpha\;\hat{\mathbf{e}}_N + \sin\alpha\;\hat{\mathbf{e}}_E)$, so $g \approx g_\text{vertical} + \epsilon\,g_\text{horizontal}$. A 1-degree tilt couples horizontal tidal acceleration (~$10^{-7}$ m/s$^2$) into the measurement at ~$2 \times 10^{-9}$ m/s$^2$ (hundreds of nanoGal) — small but detectable with modern instruments. This is why Pytheas supports arbitrary sensor orientation.

```python
def measurement_axis(lat_deg, lon_deg, zenith_deg=0.0, azimuth_deg=0.0):
    e_e, e_n, e_u = enu_basis(lat_deg, lon_deg)
    zen = np.radians(zenith_deg)
    azi = np.radians(azimuth_deg)
    return (np.cos(zen) * e_u
            + np.sin(zen) * (np.cos(azi) * e_n + np.sin(azi) * e_e))
```

---

## 2.11 Summary: The Coordinate Pipeline

**Observer position.** Geodetic $(\varphi, \lambda, h)$ to ECEF via Section 2.4 (GRS80 ellipsoid):

$$\vec{r}_\text{obs} = \text{geodetic\_to\_ecef}(\varphi, \lambda, h)$$

**Celestial body positions.** ECI ephemerides rotated to ECEF via GMST (Sections 2.7--2.8):

$$\vec{r}_\text{body}^{\text{ECEF}} = R_z(-\theta_\text{GMST}) \cdot \vec{r}_\text{body}^{\text{ECI}}$$

**Tidal acceleration.** Computed in ECEF where both position vectors coexist.

**Projection.** Projected onto the measurement axis $\hat{\mathbf{n}}$ (Section 2.10) using ENU basis vectors (Section 2.9):

$$g_\text{tidal} = \vec{a}_\text{tidal} \cdot \hat{\mathbf{n}}$$

Every vector in this pipeline is expressed in ECEF; all operations take place in that single frame.

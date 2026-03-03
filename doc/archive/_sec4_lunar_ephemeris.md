# 4. Lunar Ephemeris

Computing tidal acceleration requires the Moon's geocentric position vector at every instant. This section develops the analytical lunar ephemeris used by Pytheas.

---

## 4.1 Sensitivity to Lunar Position

The tidal acceleration from Section 6 is

$$
\mathbf{a}_{\text{tidal}} = GM_{\text{Moon}} \left[ \frac{\mathbf{R} - \mathbf{r}}{|\mathbf{R} - \mathbf{r}|^3} - \frac{\mathbf{R}}{|\mathbf{R}|^3} \right],
$$

where $\mathbf{R}$ is the Moon's geocentric position vector. Because the tidal acceleration depends on direction (through the angular position of the Moon) and scales as $a_{\text{tidal}} \propto GM/R^3$, precise knowledge of the lunar distance is essential:

| Quantity | Value |
|----------|-------|
| Mean distance | 385,001 km |
| Perigee (closest) | ~356,500 km |
| Apogee (farthest) | ~406,700 km |
| Eccentricity | ~0.055 |

The apogee-to-perigee distance ratio of ~1.14 produces a tidal force ratio of $(406{,}700 / 356{,}500)^3 \approx 1.48$ -- the tidal pull at perigee is nearly 50% stronger than at apogee. At nGal-level gravimetry, even a few hundred kilometres of distance error exceeds the target accuracy. Figure 6 illustrates this anti-correlation between lunar distance and tidal acceleration amplitude over a synodic month.

![Figure 6: Lunar distance and tidal acceleration amplitude over a synodic month, showing the inverse-cube relationship between distance variations and tidal force.](figures/fig06_lunar_distance_tides.png)

All three coordinates -- longitude, latitude, and distance -- must be known to high precision at every moment.

---

## 4.2 Keplerian Orbit Review

In the unperturbed two-body problem, the orbit is an ellipse with Earth at one focus, parameterised by the semi-major axis $a$ and eccentricity $e$ ($\approx 0.055$ for the Moon). The perigee and apogee distances are $a(1-e)$ and $a(1+e)$, respectively.

### Orbital elements

Six Keplerian elements fully specify the orbit: $a$ (semi-major axis), $e$ (eccentricity), $i$ (inclination; $\approx 5.14^\circ$ to the ecliptic), $\Omega$ (longitude of the ascending node), $\omega$ (argument of periapsis), and $M$ (mean anomaly).

### Anomalies

The mean anomaly increases uniformly in time:

$$
M = M_0 + n(t - t_0),
$$

where $n = 2\pi/P$ is the mean motion. The true anomaly $\nu$, the actual angular position measured from perigee, is obtained via the eccentric anomaly $E$ by solving **Kepler's equation**:

$$
M = E - e \sin E,
$$

iteratively, and then computing

$$
\tan\frac{\nu}{2} = \sqrt{\frac{1+e}{1-e}} \tan\frac{E}{2}, \qquad r = a(1 - e\cos E).
$$

For an isolated two-body system, these relations would suffice. The Moon, however, is far from isolated.

---

## 4.3 Perturbations and the Three-Body Problem

### Solar perturbations

The Sun's gravitational influence turns the lunar problem into a three-body system with no closed-form solution (Poincare, 1890s). The principal qualitative effects of the solar perturbation are:

**Apsidal precession.** The line of apsides rotates prograde with period **~8.85 years**.

**Nodal regression.** The line of nodes regresses (retrograde) with period **~18.61 years** -- the lunar nodal cycle that modulates tidal amplitudes.

**Eccentricity oscillation.** Solar forcing drives $e$ between ~0.044 and ~0.067 over roughly half the apsidal precession period.

**Inclination oscillation.** The orbital inclination varies between ~$4.98^\circ$ and ~$5.30^\circ$, coupled to the nodal regression.

The difficulty of lunar theory famously occupied Newton -- "his head never ached but with his studies on the Moon" -- and subsequently Euler, Clairaut, Laplace, Delaunay, Hill, and Brown, driven in large part by the practical demands of maritime navigation.

---

## 4.4 The Meeus Approach

### Analytical perturbation theory

Two strategies exist for predicting the Moon's position: numerical integration of the equations of motion (JPL DE430/DE440) and analytical perturbation series. Pytheas adopts the latter approach, following **Jean Meeus** (*Astronomical Algorithms*, Chapter 47), which distils the ELP-2000 theory of Chapront-Touze and Chapront to sub-arcminute positional accuracy and ~200 km distance accuracy.

In this framework, the Moon's ecliptic longitude, latitude, and distance are each expanded as sums of sinusoidal terms whose arguments are integer linear combinations of five fundamental angles.

### The five fundamental arguments

Each is a polynomial in $T$, the number of Julian centuries since J2000.0 ($T = (JD - 2{,}451{,}545.0)/36{,}525$):

**$L'$ -- Moon's mean longitude** (tracks mean orbital motion, ~$13.18^\circ$/day):

$$
L' = 218.3164477^\circ + 481{,}267.88123421^\circ\, T - 0.0015786^\circ\, T^2 + \frac{T^3}{538{,}841}
$$

**$D$ -- Mean elongation** (Moon-Sun angular separation; $D=0$ at new moon, $180^\circ$ at full moon, ~$12.19^\circ$/day $\to$ synodic month ~29.53 d):

$$
D = 297.8501921^\circ + 445{,}267.1114034^\circ\, T - 0.0018819^\circ\, T^2 + \frac{T^3}{545{,}868}
$$

**$M$ -- Sun's mean anomaly** (~$0.986^\circ$/day $\to$ tropical year ~365.25 d):

$$
M = 357.5291092^\circ + 35{,}999.0502909^\circ\, T - 0.0001536^\circ\, T^2
$$

**$M'$ -- Moon's mean anomaly** (position on its ellipse; anomalistic month ~27.55 d, differs from sidereal month due to apsidal precession):

$$
M' = 134.9633964^\circ + 477{,}198.8675055^\circ\, T + 0.0087414^\circ\, T^2 + \frac{T^3}{69{,}699}
$$

**$F$ -- Moon's argument of latitude** (angular distance from ascending node; draconic month ~27.21 d, shorter than sidereal month due to nodal regression):

$$
F = 93.2720950^\circ + 483{,}202.0175233^\circ\, T - 0.0036539^\circ\, T^2 - \frac{T^3}{3{,}526{,}000}
$$

Because these five arguments have distinct periods, every lunar perturbation can be decomposed into sinusoids whose arguments are integer linear combinations of $(D, M, M', F)$, with the mean longitude $L'$ serving as the reference for ecliptic longitude.

### The perturbation series structure

**Ecliptic longitude:**

$$
\lambda = L' + \frac{1}{10^6} \sum_i a_i \sin(n_1^{(i)} D + n_2^{(i)} M + n_3^{(i)} M' + n_4^{(i)} F)
$$

where each coefficient $a_i$ is given in units of $10^{-6}$ degrees and the $n_k^{(i)}$ are small integers.

**Ecliptic latitude** (no mean term, since the Moon's mean orbit lies in the ecliptic plane):

$$
\beta = \frac{1}{10^6} \sum_j c_j \sin(n_1^{(j)} D + n_2^{(j)} M + n_3^{(j)} M' + n_4^{(j)} F)
$$

**Distance** (a cosine series, since distance extrema coincide with anomaly extrema):

$$
\Delta = 385{,}000.56 \text{ km} + \frac{1}{1000} \sum_k b_k \cos(n_1^{(k)} D + n_2^{(k)} M + n_3^{(k)} M' + n_4^{(k)} F)
$$

where each $b_k$ is in metres. Pytheas retains **24 terms** for longitude, **23 for distance**, and **18 for latitude**.

### The biggest terms and what they mean

#### Longitude terms (units: $10^{-6}$ degrees)

| # | $n_D$ | $n_M$ | $n_{M'}$ | $n_F$ | Coefficient $a_i$ | Amplitude | Name / Physical origin |
|---|-------|-------|----------|-------|-------------------|-----------|----------------------|
| 1 | 0 | 0 | 1 | 0 | +6,288,774 | $6.289^\circ$ | **Equation of the centre**: elliptical motion (faster at perigee, slower at apogee). |
| 2 | 2 | 0 | $-1$ | 0 | +1,274,027 | $1.274^\circ$ | **Evection**: Sun modulates the effective eccentricity depending on apsidal-line alignment (Ptolemy, ~2nd c.). |
| 3 | 2 | 0 | 0 | 0 | +658,314 | $0.658^\circ$ | **Variation**: Moon speeds up near syzygy, slows at quadrature (Tycho Brahe, ~1590). |
| 4 | 0 | 0 | 2 | 0 | +213,618 | $0.214^\circ$ | Second harmonic of the equation of the centre. |
| 5 | 0 | 1 | 0 | 0 | $-185,116$ | $0.185^\circ$ | **Annual equation**: varying Sun-Earth distance modulates the solar perturbation over the year. |

The series converges rapidly: the 24th term, with coefficient 4,036, contributes only ~$0.004^\circ$ (~14 arcseconds).

#### Distance terms (units: metres)

| # | $n_D$ | $n_M$ | $n_{M'}$ | $n_F$ | Coefficient $b_k$ | Distance | Physical origin |
|---|-------|-------|----------|-------|-------------------|----------|----------------|
| 1 | 0 | 0 | 1 | 0 | $-20,905,355$ | 20,905 km | Equation of the centre. |
| 2 | 2 | 0 | $-1$ | 0 | $-3,699,111$ | 3,699 km | Evection in distance. |
| 3 | 2 | 0 | 0 | 0 | $-2,955,968$ | 2,956 km | Variation in distance. |
| 4 | 0 | 0 | 2 | 0 | $-569,925$ | 570 km | Second eccentricity harmonic. |
| 5 | 2 | 0 | $-2$ | 0 | +246,158 | 246 km | Mixed term. |

The dominant variation (~20,900 km) is consistent with $a \cdot e \approx 385{,}001 \times 0.055 \approx 21{,}175$ km.

#### Latitude terms (units: $10^{-6}$ degrees)

| # | $n_D$ | $n_M$ | $n_{M'}$ | $n_F$ | Coefficient $c_j$ | Amplitude | Physical origin |
|---|-------|-------|----------|-------|-------------------|-----------|----------------|
| 1 | 0 | 0 | 0 | 1 | +5,128,122 | $5.128^\circ$ | **Inclination term**: orbital tilt to ecliptic. |
| 2 | 0 | 0 | 1 | 1 | +280,602 | $0.281^\circ$ | Eccentricity-inclination coupling. |
| 3 | 0 | 0 | 1 | $-1$ | +277,693 | $0.278^\circ$ | Same coupling, opposite node crossing. |
| 4 | 2 | 0 | 0 | $-1$ | +173,237 | $0.173^\circ$ | Solar perturbation of latitude. |
| 5 | 2 | 0 | $-1$ | 1 | +55,413 | $0.055^\circ$ | Evection in latitude. |

Figure 7 ranks the longitude perturbation terms by amplitude, illustrating this rapid convergence. The first few terms capture the dominant physics, while the remaining terms refine the result to sub-arcminute accuracy.

![Figure 7: Convergence of the lunar longitude perturbation series. The largest terms (equation of the centre, evection, variation) dominate, with a rapid fall-off to smaller corrections.](figures/fig07_perturbation_terms.png)

### The eccentricity correction $E$

Terms involving the Sun's mean anomaly $M$ (i.e., those with $n_2 \neq 0$) require correction for the secular decrease of Earth's orbital eccentricity:

$$
E = 1 - 0.002516\, T - 0.0000074\, T^2.
$$

Each coefficient is multiplied by $E^{|n_2|}$. This factor accounts for Earth's orbit becoming more circular (~$-0.00004$/century in $e$), which modulates the strength of the solar perturbation on all terms that depend on $M$.

In Pytheas:

```python
E = 1.0 - 0.002516 * T - 0.0000074 * T ** 2

def _sum_series(terms, use_cos=False):
    total = 0.0
    for d, ms, mp, f, coeff in terms:
        arg = d * D_r + ms * Ms_r + mp * Mp_r + f * F_r
        e_corr = E ** abs(ms)      # eccentricity correction
        if use_cos:
            total += coeff * e_corr * np.cos(arg)
        else:
            total += coeff * e_corr * np.sin(arg)
    return total
```

### Additional corrections: $A_1$, $A_2$, $A_3$

Beyond the main perturbation series, Meeus includes three additional angular arguments:

$$
A_1 = 119.75^\circ + 131.849^\circ\, T, \qquad
A_2 = 53.09^\circ + 479{,}264.290^\circ\, T, \qquad
A_3 = 313.45^\circ + 481{,}266.484^\circ\, T.
$$

These generate additive corrections for:
- **Venus perturbation of longitude**: $+3958 \sin A_1$ and $+318 \sin A_2$.
- **Earth's oblateness effect on longitude**: $+1962 \sin(L' - F)$.
- **Venus perturbation of latitude**: $+382 \sin A_3$.
- **Further latitude corrections**: $-2235 \sin L'$, $+175 \sin(A_1 \pm F)$, $+127 \sin(L' - M')$, $-115 \sin(L' + M')$.

In Pytheas:

```python
A1 = np.radians((119.75 + 131.849 * T) % 360)
A2 = np.radians((53.09 + 479264.290 * T) % 360)
A3 = np.radians((313.45 + 481266.484 * T) % 360)

sum_l += 3958 * np.sin(A1) + 1962 * np.sin(Lp_r - F_r) + 318 * np.sin(A2)
sum_b += (-2235 * np.sin(Lp_r) + 382 * np.sin(A3)
          + 175 * np.sin(A1 - F_r) + 175 * np.sin(A1 + F_r)
          + 127 * np.sin(Lp_r - Mp_r) - 115 * np.sin(Lp_r + Mp_r))
```

### Assembling the final coordinates

$$
\lambda = L' + \frac{\Sigma_l}{10^6} \quad (\text{degrees}), \qquad
\beta = \frac{\Sigma_b}{10^6} \quad (\text{degrees}), \qquad
\Delta = 385{,}000.56 + \frac{\Sigma_r}{1000} \quad (\text{km}),
$$

where $\Sigma_l$, $\Sigma_b$, $\Sigma_r$ include eccentricity corrections and the $A_1, A_2, A_3$ terms.

---

## 4.5 Ecliptic-to-ECEF Conversion

The perturbation series yields ecliptic coordinates $(\lambda, \beta, \Delta)$, but the tidal formula requires ECEF Cartesian coordinates. The conversion proceeds through two successive rotations.

### Step 1: Ecliptic to Equatorial (ECI)

The ecliptic-equatorial tilt (obliquity) is

$$
\varepsilon = 23.439291^\circ - 0.013004^\circ\, T.
$$

First, express the Moon's position in Cartesian ecliptic coordinates:

$$
\begin{pmatrix} X_{\text{ecl}} \\ Y_{\text{ecl}} \\ Z_{\text{ecl}} \end{pmatrix}
= \Delta \begin{pmatrix} \cos\beta\,\cos\lambda \\ \cos\beta\,\sin\lambda \\ \sin\beta \end{pmatrix},
$$

Then apply the ecliptic-to-equatorial rotation about the $x$-axis (the vernal equinox direction, shared by both frames):

$$
\begin{pmatrix} X_{\text{ECI}} \\ Y_{\text{ECI}} \\ Z_{\text{ECI}} \end{pmatrix}
= R_x(\varepsilon) \begin{pmatrix} X_{\text{ecl}} \\ Y_{\text{ecl}} \\ Z_{\text{ecl}} \end{pmatrix}, \qquad
R_x(\varepsilon) = \begin{pmatrix} 1 & 0 & 0 \\ 0 & \cos\varepsilon & -\sin\varepsilon \\ 0 & \sin\varepsilon & \cos\varepsilon \end{pmatrix}.
$$

In components:

$$
X_{\text{ECI}} = \Delta\,\cos\beta\,\cos\lambda
$$

$$
Y_{\text{ECI}} = \Delta\left(\cos\beta\,\sin\lambda\,\cos\varepsilon - \sin\beta\,\sin\varepsilon\right)
$$

$$
Z_{\text{ECI}} = \Delta\left(\cos\beta\,\sin\lambda\,\sin\varepsilon + \sin\beta\,\cos\varepsilon\right)
$$

**Sanity check:** For $\beta = 0$, $\lambda = 90^\circ$: $X_{\text{ECI}} = 0$, $Y_{\text{ECI}} = \Delta\cos\varepsilon$, $Z_{\text{ECI}} = \Delta\sin\varepsilon$ -- the Moon appears above the equatorial plane by the obliquity angle, as expected.

In Pytheas:

```python
eps = np.radians(OBLIQUITY_J2000 - 0.013004 * T)
cb, sb = np.cos(beta), np.sin(beta)
cl, sl = np.cos(lam),  np.sin(lam)
ce, se = np.cos(eps),  np.sin(eps)

x_eci = dist_m * cb * cl
y_eci = dist_m * (cb * sl * ce - sb * se)
z_eci = dist_m * (cb * sl * se + sb * ce)
```

### Step 2: ECI to ECEF

The ECEF frame co-rotates with the Earth. The angle between the ECI and ECEF $x$-axes is the **Greenwich Mean Sidereal Time** (GMST) $\theta$, derived in Section 2. The transformation is a rotation about the $z$-axis:

$$
\begin{pmatrix} X_{\text{ECEF}} \\ Y_{\text{ECEF}} \\ Z_{\text{ECEF}} \end{pmatrix}
= R_z(\theta) \begin{pmatrix} X_{\text{ECI}} \\ Y_{\text{ECI}} \\ Z_{\text{ECI}} \end{pmatrix}
= \begin{pmatrix} \cos\theta & \sin\theta & 0 \\ -\sin\theta & \cos\theta & 0 \\ 0 & 0 & 1 \end{pmatrix}
\begin{pmatrix} X_{\text{ECI}} \\ Y_{\text{ECI}} \\ Z_{\text{ECI}} \end{pmatrix}.
$$

Note that the $z$-component is unchanged, since the rotation is about the polar axis.

In Pytheas:

```python
def _eci_to_ecef(x_eci, y_eci, z_eci, dt):
    theta = gmst_rad(dt)
    c, s = np.cos(theta), np.sin(theta)
    return (c * x_eci + s * y_eci,
            -s * x_eci + c * y_eci,
            z_eci)
```

### The complete pipeline

1. **Compute $T$** from the Julian Date.
2. **Evaluate the five fundamental arguments** $L'$, $D$, $M$, $M'$, $F$.
3. **Compute the eccentricity correction** $E$.
4. **Sum the perturbation series** for longitude (24 sine terms), distance (23 cosine terms), and latitude (18 sine terms), applying $E^{|n_2|}$.
5. **Apply the $A_1$, $A_2$, $A_3$ corrections.**
6. **Assemble ecliptic coordinates**: $\lambda$, $\beta$, $\Delta$.
7. **Rotate ecliptic $\to$ ECI** via obliquity $\varepsilon$.
8. **Rotate ECI $\to$ ECEF** via GMST $\theta$.

The result is $\mathbf{R}_{\text{Moon}}$ in ECEF coordinates, ready for use in the tidal acceleration formula of Section 6.

---

### Summary

The Moon's position cannot be computed from a simple Keplerian orbit, because solar perturbations continuously distort the lunar trajectory. Pytheas instead uses an analytical perturbation series (Meeus's compilation of ELP-2000): a total of 24 + 23 + 18 = 65 sinusoidal terms in five fundamental angles capture ecliptic longitude, distance, and latitude to sub-arcminute accuracy. Two rotations -- ecliptic to equatorial, then equatorial to ECEF -- place the result in Earth-fixed coordinates for the tidal calculation.

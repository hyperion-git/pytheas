# 8. General Relativity and the Lab-Frame Gravimeter Equation

Section 6 derived the tidal acceleration from Newtonian gravity and
obtained the exact formula

$$
\mathbf{a}_{\text{tidal}} = GM\left[\frac{\mathbf{R} - \mathbf{r}}{|\mathbf{R} - \mathbf{r}|^3} - \frac{\mathbf{R}}{R^3}\right],
$$

which predicts tidal signals at the $\mu$Gal level in excellent agreement with
observation.  This section places that formula on rigorous general-relativistic
footing by working entirely in the **proper reference frame of the lab** — a
non-geodesic, rotating observer on Earth's surface.  Fermi normal
coordinates make every physical effect transparent: proper
acceleration, Coriolis and centrifugal forces, tidal curvature, gravitomagnetic
cross-terms, and post-Newtonian corrections each appear as a distinct term in a
single geodesic equation.  We compute the magnitude of each term, classify it by
observational relevance, and show that the Newtonian tidal formula suffices
for surface gravimetry at the nGal level.

Every step of the derivation has been verified computationally using SymPy; the verification scripts are documented in Section 8.9.


## 8.1 The Fermi Normal Coordinate Framework

### 8.1.1 The Observer

A gravimeter on Earth's surface is a **non-geodesic** observer: the normal
force from the ground supports it against gravitational free fall.
It also **co-rotates** with Earth at angular velocity $\boldsymbol{\Omega}$.
We construct Fermi normal coordinates $(T, X^i)$ along the lab worldline
$\gamma$ following Ni & Zimmerman (1978), Li & Ni (1979), and Poisson & Will
(2014, Ch. 9).  Here $T$ is proper time along $\gamma$, and $X^i$
are spatial coordinates defined by geodesics orthogonal to $\gamma$, with an
orthonormal tetrad that is Fermi--Walker transported and then rotated to
co-rotate with the lab.

Three quantities characterize the observer:
- $a^i$: proper acceleration of the lab worldline (upward, supporting against gravity)
- $\Omega^i$: angular velocity of the spatial triad relative to Fermi--Walker transport (Earth's rotation)
- $R_{\alpha\beta\gamma\delta}$: Riemann curvature tensor evaluated on the worldline, in the Fermi tetrad basis

### 8.1.2 The Metric

Expanding the metric to quadratic order in $X^i$ about the worldline gives:

$$
g_{00} = -\!\left(1 + \frac{a_i X^i}{c^2}\right)^{\!2} + \frac{R_{0i0j}\,X^i X^j}{c^2} + \mathcal{O}(X^3)
$$

$$
g_{0i} = \frac{2}{3}\,R_{0jik}\,X^j X^k + \epsilon_{ijk}\,\Omega^j X^k + \mathcal{O}(X^3)
$$

$$
g_{ij} = \delta_{ij} - \frac{1}{3}\,R_{ikjl}\,X^k X^l + \mathcal{O}(X^3)
$$

We adopt signature $(-,+,+,+)$ with $x^0 = cT$.  The Riemann components carry
dimension $[\text{length}^{-2}]$ and are evaluated on the worldline.

**Physical content of each component:**
- $g_{00}$: the $a_i X^i$ term encodes the **gravitational redshift** (equivalence principle); the $R_{0i0j}$ term encodes **tidal curvature** from external fields.
- $g_{0i}$: the $\epsilon_{ijk}\Omega^j X^k$ term encodes **rotation** (Coriolis/centrifugal); the $R_{0jik}$ term encodes **gravitomagnetic curvature** (frame-dragging, velocity-dependent tidal effects).
- $g_{ij}$: the $R_{ikjl}$ term encodes **spatial curvature**, negligible for our purposes.

Expanding the square in $g_{00}$ produces a term
$(\Omega^2 \delta_{ij} - \Omega_i\Omega_j)X^i X^j / c^2$ from the coupling
of the $a_i X^i$ piece to the rotation-dependent part of the proper
acceleration.  This term gives rise to the centrifugal force, as verified below.


## 8.2 The Geodesic Equation in the Lab Frame

A gravimeter measures the coordinate acceleration of a freely falling test mass
in the lab frame.  For a non-relativistic test mass with $u^0 \approx c$ and
$v^i = dX^i/dT \ll c$, the spatial geodesic equation reduces to:

$$
\frac{d^2 X^i}{dT^2} = -c^2\,\Gamma^i_{00} - 2c\,\Gamma^i_{0j}\,v^j - \Gamma^i_{jk}\,v^j v^k + \cdots
$$

Computing the Christoffel symbols from the Fermi metric and collecting all
contributions yields (verified symbolically in Block 1 of the SymPy scripts):

$$
\boxed{\frac{d^2 X^i}{dT^2} = \underbrace{-a^i}_{\text{proper accel.}} \underbrace{- 2(\boldsymbol{\Omega} \times \mathbf{v})^i}_{\text{Coriolis}} \underbrace{- [\boldsymbol{\Omega} \times (\boldsymbol{\Omega} \times \mathbf{X})]^i}_{\text{centrifugal}} + \underbrace{c^2\,R^i{}_{0j0}\,X^j}_{\text{tidal}} + \underbrace{2\,R^i{}_{0jk}\,X^k\,v^j}_{\text{GM cross}}}
$$

The following subsections derive and classify each term.


## 8.3 Term-by-Term Derivation and Classification

### Term 1: Proper Acceleration $-a^i$ (Static Gravity)

The leading contribution from $\Gamma^i_{00}$ is $-a^i$, where $a^i$ denotes
the lab's four-acceleration.  By the equivalence principle, $a^i = -g^i$
(directed upward), so the test mass accelerates at $+g^i$ (downward):

$$
\frac{d^2 X^i}{dT^2}\bigg|_{\text{proper}} = g^i(\phi, h) \approx 9.81\;\text{m/s}^2 \approx 9.81 \times 10^8\;\mu\text{Gal}
$$

**Classification:** Time-independent.  This is the static normal gravity
$\gamma(\phi, h)$ derived in Section 3.

### Term 2: Centrifugal Acceleration

The coupling between the $g_{0i}$ rotation terms and the four-velocity
normalization produces the centrifugal acceleration.  Expanding $g_{00}$
shows that the $|\boldsymbol{\Omega} \times \mathbf{X}|^2 / c^2$ piece
contributes to $\Gamma^i_{00}$, giving (via the Levi-Civita identity
$\epsilon^i{}_{jk}\epsilon^k{}_{lm} = \delta^i_l\delta_{jm} - \delta^i_m\delta_{jl}$):

$$
\frac{d^2 X^i}{dT^2}\bigg|_{\text{centrifugal}} = -[\boldsymbol{\Omega} \times (\boldsymbol{\Omega} \times \mathbf{X})]^i = \Omega^2 X^i_\perp
$$

At the equator: $\omega^2 R_\oplus \approx 0.0339\;\text{m/s}^2 \approx 3.39 \times 10^6\;\mu\text{Gal}$.

**Classification:** Time-independent, absorbed into $\gamma(\phi, h)$.

**SymPy verification (Block 1):** The centrifugal term was verified symbolically
by computing $\Gamma^i_{00}$ from the full Fermi metric including the
$(\Omega^2\delta_{ij} - \Omega_i\Omega_j)X^iX^j/c^2$ contribution to $g_{00}$.
All three spatial components match $-[\boldsymbol{\Omega}\times(\boldsymbol{\Omega}\times\mathbf{X})]^i$.

### Term 3: Coriolis Force and the Eötvös Effect

The rotation term in $g_{0i}$ produces the Coriolis acceleration via
$\Gamma^i_{0j}|_\Omega = \epsilon^i{}_{jk}\,\Omega^k/c$:

$$
\frac{d^2 X^i}{dT^2}\bigg|_{\text{Coriolis}} = -2(\boldsymbol{\Omega} \times \mathbf{v})^i
$$

**Vertical projection.**  In the local ENU frame with
$\boldsymbol{\Omega} = (0,\;\omega\cos\phi,\;\omega\sin\phi)$, a test mass
falling vertically with $\mathbf{v} = (0, 0, -v_z)$ yields:

$$
-2(\boldsymbol{\Omega} \times \mathbf{v}) = 2\omega v_z \cos\phi\;\hat{E}
$$

The first-order Coriolis acceleration is **purely eastward**: it has zero
projection onto the vertical measurement axis $\hat{\mathbf{n}} = \hat{U}$.

**The Eötvös (second-order Coriolis) effect.**  The eastward velocity acquired
from Coriolis deflection feeds back through the Coriolis force to produce a
vertical second-order component, derived by perturbation theory:

*Step 1:* The eastward acceleration $a_E = 2\omega g t \cos\phi$ integrates to
$v_E(t) = \omega g t^2 \cos\phi$.

*Step 2:* The Coriolis force from $v_E\hat{E}$ has a vertical component:
$-2(\boldsymbol{\Omega} \times v_E\hat{E})_U = 2\omega\cos\phi \cdot v_E$.
Substituting $v_E$:

$$
a_U^{(2)} = 2\omega^2 g\, t^2 \cos^2\!\phi
$$

The latitude dependence is $\cos^2\!\phi$, **not** $\sin\phi\cos\phi$.  Both
factors of $\cos\phi$ trace back to the same $\Omega_N = \omega\cos\phi$
component: one governs the initial eastward deflection, the other governs
the secondary vertical deflection.  The $\sin\phi\cos\phi$ factor sometimes
seen in the literature refers to the *northward* (unmeasured) deflection
component.  At $\varphi = 45^\circ$ the two forms coincide numerically because
$\cos^2 45^\circ = \sin 45^\circ \cos 45^\circ = 1/2$.

For a drop height $h = 0.3\;\text{m}$, $t_\text{fall} = \sqrt{2h/g} \approx 0.247\;\text{s}$:

$$
\delta a_z = 2 \times (7.29 \times 10^{-5})^2 \times 9.81 \times (0.247)^2 \times \tfrac{1}{2} \approx 3.19 \times 10^{-9}\;\text{m/s}^2 \approx 319\;\text{nGal}
$$

**SymPy verification (Block 3):**  Full ODE integration using `scipy.integrate.solve_ivp`
(DOP853, rtol = $10^{-14}$) confirms the analytic formula at all five test
latitudes ($0^\circ, 30^\circ, 45^\circ, 60^\circ, 90^\circ$) with sub-nGal residuals.

**Classification:** Instrument-dependent.  Present in absolute (free-fall)
gravimeters but **absent** from superconducting gravimeters, whose test mass is
levitated rather than falling.  Routinely corrected for in absolute gravimetry.

### Term 4: Tidal Acceleration from Riemann Curvature

The $R_{0i0j}\,X^i X^j$ term in $g_{00}$ contributes to $\Gamma^i_{00}$:

$$
\frac{d^2 X^i}{dT^2}\bigg|_{\text{tidal}} = c^2\,R^i{}_{0j0}\,X^j
$$

This is the **geodesic deviation equation** -- tidal acceleration driven by
spacetime curvature.  Decomposing the Riemann tensor into Earth's self-field
(static, absorbed into $\gamma$) and time-varying external fields (Moon, Sun):

$$
R^i{}_{0j0} = R^i{}_{0j0}\big|_{\text{Earth}} + R^i{}_{0j0}\big|_{\text{ext}}
$$

Section 8.4 recovers the Newtonian tidal formula from the external piece.

### Term 5: Gravitomagnetic Curvature Cross-Terms

The $R_{0jik}$ Riemann components in $g_{0i}$ couple the test mass velocity
to the gravitomagnetic tidal field:

$$
\frac{d^2 X^i}{dT^2}\bigg|_{\text{GM}} = 2\,R^i{}_{0jk}\,X^k\,v^j
$$

The factor of 2 has two origins: the $(2/3)$ coefficient in the metric becomes
$1$ in the Christoffel symbol (via Riemann symmetries), and the overall factor
of $2$ comes from the $-2c\,\Gamma^i_{0j}\,v^j$ term in the geodesic equation.

In the weak-field limit, these components are suppressed relative to the
electric-type tidal components by $v_\text{rot}/c$, where
$v_\text{rot} = \omega R_\oplus \cos\phi$ is the lab's rotational velocity.
At $\varphi = 45^\circ$:

$$
a_\text{GM} \sim 2\,\frac{v_\text{rot}}{c}\,a_\text{tidal}^\text{Newt} \sim 2 \times \frac{329}{3 \times 10^8} \times 8 \times 10^{-7} \approx 0.17\;\text{nGal}
$$

**Classification:** Below the 10 nGal accuracy floor.  Negligible.


## 8.4 Recovering the Newtonian Tidal Formula

In the weak-field, slow-motion limit ($|\Phi|/c^2 \ll 1$, $v/c \ll 1$), the
linearized metric takes the form $g_{00} = -(1 + 2\Phi/c^2)$,
$g_{ij} = (1 - 2\Phi/c^2)\delta_{ij}$, and the Riemann tensor reduces to:

$$
c^2\,R^i{}_{0j0} = -\frac{\partial^2 \Phi_\text{ext}}{\partial X^i \partial X^j} \equiv -\mathcal{T}_{ij}
$$

where $\mathcal{T}_{ij}$ is the **Newtonian tidal tensor**.  For a point mass
$M$ at geocentric position $\mathbf{r}_M$, with the lab at $\mathbf{R}$,
the potential $\Phi = -GM/|\mathbf{r}_M - \mathbf{x}|$ gives:

$$
\mathcal{T}_{ij} = \frac{\partial^2 \Phi}{\partial x_i\,\partial x_j}\bigg|_{\mathbf{x}=\mathbf{R}} = -\frac{GM}{d^3}\left(3\hat{\mathbf{n}}_i\hat{\mathbf{n}}_j - \delta_{ij}\right)
$$

where $d = |\mathbf{r}_M - \mathbf{R}|$ and $\hat{\mathbf{n}} = (\mathbf{r}_M - \mathbf{R})/d$.
The tensor is trace-free ($\mathcal{T}_{ii} = 0$) as required by the vacuum
Laplace equation $\nabla^2\Phi = 0$.

The tidal acceleration at the lab, $a^i_\text{tidal} = -\mathcal{T}_{ij}\,R^j$,
reproduces to leading order in $R/r_M$ the exact formula from Section 6:

$$
\boxed{\mathbf{a}_{\text{tidal}} = GM\left[\frac{\mathbf{r}_M - \mathbf{R}}{|\mathbf{r}_M - \mathbf{R}|^3} - \frac{\mathbf{r}_M}{|\mathbf{r}_M|^3}\right]}
$$

**SymPy verification (Block 2):**
- All six independent components of $\mathcal{T}_{ij}$ verified symbolically
  by computing $\partial_i\partial_j\Phi$ in SymPy and simplifying against
  the $-GM/d^3(3\hat{\mathbf{n}}_i\hat{\mathbf{n}}_j - \delta_{ij})$ formula.
- Trace-free property verified: $\mathcal{T}_{ii} = 0$.
- Exact tidal formula verified symbolically by computing
  $-\nabla\Phi|_\mathbf{R} - (-\nabla\Phi|_\mathbf{0})$ and simplifying.
- Numerical cross-check with Moon parameters: tensor approximation vs. exact
  formula agree to $O(R_\oplus/d_\text{Moon}) \approx 1.7\%$.

**Magnitudes:**

| Body | $GM/d^3$ (s$^{-2}$) | Max tidal at surface ($\mu$Gal) | Moon/Sun ratio |
|:-----|:---------------------|:-------------------------------|:---------------|
| Moon | $8.63 \times 10^{-14}$ | $\sim 110$ | 2.18 |
| Sun  | $3.97 \times 10^{-14}$ | $\sim 51$ | 1.00 |

These are the $\delta \cdot [\mathbf{a}_\text{Moon} + \mathbf{a}_\text{Sun}] \cdot \hat{\mathbf{n}}$ terms in $g(t)$.


## 8.5 Post-Newtonian and Gravitomagnetic Corrections

### 8.5.1 The $g_{0i}$ Hierarchy: Rotation vs. Lense-Thirring

The $g_{0i}$ component contains two physically distinct contributions:

$$
g_{0i} = \underbrace{\epsilon_{ijk}\Omega^j X^k}_{\text{coordinate rotation}} + \underbrace{\tfrac{2}{3}\,R_{0jik}\,X^j X^k}_{\text{curvature (frame-dragging)}}
$$

Earth's rotation contributes $\Omega = \omega \approx 7.29 \times 10^{-5}\;\text{s}^{-1}$.
The Lense--Thirring gravitomagnetic field at Earth's surface is:

$$
B_g = \frac{2GJ_\oplus}{c^2 R_\oplus^3} \approx 3.4 \times 10^{-14}\;\text{s}^{-1}
$$

where $J_\oplus = 0.3307\,M_\oplus R_\oplus^2\,\omega$.  Their ratio:

$$
\frac{\omega}{B_g} \approx 2.2 \times 10^9
$$

The Fermi frame makes this nine-order-of-magnitude hierarchy manifest: both
enter as $g_{0i}$ cross-terms, but the rotation term is a coordinate effect
(removable by passing to a non-rotating frame), whereas the Lense--Thirring
term reflects genuine curvature encoded in $R_{0jik}$.

The direct Lense--Thirring acceleration on a test mass falling at
$v_\text{fall} \approx g\,t_\text{fall} \approx 2.4\;\text{m/s}$ is:

$$
a_\text{LT} \sim B_g \times v_\text{fall} \approx 3.4 \times 10^{-14} \times 2.4 \approx 0.008\;\text{nGal}
$$

### 8.5.2 First Post-Newtonian Correction to the Tidal Tensor

At first post-Newtonian (1PN) order in GR ($\gamma = \beta = 1$ in PPN), the
tidal tensor acquires corrections:

$$
\mathcal{T}^{\text{1PN}}_{ij} = \partial_i \partial_j \Phi\!\left(1 + \frac{4\Phi}{c^2}\right) + \frac{4}{c^2}\,(\partial_i \Phi)(\partial_j \Phi) - \frac{\delta_{ij}}{c^2}\,|\nabla\Phi|^2
$$

The coefficient 4 arises from the PPN combination $(2 + 2\gamma) = 4$ for
$\gamma = 1$.  The fractional correction relative to the Newtonian result is
$O(\Phi/c^2)$:

| Body | $\Phi/c^2$ | 1PN tidal correction |
|:-----|:-----------|:---------------------|
| Sun  | $9.9 \times 10^{-9}$ | $\sim 0.001\;\text{nGal}$ |
| Moon | $1.4 \times 10^{-13}$ | $\sim 3 \times 10^{-8}\;\text{nGal}$ |

### 8.5.3 Gravitomagnetic Tidal Effects: Complete Inventory

| Effect | Source | Magnitude (nGal) |
|:-------|:-------|:-----------------|
| Direct LT on test mass | Earth $J_\oplus$ | $\sim 0.008$ |
| Spin GM tidal | Sun $J_\odot$ | $\sim 2 \times 10^{-13}$ |
| Spin GM tidal | Moon $J_\text{Moon}$ | $\sim 5 \times 10^{-15}$ |
| Orbital GM tidal | Sun (via $v_\oplus$) | $\sim 10^{-9}$ |
| Orbital GM tidal | Moon (via $v_\text{Moon}$) | $\sim 10^{-10}$ |

All fall below the 10 nGal accuracy floor by orders of magnitude.


## 8.6 The Complete Hierarchy

Figure 14 displays the hierarchy of corrections graphically.  The table below
lists every contribution to the measured gravitational acceleration, ordered by
magnitude.

![Figure 14: Complete hierarchy of corrections to the gravimeter signal from the Fermi frame analysis, spanning from static normal gravity ($\sim 10^{12}$ nGal) down to spin gravitomagnetic tidal ($\sim 10^{-15}$ nGal).  Effects are classified by observational status: absorbed into $\gamma(\varphi,h)$ (gray), measured tidal signal (navy), horizontal with no $\hat{\mathbf{n}}$ projection (green), instrument-dependent (blue), and below the 10 nGal accuracy floor (orange).](figures/fig14_gr_corrections.png)

| Rank | Effect | Magnitude | Classification |
|:-----|:-------|:----------|:---------------|
| 1 | Static normal gravity $\gamma(\phi,h)$ | $9.81 \times 10^8\;\mu$Gal | (a) Time-independent |
| 2 | Centrifugal | $3.4 \times 10^6\;\mu$Gal | (a) Time-independent |
| 3 | Moon Newtonian tidal | $\sim 110\;\mu$Gal | **Measured tidal signal** |
| 4 | Earth gravity gradient (0.3 m) | $\sim 92\;\mu$Gal | (a) Instrument correction |
| 5 | Sun Newtonian tidal | $\sim 51\;\mu$Gal | **Measured tidal signal** |
| 6 | Coriolis (1st order) | Horizontal | (b) No projection onto $\hat{\mathbf{n}}$ |
| 7 | Eötvös (2nd-order Coriolis) | $\sim 319\;$nGal | (c) Instrument-dependent |
| 8 | GM curvature cross-term | $\sim 0.17\;$nGal | (d) Below 10 nGal floor |
| 9 | Direct Lense--Thirring | $\sim 0.008\;$nGal | (d) Below 10 nGal floor |
| 10 | 1PN tidal (Sun) | $\sim 0.001\;$nGal | (d) Below 10 nGal floor |
| 11 | 1PN tidal (Moon) | $\sim 3 \times 10^{-8}\;$nGal | (d) Below 10 nGal floor |
| 12 | Orbital GM tidal (Sun) | $\sim 10^{-9}\;$nGal | (d) Below 10 nGal floor |
| 13 | Orbital GM tidal (Moon) | $\sim 10^{-10}\;$nGal | (d) Below 10 nGal floor |
| 14 | Spin GM tidal (Sun) | $\sim 2 \times 10^{-13}\;$nGal | (d) Below 10 nGal floor |
| 15 | Spin GM tidal (Moon) | $\sim 5 \times 10^{-15}\;$nGal | (d) Below 10 nGal floor |

Every correction falls into exactly one of four categories:

**(a) Time-independent, absorbed into $\gamma(\phi,h)$:**
Static gravity, centrifugal force, Earth's gravity gradient.

**(b) Horizontal, zero projection onto $\hat{\mathbf{n}}$:**
First-order Coriolis on a vertically falling mass.

**(c) Instrument-dependent, corrected for:**
Eötvös effect ($\sim 319\;$nGal for 0.3 m drop at $\varphi = 45^\circ$).
Present in absolute gravimeters only; absent from superconducting gravimeters.

**(d) Below the 10 nGal accuracy floor:**
All general-relativistic corrections: gravitomagnetic cross-terms,
Lense--Thirring, 1PN tidal, and spin/orbital GM tidal effects.


## 8.7 Dimensional and Limiting-Case Checks

**Dimensional consistency** (all verified):
- Tidal tensor: $[GM/d^3] = [\text{s}^{-2}]$; multiplied by $R_\oplus\;[\text{m}]$ gives $[\text{m/s}^2]$.
- Eötvös: $[\omega^2 g t^2] = [\text{s}^{-2}][\text{m/s}^2][\text{s}^2] = [\text{m/s}^2]$.
- 1PN: $[\Phi/c^2]$ is dimensionless, multiplying $[\text{m/s}^2]$.
- Lense--Thirring: $[B_g \cdot v] = [\text{s}^{-1}][\text{m/s}] = [\text{m/s}^2]$.
- Geodesic deviation: $[c^2 R X] = [\text{m}^2/\text{s}^2][\text{m}^{-2}][\text{m}] = [\text{m/s}^2]$.

**Limiting cases:**
- $\omega \to 0$: Coriolis, centrifugal, and Eötvös terms vanish; only proper acceleration and tidal curvature survive.
- $M_\text{Moon}, M_\odot \to 0$: All tidal terms vanish; the gravimeter reads a constant $\gamma(\phi,h)$.
- $R_{\alpha\beta\gamma\delta} \to 0$: The geodesic equation reduces to that of a particle in a uniformly accelerating, rotating frame.
- $c \to \infty$: All 1PN, gravitomagnetic, and Lense--Thirring terms vanish; Newtonian rotating-frame mechanics survives.


## 8.8 Conclusion

The Fermi normal coordinate analysis provides a rigorous justification for
every term in the gravimeter equation:

$$
\boxed{g(t) = \gamma(\phi, h) + \delta \cdot [\mathbf{a}_{\text{Moon}}(t) + \mathbf{a}_{\text{Sun}}(t)] \cdot \hat{\mathbf{n}}}
$$

Starting from the full GR metric in the lab's proper frame, the geodesic
equation yields the Newtonian tidal formula as the leading time-varying
contribution.  Every correction is either (a) time-independent and absorbed
into $\gamma(\phi,h)$, (b) horizontal with zero projection onto $\hat{\mathbf{n}}$,
(c) instrument-dependent and routinely corrected for, or (d) below the 10 nGal
accuracy floor.  The largest GR correction -- the gravitomagnetic curvature
cross-term at $\sim 0.17$ nGal -- falls 60 times below Pytheas's accuracy
floor.

The Newtonian tidal formula is therefore not merely a useful approximation: it
is the *exact leading-order result* of general relativity in the weak-field
limit, with controlled, computable, and negligible corrections at every
post-Newtonian order.


## 8.9 Computational Verification

Four independent SymPy/NumPy scripts verify the derivation, each targeting a
self-contained block of the argument:

| Block | Script | Verifies |
|:------|:-------|:---------|
| 1 | `block1_fermi_metric_christoffel.py` | Fermi metric $\to$ Christoffel symbols $\to$ all 5 geodesic terms |
| 2 | `block2_riemann_to_newtonian_tidal.py` | Tidal tensor $\mathcal{T}_{ij}$, trace-free property, exact tidal formula |
| 3 | `block3_coriolis_eotvos.py` | Coriolis projection, Eötvös $\cos^2\!\phi$, ODE integration at 5 latitudes |
| 4 | `block4_numerical_bounds.py` | Complete hierarchy table (16 effects, all numerical values) |

All scripts require only SymPy, NumPy, and SciPy.

**Key computational findings:**
- Block 1: All five geodesic equation terms (proper acceleration, Coriolis $\times 3$ components, centrifugal $\times 3$ components) verified symbolically.
- Block 2: All six independent tidal tensor components match the analytic formula to SymPy `simplify()` $= 0$.  Exact tidal formula verified for all three spatial components.
- Block 3: Eötvös latitude dependence confirmed as $\cos^2\!\phi$ (vertical) and $\sin\phi\cos\phi$ (northward, unmeasured).  Full ODE integration matches the perturbation-theory prediction at all latitudes.
- Block 4: All 16 effects in the hierarchy table computed independently from physical constants; the $\omega/B_g \approx 2.17 \times 10^9$ ratio confirmed.

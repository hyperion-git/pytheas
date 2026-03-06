# Theory

From general relativity to the lab-frame equation of motion.

This document derives the equations behind Pytheas, starting from the metric in the lab's proper reference frame and ending with the three measurable quantities --- the gravity vector $\mathbf{g}$, the gradient tensor $\mathbf{T}$, and the rotation vector $\boldsymbol{\Omega}$ --- that determine the acceleration of a test mass in a terrestrial laboratory.  For what the code computes and its accuracy, see [Implementation](implementation.md).

**Conventions.**  Signature $(-,+,+,+)$, $x^0 = cT$.  Spatial indices $i,j,k$ run over 1,2,3.  ENU coordinates: East, North, Up.  $g_U < 0$ at the surface.

**References.**
Synge (1960), *Relativity: The General Theory*;
Manasse & Misner (1963), J. Math. Phys. **4**, 735;
Ni & Zimmerman (1978), Phys. Rev. D **17**, 1473;
Poisson & Will (2014), *Gravity: Newtonian, Post-Newtonian, Relativistic*, Ch. 9.

---


## 1. What Does a Gravimeter Measure?

A gravimeter on Earth's surface is a **non-geodesic**, **rotating** observer: the ground supports it against free fall, and the Earth spins.  The instrument measures the coordinate acceleration of a freely falling test mass in the laboratory frame.

Three quantities characterize this observer:

- $a^i$: **proper acceleration** of the lab worldline (directed upward, supporting against gravity),
- $\Omega^i$: **angular velocity** of the spatial triad relative to Fermi--Walker transport (Earth's rotation),
- $R_{\alpha\beta\gamma\delta}$: **Riemann curvature tensor** evaluated on the worldline (encodes tidal fields).

The task is to write down the acceleration of a test mass in terms of these three quantities, accurate to the level where the Newtonian tidal formula suffices (sub-$\mu$Gal, typically $10^2$--$10^3$ nGal in this implementation).


## 2. Fermi Normal Coordinates

Fermi normal coordinates (FNC) $(T, X^i)$ are constructed along the lab worldline $\gamma$: $T$ is proper time, $X^i$ are spatial displacements defined by geodesics orthogonal to $\gamma$, with an orthonormal tetrad that is Fermi--Walker transported and then rotated to co-rotate with the lab (Manasse & Misner 1963; Ni & Zimmerman 1978).

### The metric to second order

Expanding the metric to quadratic order in $X^i$ about the worldline:

$$g_{00} = -\!\left(1 + \frac{a_i X^i}{c^2}\right)^{\!2} + \frac{R_{0i0j}\,X^i X^j}{c^2} + \mathcal{O}(X^3) \tag{2.1}$$

$$g_{0i} = \frac{2}{3}\,R_{0jik}\,X^j X^k + \frac{1}{c}\,\epsilon_{ijk}\,\Omega^j X^k + \mathcal{O}(X^3) \tag{2.2}$$

$$g_{ij} = \delta_{ij} - \frac{1}{3}\,R_{ikjl}\,X^k X^l + \mathcal{O}(X^3) \tag{2.3}$$

The curvature coefficients ($2/3$, $1/3$) are fixed by the Riemann symmetries together with the Fermi-coordinate gauge conditions on the worldline.  The linear terms encode the observer's prescribed acceleration and rotation; because $x^0 = cT$, the rotation term in $g_{0i}$ carries an explicit factor $1/c$.

**Physical content:**

| Component | Terms | Physics |
|-----------|-------|---------|
| $g_{00}$ | $a_i X^i$ | Gravitational redshift (equivalence principle) |
| $g_{00}$ | $R_{0i0j} X^i X^j$ | Tidal curvature from external fields |
| $g_{0i}$ | $\epsilon_{ijk}\Omega^j X^k / c$ | Rotation (Coriolis, centrifugal) |
| $g_{0i}$ | $R_{0jik} X^j X^k$ | Gravitomagnetic curvature (frame-dragging) |
| $g_{ij}$ | $R_{ikjl} X^k X^l$ | Spatial curvature (negligible here) |

### The geodesic equation

For a non-relativistic test mass ($v^i \ll c$), the spatial geodesic equation reduces to:

$$\frac{d^2 X^i}{dT^2} = -c^2\,\Gamma^i_{00} - 2c\,\Gamma^i_{0j}\,v^j + \cdots \tag{2.4}$$

Computing the Christoffel symbols from Eqs. (2.1)--(2.3) and collecting terms:

$$\boxed{\frac{d^2 X^i}{dT^2} = \underbrace{-a^i}_{\text{gravity}} \underbrace{-\,2(\boldsymbol{\Omega}\times\mathbf{v})^i}_{\text{Coriolis}} \underbrace{-\,[\boldsymbol{\Omega}\times(\boldsymbol{\Omega}\times\mathbf{X})]^i}_{\text{centrifugal}} + \underbrace{c^2 R^i{}_{0j0}\,X^j}_{\text{tidal}} + \underbrace{2c\,R^i{}_{0jk}\,X^k v^j}_{\text{gravitomagnetic}}} \tag{2.5}$$

Each term is derived and classified in Section 5.


## 3. The Bitensor Framework

The covariant Taylor expansion generalizes the ordinary Taylor series to curved spacetime.  Three objects provide the machinery that underlies the FNC metric.

### Synge world function

The *world function* $\sigma(x, x')$ assigns to every pair of points connected by a unique geodesic:

$$\sigma(x, x') = \tfrac{1}{2}\,\varepsilon\,(\text{geodesic length})^2 \tag{3.1}$$

where $\varepsilon = -1$ for timelike and $+1$ for spacelike geodesics.  The covariant gradient $\sigma_a \equiv \nabla_a\sigma$ is the tangent vector to the geodesic at $x$; $-\sigma^a$ is the curved-space displacement from $x$ toward $x'$.

At coincidence ($x' \to x$): $[\sigma] = 0$, $[\sigma_{a\,b'}] = -g_{ab}$.

### Parallel propagator

The parallel propagator $g^a_{\ b'}(x,x')$ transports a vector from $x'$ to $x$ along the geodesic:

$$V^a(x) = g^a_{\ b'}(x,x')\,V^{b'}(x') \tag{3.2}$$

It satisfies $\sigma^c\nabla_c\,g^a_{\ b'} = 0$ with $[g^a_{\ b'}] = \delta^a_{\ b}$.  In terms of the world function: $g^a_{\ b'} = -g^{ac}\sigma_{cb'}$.  For rank-2 tensors, one applies one propagator per index.

### Covariant Taylor expansion

For a scalar field $f(x')$ at a nearby point:

$$f(x') = f - (\nabla_a f)\,\sigma^a + \tfrac{1}{2}(\nabla_a\nabla_b f)\,\sigma^a\sigma^b - \cdots \tag{3.3}$$

where everything on the right is evaluated at $x$.  For a vector field, the parallel propagator handles index transport.

Applying Eq. (3.3) to the metric $g_{ab}$ about the lab worldline, using the orthonormal Fermi--Walker-transported tetrad as basis, recovers the FNC metric of Eqs. (2.1)--(2.3).  The bitensor framework is thus the theoretical engine behind the FNC expansion.


## 4. The Rotating Lab Frame

### From FNC to the ADM Hamiltonian

The FNC metric decomposes into ADM form (lapse $N$, shift $N_i$, spatial metric $\gamma_{ij}$).  For a non-relativistic test particle ($|v| \ll c$, lab heights $h \ll c^2/g$), the Hamiltonian reduces to:

$$H_\text{FNC} = \frac{p^2}{2m} + m\,\Phi_\text{grav}(\mathbf{x}) \tag{4.1}$$

where the **gravitational potential in FNC** is:

$$\Phi_\text{grav}(\mathbf{x}) = a_j\,x^j - \frac{c^2}{2}\,R_{0i0j}\,x^i x^j \tag{4.2}$$

The first term is zeroth-order gravity ($g_i = -a_i$ at the origin); the second is the tidal field from spacetime curvature.

### Canonical transformation to the rotating frame

The lab co-rotates with the Earth.  The canonical transformation generated by the angular momentum $\mathbf{L} = \mathbf{x}\times\mathbf{p}$ gives:

$$H_\text{rot} = \frac{p^2}{2m} + m\,\Phi_\text{grav}(\mathbf{x}) - \boldsymbol{\Omega}\cdot(\mathbf{x}\times\mathbf{p}) \tag{4.3}$$

No centrifugal potential appears explicitly --- it is encoded in $\mathbf{p} \neq m\mathbf{v}$ in the rotating frame.

### The Lagrangian

Hamilton's equations give $\mathbf{v} = \mathbf{p}/m - \boldsymbol{\Omega}\times\mathbf{x}$, hence $\mathbf{p} = m(\mathbf{v} + \boldsymbol{\Omega}\times\mathbf{x})$.  The Legendre transform yields the rotating-frame Lagrangian (per unit mass):

$$L = \tfrac{1}{2}|\mathbf{v}|^2 + (\boldsymbol{\Omega}\times\mathbf{x})\cdot\mathbf{v} - \Phi_\text{eff}(\mathbf{x}) \tag{4.4}$$

where the effective potential is:

$$\Phi_\text{eff} = \Phi_\text{grav} - \tfrac{1}{2}|\boldsymbol{\Omega}\times\mathbf{x}|^2 \tag{4.5}$$

The three terms in Eq. (4.4) are: kinetic energy, Coriolis coupling (analogous to a magnetic vector potential $\mathbf{A} = \boldsymbol{\Omega}\times\mathbf{x}$), and the effective potential (gravity + centrifugal).

### Taylor expansion: defining g and T

Expanding $\Phi_\text{eff}$ about the lab origin to second order in $\delta\mathbf{x} = \mathbf{x} - \mathbf{x}_0$:

**The gravity vector** (including centrifugal):

$$g_i \equiv -\frac{\partial\Phi_\text{eff}}{\partial x^i}\bigg|_{\mathbf{x}_0} = -a_i \tag{4.6}$$

At the FNC origin, the gravity vector equals the negative of the proper acceleration.  In ENU coordinates: $g_i = (0, 0, -\gamma)$ where $\gamma$ is the Somigliana normal gravity (which already includes the centrifugal correction).

**The gravity gradient tensor**:

$$T_{ij} \equiv -\frac{\partial^2\Phi_\text{eff}}{\partial x^i\,\partial x^j}\bigg|_{\mathbf{x}_0} = c^2\,R_{0i0j} + \Omega^2\delta_{ij} - \Omega_i\Omega_j \tag{4.7}$$

This symmetric $3\times 3$ tensor encodes tidal stretching (from curvature $c^2 R_{0i0j}$) and the centrifugal effect ($\Omega^2\delta_{ij} - \Omega_i\Omega_j$).

**The rotation vector** $\Omega_i$ enters through the Coriolis coupling and cannot be written as the gradient of any potential.  In ENU coordinates:

$$\boldsymbol{\Omega} = (0,\;\Omega\cos\varphi,\;\Omega\sin\varphi) \tag{4.8}$$

with $\Omega = 7.292\,115 \times 10^{-5}$ rad/s.


## 5. The Equation of Motion

### Derivation

Starting from Hamilton's equations for $H_\text{rot}$ (Eq. 4.3), differentiating $\mathbf{p} = m(\mathbf{v} + \boldsymbol{\Omega}\times\mathbf{x})$ and expanding:

$$m\ddot{\mathbf{x}} = -m\nabla\Phi_\text{grav} - m\boldsymbol{\Omega}\times(\boldsymbol{\Omega}\times\mathbf{x}) - 2m\boldsymbol{\Omega}\times\mathbf{v} \tag{5.1}$$

The first two terms combine into $-m\nabla\Phi_\text{eff}$ via the identity:

$$\tfrac{1}{2}\nabla|\boldsymbol{\Omega}\times\mathbf{x}|^2 = -\boldsymbol{\Omega}\times(\boldsymbol{\Omega}\times\mathbf{x}) \tag{5.2}$$

Taylor-expanding $-\nabla_i\Phi_\text{eff} \approx g_i + T_{ij}\,\delta x^j$ yields:

$$\boxed{\ddot{\mathbf{x}} = \mathbf{g} + \mathbf{T}\cdot\delta\mathbf{x} - 2\,\boldsymbol{\Omega}\times\mathbf{v}} \tag{5.3}$$

This is the central result.  Each term:

- $\mathbf{g}$: gravity at the lab origin (centrifugal absorbed via Somigliana),
- $\mathbf{T}\cdot\delta\mathbf{x}$: gravity gradient (tidal + centrifugal),
- $-2\boldsymbol{\Omega}\times\mathbf{v}$: Coriolis (the only explicitly velocity-dependent term).

**Dimensional check.** $[\mathbf{g}] = LT^{-2}$; $[T_{ij}\delta x^j] = T^{-2}\cdot L = LT^{-2}$; $[\Omega\times\mathbf{v}] = T^{-1}\cdot LT^{-1} = LT^{-2}$.

### Centrifugal absorption

The centrifugal acceleration does **not** appear as a separate term.  It is fully absorbed:

- **Into $\mathbf{g}$:** The Somigliana formula (WGS84 normal gravity) computes the magnitude of $-\nabla\Phi_\text{eff}$, not $-\nabla\Phi_\text{grav}$ alone.  The centrifugal correction is already included.
- **Into $\mathbf{T}$:** The centrifugal contribution to the gradient tensor is $T^\text{cent}_{ij} = \Omega^2\delta_{ij} - \Omega_i\Omega_j$, which enters the Poisson trace condition (Section 5.3).

### The Poisson trace condition

The gravitational potential satisfies $\nabla^2\Phi_\text{grav} = 4\pi G\rho$, and the centrifugal potential contributes:

$$\nabla^2\!\left(-\tfrac{1}{2}|\boldsymbol{\Omega}\times\mathbf{x}|^2\right) = -2\Omega^2 \tag{5.4}$$

Since $T_{ii} = -\nabla^2\Phi_\text{eff}$:

$$\boxed{\operatorname{Tr}(\mathbf{T}) = -4\pi G\rho + 2\Omega^2} \tag{5.5}$$

In vacuum ($\rho = 0$) outside the Earth: $\operatorname{Tr}(\mathbf{T}) = 2\Omega^2 \approx 1.06 \times 10^{-8}$ s$^{-2}$ (~10.6 Eotvos).  This is the mathematically consistent acceleration-gradient convention used in Eqs. (4.7) and (5.3); the current Pytheas static Earth approximation is discussed separately in [Implementation](implementation.md).

### Tidal tracelessness

The tidal gradient tensor from an external body (Moon, Sun) is symmetric and traceless:

$$T^\text{tidal}_{ij} = T^\text{tidal}_{ji}, \qquad \operatorname{Tr}(\mathbf{T}_\text{tidal}) = 0 \tag{5.6}$$

Symmetry follows from $\partial_i\partial_j = \partial_j\partial_i$.  Tracelessness follows from the Laplace equation: the tidal potential is harmonic in the vacuum exterior, so $\nabla^2\Phi_\text{tidal} = 0$.

For a body at distance $d$ along the radial direction, the leading-order gradient tensor:

$$T^\text{tidal} = \frac{GM}{d^3}\begin{pmatrix} -1 & 0 & 0 \\ 0 & -1 & 0 \\ 0 & 0 & 2\end{pmatrix} \tag{5.7}$$

manifestly traceless ($-1-1+2=0$), reflecting the quadrupolar stretching-compression pattern.

### Limiting cases

| Limit | Behaviour |
|-------|-----------|
| $\Omega \to 0$ | $\ddot{\mathbf{x}} = \mathbf{g} + \mathbf{T}\cdot\delta\mathbf{x}$ (non-rotating geodesic deviation); $\operatorname{Tr}(\mathbf{T}) = -4\pi G\rho$ |
| $M_\text{ext} \to 0$ | All tidal terms vanish; gravimeter reads constant $\gamma(\varphi,h)$ |
| Poles ($\varphi = 90\degree$) | $\boldsymbol{\Omega} = (0,0,\Omega)$ purely vertical; centrifugal vanishes at origin |
| Flat spacetime | $R_{\alpha\beta\gamma\delta} = 0$, no tidal forces; rotating-frame mechanics only |


## 6. The Newtonian Limit

### Riemann curvature as the tidal tensor

In the weak-field, slow-motion limit ($|\Phi|/c^2 \ll 1$, $v/c \ll 1$), the linearized metric gives:

$$c^2\,R^i{}_{0j0} = -\frac{\partial^2\Phi_\text{ext}}{\partial X^i\partial X^j} \equiv -\mathcal{T}_{ij} \tag{6.1}$$

where $\mathcal{T}_{ij}$ is the Newtonian tidal tensor.  For a point mass $M$ at geocentric position $\mathbf{r}_M$, with the lab at $\mathbf{R}$ and $\mathbf{d} = \mathbf{r}_M - \mathbf{R}$:

$$\mathcal{T}_{ij} = -\frac{GM}{d^3}\left(3\hat{\mathbf{d}}_i\hat{\mathbf{d}}_j - \delta_{ij}\right) \tag{6.2}$$

This tensor is trace-free ($\mathcal{T}_{ii} = 0$) as required by the vacuum Laplace equation.
With this convention, $\mathcal{T}_{ij}$ is the Hessian-level object from Eq. (6.1); the acceleration-gradient tensor used in Pytheas is $T^\text{tidal}_{ij} = -\mathcal{T}_{ij} = \frac{GM}{d^3}\!\left(-\delta_{ij} + 3\hat{\mathbf{d}}_i\hat{\mathbf{d}}_j\right)$ (Eq. 5.7).

### The exact tidal formula

The tidal acceleration at the lab recovers the exact Newtonian expression:

$$\boxed{\mathbf{a}_\text{tidal} = GM\left[\frac{\mathbf{r}_M - \mathbf{R}}{|\mathbf{r}_M - \mathbf{R}|^3} - \frac{\mathbf{r}_M}{|\mathbf{r}_M|^3}\right]} \tag{6.3}$$

This is the difference between gravitational pull at the observer and at the Earth's center, without expanding in $R/r_M$.  The gradient (tidal tensor) approximation truncates at order $(R/r_M)^1$; for the Moon with $R/r_M \approx 1/60$, the next-order term contributes ~2.5%, or ~2.5 $\mu$Gal (~2500 nGal).  **The exact formula is not an approximation --- it is the leading-order result of GR in the weak-field limit.**

### Why "exact" is the right word

Eq. (6.3) is often called the "exact Newtonian" formula to distinguish it from the gradient approximation.  More precisely, it is the *complete* Newtonian tidal acceleration.  The only corrections to it are:

1. Post-Newtonian ($\sim \Phi/c^2 \sim 10^{-8}$) --- at most 0.001 nGal (Section 8),
2. Higher multipoles of the Moon's mass distribution --- negligible.

The formula is exact to all orders in $R/r_M$ within Newtonian gravity.

### Magnitudes

| Body | $GM/d^3$ (s$^{-2}$) | Max tidal at surface ($\mu$Gal) | Moon/Sun ratio |
|------|---------------------|-------------------------------|----------------|
| Moon | $8.63 \times 10^{-14}$ | ~110 | 2.18 |
| Sun  | $3.97 \times 10^{-14}$ | ~51 | 1.00 |


## 7. Elastic Earth

The solid Earth deforms under tidal stress, amplifying the surface gravity change.

### Love numbers

The displacement Love number $h_2$ and potential Love number $k_2$ characterize the degree-2 elastic response to a tidal potential $V_2$:

- **Radial displacement:** $\Delta r = h_2\,V_2/g$
- **Potential perturbation:** $V_\text{deform} = k_2\,V_2$

### The gravimetric factor

The measured tidal gravity change on an elastic Earth combines three effects:

| Contribution | Expression | Effect on $|\Delta g|$ |
|-------------|------------|----------------------|
| Rigid-Earth tidal gravity | $-2V_2/R$ | Baseline |
| Free-air (displacement) | $-2h_2 V_2/R$ | Amplifies (+60.8%) |
| Mass redistribution | $+3k_2 V_2/R$ | Reduces ($-44.7$%) |

Combining:

$$\Delta g_\text{elastic} = \Delta g_\text{rigid} \times \left(1 + h_2 - \tfrac{3}{2}k_2\right) \tag{7.1}$$

The gravimetric factor is:

$$\boxed{\delta = 1 + h_2 - \tfrac{3}{2}k_2} \tag{7.2}$$

With IERS 2010 values ($h_2 = 0.6078$, $k_2 = 0.2980$):

$$\delta = 1.1608 \tag{7.3}$$

The tidal signal is 16% larger than on a rigid Earth.

### Known limitations

**FCN resonance.**  Near the K1 tidal frequency (~23.93 h), the free core nutation resonance causes the Love numbers to vary rapidly.  $\delta$ can deviate from 1.1608 by ~1%, introducing a site/frequency-dependent error typically in the $10^2$--$10^3$ nGal range.  This is one of the dominant model-error terms in Pytheas.  A frequency-dependent $\delta(f)$ table would correct it.

**Higher degrees.**  The degree-3 lunar contribution is ~1.7% of degree 2; in this stripped-down model it is smaller than the dominant error terms above.


## 8. Post-Newtonian Corrections

### The complete hierarchy

Eq. (2.5) contains five terms.  The following table classifies every contribution to the measured acceleration:

| Rank | Effect | Magnitude | Classification |
|------|--------|-----------|----------------|
| 1 | Static normal gravity $\gamma(\varphi,h)$ | $9.81 \times 10^8$ $\mu$Gal | Time-independent |
| 2 | Centrifugal | $3.4 \times 10^6$ $\mu$Gal | Time-independent, absorbed into $\gamma$ |
| 3 | Lunar tidal (Newtonian) | ~110 $\mu$Gal | **Measured signal** |
| 4 | Solar tidal (Newtonian) | ~51 $\mu$Gal | **Measured signal** |
| 5 | First-order Coriolis | Horizontal | No projection onto $\hat{\mathbf{n}}$ |
| 6 | Second-order Coriolis (Eotvos) | ~319 nGal | Instrument-dependent$^*$ |
| 7 | GM curvature cross-term | ~0.17 nGal | Below practical model floor |
| 8 | Direct Lense--Thirring | ~0.008 nGal | Below practical model floor |
| 9 | 1PN tidal (Sun) | ~0.001 nGal | Below practical model floor |
| 10 | 1PN tidal (Moon) | $\sim 3 \times 10^{-8}$ nGal | Below practical model floor |

$^*$Present in absolute (free-fall) gravimeters; absent from superconducting gravimeters.

### Gravitomagnetic effects

The $R^i{}_{0jk}$ Riemann components in $g_{0i}$ couple the test mass velocity to the gravitomagnetic tidal field:

$$a^i_\text{GM} = 2c\,R^i{}_{0jk}\,X^k v^j \tag{8.1}$$

These are suppressed relative to the electric-type tidal components by $v_\text{rot}/c$, where $v_\text{rot} = \Omega R_\oplus\cos\varphi$ is the lab's rotational velocity.  At $\varphi = 45\degree$:

$$a_\text{GM} \sim 2\frac{v_\text{rot}}{c}\,a^\text{Newt}_\text{tidal} \sim 2 \times \frac{329}{3\times10^8}\times 8\times10^{-7} \approx 0.17\;\text{nGal}  \tag{8.2}$$

The Lense--Thirring gravitomagnetic field from Earth's angular momentum:

$$B_g = \frac{2GJ_\oplus}{c^2 R_\oplus^3} \approx 3.4\times10^{-14}\;\text{s}^{-1} \tag{8.3}$$

Nine orders of magnitude below Earth's rotation rate.  Direct acceleration on a falling test mass: ~0.008 nGal.

### 1PN tidal corrections

At first post-Newtonian order, the tidal tensor acquires corrections of order $\Phi/c^2$:

| Body | $\Phi/c^2$ | 1PN tidal correction |
|------|-----------|---------------------|
| Sun  | $9.9\times10^{-9}$ | ~0.001 nGal |
| Moon | $1.4\times10^{-13}$ | $\sim 3\times10^{-8}$ nGal |

### Conclusion

Every correction falls into one of four categories:

**(a) Time-independent, absorbed into $\gamma(\varphi,h)$:** Static gravity, centrifugal.

**(b) Horizontal, zero projection onto $\hat{\mathbf{n}}$:** First-order Coriolis on a vertically falling mass.

**(c) Instrument-dependent:** Eotvos effect (~319 nGal for 0.3 m drop); present in absolute gravimeters, absent from superconducting gravimeters.

**(d) Below the practical model floor (typically $10^2$--$10^3$ nGal):** All general-relativistic corrections.  The largest (GM curvature cross-term at ~0.17 nGal) remains far below that floor.

The gravimeter equation:

$$\boxed{g(t) = \gamma(\varphi,h) + \delta\left[\mathbf{a}_\text{Moon}(t) + \mathbf{a}_\text{Sun}(t)\right]\cdot\hat{\mathbf{n}}} \tag{8.4}$$

is not merely a useful approximation.  It is the *exact leading-order result* of general relativity in the weak-field limit, with controlled, computable, and negligible corrections at every post-Newtonian order.

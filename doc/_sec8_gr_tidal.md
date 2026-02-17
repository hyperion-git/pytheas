## 8.1 Why spacetime curvature?

In Section 6 we computed the tidal acceleration from straightforward Newtonian reasoning and arrived at a formula that works beautifully: take the difference between the gravitational pull at the gravimeter and the gravitational pull at the Earth's center, and the result predicts tidal signals at the $\mu$Gal level that match observation. So why invoke general relativity at all?

There are three reasons, and they go beyond mere elegance. First, GR provides the *reason* tidal forces exist: they are the observable signature of spacetime curvature. Uniform gravitational acceleration can always be removed by passing to a freely falling frame -- this is Einstein's equivalence principle, the thought experiment of the elevator in free fall. But tidal gravity, the differential acceleration across a finite baseline, *cannot* be removed. It is intrinsic geometry, as irreducible as the curvature of the Earth's surface. Second, the post-Newtonian expansion gives a controlled perturbative framework that rigorously *bounds* the error incurred by using Newton's formula. Third, and most practically, it identifies the precise conditions under which the Newtonian treatment would break down: either $v/c$ or $GM/(Rc^2)$ must grow large enough for the corrections to matter. For solar system gravimetry, as we will show, they never do.

The physical content of this section is the **geodesic deviation equation**, which makes the connection between tidal forces and curvature mathematically exact.

### 8.1.1 Geodesics and free fall

In general relativity a freely falling test particle traces a geodesic of spacetime. Let $x^\alpha(\tau)$ denote the worldline parameterized by proper time $\tau$, with four-velocity $u^\alpha = dx^\alpha / d\tau$. The geodesic equation is

$$\frac{d^2 x^\alpha}{d\tau^2} + \Gamma^\alpha_{\ \beta\gamma}\, u^\beta\, u^\gamma = 0 \,,$$

where $\Gamma^\alpha_{\ \beta\gamma}$ are the Christoffel symbols of the metric $g_{\mu\nu}$:

$$\Gamma^\alpha_{\ \beta\gamma} = \frac{1}{2}\, g^{\alpha\delta}\!\left(\partial_\beta\, g_{\gamma\delta} + \partial_\gamma\, g_{\beta\delta} - \partial_\delta\, g_{\beta\gamma}\right).$$

In the language of covariant differentiation, the geodesic equation states simply that the covariant acceleration vanishes:

$$\frac{Du^\alpha}{D\tau} \equiv u^\beta \nabla_\beta\, u^\alpha = 0 \,.$$

A single geodesic tells us about the motion of one particle. To describe tidal effects -- the *relative* acceleration between two nearby particles -- we need to study how neighboring geodesics deviate from one another.

### 8.1.2 Geodesic deviation

Consider a one-parameter family of geodesics $x^\alpha(\tau, s)$, where $s$ labels neighboring worldlines and $\tau$ is proper time along each one. Two natural vector fields emerge from this construction: the tangent vector (four-velocity) $u^\alpha = \partial x^\alpha / \partial \tau$ and the deviation vector (separation) $\xi^\alpha = \partial x^\alpha / \partial s$. Because both arise as partial derivatives of the same coordinate map, they commute:

$$u^\beta \nabla_\beta \xi^\alpha = \xi^\beta \nabla_\beta u^\alpha \,.$$

A useful geometric analogy lives on the two-sphere: great circles leaving the equator parallel to one another converge as they approach the poles. The rate of that convergence is governed by the Gaussian curvature of the sphere. In spacetime the story is the same, with the Riemann curvature tensor playing the role of the Gaussian curvature.

The relative acceleration of neighboring geodesics is the second covariant derivative of $\xi^\alpha$ along the flow:

$$A^\alpha \equiv \frac{D^2 \xi^\alpha}{D\tau^2} = u^\gamma \nabla_\gamma \!\left(u^\beta \nabla_\beta \xi^\alpha\right).$$

The derivation proceeds by three key moves. First, use the commutation relation to replace $u^\beta \nabla_\beta \xi^\alpha$ with $\xi^\beta \nabla_\beta u^\alpha$ inside the outer derivative. Second, expand by the Leibniz rule and apply the commutation relation a second time. Third -- and this is where curvature enters -- interchange the order of covariant derivatives using the Ricci identity:

$$\nabla_\gamma \nabla_\beta u^\alpha - \nabla_\beta \nabla_\gamma u^\alpha = R^\alpha_{\ \delta\gamma\beta}\, u^\delta \,,$$

where $R^\alpha_{\ \delta\gamma\beta}$ is the Riemann curvature tensor. The geodesic condition $u^\gamma \nabla_\gamma u^\alpha = 0$ kills one group of terms, and the remaining terms involving products of $\nabla u$ cancel upon relabeling dummy indices. What survives is pure curvature:

$$\frac{D^2 \xi^\alpha}{D\tau^2} = -R^\alpha_{\ \beta\gamma\delta}\, u^\beta\, \xi^\gamma\, u^\delta$$

This is the **geodesic deviation equation** (also called the Jacobi equation). It is exact -- no approximations have been made. The passage from the intermediate form $R^\alpha_{\ \delta\gamma\beta}\, u^\gamma u^\delta\, \xi^\beta$ to the standard index placement uses only the algebraic symmetries of the Riemann tensor (antisymmetry in the last pair of indices and the symmetry of $u^\gamma u^\delta$).

### 8.1.3 Physical interpretation

The physical content of the geodesic deviation equation is striking in its directness: the relative acceleration between two freely falling observers is entirely determined by the Riemann curvature tensor contracted with their four-velocity and their separation. No other quantity enters.

The message is clear: tidal gravity *is* spacetime curvature. Two test masses released side by side in Earth orbit will drift apart or together at a rate set by $R^\alpha_{\ \beta\gamma\delta}$, and no choice of coordinates or reference frame can make this effect vanish. This is the same principle that underlies gravitational-wave detection: LIGO measures the differential displacement of its test masses precisely because a passing gravitational wave modulates the Riemann tensor, producing a time-varying tidal strain across the interferometer baseline.

In Section 6 we computed the tidal acceleration from Newtonian gravity. Now we will see that formula emerge from the geometry of spacetime itself -- first through the electrogravitic tidal tensor (Section 8.2), then through the weak-field limit that recovers Newton (Section 8.3).
## 8.2 The electrogravitic tidal tensor

The geodesic deviation equation from the previous section encodes all tidal physics in a single contraction of the Riemann tensor with the observer's four-velocity and the separation vector. It is natural to package that contraction into a tensor of its own -- one that plays the role of the relativistic tidal field.

### 8.2.1 Definition and compact form of geodesic deviation

Given an observer with four-velocity $u^\alpha$, define the **electrogravitic tidal tensor** with both indices down as

$$\mathcal{E}_{\alpha\gamma} = R_{\alpha\beta\gamma\delta}\, u^\beta\, u^\delta \,.$$

To connect this with the geodesic deviation equation, which involves the mixed Riemann tensor $R^\alpha_{\ \beta\gamma\delta}$, raise the first index:

$$\mathcal{E}^\alpha_{\ \gamma} = g^{\alpha\mu}\,\mathcal{E}_{\mu\gamma} = g^{\alpha\mu}\, R_{\mu\beta\gamma\delta}\, u^\beta\, u^\delta = R^\alpha_{\ \beta\gamma\delta}\, u^\beta\, u^\delta \,.$$

With this identification, the geodesic deviation equation takes a compact and physically transparent form:

$$\frac{D^2 \xi^\alpha}{D\tau^2} = -\mathcal{E}^\alpha_{\ \gamma}\, \xi^\gamma \,.$$

This is an exact restatement -- no approximation has been made. The entire tidal physics is now captured by a single object $\mathcal{E}^\alpha_{\ \gamma}$ acting linearly on the separation vector $\xi^\gamma$. In the observer's rest frame, the spatial components $\mathcal{E}^i_{\ j}$ form a $3 \times 3$ symmetric matrix that stretches and squeezes spatial separations -- exactly the tidal pattern familiar from the Newtonian analysis of Section 6.

### 8.2.2 Symmetries and spatial character

Two key properties of $\mathcal{E}_{\alpha\gamma}$ follow directly from the algebraic symmetries of the Riemann tensor:

**Symmetry.** The pair symmetry of the Riemann tensor, $R_{\alpha\beta\gamma\delta} = R_{\gamma\delta\alpha\beta}$, implies

$$\mathcal{E}_{\alpha\gamma} = R_{\alpha\beta\gamma\delta}\, u^\beta\, u^\delta = R_{\gamma\delta\alpha\beta}\, u^\delta\, u^\beta = \mathcal{E}_{\gamma\alpha} \,.$$

The electrogravitic tensor is symmetric, just as the Newtonian tidal tensor $T_{ij} = \partial_i \partial_j \Phi$ is symmetric by equality of mixed partial derivatives.

**Spatial character.** The antisymmetry of the Riemann tensor in its last index pair, $R_{\alpha\beta\gamma\delta} = -R_{\alpha\beta\delta\gamma}$, gives

$$\mathcal{E}_{\alpha\gamma}\, u^\gamma = R_{\alpha\beta\gamma\delta}\, u^\beta u^\delta u^\gamma = 0 \,,$$

since $u^\delta u^\gamma$ is symmetric under $\delta \leftrightarrow \gamma$ while $R_{\alpha\beta\gamma\delta}$ is antisymmetric. The electrogravitic tensor is therefore purely spatial in the observer's frame: it has no components along the four-velocity direction.

These two properties -- symmetry and spatial orthogonality to $u^\alpha$ -- mean that $\mathcal{E}^i_{\ j}$ has at most six independent components in the observer's rest frame. A quick dimensional check confirms the consistency: since $R_{\alpha\beta\gamma\delta}$ has dimensions of $1/\text{length}^2$ and $u^\alpha$ has dimensions of $\text{length}/\text{time}$, the tensor $\mathcal{E}_{\alpha\gamma}$ carries dimensions of $1/\text{time}^2$ -- precisely the dimensions of the Newtonian tidal tensor $\partial^2\Phi/\partial x^i \partial x^j$.

### 8.2.3 The gravitoelectromagnetic decomposition

The name "electrogravitic" is not merely decorative. In the gravitoelectromagnetic (GEM) decomposition of the Weyl curvature tensor -- the trace-free part of Riemann -- the full tidal content splits into two pieces: the electrogravitic part $\mathcal{E}_{\alpha\gamma}$, which governs the tidal stretching and squeezing of geodesic congruences, and the magnetogravitic (or gravitomagnetic) part $\mathcal{H}_{\alpha\gamma}$, which encodes frame-dragging tidal effects from rotating sources. The electrogravitic sector is the direct relativistic generalization of the Newtonian tidal field and dominates overwhelmingly in the weak-field regime relevant to surface gravimetry. The magnetogravitic sector will make its appearance in Section 8.5, where we will find that its contributions to the tidal signal are negligible by many orders of magnitude.

With the tidal tensor in hand, our next task is to take the weak-field, slow-motion limit and show that $\mathcal{E}^i_{\ j}$ reduces to the Hessian of the Newtonian potential -- recovering the tidal formula from Section 6.
## 8.3 From curved spacetime to Newton

The electrogravitic tensor $\mathcal{E}^\alpha_{\ \gamma}$ encodes tidal gravity in fully covariant language, valid in any spacetime and for any observer. But the gravitational fields relevant to surface gravimetry -- the Moon and Sun pulling on an instrument bolted to the Earth -- are emphatically weak and slow. Our task now is to take the limit where spacetime is nearly flat and all motions are non-relativistic, and show that the full GR tidal equation collapses to the familiar Newtonian one. The payoff will be twofold: we recover the exact tidal formula from Section 6 as a special case of geodesic deviation, and we gain a precise understanding of the assumptions behind it, which will let us bound the corrections in Section 8.4.

### 8.3.1 The weak-field metric

We write the metric as a perturbation around Minkowski spacetime:

$$g_{\mu\nu} = \eta_{\mu\nu} + h_{\mu\nu} \,, \qquad |h_{\mu\nu}| \ll 1 \,,$$

where $\eta_{\mu\nu} = \mathrm{diag}(-1, +1, +1, +1)$ in the $(-,+,+,+)$ signature convention. Three assumptions define the Newtonian regime: (1) the gravitational potential is weak, $|\Phi|/c^2 \ll 1$; (2) all relevant velocities are small, $v/c \ll 1$; and (3) the field is static or slowly varying. Under these conditions, the metric in the Newtonian gauge takes the form

$$ds^2 = -\!\left(1 + \frac{2\Phi}{c^2}\right) c^2\, dt^2 + \left(1 - \frac{2\Phi}{c^2}\right)\!\left(dx^2 + dy^2 + dz^2\right),$$

where $\Phi(\mathbf{x})$ is the Newtonian gravitational potential. Reading off the perturbation components:

$$h_{00} = -\frac{2\Phi}{c^2} \,, \qquad h_{ij} = -\frac{2\Phi}{c^2}\, \delta_{ij} \,, \qquad h_{0i} = 0 \,.$$

The equal coefficients in $h_{00}$ and $h_{ij}$ are not a foregone conclusion -- they are a specific prediction of GR. In the Parameterized Post-Newtonian (PPN) framework, the spatial perturbation generalizes to $h_{ij} = -2\gamma_{\text{PPN}}\,(\Phi/c^2)\,\delta_{ij}$, where $\gamma_{\text{PPN}}$ measures the spatial curvature produced per unit rest mass. General Relativity gives $\gamma_{\text{PPN}} = 1$ exactly. Other metric theories predict different values; Brans-Dicke theory, for instance, gives $\gamma_{\text{PPN}} = (1+\omega_{\text{BD}})/(2+\omega_{\text{BD}})$. The Cassini mission's measurement of the Shapiro time delay constrains $|\gamma_{\text{PPN}} - 1| < 2.3 \times 10^{-5}$, so the distinction is negligible for surface gravimetry -- but it is worth noting that the isotropic spatial metric used here is specifically a GR prediction. As we will see in Section 8.4, the value of $\gamma_{\text{PPN}}$ directly affects the coefficient in the 1PN tidal tensor.

### 8.3.2 Recovery of Newton's law

To first order in $h_{\mu\nu}$, the Christoffel symbols reduce to

$$\Gamma^\alpha_{\ \beta\gamma} = \frac{1}{2}\, \eta^{\alpha\delta}\!\left(\partial_\beta h_{\gamma\delta} + \partial_\gamma h_{\beta\delta} - \partial_\delta h_{\beta\gamma}\right) + O(h^2) \,.$$

For the static Newtonian metric with $h_{0i} = 0$, the component that governs gravitational acceleration is

$$\Gamma^i_{\ 00} = -\frac{1}{2}\, \delta^{ij}\, \partial_j h_{00} = -\frac{1}{2}\, \partial_i\!\left(-\frac{2\Phi}{c^2}\right) = \frac{1}{c^2}\, \partial_i \Phi \,.$$

Let us check that this gives the right Newtonian limit. Using coordinates $x^\mu = (x^0, x^i)$ with $x^0 = ct$, a slow-moving particle ($v \ll c$) has four-velocity

$$u^\mu = \frac{dx^\mu}{d\tau} = \left(\frac{dx^0}{d\tau},\, \frac{dx^i}{d\tau}\right) = \left(c\,\frac{dt}{d\tau},\, \frac{dx^i}{dt}\,\frac{dt}{d\tau}\right).$$

In the weak-field, slow-motion regime, $d\tau \approx dt\,\sqrt{1 + 2\Phi/c^2} \approx dt\,(1 + \Phi/c^2) \approx dt$ to leading order, so $dt/d\tau \approx 1$ and

$$u^0 \approx c \,, \qquad u^i = v^i \ll c \,.$$

The spatial geodesic equation reads

$$\frac{d^2 x^i}{d\tau^2} + \Gamma^i_{\ \beta\gamma}\, u^\beta\, u^\gamma = 0 \,.$$

The dominant term in $\Gamma^i_{\ \beta\gamma}\, u^\beta u^\gamma$ comes from $\beta = \gamma = 0$, since $u^0 \approx c$ while $u^i \sim v \ll c$:

$$\Gamma^i_{\ \beta\gamma}\, u^\beta\, u^\gamma = \Gamma^i_{\ 00}\,(u^0)^2 + 2\,\Gamma^i_{\ 0j}\, u^0\, u^j + \Gamma^i_{\ jk}\, u^j\, u^k \approx \Gamma^i_{\ 00}\, c^2 + O(v/c) \,.$$

Since $d\tau \approx dt$, we have $d^2 x^i / d\tau^2 \approx d^2 x^i / dt^2$, and thus

$$\frac{d^2 x^i}{dt^2} \approx -c^2\, \Gamma^i_{\ 00} = -c^2 \cdot \frac{1}{c^2}\, \partial_i \Phi = -\partial_i \Phi \,,$$

recovering Newton's second law $\ddot{\mathbf{x}} = -\nabla \Phi$. This is a satisfying consistency check rather than a surprise: the weak-field metric was constructed to reproduce Newtonian gravity. The real content lies in what happens when we look not at the acceleration itself, but at its spatial variation -- the tidal field.

### 8.3.3 The linearized Riemann tensor and the Newtonian tidal tensor

At linear order in $h_{\mu\nu}$, the Riemann tensor simplifies dramatically. The quadratic $\Gamma\Gamma$ terms are second order in $h$ and drop out:

$$R^\alpha_{\ \beta\gamma\delta} = \partial_\gamma \Gamma^\alpha_{\ \delta\beta} - \partial_\delta \Gamma^\alpha_{\ \gamma\beta} + O(h^2) \,.$$

In fully covariant form:

$$R_{\alpha\beta\gamma\delta} = \frac{1}{2}\!\left(\partial_\gamma \partial_\beta\, h_{\alpha\delta} + \partial_\delta \partial_\alpha\, h_{\beta\gamma} - \partial_\gamma \partial_\alpha\, h_{\beta\delta} - \partial_\delta \partial_\beta\, h_{\alpha\gamma}\right) + O(h^2) \,.$$

The components that matter for tidal physics are the spatial ones. For a slow observer with $u^\mu \approx (c, 0, 0, 0)$, the electrogravitic tensor $\mathcal{E}^i_{\ j} = R^i_{\ 0j0}$ reduces to

$$R^i_{\ 0j0} = \partial_j \Gamma^i_{\ 00} - \partial_0 \Gamma^i_{\ j0} + O(h^2) \,.$$

For a static field the time derivative vanishes ($\partial_0 = (1/c)\,\partial_t = 0$), and substituting $\Gamma^i_{\ 00} = \partial_i \Phi / c^2$:

$$R^i_{\ 0j0} = \frac{1}{c^2}\, \partial_j \partial_i \Phi \,.$$

Now we feed this into the geodesic deviation equation. With $u^0 \approx c$ dominating and $\xi^\gamma$ purely spatial for a spatial separation:

$$\frac{D^2 \xi^i}{D\tau^2} = -R^i_{\ \beta\gamma\delta}\, u^\beta\, \xi^\gamma\, u^\delta \approx -R^i_{\ 0j0}\, (u^0)^2\, \xi^j = -c^2\, R^i_{\ 0j0}\, \xi^j \,,$$

In the slow-motion limit ($d\tau \approx dt$), the covariant derivative $D^2\xi^i/D\tau^2$ reduces to the ordinary $d^2\xi^i/dt^2$ at leading order -- the Christoffel correction terms in $D/D\tau$ contribute only at higher order in $v/c$ and $\Phi/c^2$. Thus

$$\frac{d^2 \xi^i}{dt^2} = -c^2\, R^i_{\ 0j0}\, \xi^j = -\frac{\partial^2 \Phi}{\partial x^i \partial x^j}\, \xi^j \,.$$

Defining the Newtonian tidal tensor

$$\boxed{T_{ij} = \frac{\partial^2 \Phi}{\partial x^i \partial x^j}}$$

the tidal acceleration becomes

$$a^i_{\text{tidal}} = -T_{ij}\, \xi^j \,.$$

The message could not be cleaner: the Hessian of the Newtonian potential *is* the curvature. The relativistic tidal tensor $\mathcal{E}^i_{\ j}$, projected into the observer's rest frame and stripped of its relativistic dressing by the weak-field limit, reduces to exactly the second-derivative matrix $\partial_i \partial_j \Phi$. The entire apparatus of Riemannian geometry -- covariant derivatives, Christoffel connections, the Riemann tensor -- collapses to a statement about the second spatial derivatives of a single scalar function.

The tidal tensor inherits clean properties from its geometric origin. It is symmetric, $T_{ij} = T_{ji}$, by equality of mixed partials. In vacuum, Laplace's equation $\nabla^2 \Phi = 0$ forces $T_{ii} = 0$ -- the tidal tensor is trace-free. This is the Newtonian shadow of the vacuum Einstein equation $R_{\mu\nu} = 0$, which at this order gives $R^i_{\ 0i0} = \nabla^2 \Phi / c^2 = 0$. Inside matter, $\nabla^2 \Phi = 4\pi G \rho$, so $T_{ii} = 4\pi G \rho$. These limiting cases provide a useful sanity check: the trace of the tidal tensor measures the local mass density, as it must.

### 8.3.4 Recovery of the exact tidal formula

With the tidal tensor established, we can now complete the circle and recover the exact Newtonian tidal formula that Pytheas uses. Consider a gravimeter on the Earth's surface at position $\mathbf{R}$ from the Earth's center, and an external body (Moon or Sun) of mass $M$ at position $\mathbf{r}_M$ from the Earth's center. The gravitational potential due to the external body at a point $\mathbf{x}$ near the gravimeter is

$$\Phi_M(\mathbf{x}) = -\frac{GM}{|\mathbf{r}_M - \mathbf{x}|} \,.$$

A surface gravimeter measures the acceleration of its test mass relative to its housing, which is fixed to the Earth. To leading order, the housing follows the motion of the Earth's center, which is itself in free fall in the external gravitational field. The tidal acceleration is therefore the difference between the gravitational pull at the gravimeter and the pull at the Earth's center:

$$\mathbf{g}_{\text{test}} = -\nabla \Phi_M\big|_{\mathbf{R}} = -GM\, \frac{\mathbf{R} - \mathbf{r}_M}{|\mathbf{R} - \mathbf{r}_M|^3} = GM\, \frac{\mathbf{r}_M - \mathbf{R}}{|\mathbf{r}_M - \mathbf{R}|^3} \,,$$

$$\mathbf{g}_{\text{center}} = -\nabla \Phi_M\big|_{\mathbf{0}} = -GM\, \frac{-\mathbf{r}_M}{|\mathbf{r}_M|^3} = GM\, \frac{\mathbf{r}_M}{|\mathbf{r}_M|^3} \,.$$

Their difference gives the tidal acceleration:

$$\boxed{\mathbf{a}_{\text{tidal}} = \mathbf{g}_{\text{test}} - \mathbf{g}_{\text{center}} = GM\!\left[\frac{\mathbf{r}_M - \mathbf{R}}{|\mathbf{r}_M - \mathbf{R}|^3} - \frac{\mathbf{r}_M}{|\mathbf{r}_M|^3}\right]}$$

This is exactly the formula from Section 6, derived there from elementary Newtonian reasoning. The GR derivation gives us the same answer, but now we know *why* it works -- it is the weak-field limit of geodesic deviation in curved spacetime -- and, crucially, *when* it will fail: when $|\Phi|/c^2$ or $v/c$ ceases to be negligibly small.

To verify consistency with the tidal tensor, we can expand for $|\mathbf{R}| \ll |\mathbf{r}_M|$. Define $r_M = |\mathbf{r}_M|$ and the unit vector $\hat{n} = \mathbf{r}_M / r_M$. Taylor-expanding the inverse distance:

$$\frac{1}{|\mathbf{r}_M - \mathbf{R}|} = \frac{1}{r_M} + \frac{R_i\, n_i}{r_M^2} + \frac{1}{2}\, \frac{3(R_i n_i)^2 - R^2}{r_M^3} + \cdots$$

The monopole term contributes a constant (irrelevant), and the dipole term gives the uniform acceleration of the Earth's center (which the gravimeter subtracts). The leading tidal contribution comes from the quadrupole:

$$\Phi_{\text{tidal}} = -\frac{GM}{2r_M^3}\!\left(3(\hat{n} \cdot \mathbf{R})^2 - R^2\right).$$

Taking second derivatives recovers the tidal tensor for a point-mass source:

$$T_{ij} = \frac{\partial^2 \Phi_{\text{tidal}}}{\partial R_i \partial R_j} = -\frac{GM}{r_M^3}\!\left(3\, n_i\, n_j - \delta_{ij}\right),$$

and the tidal acceleration is $a^i_{\text{tidal}} = -T_{ij}\, R_j$, which to leading order in $R/r_M$ agrees with the exact formula above. The familiar $3\,n_i\,n_j - \delta_{ij}$ structure -- stretching along the line to the perturbing body, compression perpendicular to it -- is the signature of the trace-free, quadrupolar tidal field.

The chain is now complete. Starting from the geodesic deviation equation in curved spacetime, passing through the electrogravitic tidal tensor, and taking the weak-field limit, we have arrived back at the Newtonian tidal formula that Pytheas uses to compute $g(t)$. The question that remains is how large the corrections are when we go beyond the Newtonian approximation -- which is the business of the next two sections.
## 8.4 Post-Newtonian corrections to the tidal tensor

The preceding section established that the Newtonian tidal tensor $T_{ij} = \partial_i \partial_j \Phi$ is the weak-field limit of the full Riemann curvature, and that it reproduces the exact tidal formula from Section 6. But "weak-field limit" is a controlled approximation, not an article of faith. The natural follow-up question is: how large are the corrections we neglected, and could any of them matter at Pytheas's $\sim$10 nGal accuracy floor? The post-Newtonian framework gives a systematic answer.

### 8.4.1 The 1PN metric

The first post-Newtonian (1PN) approximation extends the linearized metric by including terms of order $\epsilon \sim v^2/c^2 \sim GM/(Rc^2)$ beyond the Newtonian pieces. For a single body of mass $M$, the metric in isotropic coordinates reads (maintaining the convention $\Phi = -GM/r < 0$ from Section 8.3):

$$g_{00} = -\!\left(1 + \frac{2\Phi}{c^2} + \frac{2\Phi^2}{c^4}\right) + O(c^{-6}) \,,$$

$$g_{ij} = \left(1 - \frac{2\Phi}{c^2}\right) \delta_{ij} + O(c^{-4}) \,,$$

$$g_{0i} = -\frac{4}{c^3}\, \frac{G(\mathbf{J} \times \hat{r})_i}{r^2} + O(c^{-5}) \,,$$

where $\mathbf{J}$ is the angular momentum of the source. The off-diagonal components $g_{0i}$ encode the gravitomagnetic (Lense-Thirring) field and are treated separately in Section 8.5. Note that $g_{00}$ can equivalently be written in the PPN form $g_{00} = -(1 - 2U/c^2 + 2\beta U^2/c^4)$ with $U = GM/r > 0$ and $\beta = 1$ in GR; substituting $U = -\Phi$ recovers the expression above.

A remark on the single-body treatment is in order. In what follows we compute the 1PN tidal tensor produced by one source mass $M$ at a time, treating each external body (Sun, Moon) independently. This is justified because cross-terms between different sources in the multi-body 1PN metric enter at order $\Phi_1 \Phi_2 / c^4$. For the Sun-Moon system at the Earth, $\Phi_{\text{Moon}}/c^2 \sim 10^{-13}$ and $\Phi_{\text{Sun}}/c^2 \sim 10^{-8}$, so the cross-term correction is of order $10^{-21}$ relative to the Newtonian tidal acceleration -- suppressed by a further factor of $\sim 10^{-13}$ compared to the already small single-body 1PN correction. The independent treatment is therefore entirely adequate.

### 8.4.2 Construction of the 1PN tidal tensor

We now build up the 1PN Riemann component $R^i_{\ 0j0}$ step by step, keeping all terms through order $c^{-4}$ (one order beyond the Newtonian $c^{-2}$).

**Step 1: The 1PN Christoffel symbol.** Setting aside the gravitomagnetic $g_{0i}$ terms (treated in Section 8.5), the inverse spatial metric is $g^{ij} = (1 + 2\Phi/c^2)\,\delta^{ij} + O(c^{-4})$, and the derivative of $g_{00}$ picks up the $\Phi^2$ correction. Together these give:

$$\Gamma^i_{\ 00} = \frac{1}{c^2}\, \partial_i\Phi\!\left(1 + \frac{4\Phi}{c^2}\right) + O(c^{-6}) \,.$$

The coefficient 4 has a specific origin: one factor of $(1 + 2\Phi/c^2)$ comes from inverting $g_{ij}$ to get $g^{ij}$, and a second factor of $(1 + 2\Phi/c^2)$ comes from differentiating $g_{00} = -(1 + 2\Phi/c^2 + 2\Phi^2/c^4)$. Their product is $(1 + 2\Phi/c^2)^2 = 1 + 4\Phi/c^2 + O(c^{-4})$, which produces the combined coefficient of 4. This is a distinctive prediction of GR: in a general PPN framework the coefficient would be $2(1 + \gamma_{\text{PPN}})$, which equals 4 only when $\gamma_{\text{PPN}} = 1$.

**Step 2: Spatial derivative of the Christoffel symbol.** Taking $\partial_j$ of the expression above and applying the product rule:

$$\partial_j \Gamma^i_{\ 00} = \frac{1}{c^2}\,\partial_j\partial_i\Phi\!\left(1 + \frac{4\Phi}{c^2}\right) + \frac{4}{c^4}\,(\partial_i\Phi)(\partial_j\Phi) + O(c^{-6}) \,.$$

The second term is a genuinely nonlinear contribution -- a product of two first derivatives of the potential -- that has no Newtonian counterpart.

**Step 3: The $\Gamma\Gamma$ contribution.** At linearized order, the quadratic $\Gamma\Gamma$ terms in the Riemann tensor were negligible. At 1PN order they enter at $O(c^{-4})$ and must be included. The leading-order spatial Christoffel symbols are:

$$\Gamma^i_{\ jk} = -\frac{1}{c^2}\!\left(\delta^i_k\,\partial_j\Phi + \delta^i_j\,\partial_k\Phi - \delta_{jk}\,\partial^i\Phi\right) + O(c^{-4}) \,.$$

For a static metric, $\Gamma^i_{\ 0j} = 0$ at leading order, so the cross-term $\Gamma^i_{\ 0\lambda}\,\Gamma^\lambda_{\ j0}$ vanishes through $O(c^{-4})$. The surviving piece is:

$$\sum_l \Gamma^i_{\ jl}\, \Gamma^l_{\ 00} = -\frac{1}{c^4}\,\delta_{ij}\,|\nabla\Phi|^2 \,.$$

The algebra is worth a brief look: expanding $\Gamma^i_{\ jl}\,\Gamma^l_{\ 00}$ produces three terms in the bracket, $(\partial_j\Phi)(\partial_i\Phi) + \delta_{ij}|\nabla\Phi|^2 - (\partial_i\Phi)(\partial_j\Phi)$. The first and third cancel, leaving only the isotropic piece proportional to $\delta_{ij}$.

**Step 4: Assembly.** Combining the derivative term from Step 2 with the $\Gamma\Gamma$ term from Step 3 (and noting that $\partial_0 = 0$ for a static field eliminates the time-derivative piece in $R^i_{\ 0j0} = \partial_j \Gamma^i_{\ 00} - \partial_0 \Gamma^i_{\ j0} + \Gamma\Gamma$ terms), the full 1PN Riemann component is:

$$R^i_{\ 0j0} = \frac{1}{c^2}\,\partial_i\partial_j\Phi\!\left(1 + \frac{4\Phi}{c^2}\right) + \frac{4}{c^4}\,(\partial_i\Phi)(\partial_j\Phi) - \frac{1}{c^4}\,\delta_{ij}\,|\nabla\Phi|^2 + O(c^{-6}) \,.$$

Since the tidal acceleration is $a^i = -c^2\, R^i_{\ 0j0}\, \xi^j$, the 1PN tidal tensor (defined so that $a^i_{\text{tidal}} = -T^{\text{1PN}}_{ij}\,\xi^j$) is:

$$T_{ij}^{\text{1PN}} = \left(\partial_i\partial_j\Phi\right)\!\left(1 + \frac{4\Phi}{c^2}\right) + \frac{4}{c^2}\,(\partial_i\Phi)(\partial_j\Phi) - \frac{1}{c^2}\,\delta_{ij}\,|\nabla\Phi|^2$$

This is the central result of the section. The correction relative to the Newtonian tidal tensor $T_{ij} = \partial_i \partial_j \Phi$ consists of three terms, each of order $\Phi/c^2$:

- A multiplicative correction $(1 + 4\Phi/c^2)$ that rescales the entire Newtonian tidal tensor. Since $\Phi < 0$, this slightly *weakens* the tidal field.
- A nonlinear gradient-gradient term $\propto (\partial_i\Phi)(\partial_j\Phi)/c^2$, arising from the 1PN correction to $\Gamma^i_{\ 00}$.
- An isotropic trace contribution $\propto \delta_{ij}\,|\nabla\Phi|^2/c^2$, originating entirely from the $\Gamma\Gamma$ terms in the Riemann tensor.

Each term can be checked on dimensional grounds: $\Phi/c^2$ is dimensionless, and $(\partial_i\Phi)(\partial_j\Phi)/c^2$ has dimensions of $(\text{acceleration})^2/c^2 = 1/\text{time}^2 \cdot (\text{length}/c)^0$... which matches the dimensions of $\partial_i\partial_j\Phi$. The structure is self-consistent.

### 8.4.3 Magnitude estimates

The two small parameters controlling the post-Newtonian expansion are related by the virial theorem for gravitationally bound orbits: $v^2/c^2 \sim GM/(Rc^2)$. Let us verify this and then estimate the 1PN corrections for both the Sun and the Moon.

**The Sun.** The Earth's orbital velocity is $v_\oplus \approx 29.8$ km/s, giving:

$$\epsilon_v = \frac{v_\oplus^2}{c^2} = \left(\frac{29.8\ \text{km/s}}{3.0 \times 10^5\ \text{km/s}}\right)^2 \approx 9.9 \times 10^{-9} \,,$$

$$\epsilon_\Phi = \frac{GM_\odot}{R_{\text{ES}}\, c^2} \approx 9.9 \times 10^{-9} \,.$$

Both are of order $10^{-8}$, confirming $v^2/c^2 \sim GM/(Rc^2)$ as expected. The dominant 1PN correction to the tidal tensor is the multiplicative factor $4\Phi/c^2 = -4GM/(rc^2)$. For the Sun's tidal field at the Earth, this gives a relative correction of magnitude $4\,GM_\odot/(R_{\text{ES}}\,c^2) \approx 4.0 \times 10^{-8}$. Multiplying by the Newtonian solar tidal acceleration from Section 8.3:

$$\delta a_{\text{tidal}}^{\text{1PN, Sun}} \sim 4\epsilon_\Phi \times a_{\text{tidal}}^{\text{N, Sun}} \approx 4 \times 10^{-8} \times 50\ \mu\text{Gal} \approx 0.002\ \text{nGal} \,.$$

**The Moon.** Here one must be careful about which mass controls the expansion parameter. The tidal field is generated by the Moon's mass, so the 1PN correction to its tidal tensor is governed by the Moon's own gravitational parameter evaluated at the Earth:

$$\epsilon_\Phi^{\text{Moon}} = \frac{GM_{\text{Moon}}}{r_{\text{EM}}\, c^2} \approx 1.4 \times 10^{-13} \,.$$

Why is this so much smaller than the Sun's $\epsilon_\Phi$? Because the Moon's orbital velocity ($v_{\text{Moon}} \approx 1.02$ km/s, giving $v_{\text{Moon}}^2/c^2 \approx 1.2 \times 10^{-11}$) is *not* the correct parameter here. The virial theorem $v^2 \sim GM/r$ for the Earth-Moon system relates the Moon's orbital velocity to the *Earth's* mass: $v_{\text{Moon}}^2 \sim GM_{\text{Earth}}/r_{\text{EM}}$. But the tidal field is generated by the Moon alone, so the relevant 1PN parameter is $GM_{\text{Moon}}/(r_{\text{EM}}\,c^2)$, which is smaller by a factor of $M_{\text{Moon}}/M_{\text{Earth}} \approx 1/81$.

The resulting 1PN correction:

$$\delta a_{\text{tidal}}^{\text{1PN, Moon}} \sim 4\epsilon_\Phi^{\text{Moon}} \times a_{\text{tidal}}^{\text{N, Moon}} \approx 5.7 \times 10^{-13} \times 110\ \mu\text{Gal} \approx 6 \times 10^{-8}\ \text{nGal} \,,$$

which is negligible by many orders of magnitude.

### 8.4.4 Summary

The total 1PN tidal correction is dominated by the solar term:

$$\delta a_{\text{tidal}}^{\text{1PN, total}} \lesssim 0.002\ \text{nGal} \ll 10\ \text{nGal} \,.$$

This is more than three orders of magnitude below Pytheas's accuracy floor. Figure 14 summarizes the hierarchy of GR corrections, showing how each successive order of the post-Newtonian expansion falls further below Pytheas's accuracy floor. The 1PN corrections to the tidal tensor are real -- they are measured in satellite orbit determination, where the IERS conventions include them -- but for surface gravimetry they are irrelevant. The Newtonian tidal formula from Section 6 stands uncorrected.
## 8.5 Gravitomagnetic (Lense-Thirring) tidal effects

The 1PN corrections of Section 8.4 arise entirely from the diagonal metric components $g_{00}$ and $g_{ij}$ -- the gravitoelectric sector. But a rotating mass also generates off-diagonal components $g_{0i}$, producing a gravitomagnetic field $\mathbf{B}_g$ that has no Newtonian counterpart. The analogy to electromagnetism is immediate and instructive: just as a moving electric charge generates a magnetic field, a moving or rotating mass generates a gravitomagnetic field. The parallel is not merely heuristic -- it is encoded in the gravitoelectromagnetic (GEM) decomposition of the Weyl tensor, whose electric part $\mathcal{E}_{\alpha\gamma}$ we encountered in Section 8.2 and whose magnetic part $\mathcal{H}_{ij}$ governs the effects considered here.

There is, however, a crucial difference from electromagnetism. Gravity is mediated by a spin-2 field, not spin-1, and this shows up as a factor of 2 in the gravitomagnetic force law. We need to trace how these off-diagonal metric components propagate through the Riemann tensor and into the geodesic deviation equation, then estimate whether they matter for Pytheas.

### 8.5.1 The gravitomagnetic field

In the weak-field limit, a body with angular momentum $\mathbf{J}$ sources a gravitomagnetic vector potential

$$\mathbf{A}_g = \frac{2G}{c^2}\, \frac{\mathbf{J} \times \hat{r}}{r^2} \,,$$

and the corresponding gravitomagnetic field is its curl:

$$\mathbf{B}_g = \nabla \times \mathbf{A}_g = \frac{2G}{c^2 r^3}\!\left[3(\mathbf{J} \cdot \hat{r})\, \hat{r} - \mathbf{J}\right].$$

This has the same structure as a magnetic dipole field -- an $r^{-3}$ dipolar pattern aligned with the body's spin axis. The equation of motion for a slowly moving test mass in the combined gravitoelectric and gravitomagnetic fields is

$$\ddot{\mathbf{x}} = -\nabla \Phi - 2\, \mathbf{v} \times \mathbf{B}_g + O(c^{-4}) \,.$$

The factor of 2 in front of the gravitomagnetic term (compared to the factor of 1 in the Lorentz force law $\mathbf{F} = q\,\mathbf{v} \times \mathbf{B}$) is the hallmark of the spin-2 graviton. It is this same factor that was confirmed, at the level of frame-dragging precession, by the Gravity Probe B experiment and LAGEOS satellite laser ranging measurements.

### 8.5.2 Gravitomagnetic tidal tensor from the Riemann tensor

To understand the gravitomagnetic contribution to tidal forces, we return to the full geodesic deviation equation and expand the four-velocity sum. The spatial tidal acceleration is

$$\frac{d^2\xi^i}{dt^2} = -R^i_{\ \alpha j\beta}\, u^\alpha\, u^\beta\, \xi^j \,.$$

Expanding over $\alpha, \beta \in \{0, k\}$ with $u^0 = c$ and $u^k = v^k$:

$$\frac{d^2\xi^i}{dt^2} = -R^i_{\ 0j0}\, c^2\, \xi^j - 2\, R^i_{\ 0jk}\, c\, v^k\, \xi^j - R^i_{\ kjl}\, v^k\, v^l\, \xi^j \,.$$

The first term is the electrogravitic tidal tensor treated in Sections 8.3 and 8.4. The third term is of order $v^2/c^2$ relative to the first and is already bounded by the 1PN analysis. The middle term, proportional to $R^i_{\ 0jk}$, is the gravitomagnetic tidal contribution -- and it carries a distinctive velocity dependence that has no counterpart in the Newtonian theory.

Let us compute $R^i_{\ 0jk}$ from the metric. At linearized order, the relevant Christoffel symbols involving $g_{0i}$ are

$$\Gamma^i_{\ 0j} = \frac{1}{2}\,\delta^{il}\!\left(\partial_j g_{0l} - \partial_l g_{0j}\right) + O(h^2) \,,$$

where we used stationarity ($\partial_0 g_{jl} = 0$). Writing $g_{0i} = -(2/c)\,A_{g,i}$, this becomes

$$\Gamma^i_{\ 0j} = -\frac{1}{c}\!\left(\partial_j A_{g}^i - \partial^i A_{g,j}\right),$$

which is directly related to the curl of $\mathbf{A}_g$. In terms of $\mathbf{B}_g = \nabla \times \mathbf{A}_g$:

$$\Gamma^i_{\ 0j} = \frac{1}{c}\,\epsilon^i_{\ jk}\, B_g^k + \text{symmetric part} \,.$$

The linearized Riemann component follows by differentiation:

$$R^i_{\ 0jk} = \partial_j\,\Gamma^i_{\ 0k} - \partial_k\,\Gamma^i_{\ 0j} \,.$$

Substituting the expression for $\Gamma^i_{\ 0j}$ in terms of $\mathbf{A}_g$ and simplifying (the $\partial_j \partial_k A_g^i$ terms cancel by symmetry of mixed partials), we arrive at

$$R^i_{\ 0jk} = \frac{1}{c}\,\partial^i\!\left(\partial_j A_{g,k} - \partial_k A_{g,j}\right).$$

Recognizing that $\partial_j A_{g,k} - \partial_k A_{g,j} = \epsilon_{jkm}\, B_g^m$, this yields

$$R^i_{\ 0jk} = \frac{1}{c}\,\epsilon_{jkm}\,\partial^i B_g^m \,.$$

The gravitomagnetic tidal acceleration then follows immediately from the geodesic deviation expansion:

$$\delta a^i_{\text{gm-tidal}} = -2\,R^i_{\ 0jk}\,c\,v^k\,\xi^j = -2\,\epsilon_{jkm}\,(\partial^i B_g^m)\, v^k\, \xi^j \,.$$

The physical content is now transparent. The gravitomagnetic tidal tensor involves the spatial gradient of $\mathbf{B}_g$, not the field itself. Since $\mathbf{B}_g$ falls off as $1/r^3$ for a dipole, its gradient falls off as $1/r^4$ -- even faster than the $1/r^3$ of the Newtonian tidal tensor. Distance kills the gravitomagnetic tidal effect with an extra power of $r$.

This connects to the GEM tidal decomposition mentioned in Section 8.2. The gravitomagnetic tidal tensor $\mathcal{H}_{ij}$, defined as the magnetic part of the Weyl tensor projected onto the observer's frame, is

$$\mathcal{H}_{ij} = \frac{1}{2}\,\epsilon_{ikl}\,R^{kl}_{\ \ 0j}\, u^0 = \frac{c}{2}\,\epsilon_{ikl}\,R^{kl}_{\ \ 0j} \,.$$

This encodes the same information as the $R^i_{\ 0jk}$ components computed above, confirming that the gravitomagnetic tidal acceleration on a test mass with velocity $v^k$ scales as $\delta a^i \sim \mathcal{H}^i_{\ j}\, v^j\, \xi / c$.

### 8.5.3 Magnitude estimates

With the Riemann-tensor derivation in hand, we can bound every gravitomagnetic contribution. The scaling of the gravitomagnetic tidal acceleration for a spinning body at distance $r$ is

$$\delta a_{\text{gm-tidal}} \sim 2\, v_{\text{test}} \cdot \frac{GJ}{c^2\, r^4} \cdot R_\oplus \,.$$

**Earth's own frame-dragging.** The Earth has angular momentum $J_\oplus = I_\oplus\, \omega_\oplus \approx 5.86 \times 10^{33}$ kg m$^2$/s. The gravitomagnetic field strength at the surface is

$$B_g \sim \frac{2G J_\oplus}{c^2 R_\oplus^3} \approx 3.4 \times 10^{-14}\ \text{s}^{-1} \,.$$

For a test mass falling in a gravimeter with velocity $v_{\text{test}} \sim 1$ m/s, the direct Lense-Thirring acceleration is

$$a_{\text{LT}} \sim 2\, v_{\text{test}} \cdot B_g \sim 6.7 \times 10^{-14}\ \text{m/s}^2 \sim 0.007\ \text{nGal} \,.$$

This is not a tidal effect -- it is a direct gravitomagnetic acceleration on the test mass -- but it is the largest gravitomagnetic contribution of any kind, and it is still below the 10 nGal accuracy floor.

**Spin-induced tidal effects from external bodies.** For the Sun's angular momentum ($J_\odot \approx 1.93 \times 10^{41}$ kg m$^2$/s) at the Earth's distance, with $v_{\text{test}} \sim 1$ m/s:

$$\delta a_{\text{LT, tidal}}^{\text{Sun}} \sim 2 \times 1\ \text{m/s} \times \frac{G J_\odot\, R_\oplus}{c^2\, R_{\text{ES}}^4} \sim 4 \times 10^{-24}\ \text{m/s}^2 \sim 4 \times 10^{-13}\ \text{nGal} \,.$$

For the Moon ($J_{\text{Moon}} \approx 2.3 \times 10^{29}$ kg m$^2$/s):

$$\delta a_{\text{LT, tidal}}^{\text{Moon}} \sim 2 \times 1\ \text{m/s} \times \frac{G J_{\text{Moon}}\, R_\oplus}{c^2\, R_{\text{EM}}^4} \sim 10^{-25}\ \text{m/s}^2 \sim 10^{-14}\ \text{nGal} \,.$$

**Orbital gravitomagnetic tidal contribution.** There is also a gravitomagnetic field from the orbital motion of external bodies, not just their spin. The orbital gravitomagnetic field at the Earth due to the Moon's orbital velocity $v_{\text{Moon}} \approx 1$ km/s scales as

$$B_g^{\text{orbital}} \sim \frac{v_{\text{Moon}}}{c}\, \frac{GM_{\text{Moon}}}{c\, R_{\text{EM}}^2} \sim \frac{v_{\text{Moon}}}{c}\, \frac{\Phi_{\text{Moon}}}{c\, R_{\text{EM}}} \,.$$

The associated tidal contribution is suppressed by an additional factor of $v/c \sim 10^{-5}$ relative to the already negligible 1PN correction:

$$\delta a_{\text{LT, orbital}} \sim \frac{v_{\text{Moon}}}{c} \cdot \delta a_{\text{tidal}}^{\text{1PN, Moon}} \sim 10^{-13}\ \text{nGal} \,.$$

The following table collects all gravitomagnetic contributions:

| Source | Mechanism | Magnitude |
|--------|-----------|-----------|
| Earth spin | Direct $2\,\mathbf{v} \times \mathbf{B}_g$ on test mass | $\sim 0.007$ nGal |
| Sun spin | Gravitomagnetic tidal tensor | $\sim 4 \times 10^{-13}$ nGal |
| Moon spin | Gravitomagnetic tidal tensor | $\sim 10^{-14}$ nGal |
| Moon orbit | Orbital gravitomagnetic tidal | $\sim 10^{-13}$ nGal |

Every entry is far below the 10 nGal accuracy floor. The largest gravitomagnetic effect -- the Earth's own frame-dragging acting on the falling test mass -- contributes $\sim 0.007$ nGal, which is three orders of magnitude smaller than the Newtonian tidal signals and still a factor of $\sim 1400$ below the accuracy threshold. The spin and orbital gravitomagnetic tidal effects from external bodies are smaller still by ten or more orders of magnitude, suppressed by the $1/r^4$ gradient falloff that the Riemann-tensor derivation made explicit.

The conclusion from the gravitomagnetic sector mirrors that from the 1PN sector of Section 8.4: these are real physical effects, rigorously derived from the Riemann tensor, but they are utterly negligible for surface gravimetry. The Newtonian tidal formula remains uncorrected.
## 8.6 The hierarchy of corrections

We have now assembled every piece of the relativistic tidal picture: the geodesic deviation equation (Section 8.1), the electrogravitic tidal tensor (Section 8.2), its Newtonian limit (Section 8.3), the first post-Newtonian corrections (Section 8.4), and the gravitomagnetic contributions (Section 8.5). It is time to step back and survey the full landscape of corrections from a single vantage point.

### 8.6.1 The complete hierarchy

The following table collects all contributions to the tidal acceleration measured by a surface gravimeter, ordered by magnitude. The rightmost column expresses each effect as a fraction of Pytheas's accuracy floor of approximately 10 nGal, as established by the error budget in Section 1.

| Effect | Magnitude | Ratio to accuracy floor |
|--------|-----------|------------------------|
| Newtonian lunar tidal | $\sim 110\ \mu\text{Gal} = 1.1 \times 10^5$ nGal | $\sim 10^4$ |
| Newtonian solar tidal | $\sim 50\ \mu\text{Gal} = 5.0 \times 10^4$ nGal | $\sim 5 \times 10^3$ |
| Earth gravitomagnetic (direct) | $\sim 0.007$ nGal | $\sim 7 \times 10^{-4}$ |
| 1PN tidal correction (Sun) | $\sim 0.002$ nGal | $\sim 2 \times 10^{-4}$ |
| 1PN tidal correction (Moon) | $\sim 6 \times 10^{-8}$ nGal | $\sim 6 \times 10^{-9}$ |
| Gravitomagnetic tidal (all sources) | $< 10^{-12}$ nGal | $< 10^{-13}$ |

The pattern is unmistakable. The Newtonian tidal signals tower above the accuracy floor by three to four orders of magnitude -- these are the signals that Pytheas must model with care, and Section 6 provides the exact formula for doing so. Every relativistic correction, by contrast, falls *below* the accuracy floor by three or more orders of magnitude. There is a gap of at least six orders of magnitude separating the Newtonian tidal signals from the largest GR correction.

### 8.6.2 Where would GR corrections become relevant?

The post-Newtonian expansion gives us a precise answer to this question. The dominant correction is the factor $4\Phi/c^2$ multiplying the Newtonian tidal tensor in the 1PN expression from Section 8.4:

$$T_{ij}^{\text{1PN}} = \left(\partial_i\partial_j\Phi\right)\!\left(1 + \frac{4\Phi}{c^2}\right) + \frac{4}{c^2}\,(\partial_i\Phi)(\partial_j\Phi) - \frac{1}{c^2}\,\delta_{ij}\,|\nabla\Phi|^2$$

For the 1PN correction to reach the nGal level, one would need either $v/c$ or $GM/(Rc^2)$ to approach $\sim 10^{-4}$. No solar system body relevant to Earth-surface gravimetry comes close: for the Sun, $GM_\odot/(R_{\text{ES}} c^2) \approx 10^{-8}$, and for the Moon, $GM_{\text{Moon}}/(r_{\text{EM}} c^2) \approx 10^{-13}$. Both are four to nine orders of magnitude too small. Higher-order post-Newtonian corrections (2PN and beyond), suppressed by $(v/c)^4 \sim 10^{-16}$ relative to the Newtonian terms, are of no conceivable relevance to surface gravimetry.

It is worth noting that the IERS conventions for satellite orbit determination *do* include 1PN tidal corrections. This is not a contradiction: satellite geodesy operates in a regime where the accumulated effect of tiny accelerations over many orbital periods can become significant, and the accuracy requirements for precise orbit determination are different from those for modeling the tidal signal at a fixed surface station. For Pytheas, the distinction is academic -- the corrections are negligible by a comfortable margin.

### 8.6.3 The conceptual role of general relativity

Although Newtonian gravity is operationally sufficient for Pytheas at the nGal level, the GR derivation developed across Sections 8.1--8.5 serves three essential conceptual purposes.

First, it establishes *why* tidal forces exist. They are not merely a mathematical consequence of the $1/r^2$ law; they are the observable manifestation of spacetime curvature. The Riemann tensor governs the relative acceleration of nearby geodesics -- this is the content of the geodesic deviation equation -- and tidal gravity is intrinsic geometry that cannot be removed by any choice of reference frame. Einstein's elevator can eliminate uniform gravity through free fall, but it cannot eliminate tidal gravity. This is the fundamental reason why gravity cannot be shielded: tidal effects are woven into the fabric of spacetime itself.

Second, the post-Newtonian framework provides a *controlled perturbative expansion* that rigorously bounds the errors incurred by using the Newtonian approximation. Without GR, one could only say that Newton's theory "seems to work" for surface gravimetry. With GR, one can say precisely *how well* it works: the leading correction is of order $4GM/(Rc^2) \sim 4 \times 10^{-8}$ relative to the Newtonian tidal acceleration, contributing less than 0.002 nGal for the Sun and less than $6 \times 10^{-8}$ nGal for the Moon. The gravitomagnetic corrections, derived from the Riemann components $R^i_{\ 0jk}$, are smaller still. These are not estimates or guesses; they are rigorous bounds from a well-controlled expansion.

Third, the GR analysis identifies the *precise conditions* under which the Newtonian treatment would become insufficient. The expansion parameter is $\epsilon \sim v^2/c^2 \sim GM/(Rc^2)$, and the Newtonian approximation breaks down when $\epsilon$ approaches unity. For all bodies in the solar system relevant to Earth-surface gravimetry, $\epsilon$ ranges from $10^{-8}$ (the Sun) to $10^{-13}$ (the Moon), placing us firmly and deeply in the Newtonian regime. Figure 14 collects these bounds visually, placing each GR correction on a logarithmic scale relative to the Newtonian tidal signal and Pytheas's accuracy floor.

### 8.6.4 Closing the chain

We now have the complete chain from first principles to operational gravimetry. Coordinates (Sections 2--3) tell us where we are on the rotating, deformable Earth. Ephemerides (Sections 4--5) tell us where the Moon and Sun are at any instant. The tidal formula (Section 6) tells us how they pull:

$$\mathbf{a}_{\text{tidal}} = GM\!\left[\frac{\mathbf{r}_M - \mathbf{R}}{|\mathbf{r}_M - \mathbf{R}|^3} - \frac{\mathbf{r}_M}{|\mathbf{r}_M|^3}\right]$$

Elastic response (Section 7) amplifies and modifies the signal through Love numbers. And this section -- the GR foundation -- confirms that the Newtonian framework underlying it all is more than adequate, with all relativistic corrections bounded below $\sim 0.01$ nGal, three orders of magnitude beneath Pytheas's accuracy floor. The exact Newtonian tidal formula from Section 6 is the correct formula for Pytheas. No relativistic amendment is needed.

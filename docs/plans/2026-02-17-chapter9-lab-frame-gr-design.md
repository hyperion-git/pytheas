# Design: Chapter 9 — General Relativity in the Lab Frame

**Date:** 2026-02-17
**Status:** Approved design, pending implementation
**Replaces:** Section 8 (`doc/_sec8_gr_tidal.md`)

---

## 1. Motivation

Section 8 derives the tidal acceleration from geodesic deviation in curved spacetime,
recovering the Newtonian tidal formula and bounding post-Newtonian corrections. However,
it has a structural gap: geodesic deviation describes two **freely falling** observers,
while the gravimeter lab is **non-inertial** — supported against gravity with proper
acceleration ~9.8 m/s² and rotating with Earth at $\omega \approx 7.3 \times 10^{-5}$ s$^{-1}$.

The entire mapping from the freely-falling geodesic deviation framework to the lab frame
currently rests on one sentence (line 193):

> "To leading order, the housing follows the motion of the Earth's center, which is
> itself in free fall in the external gravitational field."

This is correct but does enormous load-bearing work without justification. Worse, Section 8
devotes ~110 lines to the Lense-Thirring gravitomagnetic field ($B_g \sim 3.4 \times
10^{-14}$ s$^{-1}$) while never mentioning the coordinate gravitomagnetic field from
Earth's rotation ($\omega \sim 7.3 \times 10^{-5}$ s$^{-1}$), which is $10^9$ times
larger and produces structurally identical $g_{0i}$ metric cross-terms.

Chapter 9 replaces Section 8 with a complete derivation in the **proper reference frame
of the lab** using Fermi normal coordinates, where all effects — static gravity, rotation,
and tidal curvature — emerge from a single metric expansion.

---

## 2. Framework: Fermi Normal Coordinates

The natural coordinate system for a non-geodesic, rotating observer in curved spacetime.
The metric expanded to quadratic order in spatial coordinates $X^i$ around the observer's
worldline encodes three distinct contributions:

1. **Proper acceleration $a^i$** — appears in $g_{00}$ at linear order in $X^i$
2. **Rotation $\Omega^i$** — appears in $g_{0i}$ at linear order in $X^i$
3. **Riemann curvature $R_{\alpha\beta\gamma\delta}$** — appears in all components at
   quadratic order in $X^i$

The metric to quadratic order (Ni & Zimmerman 1978, Li & Ni 1979):

$$g_{00} = -\left(1 + \frac{a_i X^i}{c^2}\right)^2 + R_{0i0j}\frac{X^i X^j}{c^2} + \cdots$$

$$g_{0i} = \frac{2}{3} R_{0jik} X^j X^k + \epsilon_{ijk} \Omega^j X^k + \cdots$$

$$g_{ij} = \delta_{ij} - \frac{1}{3} R_{ikjl} X^k X^l + \cdots$$

This is cited as a standard result, not re-derived.

---

## 3. Chapter Structure

### 9.1 The Lab Observer's Worldline

- **Content**: Define the lab as a non-geodesic worldline in Schwarzschild-like
  Earth spacetime. The observer has:
  - Proper acceleration $a^i = Du^i/D\tau \neq 0$ (supported against gravity)
  - Angular velocity $\Omega^i$ (co-rotating with Earth)
  - Four-velocity $u^\alpha$ tangent to the worldline
- **Key point**: Contrast with the freely-falling observers of geodesic deviation.
  The equivalence principle removes *uniform* gravity but not tidal effects — and
  the lab frame makes both visible simultaneously.
- **Recapitulate** the geodesic deviation result from Section 8 (now deleted) in
  condensed form: tidal forces = Riemann curvature contracted with separation.
  This motivates asking: what happens when we evaluate this in the lab frame?

### 9.2 The Metric in the Lab Frame

- **Pedagogical approach**: Build the metric incrementally:
  1. **Flat + acceleration**: Rindler-like metric. A test mass released at rest
     falls with acceleration $-a^i$ (the proper acceleration appears as
     "gravity" in the lab). This is the equivalence principle made quantitative.
  2. **Add rotation**: Cross-terms $g_{0i} \propto \Omega \times X$ generate
     Coriolis ($-2\mathbf{\Omega} \times \mathbf{v}$) and centrifugal
     ($-\mathbf{\Omega} \times (\mathbf{\Omega} \times \mathbf{X})$) terms in
     the equation of motion.
  3. **Add curvature**: Quadratic terms in $R_{\alpha\beta\gamma\delta} X^i X^j$
     generate tidal forces — the same Riemann curvature from geodesic deviation,
     now appearing as a correction to the lab-frame metric.
- **Present the complete Fermi metric** after the incremental buildup.
- **Note on cross-terms**: The incremental construction is pedagogical scaffolding.
  In the exact Fermi metric, acceleration, rotation, and curvature enter
  simultaneously with cross-terms. State which cross-terms are negligible and why.

### 9.3 Equation of Motion for the Test Mass

- **Content**: A freely falling test mass in the lab follows a geodesic of the full
  spacetime. In Fermi coordinates, the geodesic equation gives its coordinate
  acceleration relative to the lab:

  $$\ddot{X}^i = -a^i - 2(\mathbf{\Omega} \times \dot{\mathbf{X}})^i
  - [\mathbf{\Omega} \times (\mathbf{\Omega} \times \mathbf{X})]^i
  - c^2 R^i_{\ 0j0} X^j + \cdots$$

- **Physical identification of each term**:
  - $-a^i$: proper acceleration = "gravity" in the lab
  - $-2\Omega \times \dot{X}$: Coriolis force
  - $-\Omega \times (\Omega \times X)$: centrifugal force
  - $-c^2 R^i_{0j0} X^j$: tidal acceleration (from Riemann curvature)
  - Higher-order: gravitomagnetic curvature terms ($R^i_{0jk}$), 1PN corrections

- **What the gravimeter measures**: The housing accelerates with the lab (proper
  acceleration + centrifugal + support). The test mass is in free fall. The
  *measured* acceleration is their difference, projected onto the sensor axis $\hat{n}$.

### 9.4 Recovery of the Full g(t) Formula

- **Decompose the Riemann tensor**: $R = R_\text{Earth}(\text{static}) +
  R_\text{external}(\text{tidal})$.
  - Earth's self-gravitational curvature: contributes to $R^i_{0j0}$ but is
    time-independent. Combined with the proper acceleration and centrifugal terms,
    this produces the free-air gradient and latitude-dependent corrections already
    absorbed into $\gamma(\varphi, h)$ via Somigliana's formula (Section 3).
  - External curvature (Moon, Sun): time-varying, produces the tidal signal.

- **Proper acceleration identification**: The GR framework gives the *structure*
  $g_\text{measured} = a_\text{proper} + \text{tidal terms}$. The *value* of the
  proper acceleration is identified with Somigliana normal gravity $\gamma(\varphi, h)$
  from Section 3's geophysical model. This is a matching condition, not a GR derivation.

- **Recovery**: In the weak-field limit, the tidal curvature $c^2 R^i_{0j0} X^j$
  reduces to the Newtonian tidal tensor $T_{ij} \xi^j = \partial_i \partial_j \Phi \cdot \xi^j$,
  recovering the exact tidal formula from Section 6.

- **The elastic response factor $\delta$**: The Love number correction (Section 7)
  modifies the tidal signal through Earth deformation. This is solid-mechanics
  physics independent of the reference-frame analysis. State explicitly that
  $\delta$ enters as a multiplicative correction from Section 7, not from the
  Fermi frame derivation.

- **Final formula**:
  $$g(t) = \gamma(\varphi, h) + \delta \cdot [\mathbf{a}_\text{moon}(t) +
  \mathbf{a}_\text{sun}(t)] \cdot \hat{n}$$

### 9.5 The Complete Hierarchy of Corrections

All corrections from a single unified framework, ordered by magnitude:

| Effect | Magnitude | Source in Fermi metric | Status |
|--------|-----------|----------------------|--------|
| Static gravity $\gamma(\varphi,h)$ | ~9.8 m/s² | Proper acceleration $a^i$ | Modeled (Sec 3) |
| Centrifugal | ~3.4 Gal (equator) | $\Omega \times (\Omega \times X)$ | Time-independent, in $\gamma$ |
| Newtonian lunar tidal | ~110 $\mu$Gal | $c^2 R^i_{0j0}(\text{Moon}) X^j$ | Modeled (Sec 6) |
| Newtonian solar tidal | ~50 $\mu$Gal | $c^2 R^i_{0j0}(\text{Sun}) X^j$ | Modeled (Sec 6) |
| Coriolis (leading, horizontal) | ~15,000 $\mu$Gal | $2\Omega \times \dot{X}$ | Horizontal; no $\hat{n}$ projection |
| Second-order Coriolis (Eotvos) | ~300 nGal | $(2\omega)^2 g\, t^2$ vertical | Instrument-dependent systematic |
| GM cross-term ($v_\text{rot}$ tidal) | ~0.17 nGal | $R^i_{0jk} v^k_\text{rot} \xi^j$ | Below 10 nGal floor |
| Earth frame-dragging (Lense-Thirring) | ~0.007 nGal | $2\mathbf{v} \times \mathbf{B}_g$ (GR) | Below floor |
| 1PN tidal (Sun) | ~0.002 nGal | $4\Phi/c^2 \times T_{ij}$ | Below floor |
| 1PN tidal (Moon) | ~$6 \times 10^{-8}$ nGal | $4\Phi/c^2 \times T_{ij}$ | Below floor |
| GM tidal (spin, orbital) | $< 10^{-12}$ nGal | $R^i_{0jk}$ from $\mathbf{B}_g$ gradients | Below floor |

- **Key insight**: The coordinate gravitomagnetic field from rotation ($\omega \sim
  10^{-5}$ s$^{-1}$) is $10^9$ times larger than the Lense-Thirring field
  ($B_g \sim 10^{-14}$ s$^{-1}$). Both produce $g_{0i}$ cross-terms in the metric.
  The Fermi frame makes their relative magnitudes transparent.

- **Instrument-dependent effects**: The second-order Coriolis (Eotvos) correction
  is ~300 nGal for absolute gravimeters (free-falling retroreflector, ~0.3 m drop)
  but zero for superconducting gravimeters (levitated sphere, no macroscopic
  velocity). Note this distinction clearly.

### 9.6 The Formula from First Principles — Closing the Chain

- **Synthesis**: The complete chain from GR to the gravimeter reading:
  1. Spacetime curvature (Riemann tensor) → tidal forces (Sections 9.1-9.2)
  2. Lab frame metric (Fermi coordinates) → proper acceleration + rotation + tidal
     terms (Sections 9.2-9.3)
  3. Proper acceleration → $\gamma(\varphi, h)$ via Somigliana (Section 3)
  4. External curvature → Newtonian tidal formula (Section 6)
  5. Elastic deformation → $\delta$ factor (Section 7)
  6. Ephemerides → Moon/Sun positions (Sections 4-5)
  7. All corrections bounded below 10 nGal → Newtonian formula is sufficient (Section 9.5)

- **Coordinate matching**: The Fermi spatial coordinates of the lab are identifiable
  with ECEF coordinates to $O(\Phi_\text{Earth}/c^2) \sim 10^{-9}$, which is
  negligible. State this explicitly.

- **The GR framework's role**: Not to change the Newtonian formula, but to provide:
  (a) the reason tidal forces exist (spacetime curvature),
  (b) rigorous bounds on all corrections, and
  (c) the precise conditions under which Newton breaks down.

---

## 4. Design Decisions

### 4.1 Fermi metric: cite, don't re-derive

The Fermi metric for a non-geodesic rotating observer is a standard result in the
GR literature (Manasse & Misner 1963; Ni & Zimmerman 1978; Li & Ni 1979; Poisson
& Will 2014 ch. 9). Re-deriving it from scratch would consume 10-15 pages of tensor
algebra and is not the point of this chapter. **Cite the result, then derive physics
from it.**

### 4.2 Internal vs external curvature decomposition

The Riemann tensor in the Fermi metric includes ALL curvature — Earth's self-field
and external (Moon/Sun). Must explicitly decompose:
- $R_\text{Earth}$: static, absorbed into $\gamma(\varphi, h)$ through the
  free-air gradient
- $R_\text{external}$: time-varying, the tidal signal Pytheas models

This decomposition is justified by linearity of the Riemann tensor in the weak-field
limit and the vast timescale separation (Earth's self-field changes on geological
timescales; Moon/Sun tidal field changes on hours).

### 4.3 Proper acceleration is identified, not derived

The Fermi frame gives: "there exists a proper acceleration $a^i$."
Somigliana's formula (Section 3) gives: "$a^i = \gamma(\varphi, h)$."
These are matched, not derived from one another.

### 4.4 The $\delta$ factor comes from Section 7

The Love number correction is solid-mechanics physics (elastic deformation of
Earth's interior). It enters as a multiplicative correction applied to the tidal
signal, independent of the reference-frame analysis. Section 9 does not derive
$\delta$; it shows where $\delta$ enters the structural formula.

### 4.5 Second-order Coriolis must be discussed

The Eotvos effect on a falling test mass is ~300 nGal at mid-latitude — 30x the
accuracy floor. It is a known systematic that absolute gravimeters correct for.
Include it in the hierarchy table with a note on instrument dependence.

### 4.6 Figure 14 is retained and updated

The existing `fig14_gr_corrections.py` horizontal bar chart is updated to include
the full hierarchy from the Fermi frame analysis (centrifugal, Coriolis, Eotvos,
$v_\text{rot}$ cross-term) alongside the existing GR corrections. The accuracy
floor line remains.

---

## 5. Files Modified/Created

| File | Action |
|------|--------|
| `doc/_sec8_gr_tidal.md` | **Delete** (superseded by _sec9) |
| `doc/_sec9_gr_lab_frame.md` | **Create**: Full Chapter 9 |
| `doc/figures/fig14_gr_corrections.py` | **Edit**: Update with full hierarchy |
| `doc/figures/fig14_gr_corrections.png` | **Regenerate** |
| `doc/README.md` | **Edit**: Replace Sec 8 TOC with Sec 9, update description |

---

## 6. Implementation Plan

### Phase 1: Generate the derivation (`/alethic-derive`)

Use `/alethic-derive -p thorough` with the following problem statement:

> Derive the gravitational acceleration measured by a surface gravimeter from
> general relativity, working entirely in the proper reference frame of the lab.
>
> The lab is a non-geodesic, rotating observer on Earth's surface with proper
> acceleration a^i (supported against gravity) and angular velocity Omega^i
> (co-rotating with Earth). Use the Fermi normal coordinate metric for such an
> observer (cite Ni & Zimmerman 1978, Li & Ni 1979):
>
> g_00 = -(1 + a_i X^i / c^2)^2 + R_{0i0j} X^i X^j / c^2 + ...
> g_0i = (2/3) R_{0jik} X^j X^k + epsilon_{ijk} Omega^j X^k + ...
> g_ij = delta_ij - (1/3) R_{ikjl} X^k X^l + ...
>
> From this metric, derive the equation of motion for a freely falling test mass
> in the lab frame. Show that it contains: (1) proper acceleration (= static
> gravity), (2) Coriolis and centrifugal forces from rotation, (3) tidal forces
> from Riemann curvature of external bodies (Moon, Sun). Show that the
> Coriolis force is horizontal to leading order and does not project onto the
> vertical measurement axis. Compute the second-order Coriolis (Eotvos) effect
> and show it is ~300 nGal for a 0.3m drop — an instrument-dependent systematic.
>
> Decompose the Riemann tensor into Earth's static self-field (absorbed into
> normal gravity gamma(phi,h)) and time-varying external field (Moon/Sun tidal
> signal). In the weak-field limit, show the external tidal curvature reduces
> to the Newtonian tidal tensor T_ij = d^2 Phi / dx_i dx_j, recovering the
> exact formula a_tidal = GM[(r_M - R)/|r_M - R|^3 - r_M/|r_M|^3].
>
> Compute the gravitomagnetic cross-term from the lab's rotational velocity
> v_rot = omega R_earth cos(phi) ~ 465 m/s entering the R^i_{0jk} curvature
> components, giving ~0.17 nGal — below the 10 nGal accuracy floor. Show that
> this coordinate gravitomagnetic field (omega ~ 10^{-5} s^{-1}) is 10^9 times
> larger than the Lense-Thirring field (B_g ~ 10^{-14} s^{-1}), and that the
> Fermi frame makes this hierarchy transparent.
>
> Compute the 1PN correction to the tidal tensor and bound it at < 0.002 nGal
> (Sun) and < 10^{-7} nGal (Moon). Bound the Lense-Thirring gravitomagnetic
> tidal effects at < 0.007 nGal (direct) and < 10^{-12} nGal (tidal).
>
> Conclude with a complete hierarchy table of ALL corrections (from centrifugal
> at 3.4 Gal down to spin gravitomagnetic tidal at 10^{-14} nGal), showing that
> the Newtonian tidal formula is sufficient for surface gravimetry at the nGal
> level, and that the GR + lab-frame analysis rigorously justifies every term in
> the Pytheas formula g(t) = gamma(phi,h) + delta * [a_moon + a_sun] . n_hat.

### Phase 2: Convert to textbook format (`/alethic-textbook`)

Convert the verified derivation to textbook-quality Section 9, with:
- PhD-level presentation matching Sections 1-7 style
- Equation numbering continuing from Section 7
- References to Sections 3, 6, and 7 for the matching conditions
- Figure 14 reference for the correction hierarchy
- Code snippets showing Pytheas library usage where appropriate

Output: `doc/_sec9_gr_lab_frame.md`

### Phase 3: Update figure and cross-references

1. Edit `doc/figures/fig14_gr_corrections.py`:
   - Add centrifugal, Coriolis (horizontal), second-order Coriolis (Eotvos),
     and $v_\text{rot}$ gravitomagnetic cross-term entries
   - Distinguish categories: "absorbed into gamma" (gray), "instrument-dependent"
     (blue), "below accuracy floor" (existing colors)
   - Update axis limits and labels for the expanded range

2. Regenerate `fig14_gr_corrections.png`

3. Edit `doc/README.md`:
   - Replace Section 8 TOC entries with Section 9 entries (9.1-9.6)
   - Update intro paragraph to mention "lab frame" and Fermi coordinates
   - Keep Figure 14 entry, update description

### Phase 4: Clean up and commit

1. Delete `doc/_sec8_gr_tidal.md`
2. Verify all cross-references in Sections 1-7 that mention "Section 8" →
   update to "Section 9" (grep for "Section 8", "Sec. 8", "sec8")
3. Git commit with descriptive message

---

## 7. Estimated Effort

| Phase | Method | Est. Task calls |
|-------|--------|-----------------|
| Phase 1: `/alethic-derive -p thorough` | 8 iter, 3 best-of, 5 rev | ~30-50 |
| Phase 2: `/alethic-textbook` | Planner + writers + fidelity | ~8-12 |
| Phase 3: Figure + cross-refs | Direct edits | ~2-3 |
| Phase 4: Clean up + commit | Direct | ~1 |
| **Total** | | ~41-66 |

---

## 8. Verification Criteria

1. Section 9 renders correctly in markdown with all LaTeX
2. The complete hierarchy table includes ALL effects from centrifugal down to spin GM tidal
3. The proper acceleration ↔ Somigliana identification is stated as a matching condition
4. The $\delta$ factor is attributed to Section 7, not derived
5. The internal/external Riemann decomposition is explicit
6. The second-order Coriolis (Eotvos) is computed and flagged as instrument-dependent
7. The coordinate vs GR gravitomagnetic comparison ($10^9$ ratio) is prominent
8. Figure 14 updated with full hierarchy
9. All cross-references updated (no stale "Section 8" references)
10. No LaTeX table rendering regressions

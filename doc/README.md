# Computing $g(t)$ From First Principles

**A pedagogical derivation of the gravitational acceleration measured by a surface gravimeter**

$$g_\text{total}(t) = \gamma(\varphi, h) + \delta \times \Big[a_\text{moon}(t) + a_\text{sun}(t)\Big] \cdot \hat{n}$$

This document derives every term in the formula above from first principles, building up layer by layer from Newton's shell theorem to a complete, working tidal gravity calculator accurate to ~10 nGal at inland sites. A final section traces the Newtonian tidal model back to general relativity via the Fermi normal coordinate framework, confirming that all post-Newtonian and gravitomagnetic corrections are negligible at this accuracy level.

Each section includes worked numerical examples, code snippets from the [Pytheas](../README.md) library, and publication-quality figures.

---

## Getting Started

To run the code snippets in this document, install Pytheas and its plotting dependency:

```bash
pip install -e "/path/to/pytheas[plot]"
```

Or with conda/micromamba:

```bash
git clone https://github.com/hyperion-git/pytheas.git
cd pytheas
micromamba create -f environment.yml -y
micromamba activate pytheas
pip install -e .
```

Verify the installation:

```python
from pytheas import compute_g
from datetime import datetime

result = compute_g(datetime(2025, 3, 20, 12, 0), lat_deg=48.42, lon_deg=9.96, alt_m=620.0)
print(f"g = {result['g_total']:.10f} m/s²")  # Ulm, Eselsberg
```

**Requirements:** Python >= 3.9, NumPy >= 1.20, Matplotlib >= 3.5 (for figures).
See the [project README](../README.md) for full installation options and the CLI interface.

---

## Table of Contents

### [1. Introduction: What Does a Gravimeter Read?](_sec1_introduction.md)

- [1.1 The Question](_sec1_introduction.md#11-the-question)
- [1.2 The Formula at a Glance](_sec1_introduction.md#12-the-formula-at-a-glance)
- [1.3 A 48-Hour Portrait](_sec1_introduction.md#13-a-48-hour-portrait)
- [1.4 Roadmap](_sec1_introduction.md#14-roadmap)
- [1.5 The Error Budget at a Glance](_sec1_introduction.md#15-the-error-budget-at-a-glance)
- [1.6 Notation and Conventions](_sec1_introduction.md#16-notation-and-conventions)

### [2. Setting Up Coordinates](_sec2_coordinates.md)

- [2.1 Geodetic Coordinates: Latitude, Longitude, Altitude](_sec2_coordinates.md#21-geodetic-coordinates-latitude-longitude-altitude)
- [2.2 The GRS80 Reference Ellipsoid](_sec2_coordinates.md#22-the-grs80-reference-ellipsoid)
- [2.3 Geodetic vs. Geocentric Latitude](_sec2_coordinates.md#23-geodetic-vs-geocentric-latitude)
- [2.4 Geodetic to ECEF Conversion](_sec2_coordinates.md#24-geodetic-to-ecef-conversion)
- [2.5 The ECI Frame (Earth-Centered Inertial)](_sec2_coordinates.md#25-the-eci-frame-earth-centered-inertial)
- [2.6 The ECEF Frame (Earth-Centered Earth-Fixed)](_sec2_coordinates.md#26-the-ecef-frame-earth-centered-earth-fixed)
- [2.7 GMST: Greenwich Mean Sidereal Time](_sec2_coordinates.md#27-gmst-greenwich-mean-sidereal-time)
- [2.8 ECI to ECEF Rotation](_sec2_coordinates.md#28-eci-to-ecef-rotation)
- [2.9 ENU Basis Vectors (East-North-Up)](_sec2_coordinates.md#29-enu-basis-vectors-east-north-up)
- [2.10 The Measurement Axis](_sec2_coordinates.md#210-the-measurement-axis)
- [2.11 Summary: The Coordinate Pipeline](_sec2_coordinates.md#211-summary-the-coordinate-pipeline)

### [3. Gravity on a Quiet Earth (Normal Gravity)](_sec3_normal_gravity.md)

- [3.1 The Non-Rotating Spherical Case](_sec3_normal_gravity.md#31-the-non-rotating-spherical-case)
- [3.2 Why Rotation Matters](_sec3_normal_gravity.md#32-why-rotation-matters)
- [3.3 Why Shape Matters](_sec3_normal_gravity.md#33-why-shape-matters)
- [3.4 Somigliana's Formula](_sec3_normal_gravity.md#34-somiglianas-formula)
- [3.5 Free-Air Correction (Going Up)](_sec3_normal_gravity.md#35-free-air-correction-going-up)

### [4. Where is the Moon? (Lunar Ephemeris)](_sec4_lunar_ephemeris.md)

- [4.1 Why We Need to Know](_sec4_lunar_ephemeris.md#41-why-we-need-to-know)
- [4.2 Keplerian Orbit Review](_sec4_lunar_ephemeris.md#42-keplerian-orbit-review)
- [4.3 Why the Moon Is Hard](_sec4_lunar_ephemeris.md#43-why-the-moon-is-hard)
- [4.4 The Meeus Approach](_sec4_lunar_ephemeris.md#44-the-meeus-approach)
- [4.5 From Ecliptic to Your Location](_sec4_lunar_ephemeris.md#45-from-ecliptic-to-your-location)

### [5. Where is the Sun? (Solar Ephemeris)](_sec5_solar_ephemeris.md)

- [5.1 Why the Sun is Easier](_sec5_solar_ephemeris.md#51-why-the-sun-is-easier)
- [5.2 Mean Anomaly and the Equation of Center](_sec5_solar_ephemeris.md#52-mean-anomaly-and-the-equation-of-center)
- [5.3 Distance to the Sun](_sec5_solar_ephemeris.md#53-distance-to-the-sun)
- [5.4 From Ecliptic to ECEF](_sec5_solar_ephemeris.md#54-from-ecliptic-to-ecef)
- [5.5 Why Low Precision Suffices](_sec5_solar_ephemeris.md#55-why-low-precision-suffices)

### [6. The Tidal Acceleration](_sec6_tidal_acceleration.md)

- [6.1 Tidal Forces From First Principles](_sec6_tidal_acceleration.md#61-tidal-forces-from-first-principles)
- [6.2 How Big Is This?](_sec6_tidal_acceleration.md#62-how-big-is-this)
- [6.3 The Gradient Approximation (and Why Pytheas Doesn't Use It)](_sec6_tidal_acceleration.md#63-the-gradient-approximation-and-why-pytheas-doesnt-use-it)
- [6.4 Projecting Onto Your Sensor](_sec6_tidal_acceleration.md#64-projecting-onto-your-sensor)

### [7. The Earth Fights Back (Elastic Response)](_sec7_elastic_response.md)

- [7.1 The Earth is Not Rigid](_sec7_elastic_response.md#71-the-earth-is-not-rigid)
- [7.2 How Deformation Changes What You Measure](_sec7_elastic_response.md#72-how-deformation-changes-what-you-measure)
- [7.3 Love Numbers: Parameterizing Deformation](_sec7_elastic_response.md#73-love-numbers-parameterizing-deformation)
- [7.4 Deriving the Gravimetric Factor $\delta$](_sec7_elastic_response.md#74-deriving-the-gravimetric-factor-delta)
- [7.5 Numerical Values](_sec7_elastic_response.md#75-numerical-values)
- [7.6 What We Are Ignoring](_sec7_elastic_response.md#76-what-we-are-ignoring)

### [8. General Relativity and the Tidal Model](_sec8_gr_tidal.md)

- [8.1 The Fermi Normal Coordinate Framework](_sec8_gr_tidal.md#81-the-fermi-normal-coordinate-framework)
- [8.2 The Geodesic Equation in the Lab Frame](_sec8_gr_tidal.md#82-the-geodesic-equation-in-the-lab-frame)
- [8.3 Term-by-Term Derivation and Classification](_sec8_gr_tidal.md#83-term-by-term-derivation-and-classification)
- [8.4 Recovering the Newtonian Tidal Formula](_sec8_gr_tidal.md#84-recovering-the-newtonian-tidal-formula)
- [8.5 Post-Newtonian and Gravitomagnetic Corrections](_sec8_gr_tidal.md#85-post-newtonian-and-gravitomagnetic-corrections)
- [8.6 The Complete Hierarchy](_sec8_gr_tidal.md#86-the-complete-hierarchy)
- [8.7 Dimensional and Limiting-Case Checks](_sec8_gr_tidal.md#87-dimensional-and-limiting-case-checks)
- [8.8 Conclusion](_sec8_gr_tidal.md#88-conclusion)
- [8.9 Computational Verification](_sec8_gr_tidal.md#89-computational-verification)

---

## Figures

| # | Description | File |
|---|-------------|------|
| 1 | 48-hour $g(t)$ timeseries | [fig01](figures/fig01_g_timeseries.png) |
| 2 | ECI, ECEF, and ENU reference frames | [fig02](figures/fig02_reference_frames.png) |
| 3 | ENU basis vectors and measurement axis $\hat{n}$ | [fig03](figures/fig03_measurement_axis.png) |
| 4 | Normal gravity vs. latitude (sphere, rotation, Somigliana) | [fig04](figures/fig04_gravity_vs_latitude.png) |
| 5 | Free-air correction: 1st vs. 2nd order | [fig05](figures/fig05_free_air.png) |
| 6 | Lunar distance and tidal amplitude correlation | [fig06](figures/fig06_lunar_distance_tides.png) |
| 7 | Perturbation series convergence (top terms) | [fig07](figures/fig07_perturbation_terms.png) |
| 8 | Moon vs. Sun tidal contributions | [fig08](figures/fig08_moon_vs_sun.png) |
| 9 | Tidal geometry and quadrupolar pattern | [fig09](figures/fig09_tidal_schematic.png) |
| 10 | Exact vs. gradient-approximation tidal acceleration | [fig10](figures/fig10_exact_vs_approx.png) |
| 11 | Rigid vs. elastic Earth tidal signal | [fig11](figures/fig11_rigid_vs_elastic.png) |
| 12 | Anatomy of $g(t)$: static + lunar + solar | [fig12](figures/fig12_anatomy.png) |
| 13 | Error budget | [fig13](figures/fig13_error_budget.png) |
| 14 | Complete hierarchy of corrections from Fermi frame analysis | [fig14](figures/fig14_gr_corrections.png) |

# Computing $g(t)$ From First Principles

**A pedagogical derivation of the gravitational acceleration measured by a surface gravimeter**

$$g_\text{total}(t) = \gamma(\varphi, h) + \delta \times \Big[\mathbf{a}_\text{Moon}(t) + \mathbf{a}_\text{Sun}(t)\Big] \cdot \hat{\mathbf{n}}$$

This document derives every term in the formula above from first principles, building layer by layer from Newton's shell theorem to a tidal gravity calculator accurate to ~10 nGal at inland sites. A final section traces the Newtonian tidal model back to general relativity via the Fermi normal coordinate framework, confirming that post-Newtonian and gravitomagnetic corrections remain negligible at this accuracy level.

Each section includes worked numerical examples, code snippets using the [Pytheas](../README.md) library, and publication-quality figures.

---

## Getting Started

To run the code snippets, install Pytheas with its plotting dependency:

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

Verify the installation with a quick test:

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

### [1. Introduction](_sec1_introduction.md)

- [1.1 The Question](_sec1_introduction.md#11-the-question)
- [1.2 The Formula at a Glance](_sec1_introduction.md#12-the-formula-at-a-glance)
- [1.3 A 48-Hour Portrait](_sec1_introduction.md#13-a-48-hour-portrait)
- [1.4 Roadmap](_sec1_introduction.md#14-roadmap)
- [1.5 The Error Budget at a Glance](_sec1_introduction.md#15-the-error-budget-at-a-glance)
- [1.6 Notation and Conventions](_sec1_introduction.md#16-notation-and-conventions)

### [2. Coordinate Systems](_sec2_coordinates.md)

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

### [3. Normal Gravity](_sec3_normal_gravity.md)

- [3.1 The Non-Rotating Spherical Case](_sec3_normal_gravity.md#31-the-non-rotating-spherical-case)
- [3.2 Effect of Earth's Rotation](_sec3_normal_gravity.md#32-effect-of-earths-rotation)
- [3.3 Effect of Earth's Oblateness](_sec3_normal_gravity.md#33-effect-of-earths-oblateness)
- [3.4 Somigliana's Formula](_sec3_normal_gravity.md#34-somiglianas-formula)
- [3.5 Free-Air Correction](_sec3_normal_gravity.md#35-free-air-correction)

### [4. Lunar Ephemeris](_sec4_lunar_ephemeris.md)

- [4.1 Sensitivity to Lunar Position](_sec4_lunar_ephemeris.md#41-sensitivity-to-lunar-position)
- [4.2 Keplerian Orbit Review](_sec4_lunar_ephemeris.md#42-keplerian-orbit-review)
- [4.3 Perturbations and the Three-Body Problem](_sec4_lunar_ephemeris.md#43-perturbations-and-the-three-body-problem)
- [4.4 The Meeus Approach](_sec4_lunar_ephemeris.md#44-the-meeus-approach)
- [4.5 Ecliptic-to-ECEF Conversion](_sec4_lunar_ephemeris.md#45-ecliptic-to-ecef-conversion)

### [5. Solar Ephemeris](_sec5_solar_ephemeris.md)

- [5.1 Why a Simple Analytical Model Suffices](_sec5_solar_ephemeris.md#51-why-a-simple-analytical-model-suffices)
- [5.2 Mean Anomaly and the Equation of Center](_sec5_solar_ephemeris.md#52-mean-anomaly-and-the-equation-of-center)
- [5.3 Distance to the Sun](_sec5_solar_ephemeris.md#53-distance-to-the-sun)
- [5.4 From Ecliptic to ECEF](_sec5_solar_ephemeris.md#54-from-ecliptic-to-ecef)
- [5.5 Accuracy of the Solar Ephemeris](_sec5_solar_ephemeris.md#55-accuracy-of-the-solar-ephemeris)

### [6. Tidal Acceleration](_sec6_tidal_acceleration.md)

- [6.1 Derivation from First Principles](_sec6_tidal_acceleration.md#61-derivation-from-first-principles)
- [6.2 Magnitude Estimates](_sec6_tidal_acceleration.md#62-magnitude-estimates)
- [6.3 The Gradient Approximation](_sec6_tidal_acceleration.md#63-the-gradient-approximation)
- [6.4 Projection onto the Sensor Axis](_sec6_tidal_acceleration.md#64-projection-onto-the-sensor-axis)

### [7. Elastic Response of the Solid Earth](_sec7_elastic_response.md)

- [7.1 Body Tides and Elastic Deformation](_sec7_elastic_response.md#71-body-tides-and-elastic-deformation)
- [7.2 Two Competing Effects on Measured Gravity](_sec7_elastic_response.md#72-two-competing-effects-on-measured-gravity)
- [7.3 Love Numbers](_sec7_elastic_response.md#73-love-numbers)
- [7.4 Deriving the Gravimetric Factor $\delta$](_sec7_elastic_response.md#74-deriving-the-gravimetric-factor-delta)
- [7.5 Numerical Values](_sec7_elastic_response.md#75-numerical-values)
- [7.6 Omitted Effects and Their Magnitudes](_sec7_elastic_response.md#76-omitted-effects-and-their-magnitudes)

### [8. General Relativity and the Lab-Frame Gravimeter Equation](_sec8_gr_tidal.md)

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
| 1 | 48-hour $g(t)$ time series | [fig01](figures/fig01_g_timeseries.png) |
| 2 | ECI, ECEF, and ENU reference frames | [fig02](figures/fig02_reference_frames.png) |
| 3 | ENU basis vectors and measurement axis $\hat{n}$ | [fig03](figures/fig03_measurement_axis.png) |
| 4 | Normal gravity vs. latitude: sphere, rotation, Somigliana | [fig04](figures/fig04_gravity_vs_latitude.png) |
| 5 | Free-air correction: first- vs. second-order | [fig05](figures/fig05_free_air.png) |
| 6 | Lunar distance and tidal amplitude correlation | [fig06](figures/fig06_lunar_distance_tides.png) |
| 7 | Perturbation series convergence (dominant terms) | [fig07](figures/fig07_perturbation_terms.png) |
| 8 | Lunar vs. solar tidal contributions | [fig08](figures/fig08_moon_vs_sun.png) |
| 9 | Tidal geometry and quadrupolar pattern | [fig09](figures/fig09_tidal_schematic.png) |
| 10 | Exact vs. gradient-approximation tidal acceleration | [fig10](figures/fig10_exact_vs_approx.png) |
| 11 | Rigid vs. elastic Earth tidal signal | [fig11](figures/fig11_rigid_vs_elastic.png) |
| 12 | Anatomy of $g(t)$: static, lunar, and solar components | [fig12](figures/fig12_anatomy.png) |
| 13 | Error budget summary | [fig13](figures/fig13_error_budget.png) |
| 14 | Hierarchy of corrections from Fermi frame analysis | [fig14](figures/fig14_gr_corrections.png) |

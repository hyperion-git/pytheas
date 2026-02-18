"""
Sensor orientation sweep: all acceleration terms vs orientation.

Comprehensive 2×2 figure showing how each physical contribution varies
with sensor zenith angle θ (left) and azimuth φ (right), after
subtracting the dominant effective-gravity projection g₀ cos θ.

Top row (10⁻³ m/s²): centrifugal correction and claimed direct-pull
signals from Sun and Moon, compared to the combined tidal residual
(indistinguishable from zero). Baseline is g_grav cos θ (not g₀ cos θ).

Bottom row (10⁻⁷ m/s²): actual tidal residuals decomposed into solar
and lunar contributions — the only time-varying signals in the sweep.
"""

import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

# --- Style: publication, single-column PRX ---
plt.rcParams.update({
    "figure.figsize": (6.5, 5.0),
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
    "font.family": "sans-serif",
    "font.sans-serif": [
        "DejaVu Sans", "Arial", "Helvetica",
        "Lucida Grande", "Verdana", "sans-serif",
    ],
    "font.size": 8,
    "axes.labelsize": 9,
    "axes.titlesize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 7,
    "text.usetex": False,
    "axes.formatter.use_mathtext": True,
    "mathtext.fontset": "dejavusans",
    "axes.linewidth": 0.5,
    "xtick.direction": "in",
    "xtick.major.size": 3,
    "xtick.major.width": 0.5,
    "xtick.minor.size": 1.5,
    "xtick.minor.width": 0.5,
    "xtick.minor.visible": True,
    "xtick.top": True,
    "ytick.direction": "in",
    "ytick.major.size": 3,
    "ytick.major.width": 0.5,
    "ytick.minor.size": 1.5,
    "ytick.minor.width": 0.5,
    "ytick.minor.visible": True,
    "ytick.right": True,
    "lines.linewidth": 1.0,
    "lines.markersize": 3,
    "axes.grid": False,
    "legend.frameon": False,
    "axes.prop_cycle": plt.cycler("color", [
        "#2e2cb8",  # blue
        "#db002b",  # red
        "#1f8a70",  # teal
        "#fd7400",  # orange
        "#1c809e",  # blueberry
        "#bedb43",  # lime
        "#9b2f5c",  # pink
        "#fabd1e",  # gold
    ]),
})

# --- Physical parameters ---
g0 = 9.81  # effective gravity (m/s²)

# Solar tidal components in geodetic ENU (representative snapshot:
# Sun at elevation ~40°, azimuth ~150° SSE)
a_E_sol = 1.9e-7    # East
a_N_sol = -1.3e-7   # North
a_U_sol = 0.6e-7    # Up

# Lunar tidal components (representative: Moon at elevation ~55°,
# azimuth ~250° WSW)
a_E_lun = -7.3e-7
a_N_lun = -2.6e-7
a_U_lun = 5.6e-7

# Combined tidal
a_E = a_E_sol + a_E_lun
a_N = a_N_sol + a_N_lun
a_U = a_U_sol + a_U_lun

# Centrifugal correction (at Ulm, Germany: λ = 48.4°)
omega = 7.292e-5   # Earth angular velocity (rad/s)
r_earth = 6.371e6  # Earth radius (m)
lat = np.radians(48.4)
a_cf = omega**2 * r_earth * np.cos(lat)**2  # ≈ 1.5e-2 m/s²

# Direct pulls (canceled in reality)
a_direct_sun = 5.93e-3    # m/s²
a_direct_moon = 3.32e-5   # m/s²
phi_sun = np.radians(30)   # Sun azimuth from East (illustrative)
phi_moon = np.radians(200) # Moon azimuth from East

# --- Sweep arrays ---
theta_deg = np.linspace(0, 180, 500)
theta = np.radians(theta_deg)

phi_deg = np.linspace(0, 360, 500)
phi = np.radians(phi_deg)

# --- θ sweep at fixed φ = 0° (East) ---
phi_fixed = 0.0

# Tidal residuals: g(θ) - g₀ cos θ = -aE sinθ cosφ - aN sinθ sinφ - aU cosθ
tidal_sol_theta = (-a_E_sol * np.sin(theta) * np.cos(phi_fixed)
                   - a_N_sol * np.sin(theta) * np.sin(phi_fixed)
                   - a_U_sol * np.cos(theta))
tidal_lun_theta = (-a_E_lun * np.sin(theta) * np.cos(phi_fixed)
                   - a_N_lun * np.sin(theta) * np.sin(phi_fixed)
                   - a_U_lun * np.cos(theta))
tidal_comb_theta = tidal_sol_theta + tidal_lun_theta

# Direct pulls (if not canceled)
claimed_sun_theta = -a_direct_sun * np.sin(theta) * np.cos(phi_fixed - phi_sun)
claimed_moon_theta = -a_direct_moon * np.sin(theta) * np.cos(phi_fixed - phi_moon)

# Centrifugal correction (separated from g₀, purely cos θ)
cf_theta = -a_cf * np.cos(theta)

# --- φ sweep at fixed θ = 45° ---
theta_fixed = np.radians(45)

tidal_sol_phi = (-a_E_sol * np.sin(theta_fixed) * np.cos(phi)
                 - a_N_sol * np.sin(theta_fixed) * np.sin(phi)
                 - a_U_sol * np.cos(theta_fixed))
tidal_lun_phi = (-a_E_lun * np.sin(theta_fixed) * np.cos(phi)
                 - a_N_lun * np.sin(theta_fixed) * np.sin(phi)
                 - a_U_lun * np.cos(theta_fixed))
tidal_comb_phi = tidal_sol_phi + tidal_lun_phi

claimed_sun_phi = -a_direct_sun * np.sin(theta_fixed) * np.cos(phi - phi_sun)
claimed_moon_phi = -a_direct_moon * np.sin(theta_fixed) * np.cos(phi - phi_moon)

# Centrifugal at fixed θ = 45° (constant vs φ)
cf_phi = -a_cf * np.cos(theta_fixed) * np.ones_like(phi)

# --- Colors ---
c_sun = "#db002b"      # red — Sun direct pull
c_moon = "#fd7400"     # orange — Moon direct pull
c_comb = "#2e2cb8"     # blue — combined tidal
c_sol = "#1f8a70"      # teal — solar tidal
c_lun = "#9b2f5c"      # magenta — lunar tidal
c_cf = "#1c809e"       # blueberry — centrifugal

# --- Figure: 2×2 grid ---
fig, axes = plt.subplots(2, 2, figsize=(6.5, 5.0))
(ax_a, ax_b), (ax_c, ax_d) = axes

# ============================================================
# Top row: 10⁻³ m/s² scale — direct pulls vs combined tidal
# ============================================================

for ax, xdeg, x_sun, x_moon, x_tidal, x_cf, is_theta in [
    (ax_a, theta_deg, claimed_sun_theta, claimed_moon_theta,
     tidal_comb_theta, cf_theta, True),
    (ax_b, phi_deg, claimed_sun_phi, claimed_moon_phi,
     tidal_comb_phi, cf_phi, False),
]:
    # Centrifugal correction (real, measurable)
    ax.plot(xdeg, x_cf * 1e3, color=c_cf, linewidth=1.0)
    # Sun direct pull (dominant claimed signal)
    ax.plot(xdeg, x_sun * 1e3, color=c_sun, linewidth=0.8,
            linestyle="--")
    # Moon direct pull (smaller but still >> tidal)
    ax.plot(xdeg, x_moon * 1e3, color=c_moon, linewidth=0.8,
            linestyle="--")
    # Combined tidal (indistinguishable from zero at this scale)
    ax.plot(xdeg, x_tidal * 1e3, color=c_comb, linewidth=1.2)
    ax.axhline(0, color="#999999", linewidth=0.3, zorder=0)

    if is_theta:
        ax.set_xlabel(r"Zenith angle $\theta$ (deg)")
        ax.set_xlim(0, 180)
        ax.set_xticks([0, 45, 90, 135, 180])
    else:
        ax.set_xlabel(r"Azimuth $\varphi$ (deg)")
        ax.set_xlim(0, 360)
        ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_ylabel(r"$g - g_\mathrm{grav}\cos\theta$  ($10^{-3}$ m/s$^2$)")

# Top-row y-limits (expanded to fit centrifugal ~ ±17)
ymax_top = 20
for ax in (ax_a, ax_b):
    ax.set_ylim(-ymax_top, ymax_top)
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.xaxis.set_minor_locator(AutoMinorLocator())

# Panel labels and sweep annotations
ax_a.text(0.03, 0.95, "(a)", transform=ax_a.transAxes,
          fontsize=9, fontweight="bold", va="top")
ax_a.text(0.03, 0.82, r"$\varphi = 0\degree$ (East)",
          transform=ax_a.transAxes, fontsize=7, color="#555555")

ax_b.text(0.03, 0.95, "(b)", transform=ax_b.transAxes,
          fontsize=9, fontweight="bold", va="top")
ax_b.text(0.03, 0.82, r"$\theta = 45\degree$",
          transform=ax_b.transAxes, fontsize=7, color="#555555")

# Direct annotations on top-row curves
# Panel (a): θ sweep
ax_a.text(65, -15.0, "Centrifugal", fontsize=6.5,
          color=c_cf, ha="center")
ax_a.text(95, -10.5, "Sun direct\n(canceled)", fontsize=6.5,
          color=c_sun, ha="center", style="italic")
ax_a.text(35, 2.5, "Combined tidal", fontsize=6, color=c_comb,
          ha="center")
# Lunar direct annotation with arrow (curve is tiny at this scale)
moon_peak_a = claimed_moon_theta.min() * 1e3  # negative peak
ax_a.annotate("Lunar direct (canceled)", xy=(90, moon_peak_a),
              xytext=(150, 6.0), fontsize=6, color=c_moon,
              style="italic", ha="center",
              arrowprops=dict(arrowstyle="->", color=c_moon,
                              linewidth=0.5))

# Panel (b): φ sweep
cf_val_b = cf_phi[0] * 1e3  # constant value (~-12)
ax_b.text(50, cf_val_b - 3.0, "Centrifugal", fontsize=6.5,
          color=c_cf, ha="center")
ax_b.text(180, -7.5, "Sun direct\n(canceled)", fontsize=6.5,
          color=c_sun, ha="center", style="italic")
ax_b.text(50, 2.5, "Combined tidal", fontsize=6, color=c_comb,
          ha="center")
moon_peak_b = claimed_moon_phi.min() * 1e3
ax_b.annotate("Lunar direct (canceled)", xy=(200, moon_peak_b),
              xytext=(300, 6.0), fontsize=6, color=c_moon,
              style="italic", ha="center",
              arrowprops=dict(arrowstyle="->", color=c_moon,
                              linewidth=0.5))

# ============================================================
# Bottom row: 10⁻⁷ m/s² scale — individual tidal components
# ============================================================

for ax, xdeg, x_sol, x_lun, x_comb, is_theta in [
    (ax_c, theta_deg, tidal_sol_theta, tidal_lun_theta,
     tidal_comb_theta, True),
    (ax_d, phi_deg, tidal_sol_phi, tidal_lun_phi,
     tidal_comb_phi, False),
]:
    ax.plot(xdeg, x_sol * 1e7, color=c_sol, linewidth=1.0,
            label="Solar tidal")
    ax.plot(xdeg, x_lun * 1e7, color=c_lun, linewidth=1.0,
            label="Lunar tidal")
    ax.plot(xdeg, x_comb * 1e7, color=c_comb, linewidth=1.2,
            label="Combined")
    ax.axhline(0, color="#999999", linewidth=0.3, zorder=0)

    if is_theta:
        ax.set_xlabel(r"Zenith angle $\theta$ (deg)")
        ax.set_xlim(0, 180)
        ax.set_xticks([0, 45, 90, 135, 180])
    else:
        ax.set_xlabel(r"Azimuth $\varphi$ (deg)")
        ax.set_xlim(0, 360)
        ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_ylabel(r"$g - g_0\cos\theta$  ($10^{-7}$ m/s$^2$)")

# Bottom-row y-limits
ymax_bot = 12
for ax in (ax_c, ax_d):
    ax.set_ylim(-ymax_bot, ymax_bot)
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.xaxis.set_minor_locator(AutoMinorLocator())

# Panel labels
ax_c.text(0.03, 0.95, "(c)", transform=ax_c.transAxes,
          fontsize=9, fontweight="bold", va="top")
ax_c.text(0.03, 0.82, r"$\varphi = 0\degree$ (East)",
          transform=ax_c.transAxes, fontsize=7, color="#555555")

ax_d.text(0.03, 0.95, "(d)", transform=ax_d.transAxes,
          fontsize=9, fontweight="bold", va="top")
ax_d.text(0.03, 0.82, r"$\theta = 45\degree$",
          transform=ax_d.transAxes, fontsize=7, color="#555555")

# Direct annotations on bottom curves
# θ sweep: place labels where curves are well-separated
ax_c.text(90, -3.5, "Solar", fontsize=6.5, color=c_sol, ha="center")
ax_c.text(80, 8.5, "Lunar", fontsize=6.5, color=c_lun, ha="center")
ax_c.text(25, -9.0, "Combined", fontsize=6.5, color=c_comb, ha="center")

# φ sweep: place labels where curves are well-separated
ax_d.text(130, 1.5, "Solar", fontsize=6.5, color=c_sol, ha="center")
ax_d.text(100, 5.5, "Lunar", fontsize=6.5, color=c_lun, ha="center")
ax_d.text(55, 2.5, "Combined", fontsize=6.5, color=c_comb, ha="center")

# Off-scale annotations
offscale_text = ("Centrifugal: $10^{5}{\\times}$ off scale\n"
                 "Sun direct: $10^{4}{\\times}$ off scale\n"
                 "Moon direct: $30{\\times}$ off scale")
ax_c.text(0.97, 0.05, offscale_text,
          transform=ax_c.transAxes, fontsize=5.5, color="#888888",
          ha="right", va="bottom", linespacing=1.4)
ax_d.text(0.97, 0.95, offscale_text,
          transform=ax_d.transAxes, fontsize=5.5, color="#888888",
          ha="right", va="top", linespacing=1.4)

fig.tight_layout(h_pad=2.0, w_pad=2.5)

# Save
fig.savefig("/home/xeal/dev/pytheas/doc/figures/orientation_sweep.pdf")
fig.savefig("/home/xeal/dev/pytheas/doc/figures/orientation_sweep.png")
print("Saved orientation_sweep.pdf and .png")

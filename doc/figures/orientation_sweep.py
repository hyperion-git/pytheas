"""
Sensor orientation sweep: tidal residuals vs claimed direct-pull signal.

Shows how the measured acceleration varies when sweeping the sensor
zenith angle θ (left) and azimuth φ (right), after subtracting the
dominant effective-gravity projection g₀ cos θ.  The hypothetical
"direct solar pull" signal (6×10⁻³ m/s², dashed) is shown for
comparison — it is absent from the actual measurement, which contains
only tidal terms (~10⁻⁷ m/s²).
"""

import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

# --- Style: publication, single-column PRX ---
plt.rcParams.update({
    "figure.figsize": (5.0, 2.5),
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
g0 = 9.81          # effective gravity (m/s²)

# Tidal components (representative solar + lunar, geodetic ENU)
a_E = 3.0e-7       # East tidal (m/s²)
a_N = 2.0e-7       # North tidal (m/s²)
a_U = 5.0e-7       # Up tidal (m/s²)

# Claimed (but canceled) direct solar pull
a_direct_sun = 5.93e-3  # m/s²
phi_sun = np.radians(30)  # Sun azimuth (arbitrary, for illustration)

# --- Sweep arrays ---
theta_deg = np.linspace(0, 180, 500)
theta = np.radians(theta_deg)

phi_deg = np.linspace(0, 360, 500)
phi = np.radians(phi_deg)

# --- Panel (a): Sweep θ at fixed φ = 0° (East) ---
phi_fixed = 0.0  # radians

# Tidal residual: g(θ) - g₀ cos θ
tidal_theta = (-a_E * np.sin(theta) * np.cos(phi_fixed)
               - a_N * np.sin(theta) * np.sin(phi_fixed)
               - a_U * np.cos(theta))

# Hypothetical direct-pull signal (if not canceled)
claimed_theta = -a_direct_sun * np.sin(theta) * np.cos(phi_fixed - phi_sun)

# --- Panel (b): Sweep φ at fixed θ = 45° ---
theta_fixed = np.radians(45)

# Tidal residual: g(φ) - g₀ cos θ_fixed
tidal_phi = (-a_E * np.sin(theta_fixed) * np.cos(phi)
             - a_N * np.sin(theta_fixed) * np.sin(phi)
             - a_U * np.cos(theta_fixed))

# Hypothetical direct-pull signal
claimed_phi = -a_direct_sun * np.sin(theta_fixed) * np.cos(phi - phi_sun)

# --- Figure ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.5, 2.8))

# Panel (a): θ sweep
ax1.fill_between(theta_deg, claimed_theta * 1e3, alpha=0.08, color="#db002b",
                 linewidth=0)
ax1.plot(theta_deg, claimed_theta * 1e3, color="#db002b", linewidth=0.8,
         linestyle="--", label="Claimed direct pull")
ax1.plot(theta_deg, tidal_theta * 1e3, color="#2e2cb8", linewidth=1.2,
         label="Actual tidal residual")
ax1.axhline(0, color="#999999", linewidth=0.3, zorder=0)

ax1.set_xlabel(r"Zenith angle $\theta$ (deg)")
ax1.set_ylabel(r"$g(\theta) - g_0\cos\theta$  ($10^{-3}$ m/s$^2$)")
ax1.set_xlim(0, 180)
ax1.set_xticks([0, 45, 90, 135, 180])
ax1.text(0.03, 0.95, "(a)", transform=ax1.transAxes, fontsize=9,
         fontweight="bold", va="top")
ax1.text(0.03, 0.82, r"$\varphi = 0\degree$ (East)",
         transform=ax1.transAxes, fontsize=7, color="#555555")

# Direct annotation
ax1.text(95, -4.5, "Claimed\ndirect pull", fontsize=6.5, color="#db002b",
         ha="center", style="italic")
ax1.text(95, 0.15, "Tidal residual", fontsize=6.5, color="#2e2cb8",
         ha="center")

# Panel (b): φ sweep
ax2.fill_between(phi_deg, claimed_phi * 1e3, alpha=0.08, color="#db002b",
                 linewidth=0)
ax2.plot(phi_deg, claimed_phi * 1e3, color="#db002b", linewidth=0.8,
         linestyle="--", label="Claimed direct pull")
ax2.plot(phi_deg, tidal_phi * 1e3, color="#2e2cb8", linewidth=1.2,
         label="Actual tidal residual")
ax2.axhline(0, color="#999999", linewidth=0.3, zorder=0)

ax2.set_xlabel(r"Azimuth $\varphi$ (deg)")
ax2.set_ylabel(r"$g(\varphi) - g_0\cos\theta$  ($10^{-3}$ m/s$^2$)")
ax2.set_xlim(0, 360)
ax2.set_xticks([0, 90, 180, 270, 360])
ax2.text(0.03, 0.95, "(b)", transform=ax2.transAxes, fontsize=9,
         fontweight="bold", va="top")
ax2.text(0.03, 0.82, r"$\theta = 45\degree$",
         transform=ax2.transAxes, fontsize=7, color="#555555")

# Direct annotation
ax2.text(210, -3.3, "Claimed\ndirect pull", fontsize=6.5, color="#db002b",
         ha="center", style="italic")
ax2.text(210, 0.15, "Tidal residual", fontsize=6.5, color="#2e2cb8",
         ha="center")

# Match y-axes
ymax = max(abs(ax1.get_ylim()[0]), abs(ax1.get_ylim()[1]),
           abs(ax2.get_ylim()[0]), abs(ax2.get_ylim()[1]))
ymax = 6.5
for ax in (ax1, ax2):
    ax.set_ylim(-ymax, ymax)
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.xaxis.set_minor_locator(AutoMinorLocator())

fig.tight_layout(w_pad=2.5)

# Save
fig.savefig("/home/xeal/dev/pytheas/doc/figures/orientation_sweep.pdf")
fig.savefig("/home/xeal/dev/pytheas/doc/figures/orientation_sweep.png")
print("Saved orientation_sweep.pdf and .png")

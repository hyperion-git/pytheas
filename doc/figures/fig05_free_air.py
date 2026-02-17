#!/usr/bin/env python3
"""Figure 5: Free-air correction.

Two panels showing how gravity decreases with altitude:
  (a) Gravity reduction vs altitude (0-10 km) in mGal.
  (b) Difference between second-order and first-order corrections in uGal.
"""

import sys
sys.path.insert(0, "/home/xeal/dev/pytheas/doc/figures")

import numpy as np
import matplotlib.pyplot as plt

import pytheas
from pytheas import A_GRS80, F_GRS80, OMEGA, GM_E, B_GRS80
from style import apply_style, COLORS, add_dual_axis

apply_style()

# ---- Parameters ----
lat = 45.0  # degrees
gamma_0 = pytheas.normal_gravity(lat, 0.0)

# Somigliana parameters for free-air factor
phi = np.radians(lat)
sin2 = np.sin(phi)**2
M_RATIO = OMEGA**2 * A_GRS80**2 * B_GRS80 / GM_E
fac = 1.0 + F_GRS80 + M_RATIO - 2.0 * F_GRS80 * sin2

h = np.linspace(0, 10000, 1000)  # 0 to 10 km

# First-order free-air: Delta_gamma = gamma_0 * (-2 * fac * h / a)
delta_first = gamma_0 * (-2.0 * fac * h / A_GRS80)

# Full second-order: Delta_gamma = gamma_0 * (-2*fac*h/a + 3*h^2/a^2)
delta_second = gamma_0 * (-2.0 * fac * h / A_GRS80 + 3.0 * h**2 / A_GRS80**2)

# Verify against pytheas directly
delta_full = np.array([pytheas.normal_gravity(lat, hi) - gamma_0 for hi in h])

# Convert to mGal (1 mGal = 1e-5 m/s^2)
to_mGal = 1e5
to_uGal = 1e8  # 1 uGal = 1e-8 m/s^2

# ---- Create figure ----
fig, (ax_main, ax_inset) = plt.subplots(
    1, 2, figsize=(6.5, 3.2),
    gridspec_kw={"width_ratios": [2, 1], "wspace": 0.35},
)

# --- Panel (a): main gravity reduction ---
ax_main.plot(h / 1000, delta_first * to_mGal,
             color=COLORS["orange"], linestyle="--", linewidth=1.3,
             label="First-order")
ax_main.plot(h / 1000, delta_second * to_mGal,
             color=COLORS["navy"], linestyle="-", linewidth=1.5,
             label="Second-order")

ax_main.set_xlabel("Altitude (km)")
ax_main.set_ylabel(r"$\Delta\gamma$ (mGal)")
ax_main.legend(loc="lower left", frameon=False)
ax_main.text(0.03, 0.97, r"$\mathbf{(a)}$", transform=ax_main.transAxes,
             fontsize=10, fontweight="bold", va="top")
ax_main.set_xlim(0, 10)

# --- Panel (b): difference second - first order ---
diff = (delta_second - delta_first) * to_uGal  # in uGal

ax_inset.plot(h / 1000, diff, color=COLORS["navy"], linewidth=1.5)
ax_inset.set_xlabel("Altitude (km)")
ax_inset.set_ylabel(r"$\Delta\gamma_2 - \Delta\gamma_1$ ($\mu$Gal)")
ax_inset.text(0.05, 0.95, r"$\mathbf{(b)}$", transform=ax_inset.transAxes,
              fontsize=10, fontweight="bold", va="top")
ax_inset.set_xlim(0, 10)

# Annotate value at h=1 km
h_1km_idx = np.argmin(np.abs(h - 1000))
val_1km = diff[h_1km_idx]
ax_inset.annotate(f"{val_1km:.1f} $\\mu$Gal\nat 1 km",
                  xy=(1.0, val_1km), xytext=(3.0, val_1km + 5),
                  fontsize=7, color=COLORS["navy"],
                  arrowprops=dict(arrowstyle="->", color=COLORS["navy"],
                                  lw=0.7, shrinkB=3))

add_dual_axis(ax_main, 'mGal')
add_dual_axis(ax_inset, 'µGal')

fig.savefig("/home/xeal/dev/pytheas/doc/figures/fig05_free_air.png")
plt.close(fig)
print("Saved fig05_free_air.png")

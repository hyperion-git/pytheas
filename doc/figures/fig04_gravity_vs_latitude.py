#!/usr/bin/env python3
"""Figure 4: Normal gravity vs latitude.

Single panel comparing three models of surface gravity as a function
of geodetic latitude (0-90 degrees):
  1. Spherical approximation (constant GM/R^2)
  2. Sphere + centrifugal correction
  3. Full WGS84 Somigliana formula
"""

import sys
sys.path.insert(0, "/home/xeal/dev/pytheas/doc/figures")

import numpy as np
import matplotlib.pyplot as plt

import pytheas
from style import apply_style, COLORS, add_dual_axis

apply_style()

# ---- Parameters ----
GM_E  = 3.986004418e14   # m^3/s^2
R     = 6371000.0         # mean Earth radius (m)
OMEGA = 7.292115e-5       # rad/s

lat_deg = np.linspace(0, 90, 500)
lat_rad = np.radians(lat_deg)

# Curve 1: Spherical (constant)
g_sphere = np.full_like(lat_deg, GM_E / R**2)

# Curve 2: Sphere + centrifugal
g_rot = g_sphere - OMEGA**2 * R * np.cos(lat_rad)**2

# Curve 3: Full Somigliana (WGS84)
g_grs80 = np.array([pytheas.normal_gravity(lat, 0.0) for lat in lat_deg])

# ---- Plot ----
fig, ax = plt.subplots(figsize=(6, 3.5))

ax.plot(lat_deg, g_sphere, color="gray", linestyle="--", linewidth=1.3,
        label="Sphere")
ax.plot(lat_deg, g_rot, color=COLORS["orange"], linestyle=":", linewidth=1.5,
        label="+ Rotation")
ax.plot(lat_deg, g_grs80, color=COLORS["navy"], linestyle="-", linewidth=1.5,
        label="WGS84 (Somigliana)")

ax.set_xlabel("Geodetic latitude (degrees)")
ax.set_ylabel(r"Surface gravity (m/s$^2$)")
ax.set_xlim(0, 90)
ax.legend(loc="center right", frameon=False)

# Annotate equator and pole values
# Equator
eq_grs = g_grs80[0]
ax.annotate(f"{eq_grs:.4f}",
            xy=(0, eq_grs), xytext=(8, eq_grs - 0.008),
            fontsize=7.5, color=COLORS["navy"],
            arrowprops=dict(arrowstyle="-", color=COLORS["navy"],
                            lw=0.6, shrinkB=2))

# Pole
pole_grs = g_grs80[-1]
ax.annotate(f"{pole_grs:.4f}",
            xy=(90, pole_grs), xytext=(75, pole_grs + 0.005),
            fontsize=7.5, color=COLORS["navy"],
            arrowprops=dict(arrowstyle="-", color=COLORS["navy"],
                            lw=0.6, shrinkB=2))

add_dual_axis(ax, 'm/s²')

fig.savefig("/home/xeal/dev/pytheas/doc/figures/fig04_gravity_vs_latitude.png")
plt.close(fig)
print("Saved fig04_gravity_vs_latitude.png")

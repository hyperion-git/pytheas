#!/usr/bin/env python3
"""Figure 11: Rigid vs elastic Earth.

Single panel showing 48 hours at Munich with two curves: rigid Earth
(delta = 1) and elastic Earth (delta = 1.16), demonstrating the ~16%
amplification from solid-Earth deformation.
"""

import sys
sys.path.insert(0, "/home/xeal/dev/pytheas/doc/figures")
sys.path.insert(0, "/home/xeal/dev/pytheas")

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

import pytheas
from pytheas import (
    GM_MOON, GM_SUN, DELTA_GRAV,
    geodetic_to_ecef, moon_position_ecef, sun_position_ecef,
    tidal_acceleration, enu_basis, compute_timeseries,
)
from style import apply_style, COLORS, add_dual_axis

apply_style()

# ---- Parameters ----
lat, lon, alt = 48.14, 11.58, 500.0
start = datetime(2025, 3, 20, 0, 0, 0)
end = start + timedelta(hours=48)
interval_min = 10.0

r_obs = geodetic_to_ecef(lat, lon, alt)
_, _, e_up = enu_basis(lat, lon)

# ---- Elastic Earth (standard output) ----
ts = compute_timeseries(start, end, lat, lon, alt, interval_minutes=interval_min)
hours = np.array([(t - start).total_seconds() / 3600.0 for t in ts["times"]])
g_elastic = ts["g_tidal"]  # already includes DELTA_GRAV

# ---- Rigid Earth (delta = 1, no amplification) ----
# The elastic output has DELTA_GRAV baked in, so rigid = elastic / DELTA_GRAV
g_rigid = g_elastic / DELTA_GRAV

# Convert to um/s^2
scale = 1e6

# ---- Plot ----
fig, ax = plt.subplots(figsize=(6, 3.5))

ax.plot(hours, g_rigid * scale,
        color=COLORS["navy"], linewidth=1.3, linestyle="--",
        label=r"Rigid Earth ($\delta = 1$)")
ax.plot(hours, g_elastic * scale,
        color=COLORS["red"], linewidth=1.5,
        label=r"Elastic Earth ($\delta = 1.16$)")

# Fill between
ax.fill_between(hours, g_rigid * scale, g_elastic * scale,
                color=COLORS["orange"], alpha=0.20)

# Annotate the fill region
# Find a peak to place the annotation
idx_peak = np.argmax(np.abs(g_elastic))
x_ann = hours[idx_peak]
y_mid = 0.5 * (g_rigid[idx_peak] + g_elastic[idx_peak]) * scale

# Place annotation near a peak but offset for readability
ax.annotate(
    r"$\sim$16% amplification" "\nfrom Earth deformation",
    xy=(x_ann, y_mid),
    xytext=(x_ann + 6, y_mid * 0.5 if y_mid > 0 else y_mid * 1.5),
    fontsize=7.5,
    ha="left", va="center",
    color=COLORS["red"],
    arrowprops=dict(arrowstyle="->", color=COLORS["red"],
                    lw=0.8, connectionstyle="arc3,rad=0.2"),
)

ax.set_xlabel("Hours since 2025-03-20 00:00 UTC")
ax.set_ylabel(r"Total tidal $g_z$ ($\mu$m/s$^2$)")
ax.legend(loc="upper right", frameon=False)
ax.set_xlim(0, 48)
ax.axhline(0, color="grey", linewidth=0.4, linestyle="-")

add_dual_axis(ax, 'µm/s²')

fig.savefig("/home/xeal/dev/pytheas/doc/figures/fig11_rigid_vs_elastic.png")
plt.close(fig)
print("Saved fig11_rigid_vs_elastic.png")
print(f"DELTA_GRAV = {DELTA_GRAV:.4f}")

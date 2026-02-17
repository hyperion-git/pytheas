#!/usr/bin/env python3
"""Figure 8: Moon vs Sun tidal contribution.

Single panel showing 7 days of tidal acceleration at Munich,
with separate lunar and solar contributions highlighting
the ~2.2:1 amplitude ratio.
"""

import sys
sys.path.insert(0, "/home/xeal/dev/pytheas/doc/figures")

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

import pytheas
from style import apply_style, COLORS, add_dual_axis

apply_style()

# ---- Compute timeseries (7 days) ----
start = datetime(2025, 3, 20, 0, 0, 0)
end   = start + timedelta(days=7)
lat, lon, alt = 48.14, 11.58, 500.0

ts = pytheas.compute_timeseries(start, end, lat, lon, alt, interval_minutes=10)

hours = np.array([(t - start).total_seconds() / 3600.0 for t in ts["times"]])
days  = hours / 24.0

scale = 1e6  # m/s^2 -> um/s^2

# ---- Plot ----
fig, ax = plt.subplots(figsize=(6.5, 3.5))

ax.plot(days, ts["g_tidal"] * scale,
        color=COLORS["navy"], linewidth=1.5, label="Total tidal")
ax.plot(days, ts["g_tidal_moon"] * scale,
        color=COLORS["orange"], linewidth=1.3, linestyle="--", label="Lunar")
ax.plot(days, ts["g_tidal_sun"] * scale,
        color=COLORS["green"], linewidth=1.3, linestyle="-.", label="Solar")

ax.set_xlabel("Days since 2025-03-20 00:00 UTC")
ax.set_ylabel(r"$g_{\mathrm{tidal}}$ ($\mu$m/s$^2$)")
ax.set_xlim(0, 7)
ax.legend(loc="upper right", frameon=False)

# Annotate the amplitude ratio
lunar_amp = (np.max(ts["g_tidal_moon"]) - np.min(ts["g_tidal_moon"])) / 2 * scale
solar_amp = (np.max(ts["g_tidal_sun"]) - np.min(ts["g_tidal_sun"])) / 2 * scale
ratio = lunar_amp / solar_amp

# Add text annotation for ratio
ax.text(0.02, 0.05,
        f"Lunar/Solar amplitude ratio: {ratio:.1f}:1",
        transform=ax.transAxes, fontsize=8,
        color=COLORS["black"],
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                  edgecolor="gray", alpha=0.8))

# Look for spring/neap-like alignment: when lunar and solar peaks coincide
# or are in opposition.  Annotate if visible.
from scipy.signal import argrelextrema
total = ts["g_tidal"] * scale
maxima = argrelextrema(total, np.greater, order=10)[0]
minima = argrelextrema(total, np.less, order=10)[0]

if len(maxima) >= 2:
    # Find the largest peak-to-trough swing (spring-like)
    peak_vals = total[maxima]
    best_peak_idx = maxima[np.argmax(peak_vals)]
    best_day = days[best_peak_idx]
    best_val = total[best_peak_idx]

    ax.annotate("constructive\nalignment",
                xy=(best_day, best_val),
                xytext=(best_day + 0.5, best_val + 0.3),
                fontsize=7, color=COLORS["red"],
                arrowprops=dict(arrowstyle="->", color=COLORS["red"],
                                lw=0.7, shrinkB=3),
                ha="left")

add_dual_axis(ax, 'µm/s²')

fig.savefig("/home/xeal/dev/pytheas/doc/figures/fig08_moon_vs_sun.png")
plt.close(fig)
print("Saved fig08_moon_vs_sun.png")

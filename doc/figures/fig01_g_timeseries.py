#!/usr/bin/env python3
"""Figure 1: A day in the life of g(t).

Two vertically stacked panels (shared x-axis), showing 48 hours of
gravitational acceleration at Munich.

(a) Total g(t) in m/s^2.
(b) Tidal residual in um/s^2 with lunar, solar, and total contributions.
"""

import sys
sys.path.insert(0, "/home/xeal/dev/pytheas/doc/figures")

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

import pytheas
from style import apply_style, COLORS, add_dual_axis

apply_style()

# ---- Compute timeseries ----
start = datetime(2025, 3, 20, 0, 0, 0)
end   = start + timedelta(hours=48)
lat, lon, alt = 48.14, 11.58, 500.0

ts = pytheas.compute_timeseries(start, end, lat, lon, alt, interval_minutes=10)

hours = np.array([(t - start).total_seconds() / 3600.0 for t in ts["times"]])

# ---- Create figure ----
fig, (ax_top, ax_bot) = plt.subplots(
    2, 1, sharex=True, figsize=(6, 4.5),
    gridspec_kw={"height_ratios": [1, 2], "hspace": 0.08},
)

# --- Panel (a): Total g(t) ---
ax_top.plot(hours, ts["g_total"], color=COLORS["navy"], linewidth=1.2)
ax_top.set_ylabel(r"$g(t)$ (m/s$^2$)")
ax_top.text(0.02, 0.90, r"$\mathbf{(a)}$", transform=ax_top.transAxes,
            fontsize=10, fontweight="bold", va="top")

# Format y-axis to show enough decimal places
ax_top.ticklabel_format(axis="y", useOffset=True, style="plain")
# Show ticks with more precision
ax_top.yaxis.set_major_formatter(plt.FormatStrFormatter("%.6f"))

# --- Panel (b): Tidal residual ---
scale = 1e6  # m/s^2 -> um/s^2

ax_bot.plot(hours, ts["g_tidal"] * scale,
            color=COLORS["navy"], linewidth=1.5, label="Total tidal")
ax_bot.plot(hours, ts["g_tidal_moon"] * scale,
            color=COLORS["orange"], linewidth=1.3, linestyle="--", label="Lunar")
ax_bot.plot(hours, ts["g_tidal_sun"] * scale,
            color=COLORS["green"], linewidth=1.3, linestyle=":", label="Solar")

ax_bot.set_ylabel(r"$g_{\mathrm{tidal}}$ ($\mu$m/s$^2$)")
ax_bot.set_xlabel("Hours since 2025-03-20 00:00 UTC")
ax_bot.legend(loc="upper right", frameon=False)
ax_bot.text(0.02, 0.95, r"$\mathbf{(b)}$", transform=ax_bot.transAxes,
            fontsize=10, fontweight="bold", va="top")

# --- Annotate semidiurnal period (~12.4 h) ---
# Find two consecutive peaks in the total tidal signal
from scipy.signal import argrelextrema
tidal = ts["g_tidal"] * scale
maxima = argrelextrema(tidal, np.greater, order=10)[0]

if len(maxima) >= 2:
    h1, h2 = hours[maxima[0]], hours[maxima[1]]
    y_arrow = tidal[maxima[0]] * 0.85  # slightly below peak
    # Double-headed arrow
    ax_bot.annotate(
        "", xy=(h2, y_arrow), xytext=(h1, y_arrow),
        arrowprops=dict(arrowstyle="<->", color=COLORS["red"],
                        lw=1.2, shrinkA=0, shrinkB=0),
    )
    period_h = h2 - h1
    ax_bot.text(
        (h1 + h2) / 2, y_arrow + 0.15,
        f"~{period_h:.1f} h",
        ha="center", va="bottom", fontsize=8, color=COLORS["red"],
    )

ax_bot.set_xlim(0, 48)

add_dual_axis(ax_top, 'm/s²')
add_dual_axis(ax_bot, 'µm/s²')

fig.savefig("/home/xeal/dev/pytheas/doc/figures/fig01_g_timeseries.png")
plt.close(fig)
print("Saved fig01_g_timeseries.png")

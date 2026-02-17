#!/usr/bin/env python3
"""Figure 12: Anatomy of g(t).

Three stacked panels showing how the model components combine:
(a) Static gravity -- a horizontal line at the normal-gravity value.
(b) Lunar tidal contribution in um/s^2.
(c) Solar tidal contribution in um/s^2.
"""

import sys
sys.path.insert(0, "/home/xeal/dev/pytheas/doc/figures")
sys.path.insert(0, "/home/xeal/dev/pytheas")

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

import pytheas
from pytheas import compute_timeseries
from style import apply_style, COLORS, add_dual_axis

apply_style()

# ---- Parameters ----
lat, lon, alt = 48.14, 11.58, 500.0
start = datetime(2025, 3, 20, 0, 0, 0)
end = start + timedelta(hours=48)

ts = compute_timeseries(start, end, lat, lon, alt, interval_minutes=10)
hours = np.array([(t - start).total_seconds() / 3600.0 for t in ts["times"]])

g_static_val = ts["g_static"][0]  # constant
g_moon = ts["g_tidal_moon"]
g_sun = ts["g_tidal_sun"]

scale = 1e6  # m/s^2 -> um/s^2

# ---- Create figure ----
fig, (ax_a, ax_b, ax_c) = plt.subplots(
    3, 1, sharex=True, figsize=(6, 5.5),
    gridspec_kw={"height_ratios": [0.6, 1, 1], "hspace": 0.10},
)

# --- Panel (a): Static gravity ---
ax_a.axhline(g_static_val, color=COLORS["navy"], linewidth=1.5)
ax_a.set_ylabel(r"$g_{\mathrm{static}}$ (m/s$^2$)")
ax_a.set_ylim(g_static_val - 0.0001, g_static_val + 0.0001)
ax_a.yaxis.set_major_formatter(plt.FormatStrFormatter("%.5f"))
ax_a.text(0.02, 0.85, r"$\mathbf{(a)}$ Static gravity", transform=ax_a.transAxes,
          fontsize=9, fontweight="bold", va="top")
ax_a.annotate(
    f"{g_static_val:.6f} m/s$^2$",
    xy=(24, g_static_val), xytext=(24, g_static_val + 0.00006),
    fontsize=8, ha="center", va="bottom", color=COLORS["navy"],
    arrowprops=dict(arrowstyle="->", color=COLORS["navy"], lw=0.7),
)

# --- Panel (b): Lunar tidal ---
ax_b.plot(hours, g_moon * scale, color=COLORS["orange"], linewidth=1.3)
ax_b.fill_between(hours, 0, g_moon * scale, color=COLORS["orange"], alpha=0.12)
ax_b.axhline(0, color="grey", linewidth=0.4)
ax_b.set_ylabel(r"$g_{\mathrm{lunar}}$ ($\mu$m/s$^2$)")
ax_b.text(0.02, 0.92, r"$\mathbf{(b)}$ Lunar tidal", transform=ax_b.transAxes,
          fontsize=9, fontweight="bold", va="top")

# Annotate peak-to-peak amplitude
pp_moon = (np.max(g_moon) - np.min(g_moon)) * scale
ax_b.text(0.98, 0.92, f"peak-to-peak: {pp_moon:.1f} $\\mu$m/s$^2$",
          transform=ax_b.transAxes, fontsize=7.5, ha="right", va="top",
          color=COLORS["orange"])

# --- Panel (c): Solar tidal ---
ax_c.plot(hours, g_sun * scale, color=COLORS["green"], linewidth=1.3)
ax_c.fill_between(hours, 0, g_sun * scale, color=COLORS["green"], alpha=0.12)
ax_c.axhline(0, color="grey", linewidth=0.4)
ax_c.set_ylabel(r"$g_{\mathrm{solar}}$ ($\mu$m/s$^2$)")
ax_c.set_xlabel("Hours since 2025-03-20 00:00 UTC")
ax_c.text(0.02, 0.92, r"$\mathbf{(c)}$ Solar tidal", transform=ax_c.transAxes,
          fontsize=9, fontweight="bold", va="top")

# Annotate peak-to-peak amplitude
pp_sun = (np.max(g_sun) - np.min(g_sun)) * scale
ax_c.text(0.98, 0.92, f"peak-to-peak: {pp_sun:.1f} $\\mu$m/s$^2$",
          transform=ax_c.transAxes, fontsize=7.5, ha="right", va="top",
          color=COLORS["green"])

ax_c.set_xlim(0, 48)

add_dual_axis(ax_a, 'm/s²')
add_dual_axis(ax_b, 'µm/s²')
add_dual_axis(ax_c, 'µm/s²')

# --- Bottom annotation ---
fig.text(0.5, 0.01,
         r"Elastic amplification ($\delta$ = 1.16) already applied to tidal components",
         ha="center", fontsize=8, style="italic", color="grey")

fig.savefig("/home/xeal/dev/pytheas/doc/figures/fig12_anatomy.png")
plt.close(fig)
print("Saved fig12_anatomy.png")
print(f"g_static = {g_static_val:.6f} m/s^2")
print(f"Lunar peak-to-peak: {pp_moon:.2f} um/s^2")
print(f"Solar peak-to-peak: {pp_sun:.2f} um/s^2")

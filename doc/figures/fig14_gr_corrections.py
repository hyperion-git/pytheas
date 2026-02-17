#!/usr/bin/env python3
"""Figure 14: Complete hierarchy of corrections to the gravimeter signal.

Horizontal bar chart showing every contribution from the Fermi frame analysis,
spanning from centrifugal (10^6 µGal) down to spin GM tidal (10^{-15} nGal),
with the 10 nGal accuracy floor.
"""

import sys
sys.path.insert(0, "/home/xeal/dev/pytheas/doc/figures")
sys.path.insert(0, "/home/xeal/dev/pytheas")

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from style import apply_style, COLORS, add_dual_axis

apply_style()

# ---- Data: full hierarchy from Section 8.6 ----
# Values in nGal; categories:
#   absorbed    = time-independent, absorbed into gamma(phi,h)
#   signal      = measured tidal signal (what Pytheas models)
#   horizontal  = no projection onto measurement axis
#   instrument  = instrument-dependent systematic
#   below_floor = below 10 nGal accuracy floor
entries = [
    {"label": "Static normal gravity $\\gamma(\\varphi,h)$",
     "nGal": 9.81e11, "category": "absorbed"},
    {"label": "Centrifugal",
     "nGal": 3.39e9,  "category": "absorbed"},
    {"label": "Newtonian tidal (Moon)",
     "nGal": 110_000, "category": "signal"},
    {"label": "Newtonian tidal (Sun)",
     "nGal": 50_000,  "category": "signal"},
    {"label": "Coriolis (1st order, horizontal)",
     "nGal": 2.5e7,   "category": "horizontal"},
    {"label": "Eötvös (2nd-order Coriolis)",
     "nGal": 319,     "category": "instrument"},
    {"label": "GM curvature cross-term",
     "nGal": 0.17,    "category": "below_floor"},
    {"label": "Direct Lense–Thirring",
     "nGal": 0.008,   "category": "below_floor"},
    {"label": "1PN correction (Sun)",
     "nGal": 0.001,   "category": "below_floor"},
    {"label": "1PN correction (Moon)",
     "nGal": 3e-8,    "category": "below_floor"},
    {"label": "Orbital GM tidal (Sun)",
     "nGal": 1e-9,    "category": "below_floor"},
    {"label": "Orbital GM tidal (Moon)",
     "nGal": 1e-10,   "category": "below_floor"},
    {"label": "Spin GM tidal (Sun)",
     "nGal": 2e-13,   "category": "below_floor"},
    {"label": "Spin GM tidal (Moon)",
     "nGal": 5e-15,   "category": "below_floor"},
]

# Sort by magnitude (smallest at bottom, largest at top)
entries.sort(key=lambda e: e["nGal"])

labels = [e["label"] for e in entries]
magnitudes = np.array([e["nGal"] for e in entries])
categories = [e["category"] for e in entries]

# Color map by category
color_map = {
    "absorbed":    "0.55",           # gray — time-independent, in gamma
    "signal":      COLORS["navy"],   # the tidal signal Pytheas models
    "horizontal":  COLORS["green"],  # horizontal, no n-hat projection
    "instrument":  COLORS["blue"],   # instrument-dependent systematic
    "below_floor": COLORS["orange"], # below 10 nGal accuracy floor
}
bar_colors = [color_map[c] for c in categories]

# ---- Plot ----
fig, ax = plt.subplots(figsize=(8, 6))

y_pos = np.arange(len(entries))
bars = ax.barh(y_pos, magnitudes, color=bar_colors, height=0.6, edgecolor="none")

ax.set_yticks(y_pos)
ax.set_yticklabels(labels, fontsize=7.5)
ax.set_xlabel("Magnitude (nGal)")
ax.set_xscale("log")
ax.set_xlim(1e-16, 5e15)

# Add value labels on bars
for i, (mag, bar) in enumerate(zip(magnitudes, bars)):
    if mag >= 1:
        label = f"{mag:,.3g}"
    else:
        exp = int(np.floor(np.log10(mag)))
        coeff = mag / 10**exp
        if abs(coeff - 1.0) < 0.05:
            label = f"$10^{{{exp}}}$"
        else:
            label = f"${coeff:.0f}" + r"\times 10^{" + f"{exp}" + r"}$"
    x_text = mag * 2.5
    ax.text(x_text, i, label, va="center", fontsize=7, color="0.3")

# Accuracy floor line
ax.axvline(10, color=COLORS["red"], linestyle="--", linewidth=1.2, zorder=0)
ax.text(10, len(entries) - 0.3, "  Pytheas accuracy\n  floor (10 nGal)",
        fontsize=7.5, color=COLORS["red"], va="top")

# Legend
legend_elements = [
    Patch(facecolor="0.55",           label="Absorbed into $\\gamma(\\varphi,h)$"),
    Patch(facecolor=COLORS["navy"],   label="Measured tidal signal"),
    Patch(facecolor=COLORS["green"],  label="Horizontal (no $\\hat{n}$ projection)"),
    Patch(facecolor=COLORS["blue"],   label="Instrument-dependent"),
    Patch(facecolor=COLORS["orange"], label="Below accuracy floor"),
]
ax.legend(handles=legend_elements, loc="lower right", frameon=False,
          fontsize=7.5)

add_dual_axis(ax, 'nGal', direction='x')

fig.savefig("/home/xeal/dev/pytheas/doc/figures/fig14_gr_corrections.png")
plt.close(fig)
print("Saved fig14_gr_corrections.png")

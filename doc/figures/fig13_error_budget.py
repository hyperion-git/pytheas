#!/usr/bin/env python3
"""Figure 13: Error budget.

Horizontal bar chart showing estimated error sources and their magnitudes
in nGal, colored by modeling status.
"""

import sys
sys.path.insert(0, "/home/xeal/dev/pytheas/doc/figures")
sys.path.insert(0, "/home/xeal/dev/pytheas")

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from style import apply_style, COLORS

apply_style()

# ---- Error budget data ----
errors = [
    {"source": "Ocean loading",            "nGal": 25,   "status": "not modeled"},
    {"source": "Elastic Earth (no FCN)",    "nGal": 10,   "status": "limitation"},
    {"source": "Atmospheric pressure",      "nGal": 3,    "status": "not modeled"},
    {"source": "Lunar ephemeris (position)","nGal": 1,    "status": "modeled"},
    {"source": "Lunar ephemeris (distance)","nGal": 1,    "status": "modeled"},
    {"source": "Solar ephemeris",           "nGal": 0.5,  "status": "modeled"},
    {"source": "Normal gravity",            "nGal": 0.5,  "status": "modeled"},
    {"source": "Polar motion",              "nGal": 0.5,  "status": "not modeled"},
    {"source": "Planetary tides",           "nGal": 0.1,  "status": "not modeled"},
]

# Sort by magnitude (largest at top)
errors.sort(key=lambda e: e["nGal"])

sources = [e["source"] for e in errors]
magnitudes = np.array([e["nGal"] for e in errors])
statuses = [e["status"] for e in errors]

# Color map
color_map = {
    "modeled":     COLORS["navy"],
    "not modeled": COLORS["orange"],
    "limitation":  COLORS["pink"],
}
bar_colors = [color_map[s] for s in statuses]

# ---- Plot ----
fig, ax = plt.subplots(figsize=(6, 4))

y_pos = np.arange(len(errors))
bars = ax.barh(y_pos, magnitudes, color=bar_colors, height=0.6, edgecolor="none")

ax.set_yticks(y_pos)
ax.set_yticklabels(sources)
ax.set_xlabel("Error magnitude (nGal)")
ax.set_xscale("log")

# Set x limits to give some breathing room
ax.set_xlim(0.05, 80)

# Add value labels on bars
for i, (mag, bar) in enumerate(zip(magnitudes, bars)):
    x_text = mag * 1.3
    ax.text(x_text, i, f"{mag:g}", va="center", fontsize=7.5,
            color="0.3")

# Legend
legend_elements = [
    Patch(facecolor=COLORS["navy"],   label="Modeled by Pytheas"),
    Patch(facecolor=COLORS["orange"], label="Not modeled"),
    Patch(facecolor=COLORS["pink"],   label="Known limitation"),
]
ax.legend(handles=legend_elements, loc="lower right", frameon=False,
          fontsize=7.5)

# Accuracy annotation
ax.text(0.98, 0.98,
        "Accuracy: ~10 nGal inland, ~50 nGal coastal",
        transform=ax.transAxes, fontsize=8, ha="right", va="top",
        style="italic", color="0.4",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="0.95",
                  edgecolor="0.8", linewidth=0.5))

fig.savefig("/home/xeal/dev/pytheas/doc/figures/fig13_error_budget.png")
plt.close(fig)
print("Saved fig13_error_budget.png")

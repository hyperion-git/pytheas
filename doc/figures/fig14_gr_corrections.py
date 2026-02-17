#!/usr/bin/env python3
"""Figure 14: Hierarchy of GR corrections to the tidal signal.

Horizontal bar chart comparing Newtonian tidal magnitudes with post-Newtonian
and gravitomagnetic corrections, demonstrating their negligibility relative
to Pytheas's ~10 nGal accuracy floor.
"""

import sys
sys.path.insert(0, "/home/xeal/dev/pytheas/doc/figures")
sys.path.insert(0, "/home/xeal/dev/pytheas")

import numpy as np
import matplotlib.pyplot as plt
from style import apply_style, COLORS

apply_style()

# ---- Data: tidal signals and GR corrections ----
# Values from Section 8 derivation
entries = [
    {"label": "Newtonian tidal (Moon)",         "nGal": 110_000,    "category": "newtonian"},
    {"label": "Newtonian tidal (Sun)",          "nGal": 50_000,     "category": "newtonian"},
    {"label": "1PN correction (Sun)",           "nGal": 0.001,      "category": "1pn"},
    {"label": "Earth frame-dragging\n(direct Lense–Thirring)",
                                                "nGal": 0.008,      "category": "gm"},
    {"label": "1PN correction (Moon)",          "nGal": 3e-8,       "category": "1pn"},
    {"label": "Orbital GM tidal (Sun)",         "nGal": 1e-9,       "category": "gm"},
    {"label": "Orbital GM tidal (Moon)",        "nGal": 1e-10,      "category": "gm"},
    {"label": "Spin GM tidal (Sun)",            "nGal": 2e-13,      "category": "gm"},
    {"label": "Spin GM tidal (Moon)",           "nGal": 5e-15,      "category": "gm"},
]

# Sort by magnitude (largest at top)
entries.sort(key=lambda e: e["nGal"])

labels = [e["label"] for e in entries]
magnitudes = np.array([e["nGal"] for e in entries])
categories = [e["category"] for e in entries]

# Color map by category
color_map = {
    "newtonian": COLORS["navy"],
    "1pn":       COLORS["orange"],
    "gm":        COLORS["pink"],
}
bar_colors = [color_map[c] for c in categories]

# ---- Plot ----
fig, ax = plt.subplots(figsize=(7, 5))

y_pos = np.arange(len(entries))
bars = ax.barh(y_pos, magnitudes, color=bar_colors, height=0.6, edgecolor="none")

ax.set_yticks(y_pos)
ax.set_yticklabels(labels, fontsize=8)
ax.set_xlabel("Magnitude (nGal)")
ax.set_xscale("log")
ax.set_xlim(1e-16, 5e5)

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
    ax.text(x_text, i, label, va="center", fontsize=7.5, color="0.3")

# Accuracy floor line
ax.axvline(10, color=COLORS["red"], linestyle="--", linewidth=1.2, zorder=0)
ax.text(10, len(entries) - 0.3, "  Pytheas accuracy\n  floor (10 nGal)",
        fontsize=7.5, color=COLORS["red"], va="top")

# Legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor=COLORS["navy"],   label="Newtonian tidal"),
    Patch(facecolor=COLORS["orange"], label="1PN correction"),
    Patch(facecolor=COLORS["pink"],   label="Gravitomagnetic"),
]
ax.legend(handles=legend_elements, loc="lower right", frameon=False,
          fontsize=8)

fig.savefig("/home/xeal/dev/pytheas/doc/figures/fig14_gr_corrections.png")
plt.close(fig)
print("Saved fig14_gr_corrections.png")

#!/usr/bin/env python3
"""Figure 15: Conceptual hierarchy from bitensor theory to LabFrame code.

Vertical flowchart showing the path from the Synge world function and
parallel propagator down to the Pytheas LabFrame / GravityField API.
"""

import sys
sys.path.insert(0, "/home/xeal/dev/pytheas/doc/figures")
sys.path.insert(0, "/home/xeal/dev/pytheas")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from style import apply_style, COLORS

apply_style()

# ---- Box definitions (top to bottom) ----
boxes = [
    {
        "title": "Synge World Function $\\sigma(x,x')$\n"
                 "Parallel Propagator $g^{a}_{\\;b'}$\n"
                 "[Bitensor Framework]",
        "color": COLORS["navy"],
    },
    {
        "title": "FNC Metric Expansion\n"
                 "$g_{00},\\; g_{0i},\\; g_{ij} \\to R_{abcd}$\n"
                 "[Fermi Normal Coordinates]",
        "color": COLORS["navy"],
    },
    {
        "title": "Rotating-Frame Lagrangian\n"
                 "$L = \\frac{1}{2}|\\mathbf{v}|^2 "
                 "+ (\\boldsymbol{\\Omega} \\times \\mathbf{x})\\cdot\\mathbf{v}"
                 " - \\Phi_{eff}$\n"
                 "[Canonical Transformation]",
        "color": COLORS["navy"],
    },
    {
        "title": "Taylor Coefficients\n"
                 "$g_i,\\; T_{ij},\\; \\Omega_i$\n"
                 "[Local Expansion]",
        "color": COLORS["navy"],
    },
    {
        "title": "Equation of Motion\n"
                 "$\\mathbf{a} = \\mathbf{g} + \\mathbf{T}\\,\\delta\\mathbf{x}"
                 " - 2(\\boldsymbol{\\Omega} \\times \\mathbf{v})$\n"
                 "[Euler\u2013Lagrange]",
        "color": COLORS["orange"],
    },
    {
        "title": "LabFrame / GravityField\n"
                 ".g   .T   .omega   .eom()\n"
                 "[Pytheas Implementation]",
        "color": COLORS["green"],
    },
]

# ---- Layout parameters ----
fig_w, fig_h = 5, 8
box_w = 3.6
box_h = 0.85
gap = 0.45          # vertical gap between boxes
top_margin = 0.55   # from top of figure to top of first box

n = len(boxes)
total_h = n * box_h + (n - 1) * gap
y_start = fig_h - top_margin  # top of first box

fig, ax = plt.subplots(figsize=(fig_w, fig_h))
ax.set_xlim(0, fig_w)
ax.set_ylim(0, fig_h)
ax.axis("off")

x_center = fig_w / 2
x_left = x_center - box_w / 2

for i, box in enumerate(boxes):
    y_top = y_start - i * (box_h + gap)
    y_bot = y_top - box_h

    # Draw rounded box
    fc = box["color"]
    patch = FancyBboxPatch(
        (x_left, y_bot), box_w, box_h,
        boxstyle="round,pad=0.12",
        facecolor=fc, edgecolor="white", linewidth=1.5, alpha=0.92,
        zorder=2,
    )
    ax.add_patch(patch)

    # Text (white on dark background)
    ax.text(
        x_center, y_bot + box_h / 2, box["title"],
        ha="center", va="center", fontsize=8,
        color="white", fontweight="medium",
        linespacing=1.35, zorder=3,
    )

    # Arrow to next box
    if i < n - 1:
        arrow_top = y_bot - 0.04
        arrow_bot = y_bot - gap + 0.04
        arrow = FancyArrowPatch(
            (x_center, arrow_top), (x_center, arrow_bot),
            arrowstyle="-|>", color="0.35",
            linewidth=1.5, mutation_scale=14, zorder=1,
        )
        ax.add_patch(arrow)

fig.savefig("/home/xeal/dev/pytheas/doc/figures/fig15_lab_frame_hierarchy.png")
plt.close(fig)
print("Saved fig15_lab_frame_hierarchy.png")

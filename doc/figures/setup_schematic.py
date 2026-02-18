"""
Schematic planar solar system diagram for Section II.A.

Shows the four point masses (Sun, Earth, Moon, test mass) in an inertial
frame with all position vectors labeled.  Not to scale — the purpose is
to visualize the setup and the vector definitions used in the proof.
"""

import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle

# --- Style: publication, single-column PRX ---
plt.rcParams.update({
    "figure.figsize": (3.3, 2.5),
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
    "font.family": "sans-serif",
    "font.sans-serif": [
        "DejaVu Sans", "Arial", "Helvetica",
        "Lucida Grande", "Verdana", "sans-serif",
    ],
    "font.size": 7,
    "axes.labelsize": 9,
    "axes.titlesize": 7,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "text.usetex": False,
    "axes.formatter.use_mathtext": True,
    "mathtext.fontset": "dejavusans",
    "axes.linewidth": 0.5,
    "lines.linewidth": 1.0,
})

# --- Body positions (schematic, not to scale) ---
# Spread out so labels don't overlap; Earth-test mass gap exaggerated
O  = np.array([-0.5, -0.5])     # inertial origin (offset from center)
S  = np.array([-3.2, 1.6])      # Sun (upper left)
E  = np.array([1.2, -0.6])      # Earth (lower right)
M  = np.array([3.4, 1.6])       # Moon (upper right)

# Test mass + spring assembly: spring starts at Earth's edge, test mass
# sits at the far end.  Direction is ~15° above horizontal (East).
_spring_angle = np.radians(15 - 90)  # rotated 90° clockwise
_spring_dir = np.array([np.cos(_spring_angle), np.sin(_spring_angle)])
_spring_perp = np.array([-_spring_dir[1], _spring_dir[0]])  # perpendicular
_earth_radius = 0.14  # must match body_data below
_spring_length = 0.8  # visual length of the spring coils
_test_radius = 0.06
# Spring endpoints: from Earth's edge to test mass edge (touching both)
_total_gap = _spring_length  # gap between surfaces filled by spring
T = E + (_earth_radius + _total_gap + _test_radius) * _spring_dir
_spring_start = E + _earth_radius * _spring_dir       # Earth's surface
_spring_end   = T - _test_radius * _spring_dir         # test mass surface

# --- Colors ---
c_sun   = "#fd7400"   # orange
c_earth = "#2e2cb8"   # blue
c_moon  = "#9b2f5c"   # magenta
c_test  = "#1f8a70"   # teal
c_vec   = "#888888"   # gray for position vectors from origin
c_rel   = "#db002b"   # red for relative/geocentric vectors

fig, ax = plt.subplots(figsize=(3.5, 2.8))
ax.set_aspect("equal")
ax.axis("off")

# --- Helper: draw arrow ---
def draw_arrow(start, end, color, lw=0.8, style="-|>", zorder=2,
               shrinkA=0, shrinkB=0, ls="-"):
    arrow = FancyArrowPatch(
        start, end,
        arrowstyle=style, color=color, linewidth=lw,
        mutation_scale=8, zorder=zorder,
        shrinkA=shrinkA, shrinkB=shrinkB, linestyle=ls,
    )
    ax.add_patch(arrow)

# --- Helper: draw spring between two points ---
def draw_spring(p1, p2, n_coils=6, amplitude=0.08, color="#666666",
                lw=0.7):
    d = p2 - p1
    length = np.linalg.norm(d)
    t_hat = d / length
    n_hat = np.array([-t_hat[1], t_hat[0]])

    # Straight leads at each end
    lead = 0.12
    # Use many points for smooth zigzag
    n_pts = n_coils * 100
    ts = np.linspace(0, 1, n_pts + 1)
    pts = [p1, p1 + lead * d]
    for t in ts:
        frac = lead + (1 - 2 * lead) * t
        sign = np.sin(2 * np.pi * n_coils * t)
        pt = p1 + frac * d + sign * amplitude * n_hat
        pts.append(pt)
    pts.append(p1 + (1 - lead) * d)
    pts.append(p2)
    pts = np.array(pts)
    ax.plot(pts[:, 0], pts[:, 1], color=color, linewidth=lw,
            zorder=3, solid_capstyle="round")

# --- Draw bodies ---
body_data = [
    (S, c_sun,   0.20, "$M_S$", (-0.05, 0.30), 7),
    (E, c_earth, 0.14, "$M_E$", (0.0, -0.28), 7),
    (M, c_moon,  0.10, "$M_M$", (0.0, 0.20), 7),
    (T, c_test,  0.06, "$m$",   (0.15, 0.08), 6.5),
]

for pos, color, radius, label, offset, fs in body_data:
    circle = Circle(pos, radius, facecolor=color, edgecolor="black",
                    linewidth=0.4, zorder=5, alpha=0.85)
    ax.add_patch(circle)
    ax.text(pos[0] + offset[0], pos[1] + offset[1], label,
            fontsize=fs, ha="center", va="center", zorder=6)

# --- Inertial origin marker ---
ax.plot(*O, "+", color="black", markersize=6, markeredgewidth=0.8,
        zorder=4)
ax.text(O[0] - 0.22, O[1], r"$\mathcal{O}$", fontsize=7.5,
        ha="center", va="center", color="black")

# --- Position vectors from origin (gray, thin) ---
vec_data = [
    (S, r"$\mathbf{R}_S$", "above"),
    (E, r"$\mathbf{R}_E$", "below"),
    (M, r"$\mathbf{R}_M$", "above"),
    (T, r"$\mathbf{x}$",   "above"),
]

for pos, label, placement in vec_data:
    draw_arrow(O, pos, c_vec, lw=0.4, shrinkB=0, ls="--")
    mid = 0.45 * O + 0.55 * pos
    if placement == "above":
        # Place label above the midpoint of the vector
        d = pos - O
        n = np.array([-d[1], d[0]])
        n = n / np.linalg.norm(n) * 0.18
        lpos = mid + n
    else:
        d = pos - O
        n = np.array([-d[1], d[0]])
        n = n / np.linalg.norm(n) * 0.18
        lpos = mid - n
    ax.text(lpos[0], lpos[1], label, fontsize=6, ha="center",
            va="center", color=c_vec,
            bbox=dict(facecolor="white", edgecolor="none", pad=0.5,
                      alpha=0.8))

# --- Relative vectors (red, solid, prominent) ---
# Arrows start/end at body centers (no shrink)
# R = R_S - R_E (geocentric Sun)
draw_arrow(E, S, c_rel, lw=1.0, shrinkA=0, shrinkB=0)
mid_R = 0.5 * (E + S)
ax.text(mid_R[0] - 0.05, mid_R[1] - 0.25, r"$\mathbf{R}$",
        fontsize=7.5, ha="center", color=c_rel, fontweight="bold")

# D = R_M - R_E (geocentric Moon)
draw_arrow(E, M, c_rel, lw=1.0, shrinkA=0, shrinkB=0)
mid_D = 0.5 * (E + M)
ax.text(mid_D[0] + 0.0, mid_D[1] + 0.2, r"$\mathbf{D}$",
        fontsize=7.5, ha="center", color=c_rel, fontweight="bold")

# --- Spring between Earth and test mass ---
# Spring from Earth's surface to test mass surface (physically connected)
draw_spring(_spring_start, _spring_end, n_coils=6, amplitude=0.08,
            color="#555555", lw=0.8)

# r = x - R_E (test mass relative to Earth) — center to center
draw_arrow(E, T, c_rel, lw=0.8, shrinkA=0, shrinkB=0, zorder=4)
d_ET = T - E
n_ET_unit = np.array([-d_ET[1], d_ET[0]]) / np.linalg.norm(d_ET)
mid_r = 0.5 * (E + T) + 0.2 * n_ET_unit
ax.text(mid_r[0], mid_r[1], r"$\mathbf{r}$",
        fontsize=7.5, ha="center", va="center", color=c_rel,
        fontweight="bold")

# Set limits with padding — expanded to include test mass below Earth
ax.set_xlim(-4.0, 4.1)
ax.set_ylim(-2.2, 2.4)

fig.tight_layout()
fig.savefig("/home/xeal/dev/pytheas/doc/figures/setup_schematic.pdf")
fig.savefig("/home/xeal/dev/pytheas/doc/figures/setup_schematic.png")
print("Saved setup_schematic.pdf and .png")

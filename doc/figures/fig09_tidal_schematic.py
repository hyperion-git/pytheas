"""Figure 9: Tidal acceleration from first principles — conceptual schematic."""

import sys
sys.path.insert(0, '/home/xeal/dev/pytheas/doc/figures')

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, Circle, Arc
from style import apply_style, COLORS

apply_style()

fig, ax = plt.subplots(figsize=(7, 4))
ax.axis('off')
ax.set_xlim(-3.5, 7.5)
ax.set_ylim(-2.8, 2.8)
ax.set_aspect('equal')

# ============================================================
# Earth and Moon positions
# ============================================================
earth_center = np.array([0.0, 0.0])
earth_radius = 1.4
moon_center = np.array([5.8, 0.0])
moon_radius = 0.38

# --- Draw Earth ---
earth = Circle(earth_center, earth_radius, linewidth=1.5,
               edgecolor=COLORS['black'], facecolor='#c8dff0', alpha=0.4,
               zorder=2)
ax.add_patch(earth)
ax.text(earth_center[0], -earth_radius - 0.30, 'Earth',
        ha='center', va='top', fontsize=9, fontweight='bold',
        color=COLORS['black'])

# --- Draw Moon ---
moon = Circle(moon_center, moon_radius, linewidth=1.2,
              edgecolor=COLORS['black'], facecolor='#e0d8c8', alpha=0.5,
              zorder=2)
ax.add_patch(moon)
ax.text(moon_center[0], -moon_radius - 0.25, 'Moon',
        ha='center', va='top', fontsize=9, fontweight='bold',
        color=COLORS['black'])

# ============================================================
# Observer point (Moon-facing side of Earth)
# ============================================================
obs_angle = 0  # on the Moon-facing equator
obs_pos = earth_center + earth_radius * np.array([np.cos(obs_angle),
                                                   np.sin(obs_angle)])
ax.plot(*obs_pos, 'o', color=COLORS['black'], markersize=5, zorder=10)
ax.text(obs_pos[0] + 0.08, obs_pos[1] + 0.22, 'P',
        fontsize=9, fontweight='bold', color=COLORS['black'],
        ha='left', va='bottom')

# ============================================================
# Three key vectors
# ============================================================
arrow_scale = 1.0

# 1) g_obs: gravitational pull of Moon on observer (long arrow toward Moon)
g_obs_len = 2.2 * arrow_scale
g_obs_tip = obs_pos + np.array([g_obs_len, 0])
ax.annotate('', xy=g_obs_tip, xytext=obs_pos,
            arrowprops=dict(arrowstyle='->', color=COLORS['navy'],
                           lw=2.0, shrinkA=3, shrinkB=0),
            zorder=8)
ax.text(g_obs_tip[0] + 0.05, g_obs_tip[1] + 0.18,
        r'$\mathbf{g}_{\mathrm{obs}}$',
        color=COLORS['navy'], fontsize=10, fontweight='bold',
        ha='left', va='bottom')

# 2) g_center: gravitational pull of Moon on Earth's center (slightly shorter)
g_cen_len = 1.7 * arrow_scale
g_cen_tip = earth_center + np.array([g_cen_len, 0])
ax.annotate('', xy=g_cen_tip, xytext=earth_center,
            arrowprops=dict(arrowstyle='->', color=COLORS['orange'],
                           lw=2.0, shrinkA=3, shrinkB=0),
            zorder=8)
ax.text(g_cen_tip[0] + 0.05, g_cen_tip[1] - 0.22,
        r'$\mathbf{g}_{\mathrm{center}}$',
        color=COLORS['orange'], fontsize=10, fontweight='bold',
        ha='left', va='top')

# 3) a_tidal = g_obs - g_center (difference at observer)
# The tidal acceleration at the sub-lunar point points toward the Moon
a_tidal_len = 0.6 * arrow_scale
a_tidal_start = obs_pos + np.array([0.0, -0.55])
a_tidal_tip = a_tidal_start + np.array([a_tidal_len, 0])
ax.annotate('', xy=a_tidal_tip, xytext=a_tidal_start,
            arrowprops=dict(arrowstyle='->', color=COLORS['green'],
                           lw=3.0, shrinkA=0, shrinkB=0),
            zorder=9)
ax.text((a_tidal_start[0] + a_tidal_tip[0]) / 2, a_tidal_start[1] - 0.22,
        r'$\mathbf{a}_{\mathrm{tidal}}$',
        color=COLORS['green'], fontsize=10, fontweight='bold',
        ha='center', va='top')

# ============================================================
# Equation label
# ============================================================
ax.text(3.0, -2.2,
        r'$\mathbf{a}_{\mathrm{tidal}} = GM\!\left[\dfrac{\mathbf{R}-\mathbf{r}}'
        r'{|\mathbf{R}-\mathbf{r}|^3} - \dfrac{\mathbf{R}}{|\mathbf{R}|^3}\right]$',
        fontsize=10, ha='center', va='center',
        color=COLORS['black'],
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#f5f5dc',
                  edgecolor='gray', alpha=0.7))

# ============================================================
# Tidal bulge pattern — small arrows at several points on Earth
# ============================================================
n_arrows = 12
bulge_arrow_len = 0.35

for i in range(n_arrows):
    angle = 2 * np.pi * i / n_arrows
    # Position on Earth's surface
    pos = earth_center + earth_radius * np.array([np.cos(angle), np.sin(angle)])

    # Tidal acceleration direction at this point:
    # Toward/away from Moon at 0/pi; inward at pi/2, 3pi/2
    # The tidal field is: a_tidal ~ (2*cos(angle), -sin(angle)) in the Moon direction
    # This gives the classic quadrupole pattern
    a_r = 2 * np.cos(angle)   # radial component (along Earth-Moon line)
    a_t = -np.sin(angle)       # transverse component

    # Convert to Cartesian (Moon is in +x direction)
    a_x = a_r   # along Earth-Moon axis
    a_y = a_t   # perpendicular

    # Normalize and scale
    a_mag = np.sqrt(a_x**2 + a_y**2)
    if a_mag > 0:
        a_x /= a_mag
        a_y /= a_mag

    tip = pos + bulge_arrow_len * np.array([a_x, a_y])

    ax.annotate('', xy=tip, xytext=pos,
                arrowprops=dict(arrowstyle='->', color=COLORS['green'],
                               lw=1.2, alpha=0.7, shrinkA=0, shrinkB=0),
                zorder=7)

# ============================================================
# Labels for the tidal pattern
# ============================================================
ax.text(earth_center[0] + earth_radius + 0.55, 0.50,
        'stretch', color=COLORS['green'], fontsize=7, style='italic',
        ha='center', va='center', rotation=0)
ax.text(earth_center[0] - earth_radius - 0.55, 0.50,
        'stretch', color=COLORS['green'], fontsize=7, style='italic',
        ha='center', va='center', rotation=0)
ax.text(earth_center[0], earth_radius + 0.50,
        'squeeze', color=COLORS['green'], fontsize=7, style='italic',
        ha='center', va='center')
ax.text(earth_center[0], -earth_radius - 0.50,
        'squeeze', color=COLORS['green'], fontsize=7, style='italic',
        ha='center', va='top')

# ============================================================
# Earth center dot
# ============================================================
ax.plot(*earth_center, 'o', color=COLORS['black'], markersize=3, zorder=10)

# ============================================================
# Vector labels R and r
# ============================================================
# R vector: Earth center to Moon center
mid_R = (earth_center + moon_center) / 2
ax.annotate('', xy=moon_center, xytext=earth_center,
            arrowprops=dict(arrowstyle='->', color='gray', lw=1.0,
                           linestyle='--', shrinkA=5, shrinkB=5),
            zorder=3)
ax.text(mid_R[0], mid_R[1] + 0.25, r'$\mathbf{R}$',
        color='gray', fontsize=10, fontweight='bold',
        ha='center', va='bottom')

# r vector: Earth center to observer
mid_r = (earth_center + obs_pos) / 2
ax.text(mid_r[0], mid_r[1] + 0.18, r'$\mathbf{r}$',
        color='gray', fontsize=10, fontweight='bold',
        ha='center', va='bottom')

# ============================================================
# Title
# ============================================================
ax.set_title('Tidal Acceleration from First Principles',
             fontsize=11, fontweight='bold', pad=12)

# ============================================================
# Legend
# ============================================================
legend_items = [
    (COLORS['navy'], r'$\mathbf{g}_{\mathrm{obs}}$: Moon pull on observer'),
    (COLORS['orange'], r'$\mathbf{g}_{\mathrm{center}}$: Moon pull on Earth center'),
    (COLORS['green'], r'$\mathbf{a}_{\mathrm{tidal}} = \mathbf{g}_{\mathrm{obs}} - \mathbf{g}_{\mathrm{center}}$'),
]
lx, ly = -3.2, 2.4
for i, (color, label) in enumerate(legend_items):
    ax.plot(lx, ly - i * 0.35, 's', color=color, markersize=6)
    ax.text(lx + 0.18, ly - i * 0.35, label,
            color=COLORS['black'], fontsize=7.5, va='center', ha='left')

plt.savefig('/home/xeal/dev/pytheas/doc/figures/fig09_tidal_schematic.png',
            dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print("Figure 9 saved.")

"""Figure 2: Three reference frames (ECI, ECEF, ENU) schematic."""

import sys
sys.path.insert(0, '/home/xeal/dev/pytheas/doc/figures')

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, Arc, Ellipse
from style import apply_style, COLORS

apply_style()

fig, ax = plt.subplots(figsize=(6, 5))
ax.axis('off')
ax.set_xlim(-3.2, 3.2)
ax.set_ylim(-3.2, 3.5)
ax.set_aspect('equal')

# --- Earth as oblate ellipse (meridional cross-section) ---
earth = Ellipse((0, 0), width=3.6, height=3.4, linewidth=1.5,
                edgecolor=COLORS['black'], facecolor='#d0e8f0', alpha=0.35,
                zorder=1)
ax.add_patch(earth)

# --- 1) ECI axes (navy) ---
eci_color = COLORS['navy']
eci_len = 2.8

# ECI angle: tilted slightly to distinguish from ECEF
eci_angle_deg = 15
eci_angle = np.radians(eci_angle_deg)

# X_ECI
dx_eci = eci_len * np.cos(eci_angle)
dy_eci = eci_len * np.sin(eci_angle)
ax.annotate('', xy=(dx_eci, dy_eci), xytext=(0, 0),
            arrowprops=dict(arrowstyle='->', color=eci_color, lw=1.8),
            zorder=5)
ax.text(dx_eci + 0.12, dy_eci + 0.12, r'$X_{\mathrm{ECI}}$',
        color=eci_color, fontsize=10, fontweight='bold', ha='left', va='bottom')
# Small star symbol at tip direction
ax.text(dx_eci + 0.12, dy_eci - 0.18, r'(vernal equinox $\gamma$)',
        color=eci_color, fontsize=7, ha='left', va='top', style='italic')

# Y_ECI (perpendicular)
dx_eci_y = eci_len * np.cos(eci_angle + np.pi/2)
dy_eci_y = eci_len * np.sin(eci_angle + np.pi/2)
ax.annotate('', xy=(dx_eci_y, dy_eci_y), xytext=(0, 0),
            arrowprops=dict(arrowstyle='->', color=eci_color, lw=1.8),
            zorder=5)
ax.text(dx_eci_y - 0.15, dy_eci_y + 0.12, r'$Y_{\mathrm{ECI}}$',
        color=eci_color, fontsize=10, fontweight='bold', ha='right', va='bottom')

# --- 2) ECEF axes (orange) ---
ecef_color = COLORS['orange']
ecef_len = 2.5

# ECEF rotated by GMST from ECI
gmst_deg = 55  # angle between ECI and ECEF
ecef_angle = np.radians(eci_angle_deg + gmst_deg)

# X_ECEF
dx_ecef = ecef_len * np.cos(ecef_angle)
dy_ecef = ecef_len * np.sin(ecef_angle)
ax.annotate('', xy=(dx_ecef, dy_ecef), xytext=(0, 0),
            arrowprops=dict(arrowstyle='->', color=ecef_color, lw=1.8),
            zorder=5)
ax.text(dx_ecef + 0.08, dy_ecef + 0.12, r'$X_{\mathrm{ECEF}}$',
        color=ecef_color, fontsize=10, fontweight='bold', ha='left', va='bottom')

# Y_ECEF (perpendicular)
dx_ecef_y = ecef_len * np.cos(ecef_angle + np.pi/2)
dy_ecef_y = ecef_len * np.sin(ecef_angle + np.pi/2)
ax.annotate('', xy=(dx_ecef_y, dy_ecef_y), xytext=(0, 0),
            arrowprops=dict(arrowstyle='->', color=ecef_color, lw=1.8),
            zorder=5)
ax.text(dx_ecef_y - 0.12, dy_ecef_y + 0.08, r'$Y_{\mathrm{ECEF}}$',
        color=ecef_color, fontsize=10, fontweight='bold', ha='right', va='bottom')

# --- 3) GMST arc between X_ECI and X_ECEF ---
arc_radius = 1.2
gmst_arc = Arc((0, 0), 2*arc_radius, 2*arc_radius,
               angle=0, theta1=eci_angle_deg, theta2=eci_angle_deg + gmst_deg,
               color=COLORS['red'], lw=1.5, linestyle='--', zorder=6)
ax.add_patch(gmst_arc)

# Label the arc
mid_angle = np.radians(eci_angle_deg + gmst_deg / 2)
label_r = arc_radius + 0.25
ax.text(label_r * np.cos(mid_angle), label_r * np.sin(mid_angle),
        r'GMST ($\theta$)', color=COLORS['red'], fontsize=9,
        ha='center', va='center', fontweight='bold')

# Small arrowhead on the arc (draw a tiny arrow near the end)
arc_end_angle = np.radians(eci_angle_deg + gmst_deg - 3)
arc_tip_angle = np.radians(eci_angle_deg + gmst_deg)
ax.annotate('',
            xy=(arc_radius * np.cos(arc_tip_angle),
                arc_radius * np.sin(arc_tip_angle)),
            xytext=(arc_radius * np.cos(arc_end_angle),
                    arc_radius * np.sin(arc_end_angle)),
            arrowprops=dict(arrowstyle='->', color=COLORS['red'], lw=1.5),
            zorder=6)

# --- 4) Observer point on Earth's surface ---
# Place observer at latitude ~45 deg on the ECEF X-axis side
obs_lat_deg = 45
obs_angle_from_ecef = 20  # angular offset from X_ECEF in the ECEF frame
obs_angle_total = np.radians(eci_angle_deg + gmst_deg + obs_angle_from_ecef)

# Observer on the ellipse surface
a_earth, b_earth = 1.8, 1.7  # semi-axes
obs_r_angle = obs_angle_total
obs_x = a_earth * np.cos(obs_r_angle)
obs_y = b_earth * np.sin(obs_r_angle)

# Draw observer point
ax.plot(obs_x, obs_y, 'o', color=COLORS['black'], markersize=6, zorder=10)
ax.text(obs_x + 0.05, obs_y + 0.25, 'Observer', color=COLORS['black'],
        fontsize=9, fontweight='bold', ha='center', va='bottom')

# --- ENU axes at observer (green) ---
enu_color = COLORS['green']
enu_len = 0.7

# "Up" direction: radial outward from Earth center
up_angle = np.arctan2(obs_y, obs_x)
up_dx = enu_len * np.cos(up_angle)
up_dy = enu_len * np.sin(up_angle)

# "East" direction: perpendicular to Up, counterclockwise (tangent)
east_angle = up_angle + np.pi/2
east_dx = enu_len * 0.7 * np.cos(east_angle)
east_dy = enu_len * 0.7 * np.sin(east_angle)

# "North" direction: for a 2D meridional view, North is perpendicular
# In this equatorial-plane view, North points out of the plane (we show it
# going somewhat upward to suggest 3D)
north_angle = up_angle + np.pi/4  # compromise direction for visual clarity
north_dx = enu_len * 0.55 * np.cos(north_angle)
north_dy = enu_len * 0.55 * np.sin(north_angle)

# Draw ENU arrows
# Up
ax.annotate('', xy=(obs_x + up_dx, obs_y + up_dy), xytext=(obs_x, obs_y),
            arrowprops=dict(arrowstyle='->', color=enu_color, lw=2.0),
            zorder=10)
ax.text(obs_x + up_dx + 0.08, obs_y + up_dy + 0.08, 'U',
        color=enu_color, fontsize=9, fontweight='bold', ha='left', va='bottom')

# East
ax.annotate('', xy=(obs_x + east_dx, obs_y + east_dy), xytext=(obs_x, obs_y),
            arrowprops=dict(arrowstyle='->', color=enu_color, lw=2.0),
            zorder=10)
ax.text(obs_x + east_dx - 0.05, obs_y + east_dy + 0.12, 'E',
        color=enu_color, fontsize=9, fontweight='bold', ha='center', va='bottom')

# North (shown slightly differently to suggest out-of-plane)
ax.annotate('', xy=(obs_x + north_dx, obs_y + north_dy), xytext=(obs_x, obs_y),
            arrowprops=dict(arrowstyle='->', color=enu_color, lw=1.5,
                           linestyle='dashed'),
            zorder=10)
ax.text(obs_x + north_dx, obs_y + north_dy + 0.12, r'N ($\perp$)',
        color=enu_color, fontsize=8, fontweight='bold', ha='center', va='bottom',
        style='italic')

# --- Latitude arc (from equator/X_ECEF to observer) ---
lat_arc_r = 0.9
# Longitude: angle from X_ECEF to observer in ECEF frame
lon_angle_start = np.degrees(ecef_angle)
lon_angle_end = np.degrees(obs_r_angle)
lon_arc = Arc((0, 0), 2*lat_arc_r, 2*lat_arc_r,
              angle=0, theta1=lon_angle_start, theta2=lon_angle_end,
              color=COLORS['pink'], lw=1.2, zorder=6)
ax.add_patch(lon_arc)

# Label
lon_mid = np.radians((lon_angle_start + lon_angle_end) / 2)
ax.text(0.65 * np.cos(lon_mid), 0.65 * np.sin(lon_mid), r'$\lambda$',
        color=COLORS['pink'], fontsize=11, fontweight='bold',
        ha='center', va='center')

# --- Latitude angle indicator ---
# Draw a dashed line from center to observer
ax.plot([0, obs_x * 1.1], [0, obs_y * 1.1], '--', color='gray',
        lw=0.6, alpha=0.5, zorder=2)

# --- Small dot at Earth center ---
ax.plot(0, 0, 'o', color=COLORS['black'], markersize=3, zorder=8)

# --- Title ---
ax.set_title('Three Reference Frames: ECI, ECEF, and ENU',
             fontsize=11, fontweight='bold', pad=12)

# --- Legend-like labels ---
legend_x, legend_y = -3.0, -2.5
legend_items = [
    (COLORS['navy'], 'ECI (Earth-Centered Inertial)'),
    (COLORS['orange'], 'ECEF (Earth-Centered Earth-Fixed)'),
    (COLORS['green'], 'ENU (East-North-Up, local)'),
]
for i, (color, label) in enumerate(legend_items):
    ax.plot(legend_x, legend_y - i * 0.32, 's', color=color, markersize=7)
    ax.text(legend_x + 0.18, legend_y - i * 0.32, label,
            color=color, fontsize=8, va='center', ha='left')

plt.savefig('/home/xeal/dev/pytheas/doc/figures/fig02_reference_frames.png',
            dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print("Figure 2 saved.")

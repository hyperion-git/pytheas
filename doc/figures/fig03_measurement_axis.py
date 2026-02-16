"""Figure 3: Measurement axis geometry — ENU frame with sensor axis n-hat."""

import sys
sys.path.insert(0, '/home/xeal/dev/pytheas/doc/figures')

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Arc
from style import apply_style, COLORS

apply_style()

fig, ax = plt.subplots(figsize=(4, 4))
ax.axis('off')
ax.set_xlim(-1.6, 2.2)
ax.set_ylim(-0.8, 2.8)
ax.set_aspect('equal')

# Origin of the ENU frame
O = np.array([0.0, 0.0])

# --- Axis lengths ---
enu_len = 1.8
nhat_len = 2.0

# --- Define 3D-projected directions ---
# We project a 3D ENU frame onto 2D for the schematic.
# Up: straight up
# North: upper-right (tilted to show perspective)
# East: to the right and slightly down (foreshortened)
up_dir = np.array([0.0, 1.0])
north_dir = np.array([0.45, 0.70])
north_dir = north_dir / np.linalg.norm(north_dir)
east_dir = np.array([0.85, -0.20])
east_dir = east_dir / np.linalg.norm(east_dir)

# --- Draw ENU axes ---
arrow_kw = dict(arrowstyle='->', lw=1.8)

# East (blue)
east_tip = O + enu_len * 0.9 * east_dir
ax.annotate('', xy=east_tip, xytext=O,
            arrowprops=dict(**arrow_kw, color=COLORS['blue']), zorder=5)
ax.text(east_tip[0] + 0.08, east_tip[1] - 0.12, 'E (East)',
        color=COLORS['blue'], fontsize=9, fontweight='bold',
        ha='left', va='top')

# North (orange)
north_tip = O + enu_len * 0.85 * north_dir
ax.annotate('', xy=north_tip, xytext=O,
            arrowprops=dict(**arrow_kw, color=COLORS['orange']), zorder=5)
ax.text(north_tip[0] + 0.10, north_tip[1] + 0.05, 'N (North)',
        color=COLORS['orange'], fontsize=9, fontweight='bold',
        ha='left', va='bottom')

# Up (navy)
up_tip = O + enu_len * up_dir
ax.annotate('', xy=up_tip, xytext=O,
            arrowprops=dict(**arrow_kw, color=COLORS['navy']), zorder=5)
ax.text(up_tip[0] + 0.08, up_tip[1] + 0.02, 'U (Up)',
        color=COLORS['navy'], fontsize=9, fontweight='bold',
        ha='left', va='bottom')

# --- Sensor axis n-hat ---
# n-hat at zenith angle theta_z from Up and azimuth alpha from North
theta_z_deg = 35   # zenith angle
alpha_deg = 40     # azimuth from North toward East

# In the projected view:
# The horizontal-plane projection of n-hat combines North and East components.
# n-hat = sin(theta_z)*[sin(alpha)*East + cos(alpha)*North] + cos(theta_z)*Up
theta_z = np.radians(theta_z_deg)
alpha = np.radians(alpha_deg)

nhat_dir = (np.sin(theta_z) * (np.sin(alpha) * east_dir + np.cos(alpha) * north_dir)
            + np.cos(theta_z) * up_dir)
nhat_dir = nhat_dir / np.linalg.norm(nhat_dir)

nhat_tip = O + nhat_len * nhat_dir
ax.annotate('', xy=nhat_tip, xytext=O,
            arrowprops=dict(arrowstyle='->', color=COLORS['red'], lw=2.2),
            zorder=8)
ax.text(nhat_tip[0] + 0.10, nhat_tip[1] + 0.02,
        r'$\hat{\mathbf{n}}$ (sensor)',
        color=COLORS['red'], fontsize=10, fontweight='bold',
        ha='left', va='bottom')

# --- Projection of n-hat onto horizontal plane (dashed) ---
nhat_horiz = np.sin(theta_z) * (np.sin(alpha) * east_dir + np.cos(alpha) * north_dir)
nhat_horiz_norm = np.linalg.norm(nhat_horiz)
nhat_horiz_dir = nhat_horiz / nhat_horiz_norm if nhat_horiz_norm > 0 else north_dir

proj_len = nhat_len * nhat_horiz_norm / np.linalg.norm(
    np.sin(theta_z) * (np.sin(alpha) * east_dir + np.cos(alpha) * north_dir)
    + np.cos(theta_z) * up_dir)
proj_tip = O + proj_len * nhat_horiz_dir

# Dashed line from origin to projection
ax.plot([O[0], proj_tip[0]], [O[1], proj_tip[1]], '--',
        color=COLORS['red'], lw=1.2, alpha=0.6, zorder=4)

# Vertical dashed line from n-hat tip down to projection
ax.plot([nhat_tip[0], proj_tip[0]], [nhat_tip[1], proj_tip[1]], ':',
        color='gray', lw=1.0, alpha=0.5, zorder=3)

# --- Zenith angle arc (theta_z): between Up and n-hat ---
# Draw arc from Up direction to n-hat direction
arc_r = 0.65
up_angle_deg = 90  # Up is straight up
nhat_angle_deg = np.degrees(np.arctan2(nhat_dir[1], nhat_dir[0]))

# The arc goes from n-hat to Up (counterclockwise)
theta1 = min(nhat_angle_deg, up_angle_deg)
theta2 = max(nhat_angle_deg, up_angle_deg)
zenith_arc = Arc(O, 2*arc_r, 2*arc_r, angle=0,
                 theta1=theta1, theta2=theta2,
                 color=COLORS['navy'], lw=1.3, linestyle='-', zorder=6)
ax.add_patch(zenith_arc)

# Label theta_z
mid_zen = np.radians((theta1 + theta2) / 2)
ax.text(arc_r * 0.75 * np.cos(mid_zen) - 0.18,
        arc_r * 0.75 * np.sin(mid_zen),
        r'$\theta_z$', color=COLORS['navy'], fontsize=11,
        fontweight='bold', ha='center', va='center')

# --- Azimuth arc (alpha): from North to projection in horizontal plane ---
arc_r2 = 0.50
north_angle_deg = np.degrees(np.arctan2(north_dir[1], north_dir[0]))
proj_angle_deg = np.degrees(np.arctan2(nhat_horiz_dir[1], nhat_horiz_dir[0]))

# Ensure proper arc direction (North to projection clockwise in real world,
# but in our projected view we go from north to projection)
t1_az = min(proj_angle_deg, north_angle_deg)
t2_az = max(proj_angle_deg, north_angle_deg)
azimuth_arc = Arc(O, 2*arc_r2, 2*arc_r2, angle=0,
                  theta1=t1_az, theta2=t2_az,
                  color=COLORS['orange'], lw=1.3, linestyle='-', zorder=6)
ax.add_patch(azimuth_arc)

# Label alpha
mid_az = np.radians((t1_az + t2_az) / 2)
ax.text(arc_r2 * 0.6 * np.cos(mid_az) + 0.15,
        arc_r2 * 0.6 * np.sin(mid_az) - 0.02,
        r'$\alpha$', color=COLORS['orange'], fontsize=11,
        fontweight='bold', ha='center', va='center')

# --- Origin dot ---
ax.plot(*O, 'o', color=COLORS['black'], markersize=5, zorder=10)

# --- Title ---
ax.set_title('Measurement Axis Geometry', fontsize=10, fontweight='bold', pad=10)

plt.savefig('/home/xeal/dev/pytheas/doc/figures/fig03_measurement_axis.png',
            dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print("Figure 3 saved.")

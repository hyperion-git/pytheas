#!/usr/bin/env python3
"""Figure 6: Lunar distance and tidal amplitude correlation.

Two vertically stacked panels sharing x-axis, 30 days of data starting
2025-03-01.  Panel (a) shows Earth--Moon distance in km; panel (b) shows
the absolute value of lunar tidal acceleration at Munich in um/s^2.
Perigee and apogee are marked with vertical dashed lines.
"""

import sys
sys.path.insert(0, "/home/xeal/dev/pytheas/doc/figures")
sys.path.insert(0, "/home/xeal/dev/pytheas")

import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from style import apply_style, COLORS, add_dual_axis

apply_style()

import pytheas

# ---- Parameters ----
start = datetime(2025, 3, 1)
n_days = 30
dt_hours = 6  # 6-hour intervals
n_points = n_days * 24 // dt_hours + 1

# Observer: Munich
lat, lon, alt = 48.14, 11.58, 500.0
r_obs = pytheas.geodetic_to_ecef(lat, lon, alt)

# ---- Compute time series ----
times = [start + timedelta(hours=i * dt_hours) for i in range(n_points)]

distances_km = np.empty(n_points)
tidal_mag_ums2 = np.empty(n_points)  # microm/s^2

for i, t in enumerate(times):
    R_moon = pytheas.moon_position_ecef(t)
    distances_km[i] = np.linalg.norm(R_moon) / 1e3  # m -> km

    a_tidal = pytheas.tidal_acceleration(r_obs, R_moon, pytheas.GM_MOON)
    tidal_mag_ums2[i] = np.linalg.norm(a_tidal) * 1e6  # m/s^2 -> um/s^2

# ---- Find perigee and apogee ----
i_perigee = np.argmin(distances_km)
i_apogee = np.argmax(distances_km)

t_perigee = times[i_perigee]
t_apogee = times[i_apogee]

# ---- Plot ----
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 4.5), sharex=True,
                                gridspec_kw={'hspace': 0.08})

# Panel (a): Distance
ax1.plot(times, distances_km, color=COLORS['navy'], linewidth=1.5)
ax1.set_ylabel("Earth\u2013Moon distance (km)")
ax1.text(0.02, 0.92, "(a)", transform=ax1.transAxes,
         fontweight='bold', fontsize=10, va='top')

# Panel (b): Tidal acceleration magnitude
ax2.plot(times, tidal_mag_ums2, color=COLORS['red'], linewidth=1.5)
ax2.set_ylabel("Lunar tidal accel. (\u00b5m/s\u00b2)")
ax2.set_xlabel("Date (2025)")

ax2.text(0.02, 0.92, "(b)", transform=ax2.transAxes,
         fontweight='bold', fontsize=10, va='top')

# ---- Mark perigee and apogee ----
for ax in (ax1, ax2):
    ax.axvline(t_perigee, color=COLORS['black'], ls='--', lw=0.8, alpha=0.6)
    ax.axvline(t_apogee, color=COLORS['black'], ls='--', lw=0.8, alpha=0.6)

# Annotate at the top of panel (a)
y_top = ax1.get_ylim()[1]
ax1.annotate("perigee", xy=(t_perigee, y_top), xytext=(0, 4),
             textcoords='offset points', ha='center', va='bottom',
             fontsize=8, fontstyle='italic', annotation_clip=False)
ax1.annotate("apogee", xy=(t_apogee, y_top), xytext=(0, 4),
             textcoords='offset points', ha='center', va='bottom',
             fontsize=8, fontstyle='italic', annotation_clip=False)

# ---- Format x-axis ----
ax2.xaxis.set_major_locator(mdates.WeekdayLocator(byweekday=mdates.MO))
ax2.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
fig.autofmt_xdate(rotation=30, ha='right')

# Thousands separator on y-axis of panel (a)
ax1.yaxis.set_major_formatter(plt.FuncFormatter(
    lambda x, _: f"{x:,.0f}"))

add_dual_axis(ax2, 'µm/s²')

plt.savefig("/home/xeal/dev/pytheas/doc/figures/fig06_lunar_distance_tides.png")
plt.close()
print("Saved fig06_lunar_distance_tides.png")

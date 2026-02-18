"""CLI entry point: python -m pytheas"""

import argparse
import sys
import numpy as np
from datetime import datetime, timedelta

from . import __version__, compute_timeseries


def main():
    p = argparse.ArgumentParser(
        prog='pytheas',
        description='Compute g(t) at a point on Earth (accuracy ~10-100 nGal)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  pytheas --lat 48.14 --lon 11.58 --alt 500
  pytheas --lat 48.14 --lon 11.58 --alt 500 --start 2025-03-20 --hours 72
  pytheas --lat 48.14 --lon 11.58 --alt 500 --zenith 90 --azimuth 0
  pytheas --lat 48.14 --lon 11.58 --alt 500 --csv output.csv --plot
""")
    p.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    p.add_argument('--lat',     type=float, required=True, help='Latitude (deg)')
    p.add_argument('--lon',     type=float, required=True, help='Longitude (deg)')
    p.add_argument('--alt',     type=float, default=0.0,   help='Altitude (m, default 0)')
    p.add_argument('--zenith',  type=float, default=0.0,
                   help='Zenith angle of measurement axis (0=vertical, default 0)')
    p.add_argument('--azimuth', type=float, default=0.0,
                   help='Azimuth of measurement axis (deg from N, default 0)')
    p.add_argument('--start',   type=str,   default=None,
                   help='Start time UTC (YYYY-MM-DD or YYYY-MM-DDTHH:MM)')
    p.add_argument('--hours',   type=float, default=48.0,  help='Duration (hours, default 48)')
    p.add_argument('--interval', type=float, default=10.0, help='Cadence (minutes, default 10)')
    p.add_argument('--csv',     type=str,   default=None,  help='Output CSV file')
    p.add_argument('--plot',    action='store_true',        help='Show plot')
    args = p.parse_args()

    if args.start:
        for fmt in ('%Y-%m-%dT%H:%M', '%Y-%m-%d'):
            try:
                start = datetime.strptime(args.start, fmt)
                break
            except ValueError:
                continue
        else:
            print(f"Cannot parse start time: {args.start}", file=sys.stderr)
            sys.exit(1)
    else:
        now = datetime.utcnow()
        start = now.replace(minute=0, second=0, microsecond=0)

    end = start + timedelta(hours=args.hours)

    # Approximate local time offset from longitude (standard timezone)
    utc_offset_h = round(args.lon / 15.0)
    local_dt = timedelta(hours=utc_offset_h)
    sign = '+' if utc_offset_h >= 0 else '-'
    utc_label = f"UTC{sign}{abs(utc_offset_h):d}"

    print(f"Location : {args.lat:.4f} N, {args.lon:.4f} E, {args.alt:.1f} m")
    print(f"Axis     : zenith = {args.zenith:.1f} deg, "
          f"azimuth = {args.azimuth:.1f} deg")
    print(f"Window   : {start.isoformat()} to {end.isoformat()} UTC "
          f"({args.hours:.1f} h, {args.interval:.0f} min)")
    print(f"Local tz : {utc_label} (from longitude)")
    print()

    data = compute_timeseries(start, end, args.lat, args.lon, args.alt,
                              args.zenith, args.azimuth, args.interval)

    g_s = data['g_static'][0]
    g_t = data['g_tidal']
    print(f"Static g on axis   : {g_s:.8f} m/s^2")
    print(f"Tidal peak-to-peak : {np.ptp(g_t) * 1e6:.4f} um/s^2")
    print(f"Tidal RMS          : {np.std(g_t) * 1e6:.4f} um/s^2")
    print(f"g_total range      : [{np.min(data['g_total']):.8f}, "
          f"{np.max(data['g_total']):.8f}] m/s^2")

    # Sample table
    times = data['times']
    head = min(5, len(times))
    tail = min(3, len(times))
    rows = list(range(head))
    if len(times) > head + tail:
        rows.append(None)
        rows.extend(range(len(times) - tail, len(times)))
    elif len(times) > head:
        rows.extend(range(head, len(times)))

    print(f"\n{'Time UTC':>20s}  {'Local (' + utc_label + ')':>20s}  {'g_total':>17s}  "
          f"{'tidal (um/s2)':>14s}  {'Moon':>10s}  {'Sun':>10s}")
    print('-' * 100)
    for idx in rows:
        if idx is None:
            print(f"{'...':>20s}")
            continue
        t = times[idx]
        t_local = t + local_dt
        print(f"{t.strftime('%Y-%m-%d %H:%M'):>20s}  "
              f"{t_local.strftime('%Y-%m-%d %H:%M'):>20s}  "
              f"{data['g_total'][idx]:17.10f}  "
              f"{data['g_tidal'][idx] * 1e6:14.4f}  "
              f"{data['g_tidal_moon'][idx] * 1e6:10.4f}  "
              f"{data['g_tidal_sun'][idx] * 1e6:10.4f}")

    if args.csv:
        with open(args.csv, 'w') as f:
            f.write(f'time_utc,time_local_{utc_label},'
                    'g_total_m_s2,g_static_m_s2,g_tidal_m_s2,'
                    'g_tidal_moon_m_s2,g_tidal_sun_m_s2\n')
            for i, t in enumerate(times):
                t_local = t + local_dt
                f.write(f"{t.isoformat()},{t_local.isoformat()},"
                        f"{data['g_total'][i]:.12e},"
                        f"{data['g_static'][i]:.12e},"
                        f"{data['g_tidal'][i]:.12e},"
                        f"{data['g_tidal_moon'][i]:.12e},"
                        f"{data['g_tidal_sun'][i]:.12e}\n")
        print(f"\nSaved: {args.csv}")

    if args.plot:
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print('matplotlib not available, skipping plot', file=sys.stderr)
            return

        hours = np.array([(t - start).total_seconds() / 3600.0
                          for t in times])

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)

        ax1.plot(hours, data['g_total'], 'k-', lw=0.7)
        ax1.set_ylabel('g_total (m/s$^2$)')
        ax1.set_title(f'g(t) at ({args.lat:.2f}$\\degree$N, {args.lon:.2f}$\\degree$E, '
                      f'{args.alt:.0f} m)')
        ax1.ticklabel_format(useOffset=True, axis='y')
        ax1.grid(True, alpha=0.3)

        ax2.plot(hours, data['g_tidal_moon'] * 1e6, 'C0-', lw=0.7,
                 label='Moon', alpha=0.8)
        ax2.plot(hours, data['g_tidal_sun'] * 1e6, 'C1-', lw=0.7,
                 label='Sun', alpha=0.8)
        ax2.plot(hours, data['g_tidal'] * 1e6, 'k-', lw=0.9,
                 label='Total')
        ax2.set_ylabel('Tidal perturbation ($\\mu$m/s$^2$)')
        ax2.set_xlabel(f'Hours since {start.strftime("%Y-%m-%d %H:%M")} UTC')
        ax2.legend(loc='upper right', fontsize=8)
        ax2.axhline(0, color='gray', lw=0.5, ls='--')
        ax2.grid(True, alpha=0.3)

        fig.tight_layout()
        plt.show()


if __name__ == '__main__':
    main()

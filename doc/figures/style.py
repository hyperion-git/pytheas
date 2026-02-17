"""Shared matplotlib style for all tidal-model figures."""

import matplotlib as mpl
import matplotlib.pyplot as plt

# Okabe-Ito colorblind-safe palette
COLORS = {
    'orange':  '#E69F00',
    'blue':    '#56B4E9',
    'green':   '#009E73',
    'yellow':  '#F0E442',
    'navy':    '#0072B2',
    'red':     '#D55E00',
    'pink':    '#CC79A7',
    'black':   '#000000',
}
COLOR_LIST = list(COLORS.values())

def apply_style():
    """Apply publication-quality style for tidal-model figures."""
    mpl.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['DejaVu Sans', 'Arial', 'Helvetica'],
        'font.size': 9,
        'axes.labelsize': 10,
        'axes.titlesize': 10,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'legend.fontsize': 8,
        'figure.dpi': 150,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.prop_cycle': plt.cycler(color=COLOR_LIST),
        'lines.linewidth': 1.5,
        'axes.linewidth': 0.8,
        'xtick.major.width': 0.8,
        'ytick.major.width': 0.8,
        'figure.figsize': (6, 3.5),
    })


# Dual-unit conversion table: primary_unit -> (dual_label, factor)
# where dual_value = primary_value * factor
_DUAL_UNITS = {
    'm/s²':  ('Gal',    100),
    'µm/s²': ('µGal',   100),
    'mGal':  ('µm/s²',  10),
    'µGal':  ('nm/s²',  10),
    'nGal':  ('nm/s²',  0.01),
}


def add_dual_axis(ax, primary_unit, direction='y'):
    """Add a secondary axis showing the equivalent gravity unit.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes object that already has the primary unit.
    primary_unit : str
        One of 'm/s²', 'µm/s²', 'mGal', 'µGal', 'nGal'.
    direction : str
        'y' adds a right-hand axis; 'x' adds a top axis.
    """
    dual_label, factor = _DUAL_UNITS[primary_unit]

    if direction == 'y':
        sec = ax.secondary_yaxis(
            'right',
            functions=(lambda v: v * factor, lambda v: v / factor),
        )
        sec.set_ylabel(dual_label)
        sec.spines['right'].set_visible(True)
    else:
        sec = ax.secondary_xaxis(
            'top',
            functions=(lambda v: v * factor, lambda v: v / factor),
        )
        sec.set_xlabel(dual_label)
        sec.spines['top'].set_visible(True)

    sec.tick_params(labelsize=8)
    return sec

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

#!/usr/bin/env python3
"""Figure 7: Perturbation series convergence.

Horizontal bar chart showing the 10 largest longitude perturbation terms
from Meeus Table 47.A, sorted by absolute coefficient.  Positive
coefficients shown in navy, negative in red/orange.
"""

import sys
sys.path.insert(0, "/home/xeal/dev/pytheas/doc/figures")
sys.path.insert(0, "/home/xeal/dev/pytheas")

import numpy as np
import matplotlib.pyplot as plt
from style import apply_style, COLORS

apply_style()

# ---- Data: top 10 terms from Meeus Table 47.A ----
terms = [
    {"D": 0, "M": 0, "Mp": 1,  "F": 0, "coeff": 6_288_774,  "name": "Principal (M\u2032)"},
    {"D": 2, "M": 0, "Mp": -1, "F": 0, "coeff": 1_274_027,  "name": "Evection"},
    {"D": 2, "M": 0, "Mp": 0,  "F": 0, "coeff": 658_314,    "name": "Variation"},
    {"D": 0, "M": 0, "Mp": 2,  "F": 0, "coeff": 213_618,    "name": "2M\u2032"},
    {"D": 0, "M": 1, "Mp": 0,  "F": 0, "coeff": -185_116,   "name": "Annual eq."},
    {"D": 0, "M": 0, "Mp": 0,  "F": 2, "coeff": -114_332,   "name": "2F"},
    {"D": 2, "M": 0, "Mp": -2, "F": 0, "coeff": 58_793,     "name": "Evection 2"},
    {"D": 2, "M": -1, "Mp": -1, "F": 0, "coeff": 57_066,    "name": "(2D\u2212M\u2212M\u2032)"},
    {"D": 2, "M": 0, "Mp": 1,  "F": 0, "coeff": 53_322,     "name": "(2D+M\u2032)"},
    {"D": 2, "M": -1, "Mp": 0, "F": 0, "coeff": 45_758,     "name": "(2D\u2212M)"},
]

# Sort by absolute value (largest first -> will appear at top of horizontal bar)
terms.sort(key=lambda t: abs(t["coeff"]))

names = [t["name"] for t in terms]
coeffs = np.array([t["coeff"] for t in terms])
abs_coeffs = np.abs(coeffs)

# Colors: navy for positive, red/orange for negative
bar_colors = [COLORS['navy'] if c > 0 else COLORS['red'] for c in coeffs]

# ---- Plot ----
fig, ax = plt.subplots(figsize=(6, 3.8))

y_pos = np.arange(len(terms))
ax.barh(y_pos, abs_coeffs, color=bar_colors, height=0.65, edgecolor='none')

ax.set_yticks(y_pos)
ax.set_yticklabels(names)
ax.set_xlabel("Amplitude (\u00d710\u207b\u2076 degrees)")
ax.set_title("Convergence of the lunar longitude series", fontsize=10,
             fontweight='bold', pad=10)

# Add a legend for sign
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=COLORS['navy'], label='Positive'),
                   Patch(facecolor=COLORS['red'], label='Negative')]
ax.legend(handles=legend_elements, loc='lower right', frameon=False,
          fontsize=8)

# Format x-axis: all values in plain thousands
ax.xaxis.set_major_formatter(plt.FuncFormatter(
    lambda x, _: f"{x:,.0f}"))

plt.tight_layout()
plt.savefig("/home/xeal/dev/pytheas/doc/figures/fig07_perturbation_terms.png")
plt.close()
print("Saved fig07_perturbation_terms.png")

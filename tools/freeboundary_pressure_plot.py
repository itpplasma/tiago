#!/usr/bin/env python3
"""Plot the native synthetic pressure reference; no experimental data.

Usage: freeboundary_pressure_plot.py <pressure-reference.csv> <figure.png>
"""
import csv
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

source, output = map(Path, sys.argv[1:])
with source.open() as stream:
    rows = list(csv.DictReader(stream))
data = {key: np.array([float(row[key]) for row in rows]) for key in rows[0]}
fig, ax = plt.subplots(figsize=(7, 4), layout="constrained")
for key, label, color, style in [
    ("total_Pa", "Total pressure", "#0072B2", "-"),
    ("thermal_Pa", "Thermal pressure", "#000000", "--"),
    ("additional_Pa", "Bounded additional component", "#D55E00", "-."),
]:
    ax.plot(data["s"], data[key], label=label, color=color, linestyle=style,
            linewidth=2)
ax.set(xlabel="Normalized toroidal flux s [dimensionless]", ylabel="Pressure [Pa]",
       xlim=(0, 1), ylim=(0, None), title="Synthetic native pressure reference")
ax.text(.43, .25, "Amplitude c in [0, 2]; exponent K in [2, 8]\nLatents z_c = 0.3, z_K = -0.4",
        transform=ax.transAxes, fontsize=10)
ax.grid(alpha=.2)
ax.legend(loc="upper right", frameon=False)
fig.savefig(output, dpi=180)

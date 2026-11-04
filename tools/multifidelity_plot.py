#!/usr/bin/env python3
"""Plot native synthetic AR1 recovery against its planted truth.

Usage: multifidelity_plot.py <multifidelity-profile.csv> <figure.png>
Bars are local Gaussian 95% intervals, conditional on one planted profile.
The shared-profile lines are model constraints, not fitted old-code results.
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

fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True, layout="constrained")
for ax, value, sd, truth, old, label in [
    (axes[0], "C", "C_sd", "C_true", 1, "Temperature scale C [dimensionless]"),
    (axes[1], "delta_eV", "delta_sd_eV", "delta_true_eV", 0,
     "Discrepancy coefficient [eV]"),
]:
    ax.axhline(old, color="#D55E00", linestyle="--", linewidth=1.8,
               label="Shared-profile constraint")
    ax.plot(data["trial"], data[truth], color="black", linewidth=1.8,
            label="Planted truth")
    ax.errorbar(data["trial"], data[value], yerr=1.95996398454 * data[sd],
                fmt="o", markersize=3, linewidth=.7, color="#0072B2",
                label="AR1 fit and local 95% interval")
    ax.set(ylabel=label)
    ax.grid(axis="y", alpha=.2)
fig.legend(*axes[0].get_legend_handles_labels(), loc="outside upper center",
           ncol=3, fontsize=9)
axes[1].set(xlabel="Independent measurement-noise realization", xlim=(0, len(rows) + 1))
fig.savefig(output, dpi=180)

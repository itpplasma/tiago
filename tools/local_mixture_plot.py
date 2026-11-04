#!/usr/bin/env python3
"""Plot the synthetic abs(x) measurement from test_local_mixture.f90.

Usage: local_mixture_plot.py <test-output/two_mode> <figure.png>
All quantities are dimensionless. The reference multiplies the normalized
Gaussian prior and likelihood, then normalizes by numerical quadrature.
"""
import csv
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

root, output = Path(sys.argv[1]), Path(sys.argv[2])
with (root / "modes.csv").open() as stream:
    modes = list(csv.DictReader(stream))
assert len(modes) == 2
x = np.linspace(-8, 8, 16001)
def normal(x, mean, sd):
    return np.exp(-0.5 * ((x - mean) / sd)**2) / (sd * np.sqrt(2*np.pi))

exact = normal(np.abs(x), 6, .25) * normal(x, .5, 3)
exact /= np.trapezoid(exact, x)
components = []
for mode in modes:
    with (root / f"mode_{mode['mode']}_parameters.csv").open() as stream:
        parameter = next(csv.DictReader(stream))
    components.append((float(parameter['map']), float(parameter['sd']), float(mode['weight'])))
components.sort()
mixture = sum(weight * normal(x, mean, sd) for mean, sd, weight in components)
fig, axes = plt.subplots(1, 2, figsize=(10, 3.6), layout="constrained", width_ratios=(2, 1))
axes[0].plot(x, normal(x, components[1][0], components[1][1]), "--", color="#D55E00",
             label="One local Gaussian")
axes[0].plot(x, exact, color="black", linewidth=3, label="Analytic posterior")
axes[0].plot(x, mixture, ":", color="#0072B2", linewidth=2, label="Local mixture")
axes[0].set(xlabel="Parameter x [dimensionless]", ylabel="Posterior density [1 / x]", xlim=(-8, 8))
axes[0].legend(fontsize=8, loc="upper center")
ratio = np.exp(2*6*.5/(3**2+.25**2))
probability = np.array([1, ratio])/(1+ratio)
axes[1].bar([0, 1], [c[2] for c in components], color="#0072B2", width=.5, label="Recovered")
axes[1].plot([0, 1], probability, "_", color="black", markersize=24, markeredgewidth=2,
             label="Analytic")
axes[1].set(xticks=[0, 1], xticklabels=["Negative mode", "Positive mode"], ylim=(0, 1),
            ylabel="Posterior probability [dimensionless]")
axes[1].legend(fontsize=8)
fig.savefig(output, dpi=160)

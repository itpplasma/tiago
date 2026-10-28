#!/usr/bin/env python3
"""Plot synthetic profile fits and repeated-noise coverage in physical units."""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('baseline', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--coverage', type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    cases = [('joint', 'Joint fit'), ('fixed-calibration', 'Calibration fixed'),
             ('thomson-only', 'Thomson only')]
    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False,
                         'axes.spines.right': False})
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), sharex=True, sharey='row',
                             layout='constrained')
    for col, (name, title) in enumerate(cases):
        data = read(args.baseline/f'profiles-{name}.csv')
        for row, (profile, scale, label) in enumerate([
                ('1', 1e18, r'Electron density $n_e$ [$10^{18}$ m$^{-3}$]'),
                ('2', 1, r'Electron temperature $T_e$ [eV]')]):
            selected = [r for r in data if r['profile'] == profile]
            s, truth, mean, sd = (np.array([float(r[k]) for r in selected])
                                   for k in ['s', 'truth', 'map', 'sd'])
            ax = axes[row, col]
            ax.fill_between(s, (mean-1.95996398454*sd)/scale,
                            (mean+1.95996398454*sd)/scale, color='#0072B2', alpha=.2,
                            label='Local 95% interval')
            ax.plot(s, mean/scale, color='#0072B2', label='Posterior mode')
            ax.plot(s, truth/scale, 'k--', label='Synthetic truth')
            ax.grid(alpha=.2)
            if col == 0:
                ax.set_ylabel(label)
            if row == 0:
                ax.set_title(title)
            else:
                ax.set_xlabel('Normalised toroidal flux s [1]')
    axes[0, 0].legend(loc='best', fontsize=9)
    fig.suptitle('Synthetic stellarator: truth and local Laplace profile intervals')
    fig.savefig(args.output/'profiles.png', dpi=180)
    fig.savefig(args.output/'profiles.svg')
    plt.close(fig)
    calibration = read(args.baseline/'calibrations.csv')
    fig, axes = plt.subplots(1, 4, figsize=(12, 3.5), layout='constrained')
    labels = [('Thomson scale', 1, 'Absolute scale [1]'),
              ('Interferometer offset', 1, 'Phase offset [rad]'),
              ('ECE gain', 1e3, 'Receiver gain [mV/eV]'),
              ('Magnetic drift', 1e6, r'Integrator drift [$\mu$Wb/s]')]
    for ax, row, (title, scale, label) in zip(axes, calibration, labels):
        truth, mean, sd = (float(row[key])*scale for key in ['truth', 'map', 'sd'])
        ax.errorbar(mean, 0, xerr=1.95996398454*sd, fmt='o', capsize=4, color='#0072B2')
        ax.plot(truth, 1, 'kD')
        ax.axvline(truth, color='black', linestyle='--', alpha=.5)
        ax.set(xlabel=label, title=title, yticks=[0, 1],
               yticklabels=['Posterior mode', 'Injected truth'], ylim=(-.6, 1.6))
        ax.grid(axis='x', alpha=.2)
    fig.suptitle('Injected miscalibrations and local 95% posterior intervals')
    fig.savefig(args.output/'calibrations.png', dpi=180)
    fig.savefig(args.output/'calibrations.svg')
    plt.close(fig)
    if args.coverage:
        data = read(args.coverage)
        fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True,
                                 layout='constrained')
        for profile, ax, title in zip(['1', '2'], axes, ['Density', 'Temperature']):
            selected = [r for r in data if r['profile'] == profile]
            s = sorted({float(r['s']) for r in selected})
            rates, low, high = [], [], []
            for coordinate in s:
                values = [int(r['covered']) for r in selected
                          if float(r['s']) == coordinate]
                n = len(values)
                rate = sum(values)/n
                # Wilson binomial interval; realisations are the independent units.
                z = 1.95996398454
                centre = (rate+z*z/(2*n))/(1+z*z/n)
                width = z*np.sqrt(rate*(1-rate)/n+z*z/(4*n*n))/(1+z*z/n)
                rates.append(rate)
                low.append(centre-width)
                high.append(centre+width)
            ax.plot(s, rates, 'o-', color='#0072B2', label='Empirical coverage')
            ax.fill_between(s, low, high, color='#0072B2', alpha=.2,
                            label='95% binomial interval')
            ax.axhline(.95, color='black', linestyle='--', label='Nominal 95%')
            ax.set(xlabel='Normalised toroidal flux s [1]', title=title, ylim=(0, 1.02))
            ax.grid(alpha=.2)
        axes[0].set_ylabel('Coverage probability [1]')
        axes[0].legend(loc='best', fontsize=9)
        fig.suptitle(f'Independent Gaussian-noise realisations: N={n}')
        fig.savefig(args.output/'coverage.png', dpi=180)
        fig.savefig(args.output/'coverage.svg')
        plt.close(fig)


if __name__ == '__main__':
    main()

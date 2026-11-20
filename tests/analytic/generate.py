#!/usr/bin/env python3
"""Fixtures of tiago_geqdsk_cells: Cerfon-Freidberg equilibria (limited E2 at
65 x 65, single null E3 at 129 x 129) as G-EQDSK, a small machine, and the
plasma current and plasma part of the sensor signals from Gauss quadrature of
the exact current density (relative accuracy below 1e-11).
Usage: generate.py OUTDIR
"""
import sys
from pathlib import Path

import numpy as np
from scipy.special import ellipe, ellipk

import cerfon_freidberg as cf

MU0 = 4e-7 * np.pi
SENSORS = [  # kind, label, R, Z, angle [deg]
    ("flux", "FL1", 2.39, 0.0, 0), ("flux", "FL2", 1.65, 1.24, 0),
    ("flux", "FL3", 0.91, 0.0, 0), ("flux", "FL4", 1.65, -1.24, 0),
    ("pickup", "BP1", 2.31, 0.0, 90.0), ("pickup", "BP2", 1.65, 1.14, 0.0),
    ("pickup", "BP3", 0.99, 0.0, 90.0), ("pickup", "BP4", 1.65, -1.14, 0.0),
]


def field(R, Z, Rs, Zs, I):
    dz = Z - Zs
    d2 = (R + Rs) ** 2 + dz**2
    m = 4 * R * Rs / d2
    K, E = ellipk(m), ellipe(m)
    psi = MU0 * I * np.sqrt(R * Rs) / (2 * np.pi * np.sqrt(m)) * ((2 - m) * K - 2 * E)
    g = (Rs - R) ** 2 + dz**2
    c = MU0 * I / (2 * np.pi * np.sqrt(d2))
    bz = c * (K + (Rs**2 - R**2 - dz**2) / g * E)
    br = c * dz / R * (-K + (Rs**2 + R**2 + dz**2) / g * E)
    return psi.sum(), br.sum(), bz.sum()


def sources(eq, nt=128, nr=80):
    xr, wr = np.polynomial.legendre.leggauss(nr)
    xt, wt = np.polynomial.legendre.leggauss(nt)
    out = []
    b = eq._theta_breaks()
    for lo, hi in zip(b[:-1], b[1:]):
        for t, w in zip(0.5 * (hi - lo) * xt + 0.5 * (hi + lo), 0.5 * (hi - lo) * wt):
            rb = eq.radius(t)
            r = 0.5 * rb * (xr + 1)
            x = eq.axis[0] + r * np.cos(t)
            y = eq.axis[1] + r * np.sin(t)
            j = -(eq.psi0 / (MU0 * eq.R0**3)) * ((1 - eq.A) * x**2 + eq.A) / x
            out.append((x * eq.R0, y * eq.R0, j * w * 0.5 * rb * wr * r * eq.R0**2))
    return [np.concatenate(v) for v in zip(*out)]


def main(outdir):
    out = Path(outdir)
    with open(out / "analytic.machine", "w") as fh:
        fh.write("# sensors of tiago_geqdsk_cells, tests/analytic/generate.py\n")
        for kind, label, r, z, ang in SENSORS:
            if kind == "flux":
                fh.write(f"flux {label} {r} {z}\n")
            else:
                fh.write(f"pickup {label} {r} {z} {ang} 0 1\n")
    for name, nw in (("E2", 65), ("E3", 129)):
        eq = cf.Equilibrium(name)
        cf.write_geqdsk(eq, nw, out / f"{name.lower()}.geqdsk")
        Rs, Zs, I = sources(eq)
        with open(out / f"{name.lower()}_reference.csv", "w") as fh:
            fh.write("label,value\n")
            fh.write(f"plasma_current,{I.sum():.15e}\n")
            for kind, label, r, z, ang in SENSORS:
                psi, br, bz = field(r, z, Rs, Zs, I)
                a = np.radians(ang)
                v = 2 * np.pi * psi if kind == "flux" else br * np.cos(a) + bz * np.sin(a)
                fh.write(f"{label},{v:.15e}\n")


if __name__ == "__main__":
    main(sys.argv[1])

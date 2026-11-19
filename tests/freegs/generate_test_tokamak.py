"""Solve FreeGS's TestTokamakSensor machine and write TIAGO test fixtures.

Usage: python generate_test_tokamak.py OUTDIR

Needs FreeGS 0.9 at 5b41fe5565ff8077bf8139051c80c69416fc9b86. Writes
test_tokamak.geqdsk, test_tokamak.machine (TIAGO axisymmetric machine file)
and test_tokamak_freegs.csv (coil currents, FreeGS sensor values and
normalised flux at probe points). The equilibrium follows FreeGS's example
01-freeboundary.py with the sensor machine.
"""
import io
import sys
from pathlib import Path

import numpy as np
import freegs
import freeqdsk.geqdsk
from freegs import geqdsk
from freegs.gradshafranov import Greens, GreensBr, GreensBz
from freegs.machine import FluxLoopSensor, PoloidalFieldSensor, RogowskiSensor

out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
tokamak = freegs.machine.TestTokamakSensor()
eq = freegs.Equilibrium(tokamak=tokamak, Rmin=0.1, Rmax=2.0, Zmin=-1.0, Zmax=1.0,
                        nx=65, ny=65, boundary=freegs.boundary.freeBoundaryHagenow)
profiles = freegs.jtor.ConstrainPaxisIp(eq, 1e3, 2e5, 2.0)
constrain = freegs.control.constrain(xpoints=[(1.1, -0.6), (1.1, 0.8)],
                                     isoflux=[(1.1, -0.6, 1.1, 0.6)])
freegs.solve(eq, profiles, constrain, show=False, rtol=1e-10, atol=1e-10,
             maxits=200)
# FreeGS writes the profiles on psi_N = i/nx, without the boundary
# (freegs-plasma/freegs#163). Rewrite them on the G-EQDSK grid i/(nx - 1).
buffer = io.StringIO()
geqdsk.write(eq, buffer)
buffer.seek(0)
data = freeqdsk.geqdsk.read(buffer)
x = np.linspace(0.0, 1.0, data.nx)
data.fpol = eq.fpol(x)
data.pres = eq.pressure(x)
data.ffprime = eq.ffprime(x)
data.pprime = eq.pprime(x)
data.qpsi = np.append(np.nan, eq.q(x[1:]))
data.qpsi[0] = 2 * data.qpsi[1] - data.qpsi[2]
with open(out / "test_tokamak.geqdsk", "w") as fh:
    freeqdsk.geqdsk.write(data, fh)

lines = ["# FreeGS TestTokamakSensor: kind label fields, see docs/axisymmetric-free-boundary.md"]
for group, (label, coil) in enumerate(tokamak.coils, start=1):
    if hasattr(coil, "shape"):
        r = [p[0] for p in coil.shape]
        z = [p[1] for p in coil.shape]
        rc, zc = (min(r) + max(r)) / 2, (min(z) + max(z)) / 2
        lines.append(f"coil {label} {group} {rc:.6f} {zc:.6f} {max(r) - min(r):.6f} "
                     f"{max(z) - min(z):.6f} {coil.turns}")
    else:
        lines.append(f"coil {label} {group} {coil.R:.6f} {coil.Z:.6f} 0 0 {coil.turns}")
rows = ["kind,label,r,z,value"]
for label, coil in tokamak.coils:
    rows.append(f"current,{label},0,0,{coil.current:.17e}")
tokamak.takeMeasurements(eq)
# The flux loops and pickups take the plasma field from the Green's functions
# of FreeGS's own toroidal current on its grid. FreeGS's grid psi of the plasma
# differs from this sum by a near-constant 2.6e-4 Wb/rad at 65 and at 129
# points, which shifts its flux-loop values but not its fields.
jtor = profiles.Jtor(eq.R, eq.Z, eq.psi(), eq.psi_bndry)
mask = jtor != 0
rj, zj, ij = eq.R[mask], eq.Z[mask], jtor[mask] * eq.dR * eq.dZ


def plasma(green, r, z):
    return float(np.sum(green(rj, zj, r, z) * ij))


for s in tokamak.sensors:
    if isinstance(s, RogowskiSensor):
        pts = " ".join(f"{r:.6f} {z:.6f}" for r, z in zip(s.R, s.Z))
        lines.append(f"rogowski {s.name} {pts}")
        rows.append(f"rogowski,{s.name},0,0,{s.measurement:.17e}")
    elif isinstance(s, PoloidalFieldSensor):
        lines.append(f"pickup {s.name} {s.R:.6f} {s.Z:.6f} {np.degrees(s.theta):.12f} 0 1")
        br = tokamak.Br(s.R, s.Z) + plasma(GreensBr, s.R, s.Z)
        bz = tokamak.Bz(s.R, s.Z) + plasma(GreensBz, s.R, s.Z)
        value = br * np.cos(s.theta) + bz * np.sin(s.theta)
        rows.append(f"pickup,{s.name},{s.R},{s.Z},{value:.17e}")
    elif isinstance(s, FluxLoopSensor):
        lines.append(f"flux {s.name} {s.R:.6f} {s.Z:.6f}")
        # FreeGS measures psi per radian, TIAGO the linked flux 2 pi psi.
        value = 2 * np.pi * (tokamak.psi(s.R, s.Z) + plasma(Greens, s.R, s.Z))
        rows.append(f"flux,{s.name},{s.R},{s.Z},{value:.17e}")
rng = np.random.default_rng(20261119)
for k in range(24):
    r, z = rng.uniform(0.85, 1.45), rng.uniform(-0.55, 0.65)
    rows.append(f"psin,P{k:02d},{r},{z},{float(eq.psiNRZ(r, z)):.17e}")
# FreeGS integrates q from the first requested label, so request the axis too.
labels = np.linspace(0.0, 1.0, 101)
rhotor = eq.rhotor(eq.psi_axis + labels * (eq.psi_bndry - eq.psi_axis))
for k in (10, 30, 50, 70, 90):
    rows.append(f"rhotor,X{k},{labels[k]},0,{rhotor[k]:.17e}")
rows.append(f"current,plasma,0,0,{eq.plasmaCurrent():.17e}")
(out / "test_tokamak.machine").write_text("\n".join(lines) + "\n")
(out / "test_tokamak_freegs.csv").write_text("\n".join(rows) + "\n")

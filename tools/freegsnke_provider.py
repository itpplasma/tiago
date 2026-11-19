"""Free-boundary equilibrium provider for TIAGO on FreeGSNKE.

Usage: python freegsnke_provider.py REQUEST OUTPUT

REQUEST is a namelist-like text file with one `key = value` per line:
machine (TIAGO axisymmetric machine file), nx, ny, rmin, rmax, zmin, zmax,
tolerance (accepted relative residual), plasma_current [A], pressure_axis [Pa], fvac [T m], alpha_m,
alpha_n and current(k) [A] per coil group k. The equilibrium is a forward
static solve with fixed coil currents and FreeGSNKE's ConstrainPaxisIp
profiles, from FreeGSNKE's default initial flux. OUTPUT is the G-EQDSK and
OUTPUT.jphi the toroidal current density on its grid, R fastest [A/m^2].
OUTPUT.status gets `converged`, `residual` and `provider` lines.
"""
import contextlib
import importlib.metadata
import io
import sys

import freeqdsk.geqdsk
import numpy as np


def read_request(path):
    values = {}
    for line in open(path):
        line = line.split("#")[0].strip()
        if line:
            key, value = (s.strip() for s in line.split("=", 1))
            values[key] = value.strip("'\"")
    return values


def read_machine(path):
    coils, limiter = {}, None
    for line in open(path):
        fields = line.split()
        if not fields or fields[0].startswith("#"):
            continue
        if fields[0] == "coil":
            label, group = fields[1], int(fields[2])
            rc, zc, dr, dz, turns = map(float, fields[3:8])
            coils.setdefault(group, []).append((label, rc, zc, dr, dz, turns))
        elif fields[0] == "limiter":
            xy = np.array(fields[2:], dtype=float)
            limiter = (xy[0::2], xy[1::2])
    return coils, limiter


def active_coils(coils, n=4):
    """One FreeGSNKE coil per filament set, n x n midpoint filaments for a rectangle."""
    active = {}
    for group, members in sorted(coils.items()):
        for label, rc, zc, dr, dz, turns in members:
            if dr > 0 and dz > 0:
                x = (np.arange(n) + 0.5) / n - 0.5
                r, z = np.meshgrid(rc + dr * x, zc + dz * x, indexing="ij")
                r, z, m, size = r.ravel(), z.ravel(), turns / n**2, (dr / n, dz / n)
            else:
                r, z, m, size = np.array([rc]), np.array([zc]), turns, (0.01, 0.01)
            active[f"g{group}_{label}"] = dict(
                R=list(r), Z=list(z), dR=size[0], dZ=size[1], resistivity=1.7e-8,
                polarity=1, multiplier=m)
    return active


def write_finite_q(data, fh, label=None):
    """q is traced along field lines and fails next to an X-point. TIAGO uses
    it only for rho_tor, so the last finite value is held there."""
    q = np.asarray(data["qpsi"], dtype=float)
    finite = np.isfinite(q)
    q[~finite] = q[finite][-1]
    data["qpsi"] = q
    freeqdsk.geqdsk.write(data, fh, label=label)


def main(request_path, output):
    request = read_request(request_path)
    coils, limiter = read_machine(request["machine"])
    if limiter is None:
        sys.exit("machine file has no limiter line")
    with contextlib.redirect_stdout(io.StringIO()):
        from freegsnke import GSstaticsolver, build_machine, equilibrium_update, jtor_update
        from freegs4e import geqdsk
        lim = [dict(R=r, Z=z) for r, z in zip(*limiter)]
        tokamak = build_machine.tokamak(active_coils_data=active_coils(coils),
                                        limiter_data=lim, wall_data=lim)
        eq = equilibrium_update.Equilibrium(
            tokamak=tokamak, Rmin=float(request["rmin"]), Rmax=float(request["rmax"]),
            Zmin=float(request["zmin"]), Zmax=float(request["zmax"]),
            nx=int(request["nx"]), ny=int(request["ny"]))
        for group, members in coils.items():
            current = float(request.get(f"current({group})", 0.0))
            for label, *_ in members:
                eq.tokamak.set_coil_current(f"g{group}_{label}", current)
        profiles = jtor_update.ConstrainPaxisIp(
            eq=eq, paxis=float(request["pressure_axis"]),
            Ip=float(request["plasma_current"]), fvac=float(request["fvac"]),
            alpha_m=float(request.get("alpha_m", 1.0)),
            alpha_n=float(request.get("alpha_n", 2.0)))
        solver = GSstaticsolver.NKGSsolver(eq)
        tolerance = float(request.get("tolerance", 1e-8))
        # Aim below the accepted tolerance, since the last Newton steps stall
        # within a factor of a few of the requested one.
        solver.solve(eq=eq, profiles=profiles, constrain=None,
                     target_relative_tolerance=tolerance / 10, max_solving_iterations=400,
                     verbose=False)
        residual = float(solver.relative_change)
        with open(output, "w") as fh:
            geqdsk.write(eq, fh, fileformat=write_finite_q)
        np.savetxt(output + ".jphi", profiles.jtor.ravel(order="F"), fmt="%.17e")
    with open(output + ".status", "w") as fh:
        fh.write(f"converged {int(residual <= tolerance)}\n")
        fh.write(f"residual {residual:.6e}\n")
        fh.write(f"provider freegsnke {importlib.metadata.version('freegsnke')}\n")


if __name__ == "__main__":
    main(*sys.argv[1:3])

"""Free-boundary equilibrium provider for TIAGO on FreeGSNKE.

Usage: python freegsnke_provider.py REQUEST OUTPUT

REQUEST is a namelist-like text file with one `key = value` per line:
machine (TIAGO axisymmetric machine file), nx, ny, rmin, rmax, zmin, zmax,
tolerance (accepted relative residual), plasma_current [A], pressure_axis [Pa], fvac [T m], alpha_m,
alpha_n and current(k) [A] per coil group k. The equilibrium is a forward
static solve with fixed coil currents and FreeGSNKE's ConstrainPaxisIp
profiles, or with `profile = knots` natural cubic splines of p' and F F'
through `pprime_knots` and `ffprime_knots` at `psin_knots` (comma lists),
scaled together to the plasma current. The solve starts from the plasma flux
in the .npy file `initial` when given, otherwise from FreeGSNKE's default. OUTPUT is the G-EQDSK and
OUTPUT.jphi the toroidal current density on its grid, R fastest [A/m^2], and
OUTPUT.psi.npy the plasma flux on the grid, usable as `initial`.
OUTPUT.status gets `converged`, `residual` and `provider` lines. The solve
runs in a child process; one that runs longer than `time_limit` seconds
(default 600) is killed and reported with `converged 0` and `timed_out 1`,
without a G-EQDSK.
"""
import contextlib
import importlib.metadata
import io
import os
import subprocess
import sys
from pathlib import Path

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
    """FreeGSNKE coils: a rectangle as n x n midpoint filaments, and the single
    filaments of a group merged per number of turns, with the coil names
    returned per group."""
    active, names = {}, {}
    for group, members in sorted(coils.items()):
        merged = {}
        for label, rc, zc, dr, dz, turns in members:
            if dr > 0 and dz > 0:
                x = (np.arange(n) + 0.5) / n - 0.5
                r, z = np.meshgrid(rc + dr * x, zc + dz * x, indexing="ij")
                key = f"g{group}_{label}"
                active[key] = dict(R=list(r.ravel()), Z=list(z.ravel()), dR=dr / n, dZ=dz / n,
                                   resistivity=1.7e-8, polarity=1, multiplier=turns / n**2)
                names.setdefault(group, []).append(key)
            else:
                merged.setdefault(turns, []).append((rc, zc))
        for k, (turns, points) in enumerate(sorted(merged.items())):
            key = f"g{group}_t{k}"
            active[key] = dict(R=[p[0] for p in points], Z=[p[1] for p in points], dR=0.01, dZ=0.01,
                               resistivity=1.7e-8, polarity=1, multiplier=turns)
            names.setdefault(group, []).append(key)
    return active, names


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
        active, names = active_coils(coils)
        tokamak = build_machine.tokamak(active_coils_data=active,
                                        limiter_data=lim, wall_data=lim)
        initial = request.get("initial", "")
        psi0 = np.load(initial) if initial and Path(initial).exists() else None
        eq = equilibrium_update.Equilibrium(
            tokamak=tokamak, Rmin=float(request["rmin"]), Rmax=float(request["rmax"]),
            Zmin=float(request["zmin"]), Zmax=float(request["zmax"]),
            nx=int(request["nx"]), ny=int(request["ny"]), psi=psi0)
        for group, members in coils.items():
            current = float(request.get(f"current({group})", 0.0))
            for name in names[group]:
                eq.tokamak.set_coil_current(name, current)
        if request.get("profile", "paxis") == "knots":
            knots = [np.array(request[k].split(","), dtype=float)
                     for k in ("psin_knots", "pprime_knots", "ffprime_knots")]
            profiles = jtor_update.GeneralPprimeFFprime(
                eq=eq, Ip=float(request["plasma_current"]), fvac=float(request["fvac"]),
                psi_n=knots[0], pprime_data=knots[1], ffprime_data=knots[2],
                Ip_logic=True, interpolator="cubic_spline")
        else:
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
        limited = bool(profiles.flag_limiter)
        # FreeGS4E's writer, separatrix and q need an X-point, so a plasma
        # without one is reported as a failed solve instead of written.
        written = True
        try:
            with open(output, "w") as fh:
                geqdsk.write(eq, fh, fileformat=write_finite_q)
        except IndexError:
            written = False
        np.savetxt(output + ".jphi", profiles.jtor.ravel(order="F"), fmt="%.17e")
        np.save(output + ".psi.npy", eq.plasma_psi)
    with open(output + ".status", "w") as fh:
        fh.write(f"converged {int(residual <= tolerance and written)}\n")
        fh.write(f"limited {int(limited)}\n")
        fh.write(f"written {int(written)}\n")
        fh.write(f"residual {residual:.6e}\n")
        fh.write(f"provider freegsnke {importlib.metadata.version('freegsnke')}\n")


def supervise(request_path, output):
    # a Newton-Krylov solve can hang for hours without progress (sCluster),
    # and a signal does not interrupt FreeGSNKE's long NumPy calls
    limit = float(read_request(request_path).get("time_limit", 600))
    try:
        subprocess.run([sys.executable, "-I", __file__, request_path, output],
                       env=dict(os.environ, TIAGO_PROVIDER_CHILD="1"), check=True, timeout=limit)
    except subprocess.TimeoutExpired:
        for suffix in ("", ".jphi", ".psi.npy"):
            Path(output + suffix).unlink(missing_ok=True)
        with open(output + ".status", "w") as fh:
            fh.write("converged 0\ntimed_out 1\n")
            fh.write(f"provider freegsnke {importlib.metadata.version('freegsnke')}\n")


if __name__ == "__main__":
    if os.environ.get("TIAGO_PROVIDER_CHILD"):
        main(*sys.argv[1:3])
    else:
        supervise(*sys.argv[1:3])

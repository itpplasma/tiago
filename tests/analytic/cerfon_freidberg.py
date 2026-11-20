#!/usr/bin/env python3
"""Analytic Solov'ev and Cerfon-Freidberg equilibria, written as G-EQDSK.

Normalised Grad-Shafranov equation, x = R/R0, y = Z/R0, psi = psi0 Psi:

    x d/dx (1/x dPsi/dx) + d2Psi/dy2 = (1 - A) x^2 + A,

    Psi = x^4/8 + A (x^2 ln x / 2 - x^4/8) + sum_i c_i Psi_i,

with Psi_i the homogeneous polynomials of Cerfon and Freidberg, Phys. Plasmas
17, 032502 (2010). Then mu0 R0^4 p' = -psi0 (1 - A), R0^2 F F' = -psi0 A,
p = -(1 - A) psi0^2 Psi / (mu0 R0^4), F^2 = F0^2 - 2 A psi0^2 Psi / R0^2.
psi0 < 0 gives current along +phi.

Usage: cerfon_freidberg.py CASE NW OUT.geqdsk [OUT.json]
CASE is E1 (Solov'ev, limited), E2 (shaped, limited), E3 (single null) or E4
(double null), with R0 = 1.65 m, B0 = 2.5 T, I_P = 0.8 MA. The JSON holds the
reference values: I_P, axis, psi at axis and boundary, q and toroidal flux.
"""
import json
import sys

import numpy as np
import sympy as sp
from scipy import integrate, optimize

MU0 = 4e-7 * np.pi
NBOUNDARY = 401
x, y = sp.symbols("x y", positive=True)
L = sp.log(x)

HOMOGENEOUS = [
    sp.Integer(1),
    x**2,
    y**2 - x**2 * L,
    x**4 - 4 * x**2 * y**2,
    2 * y**4 - 9 * y**2 * x**2 + 3 * x**4 * L - 12 * x**2 * y**2 * L,
    x**6 - 12 * x**4 * y**2 + 8 * x**2 * y**4,
    8 * y**6 - 140 * y**4 * x**2 + 75 * y**2 * x**4 - 15 * x**6 * L
    + 180 * x**4 * y**2 * L - 120 * x**2 * y**4 * L,
    # up-down asymmetric terms
    y,
    y * x**2,
    y**3 - 3 * y * x**2 * L,
    3 * y * x**4 - 4 * y**3 * x**2,
    8 * y**5 - 45 * y * x**4 - 80 * y**3 * x**2 * L + 60 * y * x**4 * L,
]

CASES = {
    # R0 [m], B0 [T], I_P target [A], eps, kappa, delta, A, topology
    "E1": dict(R0=1.65, B0=2.5, Ip=0.8e6, eps=0.3, kappa=1.0, delta=0.0, A=0.0, kind="limited"),
    "E2": dict(R0=1.65, B0=2.5, Ip=0.8e6, eps=0.32, kappa=1.7, delta=0.33, A=-0.155, kind="limited"),
    "E3": dict(R0=1.65, B0=2.5, Ip=0.8e6, eps=0.32, kappa=1.7, delta=0.33, A=-0.155, kind="single_null"),
    "E4": dict(R0=1.65, B0=2.5, Ip=0.8e6, eps=0.32, kappa=1.7, delta=0.33, A=-0.155, kind="double_null"),
}


def gs_operator(f):
    return sp.simplify(x * sp.diff(sp.diff(f, x) / x, x) + sp.diff(f, y, 2))


def check_homogeneous():
    for k, f in enumerate(HOMOGENEOUS):
        assert gs_operator(f) == 0, f"Psi_{k + 1} is not homogeneous"
    for A in (0, sp.Rational(-1, 3)):
        part = x**4 / 8 + A * (x**2 * L / 2 - x**4 / 8)
        assert sp.simplify(gs_operator(part) - ((1 - A) * x**2 + A)) == 0


def solve_coefficients(c):
    eps, kappa, delta, A = c["eps"], c["kappa"], c["delta"], c["A"]
    alpha = np.arcsin(delta)
    n1 = -(1 + alpha) ** 2 / (eps * kappa**2)
    n2 = (1 - alpha) ** 2 / (eps * kappa**2)
    n3 = -kappa / (eps * np.cos(alpha) ** 2)
    outer, inner, top = (1 + eps, 0.0), (1 - eps, 0.0), (1 - delta * eps, kappa * eps)
    xsep = (1 - 1.1 * delta * eps, -1.1 * kappa * eps)
    nb = 12 if c["kind"] == "single_null" else 7
    basis = HOMOGENEOUS[:nb]
    part = x**4 / 8 + A * (x**2 * L / 2 - x**4 / 8)

    conditions = []  # (point, list of (dx, dy, weight))
    if c["kind"] == "limited":
        conditions = [
            (outer, [(0, 0, 1)]), (inner, [(0, 0, 1)]), (top, [(0, 0, 1)]),
            (top, [(1, 0, 1)]),
            (outer, [(0, 2, 1), (1, 0, n1)]), (inner, [(0, 2, 1), (1, 0, n2)]),
            (top, [(2, 0, 1), (0, 1, n3)]),
        ]
    elif c["kind"] == "double_null":
        xs = (xsep[0], -xsep[1])
        conditions = [
            (outer, [(0, 0, 1)]), (inner, [(0, 0, 1)]), (xs, [(0, 0, 1)]),
            (xs, [(1, 0, 1)]), (xs, [(0, 1, 1)]),
            (outer, [(0, 2, 1), (1, 0, n1)]), (inner, [(0, 2, 1), (1, 0, n2)]),
        ]
    else:
        conditions = [
            (outer, [(0, 0, 1)]), (inner, [(0, 0, 1)]), (top, [(0, 0, 1)]), (xsep, [(0, 0, 1)]),
            (top, [(1, 0, 1)]), (outer, [(0, 1, 1)]), (inner, [(0, 1, 1)]),
            (xsep, [(1, 0, 1)]), (xsep, [(0, 1, 1)]),
            (outer, [(0, 2, 1), (1, 0, n1)]), (inner, [(0, 2, 1), (1, 0, n2)]),
            (top, [(2, 0, 1), (0, 1, n3)]),
        ]

    def d(f, dx, dy, p):
        g = f
        if dx:
            g = sp.diff(g, x, dx)
        if dy:
            g = sp.diff(g, y, dy)
        return float(g.subs({x: p[0], y: p[1]}))

    M = np.array([[sum(w * d(f, dx, dy, p) for dx, dy, w in terms) for f in basis]
                  for p, terms in conditions])
    b = -np.array([sum(w * d(part, dx, dy, p) for dx, dy, w in terms) for p, terms in conditions])
    coef = np.linalg.solve(M, b)
    Psi = part + sum(float(ci) * f for ci, f in zip(coef, basis))
    return Psi, dict(outer=outer, inner=inner, top=top, xsep=xsep if c["kind"] != "limited" else None)


class Equilibrium:
    def __init__(self, name):
        self.c = c = CASES[name]
        self.name = name
        Psi, self.points = solve_coefficients(c)
        self.Psi = sp.lambdify((x, y), Psi, "numpy")
        self.Psi_x = sp.lambdify((x, y), sp.diff(Psi, x), "numpy")
        self.Psi_y = sp.lambdify((x, y), sp.diff(Psi, y), "numpy")
        self.R0 = c["R0"]
        self.A = c["A"]
        ax = optimize.minimize(lambda p: self.Psi(*p), [1.0, 0.0], method="Nelder-Mead",
                               options=dict(xatol=1e-14, fatol=1e-16, maxiter=20000))
        ax = optimize.root(lambda p: [self.Psi_x(*p), self.Psi_y(*p)], ax.x, tol=1e-15)
        self.axis = ax.x
        self.Psi_axis = float(self.Psi(*ax.x))
        # psi0 from the plasma current: I_P = -psi0/(mu0 R0) int [(1-A) x^2 + A]/x dx dy
        self.psi0 = -c["Ip"] * MU0 * self.R0 / self._current_integral()
        self.F0 = self.R0 * c["B0"]

    def radius(self, theta, target=0.0):
        """Distance from the axis along angle theta to the surface Psi = target
        (normalised units). Psi rises monotonically along a ray inside the
        separatrix. A ray through an X-point reaches its maximum there, so the
        march stops where Psi stops rising as well as where it reaches target."""
        ct, st = np.cos(theta), np.sin(theta)
        f = lambda r: self.Psi(self.axis[0] + r * ct, self.axis[1] + r * st) - target
        if f(0.0) >= 0:
            return 0.0
        h = 1e-3
        rs = np.arange(0.0, 1.0, h)
        fs = f(rs)
        up = np.flatnonzero(fs >= 0)
        down = np.flatnonzero(np.diff(fs) < 0)
        iu = up[0] if up.size else rs.size
        idn = down[0] + 1 if down.size else rs.size
        if iu <= idn:
            return optimize.brentq(f, rs[iu - 1], rs[iu], xtol=1e-15)
        res = optimize.minimize_scalar(lambda q: -f(q), bounds=(rs[idn - 2], rs[idn]),
                                       method="bounded", options=dict(xatol=1e-14))
        return res.x

    def _current_integral(self):
        A = self.A

        def inner(theta):
            rb = self.radius(theta)
            g = lambda r: ((1 - A) * (self.axis[0] + r * np.cos(theta)) ** 2 + A) / (
                self.axis[0] + r * np.cos(theta)) * r
            return integrate.quad(g, 0, rb, epsabs=0, epsrel=1e-13, limit=200)[0]

        breaks = self._theta_breaks()
        total = 0.0
        for a, b in zip(breaks[:-1], breaks[1:]):
            total += integrate.quad(inner, a, b, epsabs=0, epsrel=1e-12, limit=400)[0]
        return total

    def _theta_breaks(self):
        breaks = [-np.pi, np.pi]
        for p in (self.points.get("xsep"),):
            if p is not None:
                for s in (1, -1) if self.c["kind"] == "double_null" else (1,):
                    breaks.append(np.arctan2(s * p[1] - self.axis[1], p[0] - self.axis[0]))
        return sorted(breaks)

    # physical quantities, SI
    def psi(self, R, Z):
        return self.psi0 * self.Psi(R / self.R0, Z / self.R0)

    def pprime(self):
        return -self.psi0 * (1 - self.A) / (MU0 * self.R0**4)

    def ffprime(self):
        return -self.psi0 * self.A / self.R0**2

    def pressure(self, psin):
        Psi = self.Psi_axis * (1 - psin)
        return -(1 - self.A) * self.psi0**2 * Psi / (MU0 * self.R0**4)

    def fpol(self, psin):
        Psi = self.Psi_axis * (1 - psin)
        return np.sqrt(self.F0**2 - 2 * self.A * self.psi0**2 * Psi / self.R0**2)

    def surface(self, psin, n=2048):
        """Flux surface at psi_N as (R, Z) arrays by ray root finding."""
        target = self.Psi_axis * (1 - psin)
        th = np.linspace(-np.pi, np.pi, n, endpoint=False)
        out = []
        for t in th:
            r = self.radius(t, target)
            out.append((self.axis[0] + r * np.cos(t), self.axis[1] + r * np.sin(t)))
        xy = np.array(out) * self.R0
        return xy[:, 0], xy[:, 1]

    def q(self, psin):
        """q = F/(2 pi) oint dl / (R |grad psi|), as an integral over the ray angle."""
        target = self.Psi_axis * (1 - psin)
        F = self.fpol(psin)

        def integrand(t):
            ct, st = np.cos(t), np.sin(t)
            r = self.radius(t, target)
            xx, yy = self.axis[0] + r * ct, self.axis[1] + r * st
            gx, gy = self.Psi_x(xx, yy), self.Psi_y(xx, yy)
            dpsidr = gx * ct + gy * st  # d Psi / d r along the ray
            # dl / |grad Psi| = r dt / (d Psi/dr) on a ray-parameterised surface
            return r / (xx * dpsidr)

        breaks = self._theta_breaks()
        s = sum(integrate.quad(integrand, a, b, epsabs=0, epsrel=1e-11, limit=400)[0]
                for a, b in zip(breaks[:-1], breaks[1:]))
        # dl = R0 dl', |grad psi| = |psi0| |grad' Psi| / R0, R = R0 x
        return F / (2 * np.pi) * self.R0 * s / abs(self.psi0)


def write_geqdsk(eq, nw, path, psin_q_max=None):
    nh = nw
    c = eq.c
    rb, zb = eq.surface(1.0, n=NBOUNDARY)
    pad = 0.15 * c["eps"] * eq.R0
    rmin, rmax = rb.min() - pad, rb.max() + pad
    zmin, zmax = zb.min() - pad, zb.max() + pad
    R = np.linspace(rmin, rmax, nw)
    Z = np.linspace(zmin, zmax, nh)
    psi = eq.psi(R[None, :], Z[:, None])  # (nh, nw)
    psin = np.linspace(0, 1, nw)
    qmax = psin_q_max if psin_q_max is not None else (1.0 if c["kind"] == "limited" else 0.999)
    q = np.array([eq.q(min(p, qmax)) if p > 0 else eq.q(1e-6) for p in psin])
    simag = eq.psi0 * eq.Psi_axis
    sibry = 0.0
    ra, za = eq.axis * eq.R0
    fpol = eq.fpol(psin)
    pres = eq.pressure(psin)
    ffp = np.full(nw, eq.ffprime())
    pp = np.full(nw, eq.pprime())
    lim_r = np.array([rmin, rmax, rmax, rmin, rmin]) + np.array([1, -1, -1, 1, 1]) * 0.02
    lim_z = np.array([zmin, zmin, zmax, zmax, zmin]) + np.array([1, 1, -1, -1, 1]) * 0.02

    def block(values, fh):
        values = np.ravel(values)
        for i, v in enumerate(values):
            fh.write(f"{v:16.9E}")
            if i % 5 == 4:
                fh.write("\n")
        if len(values) % 5:
            fh.write("\n")

    with open(path, "w") as fh:
        fh.write(f"{'TIAGO-ANALYTIC ' + eq.name:48s}{0:4d}{nw:4d}{nh:4d}\n")
        block([rmax - rmin, zmax - zmin, eq.R0, rmin, (zmin + zmax) / 2], fh)
        block([ra, za, simag, sibry, c["B0"]], fh)
        block([c["Ip"], simag, 0.0, ra, 0.0], fh)
        block([za, 0.0, sibry, 0.0, 0.0], fh)
        block(fpol, fh)
        block(pres, fh)
        block(ffp, fh)
        block(pp, fh)
        block(psi, fh)
        block(q, fh)
        fh.write(f"{len(rb):5d}{len(lim_r):5d}\n")
        block(np.column_stack([rb, zb]), fh)
        block(np.column_stack([lim_r, lim_z]), fh)
    return dict(R=R, Z=Z, q=q, psin=psin)


def main():
    check_homogeneous()
    name, nw, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    eq = Equilibrium(name)
    grid = write_geqdsk(eq, nw, out)
    if len(sys.argv) > 4:
        # toroidal flux Phi(psi_N) = 2 pi |psi_b - psi_axis| int_0^psi_N q dpsi_N
        psin = grid["psin"]
        dpsi = abs(eq.psi0 * eq.Psi_axis)
        phi = 2 * np.pi * dpsi * integrate.cumulative_simpson(grid["q"], x=psin, initial=0.0)
        ref = dict(case=name, I_P=eq.c["Ip"], psi0=eq.psi0, R_axis=eq.axis[0] * eq.R0,
                   Z_axis=eq.axis[1] * eq.R0, psi_axis=eq.psi0 * eq.Psi_axis, psi_boundary=0.0,
                   pprime=eq.pprime(), ffprime=eq.ffprime(), psin=list(psin), q=list(grid["q"]),
                   phi=list(phi))
        with open(sys.argv[4], "w") as fh:
            json.dump(ref, fh, indent=1)


if __name__ == "__main__":
    main()

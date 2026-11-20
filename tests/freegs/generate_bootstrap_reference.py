"""Sauter bootstrap reference of TORAX on the FreeGS TestTokamakSensor equilibrium.

Usage: python generate_bootstrap_reference.py OUTDIR  (TORAX 2.0.0, freeqdsk)

Writes test_tokamak_kinetic.csv (n_e = n_i with a pedestal at psi_N = 0.95,
T_e = T_i = p / (2 e n_e), at least 10 eV, Z_eff = 1.5) and
test_tokamak_bootstrap.csv with, at nine labels, TORAX's bounce-averaged
trapped fraction (Lin-Liu and Miller), collisionalities, Sauter coefficients
and <j.B>_bs on surfaces taken as matplotlib contours of psi_N, independently
of TIAGO's surface tracing, and the equilibrium's <j.B> = F p' + F' <B^2> / mu0
there.
"""
import sys
from pathlib import Path
from types import SimpleNamespace

import freeqdsk.geqdsk
import matplotlib
import numpy as np
from scipy.interpolate import CubicSpline, RectBivariateSpline

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from torax._src.geometry import trapped_fraction  # noqa: E402
from torax._src.neoclassical.formulas import formulas, sauter  # noqa: E402
from torax._src.physics import collisions  # noqa: E402

E = 1.602176634e-19
LABELS = [0.2, 0.4, 0.6, 0.8, 0.9, 0.93, 0.95, 0.97, 0.98]
ZEFF = 1.5


def surface(d, x):
    R = np.linspace(d["rleft"], d["rleft"] + d["rdim"], d["nx"])
    Z = np.linspace(d["zmid"] - d["zdim"] / 2, d["zmid"] + d["zdim"] / 2, d["ny"])
    psin = (d["psi"] - d["simagx"]) / (d["sibdry"] - d["simagx"])
    cs = plt.contour(R, Z, psin.T, levels=[x])
    seg = min((p for p in cs.allsegs[0] if len(p) > 10),
              key=lambda p: np.hypot(p[:, 0].mean() - d["rmagx"], p[:, 1].mean() - d["zmagx"]))
    plt.close()
    if np.hypot(*(seg[-1] - seg[0])) > 1e-9:
        seg = np.vstack([seg, seg[:1]])
    s = np.concatenate([[0], np.cumsum(np.hypot(*np.diff(seg, axis=0).T))])
    t = np.linspace(0, s[-1], 4001)[:-1]
    r = CubicSpline(s, seg[:, 0], bc_type="periodic")(t)
    z = CubicSpline(s, seg[:, 1], bc_type="periodic")(t)
    spl = RectBivariateSpline(R, Z, d["psi"])
    bp = np.hypot(spl(r, z, dx=1, grid=False), spl(r, z, dy=1, grid=False)) / r
    f = CubicSpline(np.linspace(0, 1, d["nx"]), d["fpol"])(x)
    return r, np.sqrt(bp ** 2 + (f / r) ** 2), (t[1] - t[0]) / bp, f


out = Path(sys.argv[1])
with open(out / "test_tokamak.geqdsk") as fh:
    d = freeqdsk.geqdsk.read(fh)
x = np.linspace(0, 1, 101)
p = CubicSpline(np.linspace(0, 1, d["nx"]), d["pres"])(x)
ne = 5.0e17 + 1.0e18 * 0.5 * (1 - np.tanh((x - 0.95) / 0.02)) + 5.0e17 * (1 - x ** 2)
t = np.maximum(p / (2 * E * ne), 10.0)
np.savetxt(out / "test_tokamak_kinetic.csv", np.column_stack([x, ne, t, t, ne, np.full_like(x, ZEFF)]),
           delimiter=",", header="psin,ne,te,ti,ni,zeff", comments="", fmt="%.10e")
kin = CubicSpline(x, np.column_stack([ne, t, t, ne]), bc_type="natural")
q = CubicSpline(np.linspace(0, 1, d["nx"]), np.abs(d["qpsi"]))
dpsi = d["sibdry"] - d["simagx"]  # FreeGS: per radian
rows = ["psin,f_trap,nu_e,nu_i,l31,l32,l34,alpha,jdotb,jdotb_eq"]
for xl in LABELS:
    r, b, w, f = surface(d, xl)
    ft = float(trapped_fraction.calculate_bounce_averaged_trapped_fraction(b, w, n_lambda=2001))
    geo = SimpleNamespace(R_major_profile_face=(r.max() + r.min()) / 2,
                          epsilon_face=(r.max() - r.min()) / (r.max() + r.min()))
    n_e, t_e, t_i, n_i = kin(xl)
    dn_e, dt_e, dt_i, dn_i = kin(xl, 1) / dpsi
    nue = float(formulas._calculate_nu_e_star(q(xl), geo, n_e, t_e / 1e3, ZEFF,
                                              collisions.calculate_log_lambda_ei(t_e / 1e3, n_e)))
    nui = float(formulas.calculate_nu_i_star(q(xl), geo, n_i, t_i / 1e3, 1.0,
                                             collisions.calculate_log_lambda_ii(t_i / 1e3, n_i, 1.0)))
    l31, l32, l34 = (float(g(ft, nue, ZEFF)) for g in (sauter.calculate_L31, sauter.calculate_L32,
                                                       sauter.calculate_L34))
    alpha = float(sauter.calculate_alpha(ft, nui))
    pe, pi = E * n_e * t_e, E * n_i * t_i
    # TORAX's analytic bootstrap sum with d/dpsi per radian; the sign follows
    # TIAGO's convention B_p = grad psi x grad phi (docs/bootstrap.md).
    jb = f * (l31 * (pe * dn_e / n_e + pi * dn_i / n_i) + (l31 + l32) * pe * dt_e / t_e
              + (l31 + alpha * l34) * pi * dt_i / t_i)
    # FreeGS writes psi per radian with B_p = grad psi x grad phi, as TIAGO
    pp, ffp = (CubicSpline(np.linspace(0, 1, d["nx"]), d[k])(xl) for k in ("pprime", "ffprime"))
    jeq = f * pp + ffp / f * np.sum(b ** 2 * w) / np.sum(w) / (4e-7 * np.pi)
    rows.append(",".join(f"{v:.12e}" for v in (xl, ft, nue, nui, l31, l32, l34, alpha, jb, jeq)))
(out / "test_tokamak_bootstrap.csv").write_text("\n".join(rows) + "\n")

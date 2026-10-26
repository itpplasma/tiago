#!/usr/bin/env python3
"""Figures for one tiago_infer run (reads the CSV files of its output_dir).

flux_surfaces.png  surfaces rho = 0.25, 0.5, 0.75, 1 at three toroidal cuts:
                   posterior draws (thin), MAP (solid), truth (dashed), and
                   the shift of each surface's outboard/inboard extreme R
                   relative to the truth with 68/95 % bands
profiles.png       p, iota, ne, Te against s with 68 % and 95 % posterior
                   bands (NUTS draws; Laplace draws if no chain) and the truth
flux_surface_bands, profiles_bands (png + pdf)
                   truth, posterior median and 68/95 % bands of the cross
                   sections, and of p, ne, Te, Ti, iota, enclosed current
corner.png         pairwise NUTS draws of the leading latents with the truth

All quantities are computed by TIAGO (synthetic study); nothing is smoothed.
Usage: MPLBACKEND=Agg inference_plots.py <output_dir> [figure_dir]
"""
import csv
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

# Okabe-Ito (colour-vision safe); line styles and markers add redundancy.
BLUE, ORANGE, GREEN, VERMILION, GREY = "#0072B2", "#E69F00", "#009E73", "#D55E00", "#555555"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 10, "legend.fontsize": 8})


def read_csv(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    return rows


def column(rows, key, cast=float):
    return np.array([cast(r[key]) for r in rows])


def extremes(rows, cut):
    """Outboard and inboard extreme R of each (draw, rho) surface on a cut."""
    res = {}
    for r in rows:
        if int(r["cut"]) != cut:
            continue
        key = (int(r["draw"]), float(r["rho"]))
        R = float(r["R"])
        lo, hi = res.get(key, (np.inf, -np.inf))
        res[key] = (min(lo, R), max(hi, R))
    return res


def flux_surfaces(out, fig_dir):
    """Top: surfaces (draws, MAP, truth). Bottom: shift of the outboard and
    inboard extreme R of each surface relative to the truth [mm], 68/95 %
    bands over the posterior draws; the geometric extremes do not depend on
    the poloidal angle label of each equilibrium."""
    truth = read_csv(os.path.join(out, "surfaces_truth.csv"))
    fmap = read_csv(os.path.join(out, "surfaces_map.csv"))
    dpath = os.path.join(out, "surfaces_draws.csv")
    draws = read_csv(dpath) if os.path.exists(dpath) else []
    fig, axes = plt.subplots(2, 3, figsize=(10.5, 7.2), gridspec_kw={"height_ratios": [3, 2]})
    for cut in range(3):
        ax = axes[0, cut]
        for rows, style in ((draws, "draw"), (fmap, "map"), (truth, "truth")):
            sel = [r for r in rows if int(r["cut"]) == cut]
            keys = sorted({(int(r["draw"]), float(r["rho"])) for r in sel})
            for draw, rho in keys:
                pts = [r for r in sel if int(r["draw"]) == draw and float(r["rho"]) == rho]
                R = column(pts, "R")
                Z = column(pts, "Z")
                if style == "draw":
                    ax.plot(R, Z, color=BLUE, lw=0.4, alpha=0.25)
                elif style == "map":
                    ax.plot(R, Z, color=BLUE, lw=1.2)
                else:
                    ax.plot(R, Z, color=VERMILION, lw=1.2, ls="--")
        ax.set_aspect("equal")
        ax.set_xlabel("R [m]")
        ax.set_title(f"cut {cut}: phi = {cut}/3 of a half period")
        # displacement of the extremes relative to the truth
        ax = axes[1, cut]
        et = extremes(truth, cut)
        ed = extremes(draws, cut)
        em = extremes(fmap, cut)
        rhos = sorted({k[1] for k in et})
        for side, idx, marker, color in (("outboard", 1, "o", BLUE), ("inboard", 0, "s", ORANGE)):
            tr = np.array([et[(0, r)][idx] for r in rhos])
            if ed:
                ids = sorted({k[0] for k in ed})
                mat = np.array([[ed[(d, r)][idx] for r in rhos] for d in ids]) - tr
                mat *= 1e3
                ax.fill_between(rhos, np.quantile(mat, 0.025, 0), np.quantile(mat, 0.975, 0),
                                color=color, alpha=0.18, lw=0)
                ax.fill_between(rhos, np.quantile(mat, 0.16, 0), np.quantile(mat, 0.84, 0),
                                color=color, alpha=0.38, lw=0)
            mp = (np.array([em[(-1, r)][idx] for r in rhos]) - tr) * 1e3
            ax.plot(rhos, mp, marker=marker, color=color, lw=1.0, label=f"{side} (MAP, bands 68/95 %)")
        ax.axhline(0.0, color=VERMILION, ls="--", lw=1.0)
        ax.set_xlabel("rho")
        ax.grid(alpha=0.3)
    axes[0, 0].set_ylabel("Z [m]")
    axes[1, 0].set_ylabel("extreme R - truth [mm]")
    handles = [
        plt.Line2D([], [], color=BLUE, lw=0.6, alpha=0.5, label="posterior draws"),
        plt.Line2D([], [], color=BLUE, lw=1.2, label="MAP"),
        plt.Line2D([], [], color=VERMILION, lw=1.2, ls="--", label="truth"),
    ]
    axes[0, -1].legend(handles=handles, loc="upper right", frameon=False)
    axes[1, -1].legend(loc="best", frameon=False)
    fig.suptitle("Flux surfaces rho = 0.25, 0.5, 0.75, 1 (VMEC++, LI383 fixed boundary)")
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, "flux_surfaces.png"), dpi=160)
    plt.close(fig)


def bands(ax, s, q, label_prefix, color):
    ax.fill_between(s, q["q025"], q["q975"], color=color, alpha=0.18, lw=0,
                    label=f"{label_prefix} 95 %")
    ax.fill_between(s, q["q16"], q["q84"], color=color, alpha=0.38, lw=0,
                    label=f"{label_prefix} 68 %")
    ax.plot(s, q["q50"], color=color, lw=1.2, label=f"{label_prefix} median")


def profiles(out, fig_dir):
    prof = read_csv(os.path.join(out, "profiles.csv"))
    truth = read_csv(os.path.join(out, "truth_profiles.csv"))
    method = "nuts" if any(r["method"] == "nuts" for r in prof) else "laplace"
    st = column(truth, "s")
    panels = [("p", "p [kPa]", 1e-3), ("iota", "iota", 1.0), ("ne", "ne [1e19 m^-3]", 1e-19),
              ("Te", "Te [keV]", 1e-3)]
    fig, axes = plt.subplots(2, 2, figsize=(9, 6.5), sharex=True)
    for ax, (name, ylabel, scale) in zip(axes.flat, panels):
        if name == "iota":
            rows = read_csv(os.path.join(out, "iota.csv"))
            draws = sorted({int(r["draw"]) for r in rows if int(r["draw"]) > 0})
            s = column([r for r in rows if int(r["draw"]) == 0], "s")
            mat = np.array([column([r for r in rows if int(r["draw"]) == d], "iota") for d in draws])
            q = {k: np.quantile(mat, v, axis=0) for k, v in
                 (("q025", 0.025), ("q16", 0.16), ("q50", 0.5), ("q84", 0.84), ("q975", 0.975))}
            bands(ax, s, q, f"{method} ({len(draws)} solves)", BLUE)
            ax.plot(s, column([r for r in rows if int(r["draw"]) == 0], "iota"), color=VERMILION,
                    ls="--", lw=1.2, label="truth")
        else:
            sel = [r for r in prof if r["method"] == method and r["quantity"] == name]
            s = column(sel, "s")
            q = {k: column(sel, k) * scale for k in ("q025", "q16", "q50", "q84", "q975")}
            bands(ax, s, q, method, BLUE)
            ax.plot(st, column(truth, name) * scale, color=VERMILION, ls="--", lw=1.2, label="truth")
        ax.set_xlim(0.0, 1.0)
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.3)
    for ax in axes[1]:
        ax.set_xlabel("s (normalized toroidal flux)")
    axes[0, 0].legend(frameon=False, loc="upper right")
    fig.suptitle("Posterior profiles with 68 % / 95 % bands and the synthetic truth")
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, "profiles.png"), dpi=160)
    plt.close(fig)


def normal_shift(tr, pts):
    """Signed distance of the closed contour pts from the closed truth
    contour tr, measured at each truth vertex along the outward normal of
    the truth (nearest point of the densified contour pts)."""
    t = np.gradient(tr, axis=0)
    n = np.stack([t[:, 1], -t[:, 0]], 1)
    n /= np.hypot(n[:, 0], n[:, 1])[:, None]
    c = tr.mean(0)
    if np.sum(n * (tr - c)) < 0:
        n = -n
    k = np.arange(len(pts))
    fine = np.linspace(0, len(pts) - 1, 12 * len(pts))
    dense = np.stack([np.interp(fine, k, pts[:, 0]), np.interp(fine, k, pts[:, 1])], 1)
    d = dense[None, :, :] - tr[:, None, :]
    j = np.argmin(np.sum(d * d, 2), 1)
    return np.sum((dense[j] - tr) * n, 1)


def polar_radius(pts, centre, alpha):
    """Radius of a closed contour about centre on the angle grid alpha
    (star-shaped about centre; periodic interpolation)."""
    dR, dZ = pts[:, 0] - centre[0], pts[:, 1] - centre[1]
    a, r = np.arctan2(dZ, dR), np.hypot(dR, dZ)
    o = np.argsort(a)
    return np.interp(alpha, a[o], r[o], period=2 * np.pi)


def surface_bands(out, fig_dir):
    """Cross sections at three toroidal angles: truth, posterior median and
    68/95 % bands of the re-solved posterior draws (radial extent about the truth
    centre), plus the band of the displacement
    from the truth along the truth normal in mm."""
    path = os.path.join(out, "surfaces_draws.csv")
    if not os.path.exists(path):
        return
    truth, draws = read_csv(os.path.join(out, "surfaces_truth.csv")), read_csv(path)
    ids = sorted({int(r["draw"]) for r in draws})
    alpha = np.linspace(-np.pi, np.pi, 181)

    def contour(rows, draw, cut, rho):
        sel = [r for r in rows if int(r["draw"]) == draw and int(r["cut"]) == cut
               and abs(float(r["rho"]) - rho) < 1e-9]
        return np.array([[float(r["R"]), float(r["Z"])] for r in sel])

    rhos = sorted({float(r["rho"]) for r in truth})
    cuts = sorted({int(r["cut"]) for r in truth})
    shades = {rho: c for rho, c in zip(rhos, ("#9ecae1", "#6baed6", "#3182bd", "#08519c"))}
    fig, axes = plt.subplots(2, len(cuts), figsize=(3.6 * len(cuts), 7.4),
                             gridspec_kw={"height_ratios": [3, 2]}, squeeze=False)
    for ic, cut in enumerate(cuts):
        ax, bx = axes[0, ic], axes[1, ic]
        centre = contour(truth, 0, cut, rhos[0]).mean(0)
        for rho in rhos:
            tr = contour(truth, 0, cut, rho)
            rt = polar_radius(tr, centre, alpha)
            mat = np.array([polar_radius(contour(draws, d, cut, rho), centre, alpha) for d in ids])
            q = np.quantile(mat, [0.025, 0.16, 0.5, 0.84, 0.975], axis=0)
            xy = lambda r: (centre[0] + r * np.cos(alpha), centre[1] + r * np.sin(alpha))
            for lo, hi, al in ((0, 4, 0.30), (1, 3, 0.55)):
                x1, y1 = xy(q[lo])
                x2, y2 = xy(q[hi])
                ax.fill(np.r_[x2, x1[::-1]], np.r_[y2, y1[::-1]], color=BLUE, alpha=al, lw=0)
            ax.plot(*xy(q[2]), color=BLUE, lw=1.0)
            ax.plot(*xy(rt), color=VERMILION, lw=1.0, ls="--")
            if rho in (0.5, 0.75):
                col = BLUE if rho == 0.5 else ORANGE
                nm = np.array([normal_shift(tr, contour(draws, d, cut, rho)) for d in ids]) * 1e3
                qn = np.quantile(nm, [0.025, 0.16, 0.5, 0.84, 0.975], axis=0)
                x = np.linspace(0.0, 360.0, len(tr))
                bx.fill_between(x, qn[0], qn[4], color=col, alpha=0.25, lw=0)
                bx.fill_between(x, qn[1], qn[3], color=col, alpha=0.5, lw=0)
                bx.plot(x, qn[2], color=col, lw=1.1, label=f"rho = {rho:g}")
        ax.set_aspect("equal")
        ax.set_xlabel("R [m]")
        ax.set_title(f"phi = {cut}/3 of a half period")
        bx.axhline(0.0, color=VERMILION, ls="--", lw=1.0)
        bx.set_xlabel("poloidal label theta [deg]")
        bx.set_xlim(0, 360)
        bx.set_xticks(range(0, 361, 90))
        bx.grid(alpha=0.3)
    axes[0, 0].set_ylabel("Z [m]")
    axes[1, 0].set_ylabel("outward normal shift vs truth [mm]")
    axes[0, -1].legend(handles=[
        plt.Line2D([], [], color=VERMILION, ls="--", lw=1.0, label="truth"),
        plt.Line2D([], [], color=BLUE, lw=1.0, label="posterior median"),
        plt.Rectangle((0, 0), 1, 1, color=BLUE, alpha=0.55, lw=0, label="68 % band"),
        plt.Rectangle((0, 0), 1, 1, color=BLUE, alpha=0.30, lw=0, label="95 % band")],
        loc="upper right", frameon=False)
    axes[1, -1].legend(loc="best", frameon=False)
    fig.suptitle(f"Flux surfaces rho = {', '.join(f'{r:g}' for r in rhos)}: truth and posterior "
                 f"bands from {len(ids)} re-solved draws (top: bands are mm wide, see bottom row)",
                 fontsize=9)
    fig.tight_layout()
    save(fig, fig_dir, "flux_surface_bands")


def save(fig, fig_dir, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(fig_dir, f"{name}.{ext}"), dpi=200)
    plt.close(fig)


def profile_figure(out, fig_dir):
    """p, ne, Te, Ti, iota and the enclosed toroidal current against s with
    truth, posterior median and 68/95 % bands. The current is CURTOR times the
    fixed shape of VMEC's enclosed current (buco of the truth), with CURTOR
    from the NUTS chain; the bands of the others are those of tiago_infer."""
    prof = read_csv(os.path.join(out, "profiles.csv"))
    truth = read_csv(os.path.join(out, "truth_profiles.csv"))
    method = "nuts" if any(r["method"] == "nuts" for r in prof) else "laplace"
    st = column(truth, "s")
    panels = [("p", "p [kPa]", 1e-3), ("ne", "ne [1e19 m^-3]", 1e-19), ("Te", "Te [keV]", 1e-3),
              ("Ti", "Ti [keV]", 1e-3), ("iota", "iota", 1.0), ("current", "enclosed I_tor [kA]", 1.0)]
    fig, axes = plt.subplots(2, 3, figsize=(11, 6.4), sharex=True)
    qs = (("q025", 0.025), ("q16", 0.16), ("q50", 0.5), ("q84", 0.84), ("q975", 0.975))
    for ax, (name, ylabel, scale) in zip(axes.flat, panels):
        if name == "iota":
            rows = read_csv(os.path.join(out, "iota.csv"))
            draws = sorted({int(r["draw"]) for r in rows if int(r["draw"]) > 0})
            s = column([r for r in rows if int(r["draw"]) == 0], "s")
            mat = np.array([column([r for r in rows if int(r["draw"]) == d], "iota") for d in draws])
            q = {k: np.quantile(mat, v, axis=0) for k, v in qs}
            tv = column([r for r in rows if int(r["draw"]) == 0], "iota")
        elif name == "current":
            import netCDF4
            with netCDF4.Dataset(os.path.join(out, "wout_truth.nc")) as w:
                buco, ctor = np.array(w["buco"][:]), float(w["ctor"][:])
            ns = len(buco)
            s = (np.arange(ns) - 0.5) / (ns - 1)
            s[0] = 0.0
            amps = column(read_csv(os.path.join(out, "nuts_draws.csv")), "CURTOR")
            amps_true = float(next(r for r in read_csv(os.path.join(out, "parameters.csv"))
                                   if r["name"].strip() == "CURTOR")["truth"])
            shape = np.sign(ctor) * 2e-3 * np.pi * buco / (4e-7 * np.pi) / amps_true   # kA per unit of scaled CURTOR
            q = {k: np.quantile(amps, v) * shape for k, v in qs}
            tv = amps_true * shape
        else:
            sel = [r for r in prof if r["method"] == method and r["quantity"] == name]
            s = column(sel, "s")
            q = {k: column(sel, k) * scale for k in dict(qs)}
            tv = column(truth, name) * scale
            tv_s = st
        ax.fill_between(s, q["q025"], q["q975"], color=BLUE, alpha=0.25, lw=0, label="95 % band")
        ax.fill_between(s, q["q16"], q["q84"], color=BLUE, alpha=0.5, lw=0, label="68 % band")
        ax.plot(s, q["q50"], color=BLUE, lw=1.2, label="posterior median")
        ax.plot(s if name in ("iota", "current") else st, tv, color=VERMILION, ls="--", lw=1.2,
                label="truth")
        ax.set_xlim(0.0, 1.0)
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.3)
    for ax in axes[1]:
        ax.set_xlabel("s (normalized toroidal flux)")
    axes[0, 0].legend(frameon=False, loc="upper right")
    fig.suptitle("Posterior profiles with 68 % / 95 % bands against the synthetic truth")
    fig.tight_layout()
    save(fig, fig_dir, "profiles_bands")


def corner(out, fig_dir, names=("CURTOR", "PHIEDGE", "ne_z1", "Te_z1", "Ti_z1", "align_1")):
    path = os.path.join(out, "nuts_draws.csv")
    if not os.path.exists(path):
        return
    draws = read_csv(path)
    params = {r["name"]: r for r in read_csv(os.path.join(out, "parameters.csv"))}
    names = [n for n in names if n in draws[0]]
    k = len(names)
    fig, axes = plt.subplots(k, k, figsize=(1.7 * k, 1.7 * k))
    for i, yi in enumerate(names):
        for j, xj in enumerate(names):
            ax = axes[i, j]
            if j > i:
                ax.axis("off")
                continue
            x = column(draws, xj)
            if i == j:
                ax.hist(x, bins=25, color=BLUE, alpha=0.6, density=True)
                m, sd = float(params[xj]["map"]), float(params[xj]["laplace_sd"])
                g = np.linspace(x.min(), x.max(), 200)
                ax.plot(g, np.exp(-0.5 * ((g - m) / sd) ** 2) / (sd * np.sqrt(2 * np.pi)),
                        color=GREEN, lw=1.0, ls="-.")
                ax.axvline(float(params[xj]["truth"]), color=VERMILION, ls="--", lw=1.0)
                ax.set_yticks([])
            else:
                y = column(draws, yi)
                ax.plot(x, y, ".", ms=1.5, color=BLUE, alpha=0.4)
                ax.plot(float(params[xj]["truth"]), float(params[yi]["truth"]), "x",
                        color=VERMILION, ms=6, mew=1.5)
            if i == k - 1:
                ax.set_xlabel(xj)
            else:
                ax.set_xticklabels([])
            if j == 0 and i > 0:
                ax.set_ylabel(yi)
            elif j > 0:
                ax.set_yticklabels([])
            ax.tick_params(labelsize=7)
    handles = [plt.Line2D([], [], color=BLUE, marker=".", ls="", label="NUTS draws"),
               plt.Line2D([], [], color=GREEN, ls="-.", label="Laplace marginal"),
               plt.Line2D([], [], color=VERMILION, ls="--", marker="x", label="truth")]
    fig.legend(handles=handles, loc="upper right", frameon=False)
    fig.suptitle("Leading latents (scaled CURTOR, PHIEDGE; whitened profile and alignment latents)",
                 fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, "corner.png"), dpi=160)
    plt.close(fig)


def main():
    out = sys.argv[1]
    fig_dir = sys.argv[2] if len(sys.argv) > 2 else out
    os.makedirs(fig_dir, exist_ok=True)
    flux_surfaces(out, fig_dir)
    profiles(out, fig_dir)
    surface_bands(out, fig_dir)
    profile_figure(out, fig_dir)
    corner(out, fig_dir)


if __name__ == "__main__":
    main()

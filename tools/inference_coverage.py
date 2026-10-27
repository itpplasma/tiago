#!/usr/bin/env python3
"""Coverage of a tiago_infer ensemble (one output directory per seed).

For each seed the truth is drawn from the prior, so a calibrated posterior
contains the truth in its central 68 % (95 %) interval with probability 0.68
(0.95) for every latent and every derived profile value. Counted are
  latents      Laplace: |truth - MAP| < 1 sd (0.683) and < 2 sd (0.954);
               NUTS: truth in [q16, q84] (0.68) and [q2.5, q97.5] (0.95);
  profiles     ne, Te, Ti, p at 41 s values: truth inside the Laplace and NUTS
               68/95 % bands;
  iota         truth inside the 68/95 % band of the re-solved posterior draws.
Values within a seed are correlated, so the binomial intervals (Wilson, over
all counted values) understate the uncertainty; the per-seed spread is also
reported. Writes coverage.txt and coverage.png to the ensemble directory.
Usage: MPLBACKEND=Agg inference_coverage.py <ensemble_dir> [seed_glob]
"""
import csv
import glob
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

BLUE, ORANGE, VERMILION = "#0072B2", "#E69F00", "#D55E00"


def rows(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def summary(path):
    out = {}
    with open(path) as f:
        for line in f:
            if "=" in line:
                k, v = line.split("=", 1)
                out[k.strip()] = v.strip()
    return out


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def seed_counts(d):
    """Hits per category for one seed: {category: (hits68, hits95, n)}."""
    out = {}
    par = rows(os.path.join(d, "parameters.csv"))
    t = np.array([float(r["truth"]) for r in par])
    m = np.array([float(r["map"]) for r in par])
    sd = np.array([float(r["laplace_sd"]) for r in par])
    z = np.abs(t - m) / sd
    out["latents Laplace"] = (np.sum(z < 1.0), np.sum(z < 2.0), len(z))
    if os.path.exists(os.path.join(d, "nuts_draws.csv")):
        q = {k: np.array([float(r[k]) for r in par]) for k in ("nuts_q025", "nuts_q16", "nuts_q84",
                                                                   "nuts_q975")}
        out["latents NUTS"] = (np.sum((t >= q["nuts_q16"]) & (t <= q["nuts_q84"])),
                               np.sum((t >= q["nuts_q025"]) & (t <= q["nuts_q975"])), len(t))
    truth = rows(os.path.join(d, "truth_profiles.csv"))
    prof = rows(os.path.join(d, "profiles.csv"))
    for method in ("laplace", "nuts"):
        h68 = h95 = n = 0
        for qn in ("ne", "Te", "Ti", "p"):
            sel = [r for r in prof if r["method"] == method and r["quantity"] == qn]
            if not sel:
                continue
            tv = np.array([float(r[qn]) for r in truth])
            lo68, hi68 = (np.array([float(r[k]) for r in sel]) for k in ("q16", "q84"))
            lo95, hi95 = (np.array([float(r[k]) for r in sel]) for k in ("q025", "q975"))
            h68 += np.sum((tv >= lo68) & (tv <= hi68))
            h95 += np.sum((tv >= lo95) & (tv <= hi95))
            n += len(tv)
        if n:
            out[f"profiles {'Laplace' if method == 'laplace' else 'NUTS'}"] = (h68, h95, n)
    io = rows(os.path.join(d, "iota.csv"))
    draws = sorted({int(r["draw"]) for r in io if int(r["draw"]) > 0})
    if len(draws) >= 10:
        tv = np.array([float(r["iota"]) for r in io if int(r["draw"]) == 0])
        mat = np.array([[float(r["iota"]) for r in io if int(r["draw"]) == k] for k in draws])
        q = np.quantile(mat, [0.025, 0.16, 0.84, 0.975], axis=0)
        out["iota draws"] = (np.sum((tv >= q[1]) & (tv <= q[2])), np.sum((tv >= q[0]) & (tv <= q[3])),
                             len(tv))
    return out


def main():
    ens = sys.argv[1]
    pattern = sys.argv[2] if len(sys.argv) > 2 else "seed*"
    seeds = sorted(d for d in glob.glob(os.path.join(ens, pattern, "out"))
                   if os.path.exists(os.path.join(d, "summary.txt")))
    if not seeds:
        sys.exit("no finished seeds")
    total, per_seed, summ = {}, {}, []
    for d in seeds:
        c = seed_counts(d)
        for k, (a, b, n) in c.items():
            T = total.setdefault(k, [0, 0, 0])
            T[0] += a
            T[1] += b
            T[2] += n
            per_seed.setdefault(k, []).append((a / n, b / n))
        summ.append(summary(os.path.join(d, "summary.txt")))
    lines = [f"seeds: {len(seeds)} ({ens})", "",
             f"{'category':18s} {'68%':>7s} {'Wilson 95% CI':>16s} {'95%':>7s} {'Wilson 95% CI':>16s}"
             f" {'n':>6s}  per-seed sd(68/95)"]
    for k, (a, b, n) in total.items():
        w1, w2 = wilson(a, n), wilson(b, n)
        ps = np.array(per_seed[k])
        lines.append(f"{k:18s} {a / n:7.3f} [{w1[0]:.3f}, {w1[1]:.3f}] {b / n:7.3f} "
                     f"[{w2[0]:.3f}, {w2[1]:.3f}] {n:6d}  {ps[:, 0].std():.3f}/{ps[:, 1].std():.3f}")
    lines.append("nominal: 0.683/0.954 for the Laplace +-1/2 sd intervals, 0.68/0.95 otherwise")
    lines.append("")
    keys = ["psis_khat", "nuts_divergent", "nuts_mean_accept", "nuts_mean_depth", "nuts_ebfmi",
            "log_evidence", "chi2_map", "time_map_s", "time_nuts_s", "equilibrium_solves"]
    for k in keys:
        v = np.array([float(s[k]) for s in summ if k in s])
        if len(v):
            lines.append(f"{k:20s} min {v.min():12.4g} median {np.median(v):12.4g} max {v.max():12.4g}")
    kh = np.array([float(s["psis_khat"]) for s in summ if "psis_khat" in s])
    lines.append(f"psis k-hat > 0.7 in {np.sum(kh > 0.7)} of {len(kh)} seeds")
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(ens, "coverage.txt"), "w") as f:
        f.write(text + "\n")

    fig, ax = plt.subplots(figsize=(8, 4))
    cats = list(total)
    x = np.arange(len(cats))
    for off, level, color, marker, nom in ((-0.12, 0, BLUE, "o", 0.68), (0.12, 1, ORANGE, "s", 0.95)):
        frac = [total[k][level] / total[k][2] for k in cats]
        err = np.array([wilson(total[k][level], total[k][2]) for k in cats]).T
        ax.errorbar(x + off, frac, yerr=[np.array(frac) - err[0], err[1] - np.array(frac)], fmt=marker,
                    color=color, capsize=3, label=f"central {int(nom * 100)} % interval")
        ax.axhline(nom, color=color, ls="--", lw=1.0)
    ax.set_xticks(x)
    ax.set_xticklabels(cats, rotation=20, ha="right")
    ax.set_ylabel("fraction of truths inside")
    ax.set_ylim(0.4, 1.02)
    ax.set_title(f"Coverage over {len(seeds)} prior-drawn truths (dashed: nominal)")
    ax.legend(frameon=False, loc="lower left")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(ens, "coverage.png"), dpi=160)


if __name__ == "__main__":
    main()

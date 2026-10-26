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
Values within a seed are correlated (the iota values at all s of one seed
follow one dominant latent), so the binomial intervals (Wilson, over all
counted values) understate the uncertainty. The seed bootstrap interval
resamples whole seeds and is the one to read. A band made from the empirical
quantiles of N draws covers a calibrated truth with less than the nominal
probability (N = 40: 0.649/0.905, N = 400: 0.676/0.946), so the expected
coverage for the actual N is reported as well (Monte Carlo, independent
draws; thinned NUTS draws are close to that). Writes coverage.txt and
coverage.png to the ensemble directory.
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


def seed_counts(d, draws_used):
    """Hits per category for one seed: {category: (hits68, hits95, n)};
    draws_used receives the number of draws behind each draw-quantile band."""
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
    if os.path.exists(os.path.join(d, "nuts_draws.csv")):
        nd = len(rows(os.path.join(d, "nuts_draws.csv")))
        draws_used["latents NUTS"] = draws_used["profiles NUTS"] = nd
    io = rows(os.path.join(d, "iota.csv"))
    draws = sorted({int(r["draw"]) for r in io if int(r["draw"]) > 0})
    if len(draws) >= 10:
        tv = np.array([float(r["iota"]) for r in io if int(r["draw"]) == 0])
        mat = np.array([[float(r["iota"]) for r in io if int(r["draw"]) == k] for k in draws])
        q = np.quantile(mat, [0.025, 0.16, 0.84, 0.975], axis=0)
        out["iota draws"] = (np.sum((tv >= q[1]) & (tv <= q[2])), np.sum((tv >= q[0]) & (tv <= q[3])),
                             len(tv))
        draws_used["iota draws"] = len(draws)
    return out


def finite_draw_nominal(n_draws, trials=200000):
    """Expected 68/95 % coverage of the np.quantile band of n_draws
    independent draws for a truth from the same distribution."""
    rng = np.random.default_rng(12345)
    hit = np.zeros(2)
    done = 0
    while done < trials:
        m = min(20000, trials - done)
        x = rng.standard_normal((m, n_draws))
        t = rng.standard_normal(m)
        q = np.quantile(x, [0.025, 0.16, 0.84, 0.975], axis=1)
        hit += [np.sum((t >= q[1]) & (t <= q[2])), np.sum((t >= q[0]) & (t <= q[3]))]
        done += m
    return hit / trials


def seed_bootstrap(hits, level, reps=4000):
    """95 % interval of the pooled fraction when whole seeds are resampled."""
    h = np.array([x[level] for x in hits], float)
    n = np.array([x[2] for x in hits], float)
    rng = np.random.default_rng(2718)
    idx = rng.integers(0, len(h), (reps, len(h)))
    frac = h[idx].sum(1) / n[idx].sum(1)
    return np.quantile(frac, [0.025, 0.975])


def main():
    ens = sys.argv[1]
    pattern = sys.argv[2] if len(sys.argv) > 2 else "seed*"
    seeds = sorted(d for d in glob.glob(os.path.join(ens, pattern, "out"))
                   if os.path.exists(os.path.join(d, "summary.txt")))
    if not seeds:
        sys.exit("no finished seeds")
    total, per_seed, summ, draws_used = {}, {}, [], {}
    for d in seeds:
        c = seed_counts(d, draws_used)
        for k, (a, b, n) in c.items():
            T = total.setdefault(k, [0, 0, 0])
            T[0] += a
            T[1] += b
            T[2] += n
            per_seed.setdefault(k, []).append((a, b, n))
        summ.append(summary(os.path.join(d, "summary.txt")))
    lines = [f"seeds: {len(seeds)} ({ens})", "",
             f"{'category':18s} {'68%':>7s} {'seed-bootstrap CI':>16s} {'Wilson CI':>16s} {'95%':>7s} "
             f"{'seed-bootstrap CI':>16s} {'Wilson CI':>16s} {'n':>6s}  expected (draws)"]
    nominal = {}
    for k, (a, b, n) in total.items():
        w1, w2 = wilson(a, n), wilson(b, n)
        b1, b2 = seed_bootstrap(per_seed[k], 0), seed_bootstrap(per_seed[k], 1)
        if k == "latents Laplace":
            nominal[k], tag = (0.683, 0.954), "(Gaussian sd)"
        elif k in draws_used:
            nominal[k], tag = tuple(finite_draw_nominal(draws_used[k])), f"({draws_used[k]})"
        else:
            nominal[k], tag = (0.68, 0.95), ""
        lines.append(f"{k:18s} {a / n:7.3f} [{b1[0]:.3f}, {b1[1]:.3f}]   [{w1[0]:.3f}, {w1[1]:.3f}] "
                     f"{b / n:7.3f} [{b2[0]:.3f}, {b2[1]:.3f}]   [{w2[0]:.3f}, {w2[1]:.3f}] {n:6d}  "
                     f"{nominal[k][0]:.3f}/{nominal[k][1]:.3f} {tag}")
    lines.append("expected: coverage of a calibrated posterior for the interval construction and the "
                 "number of draws used (all 95 % intervals are 95 % confidence intervals)")
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
        err = np.array([seed_bootstrap(per_seed[k], level) for k in cats]).T
        ax.errorbar(x + off, frac, yerr=[np.array(frac) - err[0], err[1] - np.array(frac)], fmt=marker,
                    color=color, capsize=3, label=f"central {int(nom * 100)} % interval")
        ax.axhline(nom, color=color, ls="--", lw=1.0)
        ax.plot(x + off + 0.07, [nominal[k][level] for k in cats], "<", mfc="none", color=color, ms=6,
                ls="none", label=f"expected {int(nom * 100)} % for the draws used")
    ax.set_xticks(x)
    ax.set_xticklabels(cats, rotation=20, ha="right")
    ax.set_ylabel("fraction of truths inside")
    ax.set_ylim(0.4, 1.02)
    ax.set_title(f"Coverage over {len(seeds)} prior-drawn truths (dashed: nominal; error bars: "
                 "seed bootstrap 95 %)", fontsize=10)
    ax.legend(frameon=False, loc="lower left", fontsize=8, ncol=2)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(ens, "coverage.png"), dpi=160)
    fig.savefig(os.path.join(ens, "coverage.pdf"))
    plt.close(fig)
    parameter_coverage(ens, seeds)


def parameter_coverage(ens, seeds):
    """Coverage of each latent (Laplace and NUTS intervals) over the seeds
    with Wilson 95 % intervals (one binomial trial per seed), and the
    median Laplace sd relative to the prior sd (CURTOR 2e4 A of 174250 A and
    PHIEDGE 0.02 of 0.514386 in scaled units; the other latents are N(0,1))."""
    prior_sd = {"CURTOR": 2.0e4 / 174250.0, "PHIEDGE": 0.02 / 0.514386}
    acc = {}
    for d in seeds:
        for r in rows(os.path.join(d, "parameters.csv")):
            acc.setdefault(r["name"].strip(), []).append({k: float(v) for k, v in r.items()
                                                          if k not in ("name", None) and v.strip()})
    names = list(acc)
    res = {}
    for n in names:
        L = acc[n]
        t = np.array([x["truth"] for x in L])
        z = np.abs(t - np.array([x["map"] for x in L])) / np.array([x["laplace_sd"] for x in L])
        nu = lambda a, b: sum(x[a] <= x["truth"] <= x[b] for x in L)
        res[n] = {"Laplace": (np.sum(z < 1), np.sum(z < 2)),
                  "NUTS": (nu("nuts_q16", "nuts_q84"), nu("nuts_q025", "nuts_q975")), "n": len(L),
                  "ratio": np.median([x["laplace_sd"] for x in L]) / prior_sd.get(n, 1.0)}
    lines = [f"{'latent':9s} {'sd/prior':>8s}  {'Laplace 68 %':>18s} {'Laplace 95 %':>18s} "
             f"{'NUTS 68 %':>18s} {'NUTS 95 %':>18s}   (Wilson 95 % CI, n = seeds)"]
    for n in names:
        r = res[n]
        cells = []
        for m in ("Laplace", "NUTS"):
            for lv in (0, 1):
                lo, hi = wilson(r[m][lv], r["n"])
                cells.append(f"{r[m][lv] / r['n']:.2f} [{lo:.2f},{hi:.2f}]")
        lines.append(f"{n:9s} {r['ratio']:8.4f}  " + " ".join(f"{c:>18s}" for c in cells))
    with open(os.path.join(ens, "coverage_parameters.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(10, 6.2), sharex=True,
                                 gridspec_kw={"height_ratios": [3, 1.1]})
    x = np.arange(len(names))
    for off, lv, color, marker, nom in ((-0.15, 0, BLUE, "o", 0.683), (0.15, 1, ORANGE, "s", 0.954)):
        frac = np.array([res[n]["Laplace"][lv] / res[n]["n"] for n in names])
        ci = np.array([wilson(res[n]["Laplace"][lv], res[n]["n"]) for n in names]).T
        ax.errorbar(x + off, frac, yerr=[frac - ci[0], ci[1] - frac], fmt=marker, color=color,
                    capsize=2, ms=4, label=f"Laplace, central {int(round(nom * 100))} % interval")
        ax.axhline(nom, color=color, ls="--", lw=1.0)
    ax.set_ylabel("fraction of truths inside")
    ax.set_ylim(0.3, 1.03)
    ax.legend(frameon=False, loc="lower left", fontsize=8, ncol=2)
    ax.grid(axis="y", alpha=0.3)
    ax.set_title(f"Coverage per latent over {len(seeds)} prior-drawn truths (error bars: Wilson 95 %; "
                 "dashed: nominal)", fontsize=10)
    bx.bar(x, [res[n]["ratio"] for n in names], color=VERMILION, alpha=0.8)
    bx.set_ylabel("posterior sd / prior sd")
    bx.set_ylim(0, 1.05)
    bx.grid(axis="y", alpha=0.3)
    bx.set_xticks(x)
    bx.set_xticklabels(names, rotation=60, ha="right", fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(ens, "coverage_parameters.png"), dpi=200)
    fig.savefig(os.path.join(ens, "coverage_parameters.pdf"))
    plt.close(fig)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Thomson calibration identifiability study (tiago_infer calibration_case > 0).

<root> holds one directory per variant (e.g. c1_up, c2_dn, nolat2_up), each
with seed*/out from tiago_infer. Per variant this reports for log C (Thomson
campaign log scale): the posterior mean and sd over seeds (NUTS if present,
else Laplace), the contraction sd_post / sd_prior and the information
log2(sd_prior / sd_post) in bits, the bias (mean - truth) / sd, and the 68/95 %
coverage of the truth (Wilson intervals over seeds); for the ne profile the
fraction of (seed, s) values inside the 68/95 % bands. With a fixed injected C
the log C coverage is a frequentist check at that truth (nominal for an
identified C, near zero when the prior dominates and C is off the prior mean).
Writes calibration.txt, posterior_C.png/pdf, coverage_C.png/pdf and
profile_bands.png/pdf (ne bands of one seed for the variants given by
--bands, with and without the calibration latent) to <root>.
Usage: MPLBACKEND=Agg calibration_study.py <root> [--bands c2_up,nolat2_up]
       [--seed 101]
"""
import argparse
import csv
import glob
import math
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

COLORS = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#000000"]


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
        return (math.nan, math.nan)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def seed_result(out):
    s = summary(os.path.join(out, "summary.txt"))
    p = {r["name"].strip(): r for r in rows(os.path.join(out, "parameters.csv"))}
    res = {"truth_C": float(s["ts_cal_truth"])}
    if "logC_TS" in p:
        mu, sig = float(s["ts_cal_mean"]), float(s["ts_cal_sigma"])
        r = {k: float(v) for k, v in p["logC_TS"].items() if k not in ("name", None) and v}
        nuts = r["nuts_sd"] > 0
        z = r["nuts_mean"] if nuts else r["map"]
        sd = r["nuts_sd"] if nuts else r["laplace_sd"]
        res.update(
            mean=mu + sig * z,
            sd=sig * sd,
            prior_sd=sig,
            truth=mu + sig * r["truth"],
            in68=(r["nuts_q16"] <= r["truth"] <= r["nuts_q84"]) if nuts
            else abs(r["truth"] - r["map"]) < r["laplace_sd"],
            in95=(r["nuts_q025"] <= r["truth"] <= r["nuts_q975"]) if nuts
            else abs(r["truth"] - r["map"]) < 2 * r["laplace_sd"],
        )
    truth = {float(r["s"]): float(r["ne"]) for r in rows(os.path.join(out, "truth_profiles.csv"))}
    prof = [r for r in rows(os.path.join(out, "profiles.csv")) if r["quantity"].strip() == "ne"]
    method = "nuts" if any(r["method"].strip() == "nuts" for r in prof) else "laplace"
    k68 = k95 = n = 0
    for r in prof:
        if r["method"].strip() != method:
            continue
        t = truth[float(r["s"])]
        k68 += float(r["q16"]) <= t <= float(r["q84"])
        k95 += float(r["q025"]) <= t <= float(r["q975"])
        n += 1
    res.update(ne68=k68, ne95=k95, ne_n=n)
    return res


def variant(path):
    out = []
    for d in sorted(glob.glob(os.path.join(path, "seed*", "out"))):
        try:
            out.append(seed_result(d))
        except (OSError, KeyError, ValueError):
            pass
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("--bands", default="")
    ap.add_argument("--seed", type=int, default=101)
    a = ap.parse_args()
    names = sorted(d for d in os.listdir(a.root) if glob.glob(os.path.join(a.root, d, "seed*", "out", "summary.txt")))
    data = {v: variant(os.path.join(a.root, v)) for v in names}
    lines = ["variant  n  C_true  logC_mean  post_sd  sd/prior  bits  bias/sd  cov68 [95%CI]  cov95 [95%CI]  ne68  ne95"]
    for v in names:
        rs = data[v]
        if not rs:
            continue
        ne68 = sum(r["ne68"] for r in rs) / max(sum(r["ne_n"] for r in rs), 1)
        ne95 = sum(r["ne95"] for r in rs) / max(sum(r["ne_n"] for r in rs), 1)
        c = [r for r in rs if "mean" in r]
        if c:
            n = len(c)
            sd = np.median([r["sd"] for r in c])
            ratio = sd / c[0]["prior_sd"]
            bias = np.mean([(r["mean"] - r["truth"]) / r["sd"] for r in c])
            k68, k95 = sum(r["in68"] for r in c), sum(r["in95"] for r in c)
            lo68, hi68 = wilson(k68, n)
            lo95, hi95 = wilson(k95, n)
            lines.append(
                f"{v:9s} {n:3d} {c[0]['truth_C']:6.3f} {np.mean([r['mean'] for r in c]):9.4f} {sd:8.4f} "
                f"{ratio:8.3f} {math.log2(1 / ratio):5.2f} {bias:7.2f}  {k68 / n:.2f} [{lo68:.2f},{hi68:.2f}]"
                f"  {k95 / n:.2f} [{lo95:.2f},{hi95:.2f}]  {ne68:.3f} {ne95:.3f}")
        else:
            lines.append(f"{v:9s} {len(rs):3d} {rs[0]['truth_C']:6.3f}  (no calibration latent)"
                         f"{'':52s}{ne68:.3f} {ne95:.3f}")
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(a.root, "calibration.txt"), "w") as f:
        f.write(text + "\n")

    cal = [v for v in names if any("mean" in r for r in data[v])]
    if cal:
        fig, ax = plt.subplots(figsize=(max(6, 1.1 * len(cal)), 4), constrained_layout=True)
        for i, v in enumerate(cal):
            c = [r for r in data[v] if "mean" in r]
            x = i + np.linspace(-0.3, 0.3, len(c))
            m = np.exp([r["mean"] for r in c])
            lo = np.exp([r["mean"] - r["sd"] for r in c])
            hi = np.exp([r["mean"] + r["sd"] for r in c])
            ax.errorbar(x, m, yerr=[m - lo, hi - m], fmt="o", ms=2, lw=0.6, color=COLORS[i % len(COLORS)])
            ax.hlines(c[0]["truth_C"], i - 0.4, i + 0.4, color="k", lw=1.2)
        ax.set_yscale("log")
        ax.set_xticks(range(len(cal)), cal, rotation=30)
        ax.set_ylabel("posterior C (mean, 68 %) per seed; truth black")
        for ext in ("png", "pdf"):
            fig.savefig(os.path.join(a.root, f"posterior_C.{ext}"), dpi=150)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(max(6, 1.1 * len(cal)), 3.5), constrained_layout=True)
        for j, (key, nominal) in enumerate((("in68", 0.68), ("in95", 0.95))):
            for i, v in enumerate(cal):
                c = [r for r in data[v] if "mean" in r]
                k = sum(r[key] for r in c)
                lo, hi = wilson(k, len(c))
                ax.errorbar(i + 0.15 * (2 * j - 1), k / len(c), yerr=[[k / len(c) - lo], [hi - k / len(c)]],
                            fmt="os"[j], color=COLORS[j], label=f"{int(nominal * 100)} %" if i == 0 else None)
            ax.axhline(nominal, color=COLORS[j], ls=":", lw=1)
        ax.set_xticks(range(len(cal)), cal, rotation=30)
        ax.set_ylim(-0.02, 1.02)
        ax.set_ylabel("coverage of log C (Wilson 95 %)")
        ax.legend()
        for ext in ("png", "pdf"):
            fig.savefig(os.path.join(a.root, f"coverage_C.{ext}"), dpi=150)
        plt.close(fig)

    bands = [b for b in a.bands.split(",") if b]
    if bands:
        fig, axs = plt.subplots(1, len(bands), figsize=(4 * len(bands), 3.5), sharey=True,
                                constrained_layout=True, squeeze=False)
        for ax, v in zip(axs[0], bands):
            out = os.path.join(a.root, v, f"seed{a.seed}", "out")
            if not os.path.isdir(out):
                ax.set_title(f"{v}: seed {a.seed} missing")
                continue
            prof = [r for r in rows(os.path.join(out, "profiles.csv")) if r["quantity"].strip() == "ne"]
            method = "nuts" if any(r["method"].strip() == "nuts" for r in prof) else "laplace"
            prof = [r for r in prof if r["method"].strip() == method]
            s = np.array([float(r["s"]) for r in prof])
            q = {k: np.array([float(r[k]) for r in prof]) for k in ("q025", "q16", "q50", "q84", "q975")}
            ax.fill_between(s, q["q025"], q["q975"], color=COLORS[0], alpha=0.2, lw=0, label="95 %")
            ax.fill_between(s, q["q16"], q["q84"], color=COLORS[0], alpha=0.4, lw=0, label="68 %")
            ax.plot(s, q["q50"], color=COLORS[0], lw=1, label="median")
            t = rows(os.path.join(out, "truth_profiles.csv"))
            ax.plot([float(r["s"]) for r in t], [float(r["ne"]) for r in t], "k--", lw=1, label="truth")
            ax.set_title(f"{v} ({method}, seed {a.seed})")
            ax.set_xlabel("s")
        axs[0][0].set_ylabel("ne [m$^{-3}$]")
        axs[0][0].legend(fontsize=8)
        for ext in ("png", "pdf"):
            fig.savefig(os.path.join(a.root, f"profile_bands.{ext}"), dpi=150)
        plt.close(fig)


if __name__ == "__main__":
    main()

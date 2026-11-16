"""Compare TIAGO soft X-ray chords with STELLOPT's chisq_sxr output.

Usage: python3 compare.py <stellopt.ext> <tiago.csv> <output directory> [figure]

STELLOPT writes each signal divided by the largest one, so both codes are
compared after that normalization. Writes results.md and results.csv; with a
fourth argument also a figure (matplotlib).
"""

import csv
import statistics
import sys
from pathlib import Path


def stellopt_values(path):
    """Normalized VAL column of the last SXR block."""
    lines = Path(path).read_text().splitlines()
    starts = [i for i, line in enumerate(lines) if line.split()[:1] == ["SXR"]]
    if not starts:
        raise SystemExit(f"{path}: no SXR block")
    first = starts[-1]
    count = int(lines[first].split()[1])
    rows = lines[first + 2:first + 2 + count]
    return [float(row.split()[2]) for row in rows]


def tiago_values(path):
    with open(path, newline="") as handle:
        rows = list(csv.reader(handle))
    header = rows[0]
    labels = [row[1] for row in rows[1:]]
    columns = {name: [float(row[k]) for row in rows[1:]]
               for k, name in enumerate(header) if name.startswith("npanel_")}
    return labels, columns


def main():
    stellopt = stellopt_values(sys.argv[1])
    labels, columns = tiago_values(sys.argv[2])
    out = Path(sys.argv[3])
    finest = columns[list(columns)[-1]]
    coarse = columns[list(columns)[-2]]
    if len(finest) != len(stellopt):
        raise SystemExit("different numbers of chords")
    peak = max(finest)
    tiago = [v / peak for v in finest]
    diff = [t - s for t, s in zip(tiago, stellopt)]
    own = [abs(f - c) / peak for f, c in zip(finest, coarse)]
    with open(out / "results.csv", "w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["label", "stellopt", "tiago", "difference", "tiago_convergence"])
        for row in zip(labels, stellopt, tiago, diff, own):
            writer.writerow([row[0]] + [f"{v:.12e}" for v in row[1:]])
    cameras = sorted({label.rsplit("_", 1)[0] for label in labels},
                     key=[label.rsplit("_", 1)[0] for label in labels].index)
    lines = ["| camera | channels | max abs difference | median abs difference |",
             "| --- | --- | --- | --- |"]
    for cam in cameras:
        d = [abs(x) for x, label in zip(diff, labels) if label.rsplit("_", 1)[0] == cam]
        lines.append(f"| {cam} | {len(d)} | {max(d):.2e} | {statistics.median(d):.2e} |")
    lines.append(f"| all | {len(diff)} | {max(map(abs, diff)):.2e} | "
                 f"{statistics.median(map(abs, diff)):.2e} |")
    lines.append("")
    lines.append("Differences of signals normalized to the largest signal. TIAGO "
                 "change between its two finest rules (same normalization): "
                 f"max {max(own):.1e}.")
    (out / "results.md").write_text("\n".join(lines) + "\n")
    if len(sys.argv) > 4:
        plot(labels, stellopt, tiago, diff, sys.argv[4])


def plot(labels, stellopt, tiago, diff, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    channel = list(range(1, len(labels) + 1))
    fig, (top, bottom) = plt.subplots(2, 1, sharex=True, figsize=(7, 5),
                                      height_ratios=[2, 1])
    top.plot(channel, stellopt, "o", mfc="none", label="STELLOPT chisq_sxr")
    top.plot(channel, tiago, ".", label="TIAGO soft_xray_t")
    top.set_ylabel("signal / largest signal")
    top.legend()
    bottom.plot(channel, diff, ".", color="k")
    bottom.set_ylabel("TIAGO - STELLOPT")
    bottom.set_xlabel("channel (1-12 top, 13-24 outboard, 25-36 oblique)")
    bottom.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    fig.suptitle("Bremsstrahlung chords on LI383, ne = 1e20 (1 - s), Te = 2 keV (1 - s)")
    fig.tight_layout()
    fig.savefig(path)


if __name__ == "__main__":
    main()

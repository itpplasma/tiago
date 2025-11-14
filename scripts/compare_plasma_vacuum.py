#!/usr/bin/env python3
"""Compare TIAGO Biot–Savart plasma response against simsopt VirtualCasing output."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def load_field_csv(path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    if data.shape[1] == 3:  # Bx, By, Bz only
        points = None
        components = data
    elif data.shape[1] == 5:  # iphi, itheta, Bx, By, Bz
        points = None
        components = data[:, 2:5]
    elif data.shape[1] >= 8:  # iphi, itheta, X, Y, Z, Bx, By, Bz
        points = data[:, 2:5]
        components = data[:, -3:]
    else:
        raise ValueError(f"Unexpected column count in {path}")
    return components, points


def compare_fields(tiago_csv: Path, simsopt_total_csv: Path, simsopt_external_csv: Path,
                   output_dir: Path, label: str) -> None:
    tiago_field, _ = load_field_csv(tiago_csv)
    simsopt_total, _ = load_field_csv(simsopt_total_csv)
    simsopt_external, _ = load_field_csv(simsopt_external_csv)

    if (tiago_field.shape != simsopt_total.shape or
            simsopt_total.shape != simsopt_external.shape):
        raise ValueError("Field arrays have different shapes")

    simsopt_plasma = simsopt_total - simsopt_external

    delta = tiago_field - simsopt_plasma
    abs_err = np.linalg.norm(delta, axis=1)
    rel_err = abs_err / np.maximum(np.linalg.norm(simsopt_plasma, axis=1), 1e-12)

    summary = {
        "max_abs_error": float(np.max(abs_err)),
        "rms_abs_error": float(np.sqrt(np.mean(abs_err ** 2))),
        "max_rel_error": float(np.max(rel_err)),
        "median_rel_error": float(np.median(rel_err)),
    }

    print("[compare_plasma_vacuum] Summary:")
    for key, value in summary.items():
        print(f"  {key}: {value:.6e}")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    simsopt_mag = np.linalg.norm(simsopt_plasma, axis=1)
    tiago_mag = np.linalg.norm(tiago_field, axis=1)
    axes[0].scatter(simsopt_mag * 1e3, tiago_mag * 1e3, s=8, alpha=0.6, edgecolors="none")
    axes[0].plot([0, max(simsopt_mag.max(), tiago_mag.max()) * 1e3],
                 [0, max(simsopt_mag.max(), tiago_mag.max()) * 1e3],
                 linestyle="--", color="black", linewidth=1)
    axes[0].set_xlabel("|B_plasma| simsopt (mT)")
    axes[0].set_ylabel("|B_plasma| TIAGO (mT)")
    axes[0].set_title("Plasma Response Magnitude")
    axes[0].grid(True, linestyle=':', alpha=0.4)

    axes[1].hist(rel_err * 100.0, bins=40, color="#3498db", alpha=0.8)
    axes[1].set_xlabel("Relative Error (%)")
    axes[1].set_ylabel("Count")
    axes[1].set_title("Relative Error Distribution")
    axes[1].axvline(1.0, color="red", linestyle="--", label="1%")
    axes[1].legend()

    output_dir.mkdir(parents=True, exist_ok=True)
    plot_path = output_dir / f"plasma_response_comparison_{label}.png"
    fig.tight_layout()
    fig.savefig(plot_path, dpi=180)
    plt.close(fig)
    print(f"[compare_plasma_vacuum] Saved plot to {plot_path}")

    summary_path = output_dir / f"plasma_response_summary_{label}.txt"
    with summary_path.open("w") as handle:
        handle.write("TIAGO vs simsopt plasma response comparison\n")
        handle.write("=" * 60 + "\n")
        handle.write(f"TIAGO CSV: {tiago_csv}\n")
        handle.write(f"simsopt B_total CSV: {simsopt_total_csv}\n")
        handle.write(f"simsopt B_external CSV: {simsopt_external_csv}\n\n")
        for key, value in summary.items():
            handle.write(f"{key}: {value:.6e}\n")
    print(f"[compare_plasma_vacuum] Summary written to {summary_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare TIAGO and simsopt plasma response")
    parser.add_argument("--tiago-bext", required=True, type=Path,
                        help="CSV produced by test_tiago_with_simsopt_grid")
    parser.add_argument("--simsopt-btotal", required=True, type=Path,
                        help="simsopt B_total CSV from generate_simsopt_reference.py")
    parser.add_argument("--simsopt-bexternal", required=True, type=Path,
                        help="simsopt B_external CSV from generate_simsopt_reference.py")
    parser.add_argument("--output-dir", required=True, type=Path,
                        help="Directory for plots and summaries")
    parser.add_argument("--label", default="case", help="Label for output files")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    compare_fields(args.tiago_bext, args.simsopt_btotal, args.simsopt_bexternal,
                   args.output_dir, args.label)


if __name__ == "__main__":
    main()

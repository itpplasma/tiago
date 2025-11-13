#!/usr/bin/env python3
"""Generate side-by-side TIAGO vs simsopt B_external plots on the NCSX grid."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    from simsopt.mhd.vmec import Vmec
    from simsopt.mhd.virtual_casing import VirtualCasing
except ImportError as exc:
    print("ERROR: simsopt is required for this comparison. "
          "Install with `pip install simsopt`.", file=sys.stderr)
    raise


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot TIAGO vs simsopt B_external magnitude and differences.")
    parser.add_argument("--tiago-bin", required=True,
                        help="Path to test_tiago_with_simsopt_grid executable.")
    parser.add_argument("--gamma", required=True,
                        help="CSV file with simsopt gamma grid (X,Y,Z).")
    parser.add_argument("--b-total", required=True,
                        help="CSV file with simsopt B_total vectors (Bx,By,Bz).")
    parser.add_argument("--wout", required=True,
                        help="VMEC wout file for simsopt VirtualCasing.")
    parser.add_argument("--output-dir", default="build/tests/output",
                        help="Directory for plots + intermediate CSVs.")
    parser.add_argument("--src-nphi", type=int, default=16,
                        help="Toroidal grid samples (default: 16).")
    parser.add_argument("--src-ntheta", type=int, default=16,
                        help="Poloidal grid samples (default: 16).")
    parser.add_argument("--debug", action="store_true",
                        help="Print detailed field statistics for troubleshooting.")
    return parser.parse_args()


def run_command(cmd: list[str], verbose: bool = False) -> subprocess.CompletedProcess[str]:
    if verbose:
        print(f"[DEBUG] Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        print(result.stdout, file=sys.stderr)
        print(result.stderr, file=sys.stderr)
        raise RuntimeError(f"Command failed: {' '.join(cmd)}")
    if verbose:
        print("[DEBUG] stdout:")
        print(result.stdout.strip())
        if result.stderr.strip():
            print("[DEBUG] stderr:")
            print(result.stderr.strip())
    return result


def load_tiago_field(bext_csv: Path, nphi: int, ntheta: int) -> np.ndarray:
    df = pd.read_csv(bext_csv)
    vectors = df[["Bx", "By", "Bz"]].to_numpy()
    if vectors.shape[0] != nphi * ntheta:
        raise ValueError("TIAGO output does not match requested grid dimensions.")
    return vectors.reshape((nphi, ntheta, 3))


def compute_simsopt_field(wout_file: str, nphi: int, ntheta: int
                          ) -> tuple[np.ndarray, np.ndarray, float, int]:
    vmec = Vmec(wout_file)
    vc = VirtualCasing.from_vmec(vmec, src_nphi=nphi, src_ntheta=ntheta,
                                 use_stellsym=True, digits=6, filename=None)
    gamma = np.asarray(vc.gamma)
    bext = np.asarray(vc.B_external)
    rms = float(np.sqrt(np.mean(bext**2)))
    return (gamma.reshape(nphi, ntheta, 3),
            bext.reshape(nphi, ntheta, 3),
            rms,
            int(vmec.wout.nfp))


def log_field_stats(tag: str, field: np.ndarray, debug: bool) -> None:
    if not debug:
        return
    mag = np.linalg.norm(field, axis=2)
    print(f"[DEBUG] {tag}: shape={field.shape}, "
          f"min={mag.min():.6e} T, max={mag.max():.6e} T, "
          f"rms={np.sqrt(np.mean(mag**2)):.6e} T")


def save_vector_csv(path: Path, field: np.ndarray) -> None:
    nphi, ntheta, _ = field.shape
    with path.open("w", encoding="utf-8") as handle:
        handle.write("iphi,itheta,Bx,By,Bz\n")
        for iphi in range(nphi):
            for itheta in range(ntheta):
                bx, by, bz = field[iphi, itheta, :]
                handle.write(f"{iphi+1},{itheta+1},"
                             f"{bx:.15e},{by:.15e},{bz:.15e}\n")


def make_plots(gamma: np.ndarray, tiago_field: np.ndarray, simsopt_field: np.ndarray,
               nfp: int, output_dir: Path,
               prefix: str = "tiago_vs_simsopt") -> None:
    tiago_mag = np.linalg.norm(tiago_field, axis=2)
    simsopt_mag = np.linalg.norm(simsopt_field, axis=2)
    diff_mag = tiago_mag - simsopt_mag

    x_coords = gamma[:, :, 0]
    z_coords = gamma[:, :, 2]

    fig, axes = plt.subplots(1, 3, figsize=(18, 6), constrained_layout=True)
    plots = [
        (tiago_mag, "TIAGO |B_ext| (T)", "viridis"),
        (simsopt_mag, "simsopt |B_ext| (T)", "viridis"),
        (diff_mag, "TIAGO - simsopt (T)", "coolwarm"),
    ]

    for ax, (field, title, cmap) in zip(axes, plots):
        im = ax.contourf(x_coords, z_coords, field, levels=24, cmap=cmap)
        ax.set_xlabel("X (m)")
        ax.set_ylabel("Z (m)")
        ax.set_aspect("equal")
        ax.set_title(title)
        cbar = fig.colorbar(im, ax=ax)
        cbar.ax.set_ylabel("Tesla")

    surface_plot = output_dir / f"{prefix}_surface.png"
    fig.suptitle("TIAGO vs simsopt plasma response on the NCSX boundary", y=1.02)
    fig.savefig(surface_plot, dpi=200)
    plt.close(fig)

    phi_period = 2.0 * np.pi / max(nfp, 1)
    phi = np.linspace(0.0, phi_period, tiago_field.shape[0], endpoint=False)
    theta = np.linspace(0.0, 2.0 * np.pi, tiago_field.shape[1], endpoint=False)
    phi_grid, theta_grid = np.meshgrid(phi, theta, indexing="ij")

    fig2, axes2 = plt.subplots(1, 3, figsize=(18, 6), constrained_layout=True)
    for ax, (field, title, cmap) in zip(axes2, plots):
        im = ax.contourf(phi_grid, theta_grid, field, levels=24, cmap=cmap)
        ax.set_xlabel("Toroidal angle φ (rad)")
        ax.set_ylabel("Poloidal angle θ (rad)")
        ax.set_title(title)
        cbar = fig2.colorbar(im, ax=ax)
        cbar.ax.set_ylabel("Tesla")

    grid_plot = output_dir / f"{prefix}_grid.png"
    fig2.suptitle("TIAGO vs simsopt plasma response (φ, θ space)", y=1.02)
    fig2.savefig(grid_plot, dpi=200)
    plt.close(fig2)

    print(f"Saved plots:\n  • {surface_plot}\n  • {grid_plot}")


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    tiago_output = output_dir / "tiago_b_ext.csv"
    simsopt_output = output_dir / "simsopt_b_ext.csv"

    print("Running TIAGO reference test to capture B_external samples...")
    run_command([args.tiago_bin, args.gamma, args.b_total, str(tiago_output)],
                verbose=args.debug)
    tiago_field = load_tiago_field(tiago_output, args.src_nphi, args.src_ntheta)
    tiago_rms = float(np.sqrt(np.mean(tiago_field**2)))
    print(f"TIAGO B_external RMS: {tiago_rms:.6f} T "
          f"(samples → {tiago_output})")
    log_field_stats("TIAGO B_external", tiago_field, args.debug)
    if np.allclose(tiago_field, 0.0):
        raise RuntimeError("TIAGO B_external field is identically zero; "
                           "ensure test_tiago_with_simsopt_grid produced data.")

    print("Computing simsopt Virtual Casing reference...")
    gamma, simsopt_field, simsopt_rms, nfp = compute_simsopt_field(
        args.wout, args.src_nphi, args.src_ntheta)
    save_vector_csv(simsopt_output, simsopt_field)
    print(f"simsopt B_external RMS: {simsopt_rms:.6f} T "
          f"(samples → {simsopt_output})")
    log_field_stats("simsopt B_external", simsopt_field, args.debug)

    rel_error = abs(tiago_rms - simsopt_rms) / simsopt_rms * 100.0
    max_abs = float(np.max(np.abs(tiago_field - simsopt_field)))
    mean_abs = float(np.mean(np.abs(tiago_field - simsopt_field)))
    print(f"Relative RMS error: {rel_error:.3f}%")
    print(f"Max |ΔB|: {max_abs:.4e} T, Mean |ΔB|: {mean_abs:.4e} T")
    if args.debug:
        diff = tiago_field - simsopt_field
        print(f"[DEBUG] Difference stats: min={diff.min():.6e}, "
              f"max={diff.max():.6e}, rms={np.sqrt(np.mean(diff**2)):.6e}")
        print(f"[DEBUG] Sample TIAGO B_ext[0,0,:]={tiago_field[0,0,:]}")
        print(f"[DEBUG] Sample simsopt B_ext[0,0,:]={simsopt_field[0,0,:]}")

    make_plots(gamma, tiago_field, simsopt_field, nfp, output_dir)


if __name__ == "__main__":
    main()

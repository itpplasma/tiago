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
import time
import virtual_casing as vc_native

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
    parser.add_argument("--offset-distance", type=float, default=0.05,
                        help="Meters to displace points along outward normal for off-surface comparison "
                             "(<=0 disables off-surface evaluation).")
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


def load_tiago_field(bext_csv: Path, nphi: int | None = None,
                     ntheta: int | None = None
                     ) -> tuple[np.ndarray, np.ndarray | None, int, int]:
    df = pd.read_csv(bext_csv)
    if nphi is None:
        nphi = int(df["iphi"].max())
    if ntheta is None:
        ntheta = int(df["itheta"].max())
    vectors = df[["Bx", "By", "Bz"]].to_numpy()
    if vectors.shape[0] != nphi * ntheta:
        raise ValueError(f"TIAGO output size mismatch for {bext_csv}: "
                         f"expected {nphi*ntheta}, found {vectors.shape[0]}")
    coords = None
    if {"X", "Y", "Z"}.issubset(df.columns):
        coords = df[["X", "Y", "Z"]].to_numpy().reshape((nphi, ntheta, 3))
    field = vectors.reshape((nphi, ntheta, 3))
    return field, coords, nphi, ntheta


def compute_simsopt_field(wout_file: str, nphi: int, ntheta: int
                          ) -> tuple[np.ndarray, np.ndarray, float, int,
                                     VirtualCasing, np.ndarray]:
    vmec = Vmec(wout_file)
    vc = VirtualCasing.from_vmec(vmec, src_nphi=nphi, src_ntheta=ntheta,
                                 use_stellsym=True, digits=6, filename=None)
    gamma = np.asarray(vc.gamma)
    bext = np.asarray(vc.B_external)
    b_total = np.asarray(vc.B_total)
    rms = float(np.sqrt(np.mean(bext**2)))
    return (gamma.reshape(nphi, ntheta, 3),
            bext.reshape(nphi, ntheta, 3),
            rms,
            int(vmec.wout.nfp),
            vc,
            b_total.reshape(nphi, ntheta, 3))


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


def compute_normals(gamma: np.ndarray) -> np.ndarray:
    nphi, ntheta, _ = gamma.shape
    normals = np.zeros_like(gamma)

    for iphi in range(nphi):
        ip_next = (iphi + 1) % nphi
        ip_prev = (iphi - 1) % nphi
        for itheta in range(ntheta):
            it_next = (itheta + 1) % ntheta
            it_prev = (itheta - 1) % ntheta

            dphi = gamma[ip_next, itheta, :] - gamma[ip_prev, itheta, :]
            dtheta = gamma[iphi, it_next, :] - gamma[iphi, it_prev, :]
            normal = np.cross(dphi, dtheta)
            norm_mag = np.linalg.norm(normal)
            if norm_mag > 0:
                normals[iphi, itheta, :] = normal / norm_mag
            else:
                normals[iphi, itheta, :] = normal
    return normals


def flatten_xyz(field: np.ndarray) -> np.ndarray:
    nphi, ntheta, _ = field.shape
    flat = np.zeros(nphi * ntheta * 3)
    for comp in range(3):
        flat[comp * nphi * ntheta:(comp + 1) * nphi * ntheta] = \
            field[:, :, comp].reshape(-1, order="C")
    return flat


def unflatten_xyz(flat: np.ndarray, nphi: int, ntheta: int) -> np.ndarray:
    """Inverse of flatten_xyz for component-major virtual casing vectors."""
    vec = np.asarray(flat, dtype=float)
    expected = nphi * ntheta * 3
    if vec.size != expected:
        raise ValueError(f"unflatten_xyz size mismatch: expected {expected}, found {vec.size}")
    field = np.zeros((nphi, ntheta, 3))
    block = nphi * ntheta
    for comp in range(3):
        start = comp * block
        field[:, :, comp] = vec[start:start + block].reshape((nphi, ntheta), order="C")
    return field


def build_native_vc(gamma: np.ndarray, nfp: int, src_nphi: int, src_ntheta: int,
                    digits: int = 6, use_stellsym: bool = True,
                    trg_nphi: int | None = None, trg_ntheta: int | None = None):
    ctx = vc_native.VirtualCasing()
    gamma_flat = flatten_xyz(gamma)
    trg_nphi = trg_nphi or src_nphi
    trg_ntheta = trg_ntheta or src_ntheta
    ctx.setup(digits, nfp, use_stellsym,
              src_nphi, src_ntheta, gamma_flat,
              src_nphi, src_ntheta,
              trg_nphi, trg_ntheta)
    return ctx


def make_plots(gamma: np.ndarray, tiago_field: np.ndarray, simsopt_field: np.ndarray,
               nfp: int, output_dir: Path, prefix: str,
               title: str, timing_text: str | None = None) -> None:
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
    fig.suptitle(title, y=1.02)
    if timing_text:
        fig.text(0.02, 0.02, timing_text, fontsize=9, ha="left", va="bottom")
    fig.savefig(surface_plot, dpi=200, bbox_inches="tight")
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
    fig2.suptitle(f"{title} (φ, θ space)", y=1.02)
    if timing_text:
        fig2.text(0.02, 0.02, timing_text, fontsize=9, ha="left", va="bottom")
    fig2.savefig(grid_plot, dpi=200, bbox_inches="tight")
    plt.close(fig2)

    print(f"Saved plots:\n  • {surface_plot}\n  • {grid_plot}")


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    timings: dict[str, float] = {}

    tiago_output = output_dir / "tiago_b_ext.csv"
    tiago_offset_output = output_dir / "tiago_b_ext_offset.csv"
    simsopt_output = output_dir / "simsopt_b_ext.csv"

    print("Running TIAGO reference test to capture B_external samples...")
    t0 = time.perf_counter()
    tiago_cmd = [args.tiago_bin, args.gamma, args.b_total, str(tiago_output)]
    if args.offset_distance and args.offset_distance > 0.0:
        tiago_cmd.extend([str(tiago_offset_output),
                          f"{args.offset_distance:.9f}"])
    run_command(tiago_cmd, verbose=args.debug)
    timings["tiago_cmd"] = time.perf_counter() - t0
    tiago_field, tiago_coords, tiago_nphi, tiago_ntheta = load_tiago_field(
        tiago_output, args.src_nphi, args.src_ntheta)
    tiago_rms = float(np.sqrt(np.mean(tiago_field**2)))
    print(f"TIAGO B_external RMS: {tiago_rms:.6f} T "
          f"(samples → {tiago_output})")
    log_field_stats("TIAGO B_external", tiago_field, args.debug)
    if np.allclose(tiago_field, 0.0):
        raise RuntimeError("TIAGO B_external field is identically zero; "
                           "ensure test_tiago_with_simsopt_grid produced data.")

    tiago_offset_field = None
    tiago_offset_dims = None
    if args.offset_distance and args.offset_distance > 0.0 and tiago_offset_output.exists():
        tiago_offset_field, tiago_offset_coords, off_nphi, off_ntheta = load_tiago_field(
            tiago_offset_output)
        tiago_offset_dims = (off_nphi, off_ntheta)
        log_field_stats("TIAGO off-surface B_external", tiago_offset_field, args.debug)
    elif args.offset_distance and args.offset_distance > 0.0:
        print(f"WARNING: Expected off-surface CSV {tiago_offset_output} missing; "
              "skipping TIAGO off-surface comparison.")

    print("Computing simsopt Virtual Casing reference...")
    t0 = time.perf_counter()
    (gamma, simsopt_field, simsopt_rms, nfp,
     vc_obj, simsopt_b_total) = compute_simsopt_field(
        args.wout, args.src_nphi, args.src_ntheta)
    timings["simsopt_vc"] = time.perf_counter() - t0
    t0 = time.perf_counter()
    native_vc = build_native_vc(gamma, nfp, args.src_nphi, args.src_ntheta,
                                digits=6, use_stellsym=True)
    timings["native_setup"] = time.perf_counter() - t0
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

    t0 = time.perf_counter()
    surf_timing_text = (f"tiago_cmd={timings['tiago_cmd']:.2f}s | "
                        f"simsopt_vc={timings['simsopt_vc']:.2f}s | "
                        f"native_setup={timings['native_setup']:.2f}s")
    make_plots(gamma, tiago_field, simsopt_field, nfp, output_dir,
               prefix="tiago_vs_simsopt",
               title="TIAGO vs simsopt plasma response on the VMEC surface",
               timing_text=surf_timing_text)
    timings["plot_surface"] = time.perf_counter() - t0

    if args.offset_distance and args.offset_distance > 0.0 and tiago_offset_field is not None:
        normals = compute_normals(gamma)
        off_nphi, off_ntheta = tiago_offset_dims
        if tiago_offset_coords is None:
            raise RuntimeError("Offset CSV missing X/Y/Z columns; rebuild helper binary.")
        gamma_offset = tiago_offset_coords

        t0 = time.perf_counter()
        simsopt_offset_flat = native_vc.compute_external_B_offsurf(
            flatten_xyz(simsopt_b_total),
            flatten_xyz(gamma_offset))
        timings["native_offsurface"] = time.perf_counter() - t0
        simsopt_offset_field = unflatten_xyz(
            np.asarray(simsopt_offset_flat), off_nphi, off_ntheta)
        log_field_stats("simsopt off-surface B_external", simsopt_offset_field, args.debug)

        diff_offset = tiago_offset_field - simsopt_offset_field
        max_abs_off = float(np.max(np.abs(diff_offset)))
        mean_abs_off = float(np.mean(np.abs(diff_offset)))
        rms_tiago_off = float(np.sqrt(np.mean(tiago_offset_field**2)))
        rms_simsopt_off = float(np.sqrt(np.mean(simsopt_offset_field**2)))
        rel_error_off = abs(rms_tiago_off - rms_simsopt_off) / max(rms_simsopt_off, 1e-12) * 100.0
        print(f"Off-surface relative RMS error: {rel_error_off:.3f}% "
              f"(max |ΔB|={max_abs_off:.4e} T, mean |ΔB|={mean_abs_off:.4e} T)")

        t0 = time.perf_counter()
        off_timing_text = (f"tiago_cmd={timings['tiago_cmd']:.2f}s | "
                           f"native_offsurf={timings.get('native_offsurface', 0.0):.2f}s")
        make_plots(gamma_offset, tiago_offset_field, simsopt_offset_field,
                   nfp, output_dir,
                   prefix="tiago_vs_simsopt_offset",
                   title=f"Off-surface (+{args.offset_distance:.3f} m) TIAGO vs simsopt",
                   timing_text=off_timing_text)
        timings["plot_offsurface"] = time.perf_counter() - t0

    if timings:
        print("\n=== TIMINGS (s) ===")
        for key in sorted(timings):
            print(f"{key:>18s}: {timings[key]:6.2f}")


if __name__ == "__main__":
    main()

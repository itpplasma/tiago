#!/usr/bin/env python3
"""Generate simsopt reference data (gamma, B_total) on a VMEC surface."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Tuple

import numpy as np

try:
    from simsopt.mhd.vmec import Vmec
    from simsopt.mhd.virtual_casing import VirtualCasing
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "simsopt is required. Install with `pip install simsopt`."
    ) from exc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute simsopt Virtual Casing data on a VMEC surface."
    )
    parser.add_argument("--wout", required=True, help="Path to VMEC wout file.")
    parser.add_argument("--output-dir", required=True,
                        help="Directory for generated CSV files.")
    parser.add_argument("--src-nphi", type=int, default=16,
                        help="Number of toroidal grid points.")
    parser.add_argument("--src-ntheta", type=int, default=16,
                        help="Number of poloidal grid points.")
    parser.add_argument("--digits", type=int, default=6,
                        help="Floating-point precision for simsopt VirtualCasing.")
    return parser.parse_args()


def compute_virtual_casing(wout: str, src_nphi: int, src_ntheta: int,
                           digits: int) -> Tuple[np.ndarray, np.ndarray, float]:
    vmec = Vmec(wout)
    vc = VirtualCasing.from_vmec(
        vmec,
        src_nphi=src_nphi,
        src_ntheta=src_ntheta,
        use_stellsym=True,
        digits=digits,
        filename=None,
    )
    gamma = np.asarray(vc.gamma)
    b_total = np.asarray(vc.B_total)
    rms = float(np.sqrt(np.mean(b_total**2)))
    return gamma.reshape(src_nphi, src_ntheta, 3), b_total.reshape(src_nphi, src_ntheta, 3), rms


def write_csv(path: Path, header: Tuple[str, ...], data: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for row in data.reshape(-1, data.shape[-1]):
            writer.writerow([f"{val:.15e}" for val in row])


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    gamma_path = output_dir / "gamma.csv"
    btotal_path = output_dir / "b_total.csv"

    print(f"[generate_simsopt_reference] Loading VMEC: {args.wout}")
    gamma, b_total, rms = compute_virtual_casing(
        args.wout, args.src_nphi, args.src_ntheta, args.digits)

    write_csv(gamma_path, ("X", "Y", "Z"), gamma)
    write_csv(btotal_path, ("Bx", "By", "Bz"), b_total)

    print(f"[generate_simsopt_reference] Saved gamma -> {gamma_path}")
    print(f"[generate_simsopt_reference] Saved B_total -> {btotal_path}")
    print(f"[generate_simsopt_reference] B_total RMS: {rms:.6f} T")
    print(f"[generate_simsopt_reference] Grid: nphi={args.src_nphi}, "
          f"ntheta={args.src_ntheta}")


if __name__ == "__main__":
    main()

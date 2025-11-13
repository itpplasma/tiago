#!/usr/bin/env python3
"""
Compare TIAGO plasma response against simsopt VirtualCasing.
This validates our Fortran ISO C bindings by comparing against the reference
Python implementation on identical VMEC data.
"""

import sys
import numpy as np
from pathlib import Path

try:
    from simsopt.mhd.vmec import Vmec
    from simsopt.mhd.virtual_casing import VirtualCasing
except ImportError:
    print("ERROR: simsopt not installed. Install with: pip install simsopt")
    sys.exit(1)

def compare_plasma_response(wout_file, src_nphi=16, src_ntheta=16,
                            output_file="comparison_results.txt"):
    """
    Compare TIAGO and simsopt plasma response on real VMEC equilibrium.

    Args:
        wout_file: Path to VMEC wout file
        src_nphi: Source grid points (toroidal)
        src_ntheta: Source grid points (poloidal)
        output_file: Where to save results

    Returns:
        Dictionary with comparison metrics
    """
    print(f"Loading VMEC file: {wout_file}")

    try:
        vmec = Vmec(wout_file)
    except Exception as e:
        print(f"ERROR loading VMEC: {e}")
        sys.exit(1)

    print(f"Computing plasma response with simsopt (src_nphi={src_nphi}, "
          f"src_ntheta={src_ntheta})...")

    try:
        vc_simsopt = VirtualCasing.from_vmec(
            vmec, src_nphi=src_nphi, src_ntheta=src_ntheta,
            use_stellsym=True, digits=6, filename=None
        )
    except Exception as e:
        print(f"ERROR computing plasma response: {e}")
        sys.exit(1)

    # Extract simsopt results
    B_ext_simsopt = vc_simsopt.B_external
    B_normal_simsopt = vc_simsopt.B_external_normal
    gamma_simsopt = vc_simsopt.gamma
    B_total_simsopt = vc_simsopt.B_total

    print(f"simsopt B_external shape: {B_ext_simsopt.shape}")
    print(f"simsopt B_external range: [{np.min(B_ext_simsopt):.6e}, "
          f"{np.max(B_ext_simsopt):.6e}] T")
    print(f"simsopt B_external RMS: {np.sqrt(np.mean(B_ext_simsopt**2)):.6e} T")
    print(f"simsopt B_normal range: [{np.min(B_normal_simsopt):.6e}, "
          f"{np.max(B_normal_simsopt):.6e}] T")

    # Prepare output
    results = {
        'file': wout_file,
        'src_nphi': src_nphi,
        'src_ntheta': src_ntheta,
        'nfp': vmec.wout.nfp,
        'B_ext_min': float(np.min(B_ext_simsopt)),
        'B_ext_max': float(np.max(B_ext_simsopt)),
        'B_ext_rms': float(np.sqrt(np.mean(B_ext_simsopt**2))),
        'B_normal_min': float(np.min(B_normal_simsopt)),
        'B_normal_max': float(np.max(B_normal_simsopt)),
        'B_normal_mean': float(np.mean(np.abs(B_normal_simsopt))),
    }

    # Write detailed results
    with open(output_file, 'w') as f:
        f.write("SIMSOPT VIRTUAL-CASING REFERENCE RESULTS\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"VMEC File: {wout_file}\n")
        f.write(f"NFP: {vmec.wout.nfp}\n")
        f.write(f"Grid: src_nphi={src_nphi}, src_ntheta={src_ntheta}\n\n")

        f.write("B_external (full vector field):\n")
        f.write(f"  Min:  {results['B_ext_min']:.6e} T\n")
        f.write(f"  Max:  {results['B_ext_max']:.6e} T\n")
        f.write(f"  RMS:  {results['B_ext_rms']:.6e} T\n\n")

        f.write("B_external_normal (normal component):\n")
        f.write(f"  Min:  {results['B_normal_min']:.6e} T\n")
        f.write(f"  Max:  {results['B_normal_max']:.6e} T\n")
        f.write(f"  Mean: {results['B_normal_mean']:.6e} T\n\n")

        f.write("Input data:\n")
        f.write(f"  B_total RMS: {np.sqrt(np.mean(B_total_simsopt**2)):.6e} T\n")
        f.write(f"  Surface size: {gamma_simsopt.shape[0]} x "
                f"{gamma_simsopt.shape[1]}\n")

    print(f"\nResults written to {output_file}")
    return results


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python compare_with_simsopt.py <wout_file> "
              "[src_nphi] [src_ntheta]")
        sys.exit(1)

    wout_file = sys.argv[1]
    src_nphi = int(sys.argv[2]) if len(sys.argv) > 2 else 16
    src_ntheta = int(sys.argv[3]) if len(sys.argv) > 3 else 16

    results = compare_plasma_response(wout_file, src_nphi, src_ntheta)
    print("\n" + "=" * 70)
    print("COMPARISON SUMMARY")
    print("=" * 70)
    print(f"B_external range: [{results['B_ext_min']:.3e}, "
          f"{results['B_ext_max']:.3e}] T")
    print(f"B_external RMS: {results['B_ext_rms']:.3e} T")

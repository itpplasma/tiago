#!/usr/bin/env python3
"""
Debug script to dump ALL intermediate values from simsopt VirtualCasing.
This will help identify where TIAGO diverges from the reference.
"""

import sys
import numpy as np
from simsopt.mhd.vmec import Vmec
from simsopt.mhd.virtual_casing import VirtualCasing

def debug_virtual_casing(wout_file, src_nphi=16, src_ntheta=16):
    print("="*70)
    print("SIMSOPT VIRTUAL CASING DEBUG")
    print("="*70)

    vmec = Vmec(wout_file)
    print(f"\nVMEC loaded: {wout_file}")
    print(f"NFP: {vmec.wout.nfp}")

    vc = VirtualCasing.from_vmec(
        vmec, src_nphi=src_nphi, src_ntheta=src_ntheta,
        use_stellsym=True, digits=6, filename=None
    )

    print(f"\nGrid: src_nphi={src_nphi}, src_ntheta={src_ntheta}")
    print(f"use_stellsym: True")

    # Extract surface geometry (gamma)
    gamma = vc.gamma  # Shape: (nphi, ntheta, 3)
    print(f"\nSurface gamma shape: {gamma.shape}")
    print(f"gamma units: meters (SI)")
    print(f"gamma[0,0,:] = {gamma[0,0,:]}")
    print(f"gamma min: {np.min(gamma, axis=(0,1))}")
    print(f"gamma max: {np.max(gamma, axis=(0,1))}")
    print(f"gamma RMS: {np.sqrt(np.mean(gamma**2, axis=(0,1)))}")

    # Extract B_total
    B_total = vc.B_total  # Shape: (nphi, ntheta, 3)
    print(f"\nB_total shape: {B_total.shape}")
    print(f"B_total units: Tesla (SI)")
    print(f"B_total[0,0,:] = {B_total[0,0,:]}")
    print(f"B_total min: {np.min(B_total, axis=(0,1))}")
    print(f"B_total max: {np.max(B_total, axis=(0,1))}")
    print(f"B_total RMS: {np.sqrt(np.mean(B_total**2, axis=(0,1)))}")
    print(f"|B_total| range: [{np.min(np.linalg.norm(B_total, axis=2)):.6e}, "
          f"{np.max(np.linalg.norm(B_total, axis=2)):.6e}] T")

    # Extract B_external (result)
    B_ext = vc.B_external  # Shape: (nphi, ntheta, 3)
    print(f"\nB_external shape: {B_ext.shape}")
    print(f"B_external units: Tesla (SI)")
    print(f"B_external[0,0,:] = {B_ext[0,0,:]}")
    print(f"B_external min: {np.min(B_ext, axis=(0,1))}")
    print(f"B_external max: {np.max(B_ext, axis=(0,1))}")
    print(f"B_external RMS: {np.sqrt(np.mean(B_ext**2, axis=(0,1)))}")
    print(f"|B_external| range: [{np.min(np.linalg.norm(B_ext, axis=2)):.6e}, "
          f"{np.max(np.linalg.norm(B_ext, axis=2)):.6e}] T")

    # Save to files for comparison
    np.save('simsopt_gamma.npy', gamma)
    np.save('simsopt_B_total.npy', B_total)
    np.save('simsopt_B_external.npy', B_ext)

    print(f"\nSaved arrays to:")
    print(f"  simsopt_gamma.npy")
    print(f"  simsopt_B_total.npy")
    print(f"  simsopt_B_external.npy")

    # Print flattened indexing pattern
    print(f"\nFlattening pattern check (first 10 elements):")
    gamma_flat = gamma.flatten(order='C')  # C order: rightmost varies fastest
    print(f"gamma_flat[0:10] = {gamma_flat[0:10]}")
    print(f"This corresponds to:")
    for idx in range(min(10, len(gamma_flat))):
        k = idx // (src_nphi * src_ntheta)
        remainder = idx % (src_nphi * src_ntheta)
        iphi = remainder // src_ntheta
        itheta = remainder % src_ntheta
        print(f"  [{idx}] = gamma[{iphi},{itheta},{k}] = {gamma[iphi,itheta,k]:.6e}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python debug_simsopt_virtual_casing.py <wout_file> [nphi] [ntheta]")
        sys.exit(1)

    wout_file = sys.argv[1]
    nphi = int(sys.argv[2]) if len(sys.argv) > 2 else 16
    ntheta = int(sys.argv[3]) if len(sys.argv) > 3 else 16

    debug_virtual_casing(wout_file, nphi, ntheta)

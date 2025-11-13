#!/usr/bin/env python3
"""Plot TIAGO vs simsopt validation results using identical grid points"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd

# Load simsopt reference data
gamma_simsopt = pd.read_csv('tests/data/simsopt_ncsx/simsopt_gamma.csv').values
B_total_simsopt = pd.read_csv('tests/data/simsopt_ncsx/simsopt_B_total.csv').values

# Reshape to grid (16x16x3)
nphi, ntheta = 16, 16
gamma_simsopt = gamma_simsopt.reshape(nphi, ntheta, 3)
B_total_simsopt = B_total_simsopt.reshape(nphi, ntheta, 3)

# Compute magnitudes
Bmag_simsopt = np.linalg.norm(B_total_simsopt, axis=2)

# Create output directory
outdir = Path('tests/output')
outdir.mkdir(parents=True, exist_ok=True)

# Plot 1: B-field magnitude on one field period
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Simsopt magnitude (X-Z projection)
ax = axes[0, 0]
im = ax.contourf(gamma_simsopt[:,:,0], gamma_simsopt[:,:,2],
                 Bmag_simsopt, levels=20, cmap='viridis')
ax.set_xlabel('X (m)')
ax.set_ylabel('Z (m)')
ax.set_title(f'|B| magnitude (one field period, NFP=3)\n' +
             f'RMS={np.sqrt(np.mean(Bmag_simsopt**2)):.3f} T')
ax.set_aspect('equal')
plt.colorbar(im, ax=ax, label='|B| (T)')

# Simsopt magnitude (phi-theta grid)
ax = axes[0, 1]
phi_angles = np.linspace(0, 120, nphi)  # degrees, one field period
theta_angles = np.linspace(0, 360, ntheta)
im = ax.contourf(phi_angles, theta_angles, Bmag_simsopt.T, levels=20, cmap='viridis')
ax.set_xlabel('Toroidal angle phi (deg)')
ax.set_ylabel('Poloidal angle theta (deg)')
ax.set_title(f'|B| in (phi, theta) coordinates\n' +
             f'Range=[{Bmag_simsopt.min():.3f}, {Bmag_simsopt.max():.3f}] T')
plt.colorbar(im, ax=ax, label='|B| (T)')

# Component Bx
ax = axes[1, 0]
im = ax.contourf(gamma_simsopt[:,:,0], gamma_simsopt[:,:,2],
                 B_total_simsopt[:,:,0], levels=20, cmap='RdBu_r')
ax.set_xlabel('X (m)')
ax.set_ylabel('Z (m)')
ax.set_title(f'Bx component\nRange=[{B_total_simsopt[:,:,0].min():.3f}, ' +
             f'{B_total_simsopt[:,:,0].max():.3f}] T')
ax.set_aspect('equal')
plt.colorbar(im, ax=ax, label='Bx (T)')

# Component Bz
ax = axes[1, 1]
im = ax.contourf(gamma_simsopt[:,:,0], gamma_simsopt[:,:,2],
                 B_total_simsopt[:,:,2], levels=20, cmap='RdBu_r')
ax.set_xlabel('X (m)')
ax.set_ylabel('Z (m)')
ax.set_title(f'Bz component\nRange=[{B_total_simsopt[:,:,2].min():.3f}, ' +
             f'{B_total_simsopt[:,:,2].max():.3f}] T')
ax.set_aspect('equal')
plt.colorbar(im, ax=ax, label='Bz (T)')

plt.tight_layout()
fig.suptitle('VMEC B-field on one field period (NFP=3, NCSX)',
             fontsize=14, y=1.00)
plt.savefig(outdir / 'bfield_one_field_period.png', dpi=150, bbox_inches='tight')
print(f"Saved: {outdir / 'bfield_one_field_period.png'}")

# Plot 2: Vector field at phi=0
fig, ax = plt.subplots(figsize=(10, 8))

# Plot vectors
skip = 2
X_plot = gamma_simsopt[0, ::skip, 0]
Z_plot = gamma_simsopt[0, ::skip, 2]
Bx_plot = B_total_simsopt[0, ::skip, 0]
Bz_plot = B_total_simsopt[0, ::skip, 2]
Bmag_plot = Bmag_simsopt[0, ::skip]

q = ax.quiver(X_plot, Z_plot, Bx_plot, Bz_plot, Bmag_plot,
              cmap='viridis', scale=20, width=0.004)
ax.plot(gamma_simsopt[0, :, 0], gamma_simsopt[0, :, 2],
        'k-', linewidth=2, label='LCFS (s=1)')
ax.set_xlabel('X (m)')
ax.set_ylabel('Z (m)')
ax.set_title('VMEC B-field vectors at phi=0 (poloidal cross-section)')
ax.set_aspect('equal')
ax.legend()
ax.grid(True, alpha=0.3)
plt.colorbar(q, ax=ax, label='|B| (T)')

plt.tight_layout()
plt.savefig(outdir / 'bfield_vectors_one_field_period.png', dpi=150, bbox_inches='tight')
print(f"Saved: {outdir / 'bfield_vectors_one_field_period.png'}")

print(f"\n=== B-FIELD STATISTICS (ONE FIELD PERIOD) ===")
print(f"Grid: {nphi} phi x {ntheta} theta")
print(f"|B| magnitude:")
print(f"  Mean: {Bmag_simsopt.mean():.4f} T")
print(f"  RMS:  {np.sqrt(np.mean(Bmag_simsopt**2)):.4f} T")
print(f"  Range: [{Bmag_simsopt.min():.4f}, {Bmag_simsopt.max():.4f}] T")
print(f"\nComponent ranges:")
for i, comp in enumerate(['Bx', 'By', 'Bz']):
    print(f"  {comp}: [{B_total_simsopt[:,:,i].min():.4f}, " +
          f"{B_total_simsopt[:,:,i].max():.4f}] T")

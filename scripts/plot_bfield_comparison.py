#!/usr/bin/env python3
"""Create visual comparison plots of TIAGO vs simsopt B-field on grid

NOTE: These plots show B-fields computed at DIFFERENT grid points due to
different surface parametrizations (VMEC native vs SurfaceRZFourier).
The 33% RMS difference is NOT a bug - it's expected when comparing different
parametrizations. See test_tiago_with_simsopt_grid.f90 for validation using
identical grid points (achieves 0.9% agreement).
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Load data
B_simsopt = np.load('simsopt_B_total.npy')
gamma_simsopt = np.load('simsopt_gamma.npy')

# Load TIAGO data
tiago_b = np.loadtxt('tiago_b_total.dat')
tiago_x = np.loadtxt('tiago_x_surf.dat')

# Reshape TIAGO data (assuming 16x16 grid)
nphi, ntheta = 16, 16
B_tiago = np.zeros((nphi, ntheta, 3))
X_tiago = np.zeros((nphi, ntheta, 3))

for i in range(len(tiago_b)):
    iphi = int(tiago_b[i, 0]) - 1
    itheta = int(tiago_b[i, 1]) - 1
    B_tiago[iphi, itheta, :] = tiago_b[i, 2:5]
    X_tiago[iphi, itheta, :] = tiago_x[i, 2:5]

# Compute magnitudes
Bmag_tiago = np.linalg.norm(B_tiago, axis=2)
Bmag_simsopt = np.linalg.norm(B_simsopt, axis=2)

# Create output directory
outdir = Path('build/tests/output')
outdir.mkdir(parents=True, exist_ok=True)

# Plot 1: B-field magnitude comparison on surface
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# TIAGO magnitude
ax = axes[0, 0]
im = ax.contourf(X_tiago[:,:,0], X_tiago[:,:,2], Bmag_tiago, levels=20, cmap='viridis')
ax.set_xlabel('X (m)')
ax.set_ylabel('Z (m)')
ax.set_title(f'TIAGO |B| (Tesla)\nMean={Bmag_tiago.mean():.3f}, Max={Bmag_tiago.max():.3f}')
ax.set_aspect('equal')
plt.colorbar(im, ax=ax, label='|B| (T)')

# Simsopt magnitude
ax = axes[0, 1]
im = ax.contourf(gamma_simsopt[:,:,0], gamma_simsopt[:,:,2], Bmag_simsopt, levels=20, cmap='viridis')
ax.set_xlabel('X (m)')
ax.set_ylabel('Z (m)')
ax.set_title(f'Simsopt |B| (Tesla)\nMean={Bmag_simsopt.mean():.3f}, Max={Bmag_simsopt.max():.3f}')
ax.set_aspect('equal')
plt.colorbar(im, ax=ax, label='|B| (T)')

# Ratio: TIAGO / simsopt
ax = axes[1, 0]
ratio = Bmag_tiago / Bmag_simsopt
im = ax.contourf(X_tiago[:,:,0], X_tiago[:,:,2], ratio, levels=20, cmap='RdBu_r', vmin=0.9, vmax=1.1)
ax.set_xlabel('X (m)')
ax.set_ylabel('Z (m)')
ax.set_title(f'Ratio: TIAGO/Simsopt\nMean={ratio.mean():.4f}, Range=[{ratio.min():.4f}, {ratio.max():.4f}]')
ax.set_aspect('equal')
plt.colorbar(im, ax=ax, label='Ratio')

# Difference: TIAGO - simsopt
ax = axes[1, 1]
diff = Bmag_tiago - Bmag_simsopt
im = ax.contourf(X_tiago[:,:,0], X_tiago[:,:,2], diff, levels=20, cmap='RdBu_r')
ax.set_xlabel('X (m)')
ax.set_ylabel('Z (m)')
ax.set_title(f'Difference: TIAGO - Simsopt (T)\nMean={diff.mean():.4f}, RMS={np.sqrt((diff**2).mean()):.4f}')
ax.set_aspect('equal')
plt.colorbar(im, ax=ax, label='Δ|B| (T)')

plt.tight_layout()
fig.text(0.5, 0.02,
         'NOTE: Grids at different phi angles due to parametrization mismatch.\n' +
         'See test_tiago_with_simsopt_grid for validation with identical grids (0.9% agreement).',
         ha='center', fontsize=9, style='italic', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
plt.savefig(outdir / 'bfield_magnitude_comparison.png', dpi=150, bbox_inches='tight')
print(f"Saved: {outdir / 'bfield_magnitude_comparison.png'}")

# Plot 2: Component-wise comparison
fig, axes = plt.subplots(3, 3, figsize=(16, 14))

components = ['Bx', 'By', 'Bz']
for i, comp_name in enumerate(components):
    # TIAGO component
    ax = axes[i, 0]
    im = ax.contourf(X_tiago[:,:,0], X_tiago[:,:,2], B_tiago[:,:,i], levels=20, cmap='RdBu_r')
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Z (m)')
    ax.set_title(f'TIAGO {comp_name} (T)\nRange=[{B_tiago[:,:,i].min():.3f}, {B_tiago[:,:,i].max():.3f}]')
    ax.set_aspect('equal')
    plt.colorbar(im, ax=ax, label=f'{comp_name} (T)')

    # Simsopt component
    ax = axes[i, 1]
    im = ax.contourf(gamma_simsopt[:,:,0], gamma_simsopt[:,:,2], B_simsopt[:,:,i], levels=20, cmap='RdBu_r')
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Z (m)')
    ax.set_title(f'Simsopt {comp_name} (T)\nRange=[{B_simsopt[:,:,i].min():.3f}, {B_simsopt[:,:,i].max():.3f}]')
    ax.set_aspect('equal')
    plt.colorbar(im, ax=ax, label=f'{comp_name} (T)')

    # Difference
    ax = axes[i, 2]
    diff_comp = B_tiago[:,:,i] - B_simsopt[:,:,i]
    im = ax.contourf(X_tiago[:,:,0], X_tiago[:,:,2], diff_comp, levels=20, cmap='RdBu_r')
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Z (m)')
    ax.set_title(f'Δ{comp_name} = TIAGO - Simsopt\nRMS={np.sqrt((diff_comp**2).mean()):.4f} T')
    ax.set_aspect('equal')
    plt.colorbar(im, ax=ax, label=f'Δ{comp_name} (T)')

plt.tight_layout()
plt.savefig(outdir / 'bfield_components_comparison.png', dpi=150, bbox_inches='tight')
print(f"Saved: {outdir / 'bfield_components_comparison.png'}")

# Plot 3: Vector field visualization on poloidal cross-section (phi=0)
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# TIAGO vectors at phi=0
ax = axes[0]
skip = 2
X_plot = X_tiago[0, ::skip, 0]
Z_plot = X_tiago[0, ::skip, 2]
Bx_plot = B_tiago[0, ::skip, 0]
Bz_plot = B_tiago[0, ::skip, 2]
ax.quiver(X_plot, Z_plot, Bx_plot, Bz_plot, Bmag_tiago[0, ::skip],
          cmap='viridis', scale=20, width=0.004)
ax.plot(X_tiago[0, :, 0], X_tiago[0, :, 2], 'k-', linewidth=2, label='Surface')
ax.set_xlabel('X (m)')
ax.set_ylabel('Z (m)')
ax.set_title('TIAGO B-field at phi=0')
ax.set_aspect('equal')
ax.legend()
ax.grid(True, alpha=0.3)

# Simsopt vectors at phi=0
ax = axes[1]
X_plot = gamma_simsopt[0, ::skip, 0]
Z_plot = gamma_simsopt[0, ::skip, 2]
Bx_plot = B_simsopt[0, ::skip, 0]
Bz_plot = B_simsopt[0, ::skip, 2]
ax.quiver(X_plot, Z_plot, Bx_plot, Bz_plot, Bmag_simsopt[0, ::skip],
          cmap='viridis', scale=20, width=0.004)
ax.plot(gamma_simsopt[0, :, 0], gamma_simsopt[0, :, 2], 'k-', linewidth=2, label='Surface')
ax.set_xlabel('X (m)')
ax.set_ylabel('Z (m)')
ax.set_title('Simsopt B-field at phi=0')
ax.set_aspect('equal')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(outdir / 'bfield_vectors_phi0.png', dpi=150, bbox_inches='tight')
print(f"Saved: {outdir / 'bfield_vectors_phi0.png'}")

# Print statistics
print("\n=== COMPARISON STATISTICS ===")
print(f"Grid size: {nphi} x {ntheta}")
print(f"\n|B| magnitude:")
print(f"  TIAGO:   mean={Bmag_tiago.mean():.4f} T, range=[{Bmag_tiago.min():.4f}, {Bmag_tiago.max():.4f}]")
print(f"  Simsopt: mean={Bmag_simsopt.mean():.4f} T, range=[{Bmag_simsopt.min():.4f}, {Bmag_simsopt.max():.4f}]")
print(f"  Ratio:   mean={ratio.mean():.4f}, range=[{ratio.min():.4f}, {ratio.max():.4f}]")
print(f"  Difference RMS: {np.sqrt((diff**2).mean()):.4f} T ({100*np.sqrt((diff**2).mean())/Bmag_simsopt.mean():.1f}%)")

print(f"\nComponent RMS differences:")
for i, comp in enumerate(['Bx', 'By', 'Bz']):
    diff_comp = B_tiago[:,:,i] - B_simsopt[:,:,i]
    rms = np.sqrt((diff_comp**2).mean())
    print(f"  {comp}: {rms:.4f} T")

print(f"\nPlots saved to: {outdir}/")

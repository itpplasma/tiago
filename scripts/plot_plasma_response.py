#!/usr/bin/env python3
"""
Visualize plasma response B_external field from Tiago computation
Generates publication-quality plots comparing B_field components
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import sys
from pathlib import Path

def load_bext_data(filename):
    """Load B_external data from Fortran output file"""
    try:
        data = np.loadtxt(filename, skiprows=3)
        iphi = data[:, 0].astype(int)
        itheta = data[:, 1].astype(int)
        bx = data[:, 2]
        by = data[:, 3]
        bz = data[:, 4]
        bmag = data[:, 5]

        # Get grid dimensions
        nphi = iphi.max()
        ntheta = itheta.max()

        # Reshape to grid
        bx_grid = bx.reshape(nphi, ntheta)
        by_grid = by.reshape(nphi, ntheta)
        bz_grid = bz.reshape(nphi, ntheta)
        bmag_grid = bmag.reshape(nphi, ntheta)

        return {
            'nphi': nphi, 'ntheta': ntheta,
            'bx': bx_grid, 'by': by_grid, 'bz': bz_grid, 'bmag': bmag_grid
        }
    except Exception as e:
        print(f"Error loading {filename}: {e}", file=sys.stderr)
        return None

def create_visualization(data, output_file):
    """Create comprehensive plasma response visualization"""
    if data is None:
        return False

    fig = plt.figure(figsize=(16, 12))
    gs = GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.3)

    nphi = data['nphi']
    ntheta = data['ntheta']

    # Create coordinate grids
    phi_grid = np.linspace(0, 2*np.pi, nphi+1)
    theta_grid = np.linspace(0, 2*np.pi, ntheta+1)

    # Plot 1: Magnitude of B_external
    ax1 = fig.add_subplot(gs[0, 0])
    im1 = ax1.contourf(phi_grid[:-1], theta_grid[:-1], data['bmag'].T,
                       levels=20, cmap='viridis')
    ax1.set_xlabel('Toroidal angle φ (rad)', fontsize=11)
    ax1.set_ylabel('Poloidal angle θ (rad)', fontsize=11)
    ax1.set_title('|B_external| - Magnitude of Plasma Response Field', fontsize=12, fontweight='bold')
    cbar1 = plt.colorbar(im1, ax=ax1)
    cbar1.set_label('|B| (Tesla)', fontsize=10)

    # Plot 2: B_x component
    ax2 = fig.add_subplot(gs[0, 1])
    im2 = ax2.contourf(phi_grid[:-1], theta_grid[:-1], data['bx'].T,
                       levels=20, cmap='RdBu_r')
    ax2.set_xlabel('Toroidal angle φ (rad)', fontsize=11)
    ax2.set_ylabel('Poloidal angle θ (rad)', fontsize=11)
    ax2.set_title('B_x - Radial Component', fontsize=12, fontweight='bold')
    cbar2 = plt.colorbar(im2, ax=ax2)
    cbar2.set_label('B_x (Tesla)', fontsize=10)

    # Plot 3: B_z component
    ax3 = fig.add_subplot(gs[1, 0])
    im3 = ax3.contourf(phi_grid[:-1], theta_grid[:-1], data['bz'].T,
                       levels=20, cmap='RdBu_r')
    ax3.set_xlabel('Toroidal angle φ (rad)', fontsize=11)
    ax3.set_ylabel('Poloidal angle θ (rad)', fontsize=11)
    ax3.set_title('B_z - Vertical Component', fontsize=12, fontweight='bold')
    cbar3 = plt.colorbar(im3, ax=ax3)
    cbar3.set_label('B_z (Tesla)', fontsize=10)

    # Plot 4: Statistics and info
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.axis('off')

    # Compute statistics
    bmag_mean = np.mean(data['bmag'])
    bmag_max = np.max(data['bmag'])
    bmag_min = np.min(data['bmag'])
    bx_rms = np.sqrt(np.mean(data['bx']**2))
    by_rms = np.sqrt(np.mean(data['by']**2))
    bz_rms = np.sqrt(np.mean(data['bz']**2))

    stats_text = f"""PLASMA RESPONSE FIELD STATISTICS
{'='*45}

Grid Resolution:
  • Toroidal (φ): {nphi} points
  • Poloidal (θ): {ntheta} points
  • Total surface points: {nphi*ntheta}

Magnitude |B_external|:
  • Mean: {bmag_mean:.4e} T
  • Max: {bmag_max:.4e} T
  • Min: {bmag_min:.4e} T

RMS Component Values:
  • B_x (radial): {bx_rms:.4e} T
  • B_y (toroidal): {by_rms:.4e} T
  • B_z (vertical): {bz_rms:.4e} T

Algorithm:
  • Virtual Casing Principle
  • C++ library: HiddenSymmetries/virtual-casing
  • Method: FFT-based surface integral
  • Geometry: NCSX stellarator boundary
"""

    ax4.text(0.05, 0.95, stats_text, transform=ax4.transAxes,
             fontsize=10, verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    fig.suptitle('TIAGO Plasma Response Field - B_external from Virtual Casing',
                 fontsize=14, fontweight='bold', y=0.995)

    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"✓ Saved visualization to {output_file}")
    return True

def create_field_comparison_plot(data, output_file):
    """Create field component comparison plot"""
    if data is None:
        return False

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    phi_grid = np.linspace(0, 2*np.pi, data['nphi']+1)
    theta_grid = np.linspace(0, 2*np.pi, data['ntheta']+1)

    components = [
        (data['bx'], 'B_x (Radial)', 'RdBu_r', axes[0]),
        (data['by'], 'B_y (Toroidal)', 'RdBu_r', axes[1]),
        (data['bz'], 'B_z (Vertical)', 'RdBu_r', axes[2])
    ]

    for component, title, cmap, ax in components:
        im = ax.contourf(phi_grid[:-1], theta_grid[:-1], component.T,
                        levels=20, cmap=cmap)
        ax.set_xlabel('Toroidal angle φ (rad)', fontsize=11)
        ax.set_ylabel('Poloidal angle θ (rad)', fontsize=11)
        ax.set_title(title, fontsize=12, fontweight='bold')
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('B (Tesla)', fontsize=10)

    fig.suptitle('TIAGO Plasma Response - Field Components',
                fontsize=14, fontweight='bold')
    plt.tight_layout()

    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"✓ Saved field comparison to {output_file}")
    return True

if __name__ == '__main__':
    # Find output files from tests
    build_dir = Path('/home/ert/code/tiago/build')
    output_dir = build_dir / 'tests' / 'output'
    plot_dir = Path('/home/ert/code/tiago/plots')
    plot_dir.mkdir(exist_ok=True)

    # Look for bext output files
    bext_file = Path('bext_ncsx.txt')

    if bext_file.exists():
        print(f"Loading B_external data from {bext_file}")
        data = load_bext_data(str(bext_file))

        if data:
            # Create visualizations
            create_visualization(data, str(plot_dir / 'plasma_response_ncsx.png'))
            create_field_comparison_plot(data, str(plot_dir / 'plasma_response_components.png'))

            print(f"\n{'='*60}")
            print("PLASMA RESPONSE VISUALIZATION COMPLETE")
            print(f"{'='*60}")
            print(f"Output plots:")
            print(f"  • {plot_dir / 'plasma_response_ncsx.png'}")
            print(f"  • {plot_dir / 'plasma_response_components.png'}")
            print(f"{'='*60}")
        else:
            print("Failed to load data", file=sys.stderr)
            sys.exit(1)
    else:
        print(f"ERROR: {bext_file} not found", file=sys.stderr)
        sys.exit(1)

#!/usr/bin/env python3
"""Debug script to show how simsopt extracts B-field from VMEC"""

import numpy as np
from simsopt.mhd import Vmec

# Load VMEC
wout_file = "tests/cases/ncsx_nfp3/wout_ncsx.nc"
vmec = Vmec(wout_file)

# Get first grid point
theta = 0.0
phi = 0.0

# Get the surface directly from simsopt
s_full = 1.0
r, z = vmec.vmec.cfunct(theta, phi, s_full, 0)
drds, dzds = vmec.vmec.cfunct(theta, phi, s_full, 1)
drdt, dzdt = vmec.vmec.cfunct(theta, phi, s_full, 2)
drdp, dzdp = vmec.vmec.cfunct(theta, phi, s_full, 3)

print("=== VMEC GEOMETRY at (theta=0, phi=0, s=1) ===")
print(f"r = {r:.10f} cm")
print(f"z = {z:.10f} cm")
print(f"drds = {drds:.10e}")
print(f"dzds = {dzds:.10e}")
print(f"drdt = {drdt:.10e}")
print(f"dzdt = {dzdt:.10e}")
print(f"drdp = {drdp:.10e}")
print(f"dzdp = {dzdp:.10e}")
print()

# Get lambda values
lambd, dlamds, dlamdt, dlamdp = vmec.vmec.lfunct(theta, phi, s_full)
print(f"lambda = {lambd:.10e}")
print(f"dlambda/ds = {dlamds:.10e}")
print(f"dlambda/dtheta = {dlamdt:.10e}")
print(f"dlambda/dphi = {dlamdp:.10e}")
print()

# Get iota
iota = vmec.vmec.iota_func(s_full)
print(f"iota = {iota:.10e}")
print()

# Compute basis vectors (in cm units as VMEC stores them)
cos_phi = np.cos(phi)
sin_phi = np.sin(phi)

e_s = np.array([drds * cos_phi, drds * sin_phi, dzds])
e_theta = np.array([drdt * cos_phi, drdt * sin_phi, dzdt])
e_phi = np.array([(drdp * cos_phi - r * sin_phi),
                   (drdp * sin_phi + r * cos_phi),
                   dzdp])

print("=== BASIS VECTORS (CGS, cm) ===")
print(f"e_s = {e_s}")
print(f"|e_s| = {np.linalg.norm(e_s):.10f}")
print(f"e_theta = {e_theta}")
print(f"|e_theta| = {np.linalg.norm(e_theta):.10f}")
print(f"e_phi = {e_phi}")
print(f"|e_phi| = {np.linalg.norm(e_phi):.10f}")
print()

# Compute contravariant jacobian
cjac = 1.0 / (1.0 + dlamdt)
print(f"cjac = 1/(1 + dlambda/dtheta) = {cjac:.10e}")
print()

# Compute VMEC contravariant basis vectors
e_vartheta = cjac * (e_theta - dlamds * e_s - dlamdp * e_phi)
e_varphi = e_phi - dlamdp * cjac * e_theta

print("=== CONTRAVARIANT BASIS VECTORS (CGS, cm) ===")
print(f"e_vartheta = {e_vartheta}")
print(f"|e_vartheta| = {np.linalg.norm(e_vartheta):.10f}")
print(f"e_varphi = {e_varphi}")
print(f"|e_varphi| = {np.linalg.norm(e_varphi):.10f}")
print()

# Get contravariant B-field components from VMEC
bsupu, bsupv = vmec.vmec.bfunct(theta, phi, s_full)
print(f"=== CONTRAVARIANT B COMPONENTS (CGS, Gauss) ===")
print(f"B^vartheta (bsupu) = {bsupu:.10e} Gauss")
print(f"B^varphi (bsupv) = {bsupv:.10e} Gauss")
print()

# Compute Cartesian B-field
b_vec_cgs = bsupu * e_vartheta + bsupv * e_varphi
print(f"=== B-FIELD IN CARTESIAN (CGS) ===")
print(f"B_vec (Gauss·cm basis) = {b_vec_cgs}")
print(f"|B| (Gauss) = {np.linalg.norm(b_vec_cgs):.10f} Gauss")
print()

# Convert to SI
cm_to_m = 1e-2
gauss_to_tesla = 1e-4
b_vec_si = b_vec_cgs * gauss_to_tesla
print(f"=== B-FIELD IN SI ===")
print(f"B_vec (Tesla) = {b_vec_si}")
print(f"|B| (Tesla) = {np.linalg.norm(b_vec_si):.10f} Tesla")
print()

# Compare with simsopt surface extraction
from simsopt.geo import SurfaceRZFourier
surf = SurfaceRZFourier.from_vmec_input(wout_file, range="full torus", nphi=16, ntheta=16)
gamma = surf.gamma().reshape((16,16,3))
B_total = np.zeros((16, 16, 3))
vmec.boundary = surf
for i in range(16):
    for j in range(16):
        B_total[i,j,:] = vmec.B_external_normal(gamma[i,j,:])

print(f"=== SIMSOPT DIRECT EXTRACTION ===")
print(f"gamma[0,0,:] = {gamma[0,0,:]} meters")
print(f"B_total[0,0,:] = {B_total[0,0,:]} Tesla")
print(f"|B_total[0,0]| = {np.linalg.norm(B_total[0,0,:]):.10f} Tesla")

# Plasma Response Integration with VMEC Equilibria

## Overview

Tiago now supports plasma response calculations for VMEC equilibria using the virtual casing principle. This document explains how to integrate VMEC data with the plasma response module to compute the magnetic field contribution from plasma currents.

## Architecture

### Data Flow

```
VMEC wout file
    ↓
[libneo VMEC reader] → Fourier series coefficients (R_mn, Z_mn, λ_mn)
    ↓
[magfie_vmec] → Evaluate B-field at arbitrary points
    ↓
[plasma_response_t] → Extract B_external (plasma contribution)
```

### Key libneo Components

1. **vmecinm_sub**: Reads netCDF wout files
   - Extracts Fourier coefficients from VMEC equilibrium
   - Scales to SI units

2. **spline_vmec_data**: Provides field evaluation
   - `vmec_field()`: Evaluates B-field components at given (s, θ, φ) coordinates
   - Handles Fourier series interpolation

3. **magfie_vmec**: High-level interface
   - Wraps vmec_field for specific calculations
   - Computes B-module, derivatives, metrics

## Implementation Strategy

### Step 1: Load VMEC Equilibrium

```fortran
use new_vmec_stuff_mod, only: netcdffile, nper
use vmecin_sub, only: vmecin

! Set VMEC file path
netcdffile = 'path/to/wout.nc'

! Read Fourier coefficients from wout file
call vmecin(rmnc, zmns, almns, rmns, zmnc, almnc, aiota, &
            phi, sps, axm, axn, s, nsurfm, nstrm, kpar, flux)
```

### Step 2: Evaluate VMEC Surface Geometry

```fortran
! On the VMEC surface (s = 1.0), evaluate position
do iphi = 1, nphi
    varphi = 2.0_dp * pi * (iphi - 1) / nphi
    do itheta = 1, ntheta
        theta = 2.0_dp * pi * (itheta - 1) / ntheta

        ! Call VMEC field evaluator
        call vmec_field(s=1.0_dp, theta, varphi, &
                       A_theta, A_phi, dA_theta_ds, dA_phi_ds, aiota,&
                       sqrtg, lambda, dl_ds, dl_dt, dl_dp,&
                       B_theta, B_phi, B_r, ...)

        ! Extract R, Z coordinates from geometry
        x_surf(iphi, itheta, 1:3) = [R, phi_cart, Z]
    end do
end do
```

### Step 3: Compute Total B-Field on Surface

The total B-field on the VMEC surface consists of:
- B_vacuum (from external coils via Biot-Savart)
- B_plasma (from plasma currents, which we want to extract)

```fortran
! Evaluate B-field at surface points
do iphi = 1, nphi
    do itheta = 1, ntheta
        ! Get B from VMEC (this is B_total)
        call vmec_field(s=1.0_dp, theta(itheta), phi(iphi), ...)

        ! b_total(:, iphi, itheta) contains components in Cartesian or field-aligned coords
    end do
end do
```

### Step 4: Extract Plasma Response

```fortran
use tiago_plasma_response, only: plasma_response_t

type(plasma_response_t) :: plasma_response

! Initialize virtual casing with VMEC surface and total B-field
call plasma_response%init(nfp=nper, x_surf=x_surf, b_total=b_total,&
                          src_nphi=nphi, src_ntheta=ntheta)

! Compute plasma response
call plasma_response%compute_bext(b_total, b_external)

! b_external now contains B_plasma (the plasma contribution)
! B_vacuum = B_total - B_plasma
```

## Integration with Existing Tiago Code

### For Vacuum-Plus-Plasma Diagnostics

The envisioned workflow:

```fortran
! 1. Load coils and compute B_vacuum
type(vacuum_solver_t) :: vac_solver
call vac_solver%init(coil_file)
call vac_solver%flux_loops(loops, flux_vac)

! 2. Load VMEC equilibrium
! ... (call vmecin and surface evaluation as above)

! 3. Compute total B-field on VMEC surface from VMEC
! ... (call vmec_field repeatedly)

! 4. Extract plasma response
type(plasma_response_t) :: plasma_response
call plasma_response%init(nfp, x_surf, b_total, ...)
call plasma_response%compute_bext(b_total, b_external)

! 5. Combine for diagnostics
! B_combined = B_vacuum + B_external
! Compute diagnostic signals from combined field
```

## Dependencies

### Required Components
- **libneo**: For VMEC field evaluation
- **virtual-casing library**: For plasma response computation
- **FFTW, BLAS, LAPACK**: For linear algebra

### Optional
- **NetCDF**: Already required by libneo for reading wout files

## Validation Against DIAGNO

DIAGNO (STELLOPT's diagnostic code) also computes plasma response via virtual casing. To validate:

1. Run DIAGNO on same VMEC+coil configuration
2. Extract VMEC surface geometry in Tiago
3. Compute B_external using plasma_response_t
4. Compare plasma response components
5. Generate diagnostic signals and compare

## Performance Considerations

- **Virtual casing accuracy**: Control with `digits` parameter (default 6)
- **Surface resolution**: Balance accuracy vs. memory/computation
- **Target grid resolution**: Can be different from source (trg_nphi, trg_ntheta)
- **Stellarator symmetry**: Use `use_stellsym=true` for symmetric geometries

## Known Limitations

1. Virtual casing requires:
   - Complete surface field data (cannot interpolate missing points)
   - Field defined on full toroidal extent

2. VMEC interface:
   - Currently evaluates field at points in VMEC coordinates
   - Fourier series evaluation limited to band-limited resolution

3. Plasma response module:
   - Assumes surface entirely encloses plasma
   - Works best with closed geometries (stellarators, tokamaks)

## Next Steps

Future enhancements:
1. Integrate into `plasma_forward_t` solver for diagnostic calculations
2. Add VMEC field to libneo's field_t interface
3. Create Python bindings for VMEC+plasma response workflow
4. Optimize array layouts for cache locality
5. Add MPI support for large grids

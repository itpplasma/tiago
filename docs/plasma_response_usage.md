# Plasma Response via Virtual Casing

## Overview

Tiago now supports plasma response calculation using the virtual casing principle via the [HiddenSymmetries/virtual-casing](https://github.com/hiddenSymmetries/virtual-casing) library. This allows computation of the magnetic field contribution from plasma currents outside the plasma boundary.

## Building with Plasma Response Support

By default, Tiago builds with plasma response enabled (`TIAGO_ENABLE_PLASMA=ON`). To disable it:

```bash
cmake -S . -B build -DTIAGO_ENABLE_PLASMA=OFF
cmake --build build -j
```

## Using the Plasma Response Module

The plasma response functionality is encapsulated in the `tiago_plasma_response` module, which provides the `plasma_response_t` derived type.

### Basic Usage

```fortran
use tiago_plasma_response, only: plasma_response_t

type(plasma_response_t) :: plasma_response
real(dp), allocatable :: x_surf(:,:,:)     ! Surface geometry (nphi, ntheta, 3)
real(dp), allocatable :: b_total(:,:,:)    ! Total B-field (nphi, ntheta, 3)
real(dp), allocatable :: b_external(:,:,:) ! Plasma contribution (nphi, ntheta, 3)

! Initialize from VMEC surface and total B-field
! nfp: number of field periods
! x_surf: (nphi, ntheta, 3) surface coordinates [m]
! b_total: (nphi, ntheta, 3) total magnetic field [T]
! src_nphi, src_ntheta: resolution of input B-field grid
call plasma_response%init(nfp=3, x_surf=x_surf, b_total=b_total, &
                          src_nphi=size(b_total,1), src_ntheta=size(b_total,2))

! Compute plasma response (B_external)
call plasma_response%compute_bext(b_total, b_external)

! Check if initialized successfully
if (plasma_response%is_initialized()) then
    print *, 'Plasma response computed successfully'
end if

! Clean up resources
call plasma_response%finalize()
```

### Optional Parameters

```fortran
call plasma_response%init(nfp, x_surf, b_total, src_nphi, src_ntheta, &
                          trg_nphi=64,           ! Target resolution (default: same as source)
                          trg_ntheta=64,         ! Target resolution (default: same as source)
                          use_stellsym=.true.,   ! Use stellarator symmetry (default: true)
                          digits=6)              ! Accuracy digits (default: 6)
```

## Implementation Details

### C Bindings

The virtual-casing C++ library exposes four key functions via C bindings:

- `VirtualCasingCreateContextD()` - Create opaque context handle
- `VirtualCasingSetupD()` - Initialize from surface geometry and B-field
- `VirtualCasingComputeBextD()` - Compute plasma contribution B_external
- `VirtualCasingDestroyContextD()` - Free context resources

The Fortran module wraps these in an object-oriented interface.

### Array Layout

Fortran arrays are flattened to C ordering for interoperability:

- Surface coordinates: `{x11, x12, ..., x1Np, x21, ..., xNtNp, y11, ..., z11, ...}`
- B-field arrays: Same flattening order

The module handles this conversion transparently.

### Virtual Casing Principle

The virtual casing principle solves for B_external (plasma contribution) as:

```
B_external = B/2 + gradG[B·n] + BiotSavart[n×B]
```

where:
- B is the total magnetic field
- n is the surface normal
- G is the Green's function for the domain outside the plasma

## Integration with Diagnostics

To compute diagnostic signals including plasma response:

1. Get VMEC surface geometry and total B-field from libneo
2. Use `plasma_response_t` to compute B_external
3. Pass (B_coil + B_external) to diagnostic solvers for flux loop and Rogowski signals

## Dependencies

Virtual casing requires:
- FFTW3 (double precision)
- BLAS/LAPACK
- OpenMP

These are detected automatically during CMake configuration.

## References

- Virtual Casing Library: https://github.com/hiddenSymmetries/virtual-casing
- Theory: Landreman et al., Plasma Physics and Controlled Fusion (2017)

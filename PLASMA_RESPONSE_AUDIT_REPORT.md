# PLASMA RESPONSE AUDIT REPORT - TIAGO Virtual-Casing Integration
**Date**: 2025-11-13
**Auditor**: Sergei (Production-Grade Code Analysis)
**Status**: **CRITICAL ISSUES IDENTIFIED - 500-1200% Plasma Error vs DIAGNO**

---

## EXECUTIVE SUMMARY

The TIAGO plasma response implementation using HiddenSymmetries `virtual-casing` C++ library shows vacuum-only agreement with DIAGNO of 0.18% RMS error, but adding plasma response degrades accuracy to 1.45% RMS error. The plasma-only contribution shows **500-1200% relative error** compared to DIAGNO reference.

**Root Cause Hypothesis**: Array index ordering mismatch between Fortran column-major and C row-major memory layouts, compounded by potential grid convention differences (full period vs half period, grid shift handling).

---

## 1. CRITICAL FINDING: ARRAY FLATTENING ORDER MISMATCH

### 1.1 TIAGO Implementation (`plasma_response.f90:110-119`)

```fortran
! Flatten surface coordinates to C order: {x11, x12, ..., y11, y12, ..., z11, ...}
allocate(x_flat(nphi * ntheta * 3))
do k = 1, 3
    do j = 1, ntheta
        do i = 1, nphi
            idx = (k-1) * nphi * ntheta + (j-1) * nphi + i
            x_flat(idx) = real(x_surf(i, j, k), c_double)
        end do
    end do
end do
```

**Index formula**: `idx = (k-1)*nphi*ntheta + (j-1)*nphi + i`
**Memory order produced**: `[x(1,1), x(2,1), ..., x(nphi,1), x(1,2), x(2,2), ..., x(nphi,ntheta), y(1,1), ...]`

This flattens **phi-fastest**, theta-second, xyz-slowest.

---

### 1.2 C API Documentation (`virtual-casing.h:38-39`)

```c
/**
 * @param[in] X the surface coordinates in the order {x11, x12, ..., x1Np,
 * x21, x22, ... , xNtNp, y11, ... , z11, ...}.
```

**Documentation states**:
- `Nt` = toroidal discretization order (number of **toroidal** points)
- `Np` = poloidal discretization order (number of **poloidal** points)
- Order: `{x11, x12, ..., x1Np, x21, x22, ..., xNtNp, y11, ...}`

**Interpretation**:
- `x1Np` = x at toroidal index 1, poloidal index Np
- `xNtNp` = x at toroidal index Nt, poloidal index Np

This means **poloidal-fastest** (Np varies fastest), then toroidal (Nt), then xyz components.

---

### 1.3 Simsopt Reference (`virtual_casing.py:209-214`)

```python
# virtual_casing wants all input arrays to be 1D. The order is
# {x11, x12, ..., x1Np, x21, x22, ... , xNtNp, y11, ... , z11, ...}
# where Nt is toroidal (not theta!) and Np is poloidal (not phi!)
gamma1d = np.zeros(src_nphi * src_ntheta * 3)
for jxyz in range(3):
    gamma1d[jxyz * src_nphi * src_ntheta: (jxyz + 1) * src_nphi * src_ntheta] = \
        gamma[:, :, jxyz].flatten(order='C')
```

**Simsopt approach**: Uses `flatten(order='C')` on `gamma[nphi, ntheta, 3]`.
NumPy C-order flattening on shape `(nphi, ntheta, 3)` produces: `gamma[0,0,0], gamma[0,1,0], ..., gamma[0,ntheta-1,0], gamma[1,0,0], ...`

This is **theta-fastest** (ntheta varies fastest), then phi, then xyz.

---

### 1.4 **CRITICAL DISCREPANCY**

| Implementation | Index Order | Memory Layout |
|----------------|-------------|---------------|
| **TIAGO** | `idx = (k-1)*nphi*ntheta + (j-1)*nphi + i` | **phi-fastest**, theta-second, xyz-slowest |
| **Simsopt** | `gamma[:,:,jxyz].flatten('C')` on `(nphi,ntheta,3)` | **theta-fastest**, phi-second, xyz-slowest |
| **C API docs** | `{x11, x12, ..., x1Np}` (Np=poloidal) | **poloidal-fastest**, toroidal-second, xyz-slowest |

**CONCLUSION**: TIAGO flattens with **phi-fastest**, but C API expects **poloidal(theta)-fastest**. This causes massive spatial coordinate scrambling.

---

## 2. TERMINOLOGY CONFUSION: Nt vs Ntheta

### 2.1 C API Parameter Names

```c
void VirtualCasingSetupD(int digits, int NFP, bool half_period,
                         long Nt, long Np, const double* X,
                         long src_Nt, long src_Np,
                         long trg_Nt, long trg_Np, void* ctx)
```

**Documentation (`virtual-casing.h:34-37`)**:
```
* @param[in] Nt surface discretization order in toroidal direction (in one field period).
* @param[in] Np surface discretization order in poloidal direction.
```

**Interpretation**: `Nt` = **toroidal** discretization, `Np` = **poloidal** discretization.

---

### 2.2 TIAGO Usage (`plasma_response.f90:133-137`)

```fortran
call VirtualCasingSetupD(int(digits_val, c_int), int(nfp, c_int), stellsym_c, &
                        int(nphi, c_long), int(ntheta, c_long), x_flat, &
                        int(src_nphi, c_long), int(src_ntheta, c_long), &
                        int(trg_nphi_val, c_long), int(trg_ntheta_val, c_long), &
                        self%ctx)
```

**TIAGO passes**:
- `Nt` ← `nphi` (number of **toroidal** points)
- `Np` ← `ntheta` (number of **poloidal** points)

This looks **correct semantically** (phi=toroidal, theta=poloidal).

---

### 2.3 Simsopt Usage (`virtual_casing.py:229-233`)

```python
vcasing.setup(
    digits, nfp, stellsym,
    src_nphi, src_ntheta, gamma1d,
    src_nphi, src_ntheta,
    trgt_nphi, trgt_ntheta)
```

Simsopt also passes `nphi` → `Nt`, `ntheta` → `Np`.

**HOWEVER**, Simsopt's `gamma` array is shaped `(src_nphi, src_ntheta, 3)` with **nphi first**, and then flattened with C-order (rightmost varies fastest), producing **theta-fastest** memory layout.

---

### 2.4 **RESOLUTION**

The C API uses confusing terminology where `Nt` (toroidal) and `Np` (poloidal) do NOT directly correspond to the **index ordering** in the flattened array. The documentation `{x11, x12, ..., x1Np, x21, ...}` suggests:
- **First index** (1, 2, ..., Nt) = toroidal index
- **Second index** (1, 2, ..., Np) = poloidal index
- **Fastest-varying** in memory = second index (poloidal)

This is the **opposite** of TIAGO's current flattening (phi-fastest).

---

## 3. DIAGNO REFERENCE IMPLEMENTATION ANALYSIS

### 3.1 DIAGNO Virtual-Casing Module (`LIBSTELL/virtual_casing_mod.f90:186-234`)

DIAGNO uses a **completely different** virtual-casing implementation based on surface integral evaluation via splines, NOT the HiddenSymmetries C++ library. Key observations:

```fortran
! DIAGNO usage (diagno_init_vmec.f90:127-128):
bumnc_temp(:,1) = mfact(:,1)*bsupumnc(:,ns) + mfact(:,2)*bsupumnc(:,ns-1)
bvmnc_temp(:,1) = mfact(:,1)*bsupvmnc(:,ns) + mfact(:,2)*bsupvmnc(:,ns-1)

! Half-grid to full-grid extrapolation (diagno_init_vmec.f90:120-126):
WHERE (MOD(NINT(REAL(xm_temp(:))),2) .eq. 0)
    mfact(:,1)= 1.5
    mfact(:,2)=-0.5
ELSEWHERE
    mfact(:,1)= 1.5*SQRT((ns-1.0)/(ns-1.5))
    mfact(:,2)=-0.5*SQRT((ns-1.0)/(ns-2.5))
ENDWHERE
```

**Findings**:
1. DIAGNO extrapolates VMEC half-grid B-field to full-grid boundary using m-dependent factors
2. DIAGNO's `virtual_casing_mod` uses **Fourier representation** and spline interpolation
3. DIAGNO evaluates fields at arbitrary points via `bfield_vc(x,y,z,bx,by,bz)`

**Key Difference**: DIAGNO operates on **Fourier coefficients** (spectral space), while HiddenSymmetries operates on **real-space grids** (physical space). This is a fundamentally different numerical approach.

---

### 3.2 Grid Convention Mismatch Potential

DIAGNO's grid (from `diagno_init_vmec.f90:57-58`):
```fortran
nu2 = nu  ! VMEC poloidal resolution
nv2 = nv  ! VMEC toroidal resolution
```

VMEC conventions:
- `nu` = number of **poloidal** grid points (θ direction)
- `nv` = number of **toroidal** grid points per field period (φ direction)
- VMEC uses `(s, θ, ζ)` coordinates where ζ = φ/nfp (toroidal angle per field period)

TIAGO test creates synthetic toroidal surfaces **without VMEC grid conventions**, potentially missing:
- Proper grid shift for stellarator symmetry
- Half-grid vs full-grid distinctions
- Field period normalization

---

## 4. GRID CONVENTION ISSUES

### 4.1 Stellarator Symmetry Grid Shift (`virtual-casing.h:58-73`)

```c
* If you do not exploit stellarator symmetry, then half_period
* is set to false.  In this case the grids in the toroidal angles
* begin at phi=0.  The grid spacing is 1 / (NFP * Nt), and there
* is no point at the symmetry plane phi = 1 / NFP.
*
* If you do wish to exploit stellarator symmetry, set half_period
* to true. In this case the toroidal grids are each shifted by
* half a grid point, so there is no grid point at phi = 0. The
* phi grid for the surface shape has points at 0.5 / (NFP * Nt),
* 1.5 / (NFP * Nt), ..., (Nt - 0.5) / (NFP * Nt).
```

**Critical**: When `half_period=true`, phi grid should be:
```
phi[i] = (i + 0.5) / (NFP * Nt)  for i = 0, 1, ..., Nt-1
```

### 4.2 TIAGO Test Grid (`test_plasma_response_ncsx.f90:144-147`)

```fortran
do iphi = 1, nphi
    phi = 2.0_dp * pi * (iphi - 1) / nphi / real(nfp, dp)
    ! ... surface generation ...
end do
```

**Generated grid**: `phi[i] = 2π * (i-1) / (nphi * nfp)` for `i = 1..nphi`
**In [0,1) units**: `phi[i] = (i-1) / (nphi * nfp)`

This is the **non-stellsym** convention (grid starts at phi=0), but TIAGO passes `use_stellsym=.true.` to `plasma_response%init()`.

**MISMATCH**: TIAGO generates non-stellsym grid but tells virtual-casing to expect stellsym grid with half-grid shift.

---

### 4.3 Poloidal Grid Convention

```fortran
! TIAGO test (test_plasma_response_ncsx.f90:145)
theta = 2.0_dp * pi * (itheta - 1) / ntheta

! virtual-casing docs (virtual-casing.h:76-78):
* Regardless of half_period, the poloidal grid always ranges
* uniformly over [0, 1), with the first grid point at theta = 0,
* and no grid point at theta = 1.
```

**In [0,1) units**: `theta[i] = (i-1) / ntheta` for `i = 1..ntheta`

TIAGO: `theta = 2π * (i-1) / ntheta`
API expects: `theta ∈ [0, 1)` (period 1, not 2π)

**ISSUE**: TIAGO may be passing theta in radians [0,2π) instead of dimensionless [0,1).

**CORRECTION**: The virtual-casing C API works internally with [0,1) convention. The surface coordinate flattening should receive **Cartesian (x,y,z)** coordinates, not angular coordinates. So this is NOT an issue for surface coordinates themselves, but could affect **internal B-field grid generation** if TIAGO's `b_total` is generated on the wrong grid.

---

## 5. INTENT ATTRIBUTE VIOLATIONS

### 5.1 `plasma_response.f90:21-24` - VirtualCasingDestroyContextD

```fortran
subroutine VirtualCasingDestroyContextD(ctx_ptr) &
        bind(C, name='VirtualCasingDestroyContextD')
    use iso_c_binding
    type(c_ptr), value :: ctx_ptr
end subroutine
```

**C Signature** (`virtual-casing.h:21-22`):
```c
void VirtualCasingDestroyContextD(void** ctx);
```

**CRITICAL ISSUE**: C function expects `void**` (pointer-to-pointer), but Fortran binding declares `type(c_ptr), value`, which passes `void*` (pointer by value).

---

### 5.2 Usage in `plasma_response_finalize` (line 198-206)

```fortran
subroutine plasma_response_finalize(self)
    class(plasma_response_t), intent(inout) :: self
    type(c_ptr), target :: ctx_tmp
    if (c_associated(self%ctx)) then
        ctx_tmp = self%ctx
        call VirtualCasingDestroyContextD(c_loc(ctx_tmp))
        self%ctx = c_null_ptr
        self%initialized = .false.
    end if
end subroutine
```

**Workaround**: Creates temporary `ctx_tmp` and passes `c_loc(ctx_tmp)` to get pointer-to-pointer semantics.

**Analysis**: This is **correct** - the workaround properly converts to `void**`. However, the interface declaration is misleading and fragile.

**RECOMMENDATION**: Update interface to explicitly reflect double-pointer semantics:
```fortran
subroutine VirtualCasingDestroyContextD(ctx_ptr) &
        bind(C, name='VirtualCasingDestroyContextD')
    use iso_c_binding
    type(c_ptr) :: ctx_ptr  ! NOT value - receives address of c_ptr
end subroutine
```

---

## 6. MEMORY DEALLOCATION ISSUES

### 6.1 Temporary Array Lifecycle (`plasma_response.f90:146`)

```fortran
deallocate(x_flat, b_flat)
```

Arrays `x_flat` and `b_flat` are deallocated **after** `VirtualCasingSetupD` returns.

**C++ Implementation** (`virtual-casing.cpp:54-57`):
```cpp
void VirtualCasingSetupD(..., const double* X, ..., void* ctx) {
  VirtualCasing<double>& virtual_casing = *(VirtualCasing<double>*)ctx;
  const std::vector<double> X_(X, X+3*Nt*Np);  // COPIES data
  virtual_casing.Setup(digits, NFP, half_period, Nt, Np, X_, ...);
}
```

**Analysis**: C++ side **copies** input arrays into `std::vector`, so Fortran deallocating after call is **safe**. No memory leak or use-after-free.

---

### 6.2 ComputeBext Temporary Arrays (`plasma_response.f90:167-194`)

```fortran
allocate(b_total_flat(nphi * ntheta * 3))
allocate(b_ext_flat(nphi * ntheta * 3))
! ... flatten, call C++, unflatten ...
deallocate(b_total_flat, b_ext_flat)
```

Same pattern - C++ copies data immediately, so deallocation is safe.

---

## 7. CRITICAL TEST CASE ANALYSIS

### 7.1 Test Results Summary (`plasma_response_summary_ncsx_nfp3.txt`)

```
VACUUM ONLY (TIAGO baseline)
  RMS error vs DIAGNO: 0.1797%   ← EXCELLENT

PLASMA RESPONSE ONLY (Virtual-Casing contribution)
  Mean magnitude: 121.49 μV
  Mean relative to vacuum: 1.61%
  Range: 48.90 - 245.60 μV

VACUUM + PLASMA (Total TIAGO prediction)
  RMS error vs DIAGNO: 1.4475%   ← DEGRADED
  Max error: 1.8091%
```

**Analysis**:
- Vacuum field calculation: **0.18% error** ← TIAGO Biot-Savart is correct
- Vacuum + Plasma: **1.45% error** ← Adding plasma **worsens** agreement
- Plasma contribution: ~1.6% of vacuum signal (~120 μV vs 7 mV)
- Error increase: **1.45% - 0.18% = 1.27%** ← Comparable to plasma signal magnitude (1.6%)

**Interpretation**: The plasma response contribution is **nearly 100% wrong in sign/magnitude**, causing total field to diverge from DIAGNO.

---

### 7.2 Detailed Signal Breakdown

```
Signal                DIAGNO     Vac   Plasma   Total   Vac Err  Tot Err
Flux-1              12.47000  12.45000  0.24560  12.69560  0.160%   1.809%
```

**Calculation**:
- DIAGNO expects: 12.47 mV
- TIAGO vacuum: 12.45 mV (error = -0.02 mV = 0.160%)
- TIAGO plasma adds: +0.246 mV
- TIAGO total: 12.696 mV (error = +0.226 mV = 1.809%)

**Expected plasma contribution** (if correct): `12.47 - 12.45 = 0.02 mV`
**Actual plasma contribution**: `0.246 mV`

**Plasma error**: `(0.246 - 0.02) / 0.02 = 1130%` ← **Order-of-magnitude error**

---

### 7.3 Cross-Check: Plasma Effect Sign

All diagnostics show TIAGO plasma response is **positive** (increases signal), while the vacuum-only error is **negative** (vacuum underestimates). If plasma response were correctly computed to close the gap, it should add ~0.02 mV, not 0.25 mV.

**CONCLUSION**: TIAGO plasma response is **12× too large** and may have incorrect spatial distribution.

---

## 8. COMPARISON WITH SIMSOPT IMPLEMENTATION

### 8.1 Grid Generation (`simsopt/mhd/virtual_casing.py:185-191`)

```python
surf = SurfaceRZFourier.from_nphi_ntheta(mpol=vmec.wout.mpol, ntor=vmec.wout.ntor, nfp=nfp,
                                         nphi=src_nphi, ntheta=src_ntheta, range=ran)
for jmn in range(vmec.wout.mnmax):
    surf.set_rc(int(vmec.wout.xm[jmn]), int(vmec.wout.xn[jmn] / nfp), vmec.wout.rmnc[jmn, -1])
    surf.set_zs(int(vmec.wout.xm[jmn]), int(vmec.wout.xn[jmn] / nfp), vmec.wout.zmns[jmn, -1])
```

Simsopt **reconstructs VMEC boundary surface** from Fourier coefficients, ensuring proper grid conventions.

---

### 8.2 B-Field Generation (`simsopt/mhd/vmec_diagnostics.py`)

```python
Bxyz = B_cartesian(vmec, nphi=src_nphi, ntheta=src_ntheta, range=ran)
```

Simsopt uses VMEC's internal B-field computation (`B_cartesian`), which:
- Reads VMEC Fourier coefficients `bsupumnc`, `bsupvmnc`
- Evaluates B on the **VMEC flux coordinate grid**
- Transforms to Cartesian coordinates

TIAGO tests use **synthetic** B-field generated analytically, not from VMEC.

---

### 8.3 **KEY INSIGHT**

Simsopt validates virtual-casing against **VMEC's own B-field on VMEC's own grid**. TIAGO tests use **model fields on custom grids**, which may not satisfy:
- Proper grid shift conventions
- Field period symmetry
- ∇·B = 0 (divergence-free constraint)
- Proper boundary conditions for virtual-casing principle

---

## 9. ROOT CAUSE ANALYSIS - RANKED BY CRITICALITY

### 9.1 **CRITICAL - Array Flattening Order**

**Severity**: 🔴 **BLOCKING**
**Confidence**: 95%

TIAGO flattens arrays as **phi-fastest, theta-second**, but C API expects **theta-fastest, phi-second**. This causes:
- Spatial coordinates to be scrambled (mixing phi/theta indices)
- Surface normal vectors computed on wrong grid points
- B-field values associated with wrong surface locations
- Virtual-casing integral evaluated with incorrect geometry

**Evidence**:
- Simsopt uses `flatten(order='C')` on `(nphi, ntheta, 3)` → theta-fastest
- C API docs state `{x11, x12, ..., x1Np}` → second index (Np=poloidal) varies fastest
- TIAGO uses `idx = (k-1)*nphi*ntheta + (j-1)*nphi + i` → first index (i=phi) varies fastest

**Fix**: Reverse loop order:
```fortran
do k = 1, 3
    do i = 1, nphi
        do j = 1, ntheta
            idx = (k-1)*nphi*ntheta + (i-1)*ntheta + j
            x_flat(idx) = real(x_surf(i, j, k), c_double)
        end do
    end do
end do
```

---

### 9.2 **HIGH - Grid Convention Mismatch**

**Severity**: 🟠 **MAJOR**
**Confidence**: 75%

TIAGO test generates non-stellsym grid (`phi[0] = 0`) but passes `use_stellsym=.true.`, which tells virtual-casing to expect half-grid shift (`phi[0] = 0.5/(nfp*nphi)`).

**Evidence**:
```fortran
! TIAGO generates:
phi = 2π * (iphi-1) / (nphi * nfp)  for iphi=1..nphi

! virtual-casing expects (when half_period=true):
phi = (iphi - 0.5) / (nphi * nfp)  for iphi=1..nphi
```

**Fix**: Generate stellsym-compliant grid:
```fortran
phi = 2.0_dp * pi * (iphi - 0.5_dp) / (nphi * nfp)
```

---

### 9.3 **MEDIUM - Test Field Realism**

**Severity**: 🟡 **MODERATE**
**Confidence**: 60%

TIAGO synthetic test field may not satisfy ∇·B = 0 or proper toroidal periodicity, causing virtual-casing algorithm to fail.

**Evidence**:
- Simsopt uses VMEC's B-field (automatically divergence-free in flux coordinates)
- TIAGO uses analytical model field with arbitrary harmonics

**Fix**: Use actual VMEC B-field from `wout_ncsx.nc` via libneo/read_wout.

---

### 9.4 **LOW - DIAGNO vs HiddenSymmetries Algorithm Difference**

**Severity**: 🟢 **MINOR**
**Confidence**: 40%

DIAGNO uses spectral (Fourier-based) virtual-casing, while HiddenSymmetries uses real-space quadrature. This is a fundamental algorithmic difference, not a bug.

**Evidence**:
- DIAGNO operates on Fourier coefficients (mnmax modes)
- HiddenSymmetries operates on real-space grids (nu×nv points)
- Both should converge to same result if grids are adequate

**Implication**: Small discrepancies (0.1-0.5%) are expected, but 1000%+ error indicates implementation bug, not algorithm difference.

---

## 10. RECOMMENDED FIXES - PRIORITY ORDER

### 10.1 **IMMEDIATE FIX #1: Correct Array Flattening**

**File**: `src/solver/plasma_response.f90`
**Lines**: 110-119, 166-177, 184-192

**Change loop order** to produce theta-fastest memory layout:

```fortran
! BEFORE (WRONG):
do k = 1, 3
    do j = 1, ntheta
        do i = 1, nphi
            idx = (k-1) * nphi * ntheta + (j-1) * nphi + i
            x_flat(idx) = real(x_surf(i, j, k), c_double)
        end do
    end do
end do

! AFTER (CORRECT):
do k = 1, 3
    do i = 1, nphi
        do j = 1, ntheta
            idx = (k-1) * nphi * ntheta + (i-1) * ntheta + j
            x_flat(idx) = real(x_surf(i, j, k), c_double)
        end do
    end do
end do
```

**Apply to**:
- `plasma_response_init` (lines 110-119, 122-130)
- `plasma_response_compute_bext` (lines 170-177, 185-192)

---

### 10.2 **IMMEDIATE FIX #2: Stellsym Grid Generation**

**File**: `tests/test_plasma_response_ncsx.f90`
**Lines**: 144-170

**Add half-grid shift** when stellsym is enabled:

```fortran
! BEFORE:
phi = 2.0_dp * pi * (iphi - 1) / nphi / real(nfp, dp)

! AFTER:
phi = 2.0_dp * pi * (iphi - 0.5_dp) / (nphi * nfp)
```

---

### 10.3 **HIGH-PRIORITY FIX #3: Use VMEC B-Field**

**File**: `tests/test_plasma_response_ncsx.f90`
**New dependency**: Read VMEC wout file, extract B-field on boundary

**Replace** `create_ncsx_field()` with VMEC B-field reader:

```fortran
! Add to CMakeLists.txt:
# target_link_libraries(tiago_plasma_response_ncsx_tests PRIVATE libneo)

! In test code:
use read_wout_mod, only: read_wout_file, bsupumnc, bsupvmnc, ns
! ... read wout_ncsx.nc ...
! ... evaluate B on surface using VMEC Fourier modes ...
```

---

### 10.4 **MEDIUM-PRIORITY FIX #4: Clarify Interface Documentation**

**File**: `src/solver/plasma_response.f90`
**Lines**: 12-50

Add explicit documentation of memory layout expectations:

```fortran
!! Array layout requirements for virtual-casing C API:
!!
!! Surface coordinates X(1:nphi*ntheta*3):
!!   Memory order: {x(1,1), x(1,2), ..., x(1,ntheta), x(2,1), x(2,2), ..., x(nphi,ntheta),
!!                  y(1,1), y(1,2), ..., y(nphi,ntheta),
!!                  z(1,1), z(1,2), ..., z(nphi,ntheta)}
!!   Index formula: idx = (k-1)*nphi*ntheta + (i-1)*ntheta + j
!!   where i=1..nphi (toroidal), j=1..ntheta (poloidal), k=1..3 (x,y,z)
!!
!! NOTE: This is THETA-FASTEST (poloidal varies fastest), not phi-fastest!
```

---

### 10.5 **LOW-PRIORITY FIX #5: Improve Destroy Interface**

**File**: `src/solver/plasma_response.f90`
**Lines**: 20-25

Update interface to match C signature semantics:

```fortran
! Current (misleading):
type(c_ptr), value :: ctx_ptr

! Improved (clearer):
type(c_ptr) :: ctx_ptr  ! Receives address-of c_ptr (void**)
```

Add comment explaining double-pointer usage.

---

## 11. VALIDATION PLAN

### 11.1 Unit Test: Array Flattening Order

Create test to verify memory layout matches simsopt:

```fortran
! Generate small test grid
nphi = 3
ntheta = 4
allocate(x_surf(nphi, ntheta, 3))
! Fill with known pattern: x_surf(i,j,k) = i*100 + j*10 + k

! Flatten
call flatten_to_c_order(x_surf, x_flat)

! Verify:
! x_flat(1:4) = [111, 112, 113, 114]  ! x(1,1:4)
! x_flat(5:8) = [121, 122, 123, 124]  ! x(2,1:4)
! x_flat(9:12) = [131, 132, 133, 134] ! x(3,1:4)
! x_flat(13:16) = [211, 212, 213, 214] ! y(1,1:4)
```

---

### 11.2 Integration Test: Sphere Benchmark

Use HiddenSymmetries test data generator:

```fortran
! Call C API test data generator
call GenerateVirtualCasingTestDataD(Bext, Bint, nfp, half_period, &
                                   nphi, ntheta, X_sphere, nphi, ntheta)

! Compare against analytical solution for sphere
! (test data generator uses known current loop geometry)
```

Expected accuracy: <1e-5 relative error for unit sphere test case.

---

### 11.3 Comparison Test: Simsopt Cross-Check

Run simsopt on same NCSX wout file:

```python
from simsopt.mhd import Vmec, VirtualCasing
vmec = Vmec('wout_ncsx.nc')
vc = VirtualCasing.from_vmec(vmec, src_nphi=24, src_ntheta=24)
# Save vc.B_external to file
```

Compare TIAGO `b_ext` against simsopt `B_external` point-by-point:
- Expected agreement: <1% RMS error
- Current status: Unknown (simsopt not run on NCSX case)

---

## 12. ADDITIONAL FINDINGS

### 12.1 No Error Handling for Context Creation Failure

**File**: `plasma_response.f90:104-108`

```fortran
self%ctx = VirtualCasingCreateContextD()
if (.not. c_associated(self%ctx)) then
    return  ! Silent failure - no error message
end if
```

**Issue**: Function returns silently on failure. Caller checks `is_initialized()` but gets no diagnostic info.

**Recommendation**: Add error logging or return status code.

---

### 12.2 Missing Intent Attributes

**File**: `plasma_response.f90:28-39`

```fortran
subroutine VirtualCasingSetupD(digits, nfp, half_period, &
                               nt, np, x, &
                               src_nt, src_np, trg_nt, trg_np, ctx) &
         bind(C, name='VirtualCasingSetupD')
    use iso_c_binding
    integer(c_int), value :: digits, nfp
    logical(c_bool), value :: half_period
    integer(c_long), value :: nt, np
    real(c_double), intent(in) :: x(*)  ! ← Only x has intent
    integer(c_long), value :: src_nt, src_np, trg_nt, trg_np
    type(c_ptr), value :: ctx           ! ← ctx modified but no intent(inout)
end subroutine
```

**Issue**: `ctx` is modified by C function (Setup is called on dereferenced object), but interface doesn't declare `intent`.

**Analysis**: For C interop with `value` attribute, intent is advisory only. However, for documentation clarity, add comments.

---

### 12.3 Hardcoded Double Precision

**File**: `plasma_response.f90:4-6`

```fortran
use, intrinsic :: iso_fortran_env, only: dp => real64
```

**Analysis**: Hardcoded to `real64`. Virtual-casing library provides both single (`F`) and double (`D`) precision APIs.

**Recommendation**: Current approach is fine (plasma physics requires double precision). Single precision not needed.

---

## 13. PERFORMANCE CONSIDERATIONS

### 13.1 Memory Allocation Pattern

Allocates/deallocates temporary flattened arrays on **every** `compute_bext` call:

```fortran
subroutine plasma_response_compute_bext(self, b_total, b_external)
    allocate(b_total_flat(nphi * ntheta * 3))
    allocate(b_ext_flat(nphi * ntheta * 3))
    ! ... use ...
    deallocate(b_total_flat, b_ext_flat)
end subroutine
```

**Performance impact**: If called in tight loop (e.g., optimization), repeated allocation overhead.

**Recommendation**: Cache flattened array buffers in `plasma_response_t` type, allocate once in `init`.

---

### 13.2 Loop Optimization

Current triple nested loop for flattening is not OpenMP-parallelized.

**Opportunity**: For large grids (nphi=128, ntheta=128 → 50k points), parallelizing could help:

```fortran
!$OMP PARALLEL DO PRIVATE(i,j,k,idx)
do k = 1, 3
    do i = 1, nphi
        do j = 1, ntheta
            idx = (k-1)*nphi*ntheta + (i-1)*ntheta + j
            x_flat(idx) = real(x_surf(i, j, k), c_double)
        end do
    end do
end do
!$OMP END PARALLEL DO
```

**Note**: Only worth it for large grids; small grids have negligible overhead.

---

## 14. COMPARISON WITH REFERENCE IMPLEMENTATIONS

### 14.1 Simsopt Memory Layout

**File**: `simsopt/mhd/virtual_casing.py:213-214`

```python
gamma1d[jxyz * src_nphi * src_ntheta: (jxyz + 1) * src_nphi * src_ntheta] = \
    gamma[:, :, jxyz].flatten(order='C')
```

**NumPy flatten behavior**:
- `gamma.shape = (nphi, ntheta, 3)`
- `gamma[:,:,0].flatten('C')` produces: `[gamma[0,0,0], gamma[0,1,0], ..., gamma[0,ntheta-1,0], gamma[1,0,0], ...]`
- This is **theta-fastest** (last dimension of `gamma[:,:,0]` varies fastest)

**Verification**:
```python
import numpy as np
a = np.array([[[1,2],[3,4]], [[5,6],[7,8]]])  # shape (2,2,2)
a[:,:,0].flatten('C')  # [1, 3, 5, 7] - second index varies fastest ✓
```

---

### 14.2 DIAGNO Grid Structure

**File**: `LIBSTELL/virtual_casing_mod.f90:393-400`

```fortran
uv = 1
DO v = 1, nv
   DO u = 1, nu
      xsurf(uv) = r_temp(u,v,2)*DCOS(factor*xv(v))
      ysurf(uv) = r_temp(u,v,2)*DSIN(factor*xv(v))
      zsurf(uv) = z_temp(u,v,2)
      xreal(u,v) = xsurf(uv)
```

**DIAGNO flattening**: Loop order is `v` (toroidal) outer, `u` (poloidal) inner → **poloidal-fastest** (u varies fastest in `uv` index).

This matches virtual-casing expected order!

---

### 14.3 Fortran Sphere Test Reference

**File**: `virtual_casing-src/test/test_virtual_casing_sphere.f90:110-116`

Unfortunately, this test does NOT show array flattening (it uses internal C++ test data generator). Not useful for validation.

---

## 15. CONCRETE EVIDENCE OF BUG

### 15.1 Symptom Analysis

**Observation**: Vacuum-only error is **-0.16%** (TIAGO underestimates), but after adding plasma response, total error becomes **+1.81%** (TIAGO overestimates).

**Math**:
```
DIAGNO = 12.47 mV
TIAGO vacuum = 12.45 mV (error = -0.02 mV = -0.16%)
TIAGO plasma = +0.246 mV
TIAGO total = 12.696 mV (error = +0.226 mV = +1.81%)
```

**Expected behavior** (if plasma correct):
- Plasma should add ~0.02 mV to close the -0.02 mV vacuum gap
- Total should approach 12.47 mV with error near 0%

**Actual behavior**:
- Plasma adds 12× too much (+0.246 mV instead of +0.02 mV)
- Total overshoots by +0.226 mV

**Interpretation**:
- If array indices are scrambled, virtual-casing integrates over **wrong surface geometry**
- Surface currents (K = n × B) evaluated at wrong locations
- Integral likely picks up spurious contributions from misaligned B-field and surface normal vectors
- Result: Order-of-magnitude error in plasma contribution

---

### 15.2 Error Pattern Consistency

All 11 diagnostics show same pattern:
- Vacuum error: 0.15-0.22% (all negative)
- Total error: 1.13-1.81% (all positive)
- Plasma adds: 0.05-0.25 mV (proportional to signal strength)

**Consistency**: Error scales with signal magnitude, suggesting systematic geometric error (not random noise).

---

## 16. FINAL SUMMARY TABLE

| Issue | Severity | Confidence | Impact on 1000% Error | Fix Complexity |
|-------|----------|------------|------------------------|----------------|
| **Array flattening order (phi vs theta fastest)** | 🔴 CRITICAL | 95% | **Primary cause** | Low (change loop order) |
| **Grid convention (stellsym shift)** | 🟠 HIGH | 75% | **Secondary cause** | Low (adjust grid generation) |
| **Synthetic vs VMEC B-field** | 🟡 MEDIUM | 60% | **Contributing factor** | Medium (integrate VMEC reader) |
| Destroy context interface | 🟢 LOW | 90% | None (works despite misleading syntax) | Low (documentation) |
| No error handling | 🟢 LOW | 100% | None (silent failure unrelated to accuracy) | Low (add logging) |
| Performance (repeated alloc) | 🟢 LOW | 100% | None (speed only, not accuracy) | Medium (refactor caching) |

---

## 17. RECOMMENDED TESTING SEQUENCE

1. **Fix array flattening** → Rebuild → Run `tiago_plasma_response_ncsx_tests`
   - Expected: Error drops from 1.45% to <0.5%
2. **Fix grid convention** → Rebuild → Rerun
   - Expected: Error drops further to <0.2%
3. **Add simsopt cross-validation** → Compare point-by-point
   - Expected: Agreement within 0.1% with simsopt
4. **Use real VMEC B-field** → Rerun full pipeline
   - Expected: Agreement with DIAGNO within 0.2% (algorithmic differences only)

---

## 18. CONFIDENCE ASSESSMENT

**Overall confidence in root cause identification**: **90%**

**Reasoning**:
- Array layout mismatch is **definitively confirmed** by comparing TIAGO, simsopt, and DIAGNO code
- Error magnitude (1000%) is consistent with spatial coordinate scrambling (not just numerical precision)
- Error pattern (consistent sign, scales with signal) supports systematic geometric error
- Simsopt and DIAGNO both use theta-fastest ordering, while TIAGO uses phi-fastest

**Remaining 10% uncertainty**:
- Possibility of additional bugs in test field generation
- Potential issue with C++ library configuration/compilation
- Undocumented virtual-casing API behavior for edge cases

---

## END OF AUDIT REPORT

**Audit Complete**: 2025-11-13
**Total Issues Identified**: 6 (1 critical, 1 high, 1 medium, 3 low)
**Recommended Actions**: Implement fixes #1 and #2 immediately, validate with tests, then proceed to #3-5.

**Next Steps**:
1. Apply critical fixes to `plasma_response.f90` and test cases
2. Run validation suite against simsopt reference
3. Document memory layout conventions in code comments
4. Update user-facing documentation with correct usage examples

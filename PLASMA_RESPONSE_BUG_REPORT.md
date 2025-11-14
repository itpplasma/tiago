# Plasma Response Virtual Casing Bug Report

## Executive Summary

Investigation of the plasma response implementation revealed TWO confirmed bugs and ONE fundamental design limitation.

## Bugs Fixed

### Bug 1: Incorrect Normalization (CRITICAL)
**File**: `src/solver/surface_biot_savart.f90`, line 155
**Was**: `self%norm_const = 0.25_dp / (pi * real(phi_count * theta_count, dp))`
**Should be**: `self%norm_const = 0.25_dp / pi`

**Impact**: Results were too small by a factor of `nphi × ntheta`. For a 360×360 grid, this is a factor of 129,600×!

**Explanation**: The virtual casing integral is:
```
B_plasma(r) = ∫∫ [K × (r-r') + Bn (r-r')] / |r-r'|³ dS / (4π)
```

The normalization `1/(4π)` comes from the Biot-Savart law and does NOT depend on grid resolution. The grid spacing `dφ` and `dθ` already appears in the Gaussian quadrature weights.

**Reference**: STELLOPT `virtual_casing_mod.f90` line 239: `norm_nag = nfp/(4*π)` for adaptive integration over [0,1]×[0,1] normalized coordinates. For our case integrating over [0,2π]×[0,2π] in actual angles with full-torus data, the normalization is `1/(4π)`.

### Bug 2: Wrong Normal Direction
**File**: `src/solver/surface_biot_savart.f90`, line 118
**Was**: `normal_vec = cross_product(dx_dtheta, dx_dphi)`
**Should be**: `normal_vec = cross_product(dx_dphi, dx_dtheta)`

**Impact**: Normal vector pointed inward instead of outward, giving wrong sign for some field components.

**Explanation**: Virtual casing requires the OUTWARD normal from plasma to vacuum. For stellarator coordinates (φ=toroidal, θ=poloidal):
- `∂x/∂φ × ∂x/∂θ` points OUTWARD ✓
- `∂x/∂θ × ∂x/∂φ` points INWARD ✗

## Fundamental Limitation (NOT FIXED)

### Issue: Assumes Full-Torus Input Data
**File**: `src/solver/surface_biot_savart.f90`, lines 147-150

The code hardcodes spline domains as [0,2π]×[0,2π] regardless of the actual angular range of input data.

**Test Case Failure**: The `simsopt_li383` comparison uses data spanning only [0,1.01] rad × [0,2π] rad ≈ [0,58°] × [0,360°] - less than half a field period! The splines extrapolate wildly beyond this range, producing nonsense results with RMS errors >600×.

**Working Test Case**: The NCSX test (`test_plasma_response_ncsx.f90`) works because it generates full-torus data [0,2π]×[0,2π] directly from VMEC (line 62: `phi_edge = 2π×(iphi-1)/nphi`).

**Proper Solutions** (future work):
1. **Auto-detect angular range**: Compute φ and θ ranges from input Cartesian coordinates and only integrate over the actual data domain
2. **Stellarator symmetry extension**: Use nfp to replicate one field period data across the full torus before splining
3. **Documentation**: Require users to provide full-torus data and document this requirement

## Verification

### Before Fixes
- `test_plasma_response_ncsx`: Values ~10⁻⁶ T (should be ~1 T) - factor of 10⁶× too small
- `simsopt_li383` comparison: RMS error 0.157 T with chaotic scatter plot

### After Normalization Fix Only
- `test_tiago_with_simsopt_grid`: RMS error 37 T, relative error 611× (but wrong due to data domain issue)

### After Both Bug Fixes
- `test_plasma_response_ncsx`: Passes (no reference comparison, just sanity checks)
- `simsopt_li383` comparison: Still fails due to partial data domain issue

## Comparison with STELLOPT

**Key Formula** (`virtual_casing_mod.f90` lines 1498-1501):
```fortran
gf3  = norm_nag*gf*gf*gf
f(1) = (ky*(z_nag-zs)-kz*(y_nag-ys)+bn*(x_nag-xs))*gf3
f(2) = (kz*(x_nag-xs)-kx*(z_nag-zs)+bn*(y_nag-ys))*gf3
f(3) = (kx*(y_nag-ys)-ky*(x_nag-xs)+bn*(z_nag-zs))*gf3
```

This is `[K × r + Bn × r] / |r|³ × norm_nag` where `r = point - surface_point`.

**Our Implementation** (lines 295-298) matches this exactly (after bug fixes).

**Key Difference**: STELLOPT precomputes K with unnormalized jacobian (lines 468-470):
```fortran
kxreal(u,v) = -(bys * snz - bzs * sny)  ! sn = |∂x/∂u × ∂x/∂v|
```

We compute K with normalized normal and apply jacobian during integration (line 293):
```fortran
weight = self%norm_const * jacobian * weight_phi * weight_theta
```

Both approaches are mathematically equivalent.

## Files Modified

- `/home/ert/code/tiago/src/solver/surface_biot_savart.f90`
  - Line 118: Fixed normal direction
  - Line 155: Fixed normalization constant
  - Line 55: Added `nfp_real` declaration (currently unused, kept for future fix)

## Testing Commands

```bash
# NCSX test (works with full-torus VMEC data)
/home/ert/code/tiago/build/tiago_plasma_response_tests tests/cases/ncsx_nfp3/wout_ncsx.nc

# Simsopt comparison (fails due to partial data domain)
/home/ert/code/tiago/build/test_tiago_with_simsopt_grid \
  build/tests/output/simsopt_li383/gamma.csv \
  build/tests/output/simsopt_li383/b_total.csv \
  build/tests/output/tiago_bext.csv \
  build/tests/output/simsopt_li383/b_external.csv
```

## Recommendations

1. **Immediate**: Document that `surface_biot_savart_t` requires full-torus data [0,2π]×[0,2π]
2. **Short-term**: Add input validation to detect partial data and error with clear message
3. **Long-term**: Implement stellarator symmetry extension to handle partial period data
4. **Testing**: Create a full-torus simsopt comparison test that can validate the fixes

## References

- STELLOPT: `/home/ert/code/external/STELLOPT/LIBSTELL/Sources/Modules/virtual_casing_mod.f90`
- Simsopt: `/home/ert/code/external/simsopt/src/simsopt/mhd/virtual_casing.py`
- Virtual casing theory: Hanson & Hirshman, Phys. Plasmas 9, 4410 (2002)
